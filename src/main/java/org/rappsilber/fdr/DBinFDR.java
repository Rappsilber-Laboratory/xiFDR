/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr;

import java.awt.GraphicsEnvironment;
import java.io.File;
import org.rappsilber.fdr.result.FDRResult;
import java.io.FileNotFoundException;
import java.lang.invoke.MethodHandles;
import java.sql.Array;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.text.DecimalFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Timer;
import java.util.TimerTask;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;
import org.rappsilber.config.LocalProperties;
import org.rappsilber.fdr.dataimport.Xi2Xi1Config;
import rappsilber.config.DBRunConfig;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.fdr.entities.PeptidePair;
import org.rappsilber.fdr.entities.Protein;
import org.rappsilber.fdr.entities.ProteinGroupLink;
import org.rappsilber.fdr.entities.ProteinGroupPair;
import org.rappsilber.fdr.entities.ProteinGroup;
import org.rappsilber.fdr.gui.components.MZIdentMLOwnerGUI;
import org.rappsilber.fdr.utils.CalculateWriteUpdate;
import org.rappsilber.fdr.utils.MZIdentMLExport;
import org.rappsilber.fdr.utils.MZIdentMLOwner;
import org.rappsilber.fdr.utils.MaximisingStatus;
import rappsilber.ms.crosslinker.CrossLinker;
import rappsilber.ms.lookup.peptides.PeptideTree;
import rappsilber.ms.sequence.AminoAcid;
import rappsilber.ms.sequence.Peptide;
import rappsilber.ms.sequence.Sequence;
import rappsilber.ms.sequence.digest.Digestion;
import rappsilber.utils.CountOccurence;
import org.rappsilber.utils.IntArrayList;
import org.rappsilber.utils.RArrayUtils;
import org.rappsilber.utils.SelfAddHashSet;
import org.rappsilber.utils.Version;
import rappsilber.config.DBConnectionConfig;
import rappsilber.config.RunConfig;
import rappsilber.gui.components.db.DatabaseProvider;
import rappsilber.ms.sequence.AminoModification;
import rappsilber.ms.statistics.utils.UpdateableLong;

/**
 *
 * @author lfischer
 */
public class DBinFDR extends org.rappsilber.fdr.OfflineFDR implements XiInFDR {

    public static final String subscoreroperty = "xiFDR.FILTER_SUBSCORES";
    private Connection m_db_connection;
    private UpdateableLong       m_db_last_used = new UpdateableLong(0);
    private AtomicBoolean       keepConnection  = new AtomicBoolean(false);
    private Timer      m_db_autoclosetimer;
    private PreparedStatement updateValidateOverWrite;
    private PreparedStatement updateValidateNonOverWrite;
//    private PreparedStatement updateValidateOverWriteSetFDR;
//    private PreparedStatement updateValidateNonOverWriteSetFDR;
//    private PreparedStatement updateSetFDR;
    private PreparedStatement updateSetConfidence;
//    private PreparedStatement updateValidateOverWriteBatch;
//    private PreparedStatement updateValidateNonOverWriteBatch;
//    private PreparedStatement updateValidateOverWriteSetFDRBatch;
//    private PreparedStatement updateValidateNonOverWriteSetFDRBatch;
//    private PreparedStatement updateSetFDRBatch;
//    private PreparedStatement updateSetConfidenceBatch;
    private ArrayList<String> m_search_ids = new ArrayList<>();
    private String m_connectionString = null;
    private String m_dbuser = null;
    private String m_dbpass = null;
    private HashMap<Long, Sequence> m_proteinSequences = new HashMap<Long, Sequence>();
    private Sequence m_noSequence = new Sequence(new AminoAcid[0]);
    private DBRunConfig m_conf;
    private HashMap<String,RunConfig> m_configs;
    private boolean m_writePrePostAA = false;
    public static DecimalFormat sixDigits = new DecimalFormat("###0.000000");

    String[] searchIDSetting;
    String filterSetting = "";
    private ArrayList<String> sequenceDBs;

    private DatabaseProvider databaseProvider;
    private static boolean versionSet = false;
    private boolean flagAutoValidated = false;
    private boolean topOnly = true;
    /** should variable modifications (if detected) be marked out */
    private boolean markModifications;
    private boolean markSearchID;
    private boolean markRun;
    private HashSet<String> additionalInfoNames = new HashSet<>();
    
    private String command_line_auto_validate = null;
    
    private Pattern subScoresToForward;
    
    private HashMap<String, Version> m_xi_versions = new HashMap<>();


    // reuse the last mzidentml owner infos
    boolean lastMzIDowner = false;

    private Sequence loadSequence(long id) throws SQLException {
        Connection con = getDBConnection();
        Statement st = con.createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
        ResultSet rs = st.executeQuery("SELECT "
                + " sequence, "
                + " CASE "
                + "    WHEN header IS NULL AND accession_number IS NULL THEN name "
                + "    WHEN header IS NULL THEN accession_number ELSE header "
                + " END as header, "
                + " accession_number "
                + " FROM protein where ID = " + id);
        if (rs.next()) {
            String seqStr = rs.getString(1);
            if (seqStr.matches("'.*'")) {
                seqStr = seqStr.substring(1, seqStr.length() - 2);
            }

            Sequence seq = new Sequence(seqStr, getConfig());
            String headerStr = rs.getString(2);
            if (headerStr.matches("'.*'")) {
                headerStr = headerStr.substring(1, headerStr.length() - 2);
            }

            seq.setFastaHeader(rs.getString(2));

            m_proteinSequences.put(id, seq);
            rs.close();
            return seq;
        } else {
            m_proteinSequences.put(id, m_noSequence);
            rs.close();
            return m_noSequence;
        }
    }

    private Sequence getSequence(Protein p) {
        try {
            long id = p.getId();
            Sequence seq = m_proteinSequences.get(id);

            if (seq == null) {
                return loadSequence(id);
            }
            return seq;
        } catch (SQLException ex) {
            Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, null, ex);
        }
        return m_noSequence;
    }

    public String getLinkWindow(PeptidePair pp, int proteingroup, int window) {
//        return getLinkWindow(pp.getLink(), proteingroup, window);
        try {
            int pepLink = 0;
            org.rappsilber.fdr.entities.Peptide pep;
            if (proteingroup == 0) {
                pep = pp.getPeptide1();
                pepLink = pp.getPeptideLinkSite1() - 1;
            } else {
                pep = pp.getPeptide2();
                pepLink = pp.getPeptideLinkSite2() - 1;
            }
            ProteinGroup pg = pep.getProteinGroup();
            HashMap<Protein, HashSet<Integer>> positions = pep.getPositions();

            Sequence pepSeq;
            if (pep.getSequence().matches("X[0-9]+(\\.[0-9]+)?")) {
                pepSeq = new Sequence("X", getConfig());
            } else {
                pepSeq = new Sequence(pep.getSequence(), getConfig());
            }

            int pepLen = pepSeq.length();

            // prepare the return window 
            ArrayList<CountOccurence<AminoAcid>> sequenceWindow = new ArrayList<CountOccurence<AminoAcid>>(2 * window + 1);

            int from = window - pepLink;
            int to = window + (pepLen - pepLink - 1);

            int pepPos = 0;
            for (int i = 0; i <= 2 * window; i++) {
                CountOccurence<AminoAcid> o = new CountOccurence<AminoAcid>();
                sequenceWindow.add(o);
                if (i >= from && i <= to) {
                    o.add(pepSeq.aminoAcidAt(pepPos++));
                } else {
                    o.add(null);
                }
            }

            // add everything from the surounding proteins into the window
            // load in all sequence windows
            int pid = 0;
            for (Protein p : pg) {
                Sequence seq = getSequence(p);

                HashSet<Integer> sites = positions.get(p);
                for (int protPepFrom : sites) {
                    protPepFrom--;
                    int protPepTo = protPepFrom + pepLen;

                    // protein link
                    int lp = pepLink + protPepFrom;

                    int protFrom = lp - window;
                    int protTo = lp + window;

                    int f = Math.max(protFrom, 0);
                    int swp = f - protFrom;

                    for (int aap = f; aap < protPepFrom; aap++) {
                        sequenceWindow.get(swp).add(seq.aminoAcidAt(aap));
                        swp++;
                    }
                    swp = protPepTo - protFrom;
                    int t = Math.min(seq.length() - 1, protTo);
                    for (int aap = protPepTo; aap <= t; aap++) {
                        sequenceWindow.get(swp).add(seq.aminoAcidAt(aap));
                        swp++;
                    }
                    // next sequence
                    pid++;
                }
            }

            StringBuffer ret = new StringBuffer(2 * window + 1);
            for (CountOccurence<AminoAcid> o : sequenceWindow) {
                ArrayList<AminoAcid> so = o.getSortedList(true);
                AminoAcid aa = so.get(0);
                if (aa == null) {
                    if (so.size() == 1) {
                        ret.append(".");
                    } else {
                        ret.append(so.get(1).SequenceID);
                    }
                } else {
                    ret.append(aa.SequenceID);
                }
            }

            return ret.toString();
        } catch (Exception e) {
            return e.getMessage();
        }
    }

    public String getLinkWindow(ProteinGroupLink l, int proteingroup, int window) {
        try {
            ProteinGroup pg = null;
            HashMap<Protein, IntArrayList> positions = null;
            if (proteingroup == 0) {
                pg = l.getProteinGroup1();
                positions = l.getPosition1();
            } else {
                pg = l.getProteinGroup2();
                positions = l.getPosition2();
            }

            ArrayList<CountOccurence<AminoAcid>> sequenceWindow = new ArrayList<CountOccurence<AminoAcid>>(2 * window + 1);
            for (int i = 2 * window; i >= 0; i--) {
                sequenceWindow.add(new CountOccurence<AminoAcid>());
            }

            // load in all sequence windows
            int pid = 0;
            for (Protein p : pg) {
                Sequence seq = getSequence(p);

                IntArrayList sites = positions.get(p);
                for (int lp : sites) {
                    // make it start at 0
                    lp--;

                    int swp = 0;
                    for (int aap = lp - window; aap <= lp + window; aap++) {
                        if (aap < 0 || aap >= seq.length()) {
                            sequenceWindow.get(swp).add(null);
                        } else {
                            sequenceWindow.get(swp).add(seq.baseAminoAcidAt(aap));
                        }
                        swp++;
                    }
                    // next sequence
                    pid++;
                }
            }

            StringBuffer ret = new StringBuffer(2 * window + 1);
            for (CountOccurence<AminoAcid> o : sequenceWindow) {
                ArrayList<AminoAcid> so = o.getSortedList(true);
                AminoAcid aa = so.get(0);
                if (aa == null) {
                    if (so.size() == 1) {
                        ret.append(".");
                    } else {
                        ret.append(so.get(1).SequenceID);
                    }
                } else {
                    ret.append(aa.SequenceID);
                }
            }

            return ret.toString();
        } catch (Exception e) {
            return e.getMessage();
        }

    }
    
    public void flagDBUsage() {
        synchronized (m_db_last_used) {
            m_db_last_used.value = Calendar.getInstance().getTimeInMillis();
        }
    }
    
    /**
     * @return the m_db_connection
     */
    public Connection getDBConnection() {
        flagDBUsage();
        if (m_db_autoclosetimer == null) {
            m_db_autoclosetimer = new Timer("Timer - auto close db", true);
            m_db_autoclosetimer.schedule(new TimerTask() {
                boolean running = false;
                @Override
                public void run() {
                    if (!running) {
                        running = true;
                        Calendar n = Calendar.getInstance();
                        // if the database connection was not requested for 30 minutes close it.
                        long timeUnused;
                        synchronized (m_db_last_used) {
                            timeUnused = n.getTimeInMillis() - m_db_last_used.value;
                        }
                        if (!keepConnection.get()) {
                            if (timeUnused > 3600000) {
                                closeConnection();
                            }
                        }
                        running = false;
                    }
                }
            }, 60000, 60000);
        }
        try {
            if (!ensureConnection()) {
                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "The database connection seem to be closed and can't be opend now");
                return null;
            }
        } catch (SQLException ex) {
            return null;
        }
        return m_db_connection;
    }

    public void setDBConnection(String connection, String user, String password) throws SQLException {
        m_connectionString = connection;
        m_dbuser = user;
        m_dbpass = password;
        this.m_db_connection = DriverManager.getConnection(connection, user, password);
        setupPreparedStatements();
    }

    /**
     * @param m_db_connection the m_db_connection to set
     */
    public void setDBConnection(Connection connection) throws SQLException {
        this.m_db_connection = connection;
        setupPreparedStatements();
    }

    public void setDatabaseProvider(DatabaseProvider getSearch) throws SQLException {
        databaseProvider = getSearch;
        setDBConnection(databaseProvider.getConnection());
    }

    public synchronized void closeConnection()  {
        try {
            if (m_db_connection!=null) {
                m_db_connection.close();
                m_db_connection=null;
            }
        } catch (SQLException ex) {
            Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public synchronized boolean ensureConnection() throws SQLException {
        flagDBUsage();
        boolean isconnected = false;
//        try {
//            isconnected = m_db_connection.isValid(30);
//            
//        } catch (SQLException ex) {
//            Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, null, ex);
//        } catch (Exception e) {
        try {
            Statement st = m_db_connection.createStatement();
            ResultSet rs = st.executeQuery("select 1");
            rs.close();
            st.close();
            isconnected = true;
        } catch (Exception sex) {
            //Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, "Database connection is closed/ non-functioning. Will try to reopen", sex);
            try {
                m_db_connection.close();
            } catch (Exception e) {
            }

        }
//        }
//        }

        if (!isconnected) {
            Logger.getLogger(this.getClass().getName()).log(Level.WARNING, "Lost connection to database. Will try to reopen");

            try {
                if (m_connectionString != null && m_dbuser != null && m_dbpass != null) {
                    this.m_db_connection = DriverManager.getConnection(m_connectionString, m_dbuser, m_dbpass);
                } else if (databaseProvider != null) {
                    setDBConnection(databaseProvider.getConnection());
                } else {
                    Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, "Currently I don't have a way to reopen the connection");
                    return false;
                }

                isconnected = true;
                this.m_db_connection.setAutoCommit(true);

                setupPreparedStatements();

                if (!GraphicsEnvironment.isHeadless()) {
                    SwingUtilities.invokeLater(new Runnable() {

                        public void run() {
                            try {
                                JOptionPane.showMessageDialog(null, "Reopened the database connection");
                            } catch (Exception e) {}
                        }
                    });
                }

            } catch (SQLException ex) {
                Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Could not reopen the connection to the database", ex);
                throw ex;
            }
        }
        return isconnected;

    }

    protected void setupPreparedStatements() throws SQLException {
        updateValidateOverWrite = this.m_db_connection.prepareStatement("UPDATE spectrum_match set validated = ? WHERE id = ?");
        updateValidateNonOverWrite = this.m_db_connection.prepareStatement("UPDATE spectrum_match set validated = ? WHERE id = ? AND validated is null");
//        updateValidateNonOverWriteSetFDR = this.m_db_connection.prepareStatement("UPDATE spectrum_match SET validated = CASE WHEN validated is null THEN ? ELSE validated END, set peptide_fdr = ?, set link_fdr = ?, ppi_fdr = ? SET WHERE id = ?");
//        updateValidateOverWriteSetFDR = this.m_db_connection.prepareStatement("UPDATE spectrum_match SET validated = ?, set peptide_fdr = ?, set link_fdr = ?, ppi_fdr = ? SET WHERE id = ?");
//        updateSetFDR = this.m_db_connection.prepareStatement("UPDATE spectrum_match SET peptide_fdr = ?, set link_fdr = ?, ppi_fdr = ? SET WHERE id = ?");
        updateSetConfidence = this.m_db_connection.prepareStatement("UPDATE spectrum_match rescored = ? SET WHERE id = ?");

//        updateValidateOverWriteBatch = this.m_db_connection.prepareStatement("UPDATE spectrum_match set validated = ? WHERE id = any(?)");
//        updateValidateNonOverWriteBatch = this.m_db_connection.prepareStatement("UPDATE spectrum_match set validated = ? WHERE id = any(?) AND validated is null");
//        updateValidateNonOverWriteSetFDRBatch = this.m_db_connection.prepareStatement("UPDATE spectrum_match SET validated = CASE WHEN validated is null THEN ? ELSE validated END, set peptide_fdr = ?, set link_fdr = ?, ppi_fdr = ? SET WHERE id = any(?)");
//        updateValidateOverWriteSetFDRBatch = this.m_db_connection.prepareStatement("UPDATE spectrum_match SET validated = ?, set peptide_fdr = ?, set link_fdr = ?, ppi_fdr = ? SET WHERE id = any(?)");
//        updateSetFDRBatch = this.m_db_connection.prepareStatement("UPDATE spectrum_match SET peptide_fdr = ?, set link_fdr = ?, ppi_fdr = ? SET WHERE id = any(?)");
//        updateSetConfidenceBatch = this.m_db_connection.prepareStatement("UPDATE spectrum_match rescored = ? SET WHERE id = any(?)");
    }

    @Override
    public String getSource() {
        if (m_connectionString != null && !m_connectionString.isEmpty()) {
            return "\n\tDB:" + m_connectionString + "\n\tIDs:" + RArrayUtils.toString(m_search_ids, ";") + (filterSetting == null || filterSetting.isEmpty() ? "" : "\n\tFilter: " + filterSetting);
        } else {
            return "\n\tIDs:" + RArrayUtils.toString(m_search_ids, ",") + (filterSetting == null || filterSetting.isEmpty() ? "" : "\n\tFilter: " + filterSetting);
        }
    }

    public ArrayList getSearchIDs() {
        return m_search_ids;
    }

    /**
     * @return the m_conf
     */
    @Override
    public DBRunConfig getConfig() {
        return m_conf;
    }

    @Override
    public RunConfig getConfig(String searchid) {
        if (m_search_ids.size() == 1 && searchid.contentEquals(m_search_ids.get(0))) {
            return getConfig();
        }
        if (m_configs == null) {
            m_configs = new HashMap<>();
        }
        RunConfig ret = m_configs.get(searchid);
        if (ret == null) {
            try {
                ret = new DBRunConfig(getDBConnection());
                ((DBRunConfig)ret).readConfig(searchid);
                m_configs.put(searchid, ret);
            } catch (SQLException ex) {
                Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return ret;
    }

    /**
     * @param m_conf the m_conf to set
     */
    public void setConfig(RunConfig m_conf) {
        if (m_conf instanceof DBRunConfig)
            this.m_conf = (DBRunConfig)m_conf;
        else
            throw new UnsupportedOperationException("Can't assigne a non DBRunConfig here");
    }

    /**
     * @return the writePrePostAA
     */
    public boolean isWritePrePostAA() {
        return m_writePrePostAA;
    }

    /**
     * @param writePrePostAA the writePrePostAA to set
     */
    public void setWritePrePostAA(boolean writePrePostAA) {
        this.m_writePrePostAA = writePrePostAA;
    }

    /**
     * @return the sequenceDBs
     */
    public ArrayList<String> getSequenceDBs() {
        return sequenceDBs;
    }

    public void setMarkSearchID(boolean b) {
        this.markSearchID = b;
    }

    public void setMarkRun(boolean b) {
        this.markRun = b;
    }

    public boolean markSearchID() {
        return this.markSearchID;
    }

    public boolean markRun() {
        return this.markRun;
    }

    private class hashableXiPeptide {

        public rappsilber.ms.sequence.Peptide peptide;
        int hashcode;

        public hashableXiPeptide(Peptide peptide) {
            this.peptide = peptide;
            hashcode = peptide.toString().hashCode();
        }

        public int hashCode() {
            return hashcode;
        }

        public boolean equals(Object o) {
            return ((hashableXiPeptide) o).peptide.equalSequence(peptide);
        }

    }

    public DBinFDR() {
    }

    public DBinFDR(Connection connection) {
        this.m_db_connection = connection;
    }

    public DBinFDR(Connection connection, int[] peptideLengthGroups) {
        super(peptideLengthGroups);
        this.m_db_connection = connection;
    }

    public void readDB() throws SQLException {
        readDB(searchIDSetting, filterSetting, true);
    }

    public void readDB(String searchId, String filter, boolean topOnly) throws SQLException {
        readDB(new String[]{searchId}, filter, topOnly);
    }

    public void readDB(String[] searchIds, String filter, boolean topOnly) throws SQLException {
        this.readDBSteps(searchIds, filter, topOnly, new HashMap<String, HashSet<Long>>(), new HashSet<Long>(), 0, new HashMap<String, Double>());
    }

    public void readDBSteps(String[] searchIds, String filter, boolean topOnly, HashMap<String,HashSet<Long>> allProteinIds, HashSet<Long> skip, int tries, HashMap<String,Double>  lastScore) throws SQLException {

        if (!ensureConnection()) {
            return;
        }

        setConfig(new DBRunConfig(getDBConnection()));
        getConfig().readConfig(searchIds);
        boolean isTargted = false;

        for (CrossLinker xl : getConfig().getCrossLinker()) {
            if (xl.getName().toLowerCase().contains("targetmodification")) {
                isTargted = true;
            }
        }

        String dbNameQuerry = "Select id,name from search_sequencedb ss inner join sequence_file sf ON ss.search_id in (" + RArrayUtils.toString(searchIds, ",") + ") and ss.seqdb_id = sf.id";

        Statement dbnSt = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
        ResultSet rsDBn = dbnSt.executeQuery(dbNameQuerry);
        StringBuilder sb = new StringBuilder();
        sequenceDBs = new ArrayList<String>(1);
        HashSet<Integer> dbIds = new HashSet<Integer>();

        while (rsDBn.next()) {
            int id = rsDBn.getInt(1);
            if (!dbIds.contains(id)) {
                sequenceDBs.add(rsDBn.getString(2));
                dbIds.add(id);
            }
        }

        PeptidePair.ISTARGETED = isTargted;
        boolean xi3db = false;
        try {
            Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
            stm.setFetchSize(100);
            ResultSet rs = stm.executeQuery("SELECT * FROM spectrum_source limit 0;");
            rs.close();
            xi3db = true;
        } catch (SQLException ex) {
            xi3db = false;
        }

        if (!m_search_ids.isEmpty()) {
            if ((filter == null && !filterSetting.isEmpty()) || (!filter.contentEquals(filterSetting))) {
                filterSetting = filterSetting + "\n Search IDS:" + RArrayUtils.toString(searchIds, ",") +":" + filter;
            }
        } else if (filter != null && !filter.isEmpty()) {
            filterSetting = filter;
        }
        boolean shownMinPepWarning =false;

        for (int currentsearch = 0 ; currentsearch<searchIds.length;currentsearch++){
            ensureConnection();
            String searchId = searchIds[currentsearch];
            String searchfilter = filter;
            // wee read that one already
            if (m_search_ids.contains(searchId)) {
                continue;
            }
            HashSet<Long> proteinIds = allProteinIds.get(searchId);
            if (proteinIds == null) {
                proteinIds = new HashSet<>();
                allProteinIds.put(searchId, proteinIds);
            }
                        
            String sPepCoverage1 = "peptide1 unique matched non lossy coverage";
            String sPepCoverage2 = "peptide2 unique matched non lossy coverage";
            String sPepStubs = "fragment CCPepFragment";
            String sPepDoublets = "fragment CCPepDoubletFound";
            int cPepCoverage1 = -1;
            int cPepCoverage2 = -1;
            Integer cPepStubs = -1;
            Integer cPepDoublets = -1;

            ArrayList<Integer> peaks = new ArrayList<>();
            ArrayList<String> scorenames = new ArrayList<>();
            ArrayList<Integer> scoresForwarded = new ArrayList<>();
            
            if (xi3db) {
                // retrive the subscore names
                String scoreNameQuerry = 
                        "SELECT scorenames FROM search WHERE id = " + searchId +";";
                Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
                stm.setFetchSize(500);
                ResultSet rs = stm.executeQuery(scoreNameQuerry);
                if (rs.next()) {
                    java.sql.Array sn = rs.getArray(1);
                    if (sn != null) {
                        for (String s : (String[])sn.getArray()) {
                            if (s.contains("recoursor"))
                                s = s.replace("recoursor", "recursor");
                            additionalInfoNames.add(s);

                            if (s.contentEquals(sPepCoverage1)) {
                                cPepCoverage1 = scorenames.size();
                            } else if (s.contentEquals(sPepCoverage2)) {
                                cPepCoverage2 = scorenames.size();
                            } if (s.startsWith("peak_")) {
                                peaks.add(scorenames.size());
                            }
                            if (subScoresToForward != null && subScoresToForward.matcher(s).matches()) {
                                scoresForwarded.add(scorenames.size());
                            }

                            scorenames.add(s);
                        }
                    }
                }
                cPepCoverage1 = scorenames.indexOf(sPepCoverage1);
                cPepCoverage2 = scorenames.indexOf(sPepCoverage2);
                cPepStubs = scorenames.indexOf(sPepStubs);
                cPepDoublets = scorenames.indexOf(sPepDoublets);
            }

            if (cPepCoverage2 <0 && !shownMinPepWarning) {
                shownMinPepWarning = true;
                Logger.getLogger(this.getClass().getName()).log(Level.WARNING, "Warning - no relative peptide coverage for peptide 2 - bossting on minimum peptide coverage likely not helpfull");
            }

            // if any custom scores where filtered - adapt these
            Pattern subscorepattern = Pattern.compile("\\[%([^%]*)%\\]");
            Matcher subscorematcher = subscorepattern.matcher(searchfilter);
            ArrayList<String> scorefilters_used = new ArrayList<>();
            while(subscorematcher.find()) {
                String scorename = subscorematcher.group(1).replace("oursor", "ursor");
                
                int scoreid=-1;
                for (int i = 0;i<scorenames.size();i++) {
                    if (scorenames.get(i).toLowerCase().equals(scorename.trim().toLowerCase())) {
                        scoreid = i;
                    }
                }
                if (scoreid == -1) {
                    throw new SQLException("subscore [" + scorename + "] not found");
                }
                
                scorefilters_used.add(subscorematcher.group(1));
                
                searchfilter = searchfilter.substring(0, subscorematcher.start()) 
                        +"subscores[" + (scoreid+1) + "]" + searchfilter.substring(subscorematcher.end());
                subscorematcher = subscorepattern.matcher(searchfilter);
                if (!scoresForwarded.contains(scoreid)) {
                    scoresForwarded.add(scoreid);
                }
            }   
            // make sure the gui shows previously used subscores higher up in the list
            if (scorefilters_used.size()>0) {
                String prev = LocalProperties.getProperty(DBinFDR.subscoreroperty, "");
                if (prev.length() > 0) {
                    for (String s : prev.split(";")) {
                        if (!scorefilters_used.contains(s))
                            scorefilters_used.add(s);
                    }
                    
                }
                LocalProperties.setProperty(DBinFDR.subscoreroperty, RArrayUtils.toString(scorefilters_used,";"));
            }
            
            Double searchLastScore = lastScore.get(searchId);
            String matchQuerry;
            matchQuerry
                    = "SELECT * FROM (SELECT sm.id AS psmID, \n"
                    + "p1.sequence AS pepSeq1, \n"
                    + "p2.sequence AS pepSeq2, \n"
                    + "p1.peptide_length as peplen1, \n"
                    + "p2.peptide_length as peplen2, \n"
                    + "mp1.link_position + 1 as site1, \n"
                    + "mp1.link_site_score as link_site_score1, \n"
                    + "mp2.link_position + 1 as site2, \n"
                    + "mp2.link_site_score as link_site_score2, \n"
                    + "pr1.is_decoy AS isDecoy1, \n"
                    + "pr2.is_decoy AS isDecoy2, \n"
                    + "sm.precursor_charge AS calc_charge, \n"
                    + "sm.score, \n"
                    + "array_agg(hp1.peptide_position + 1) AS pepPosition1,  \n"
                    + "array_agg(hp2.peptide_position + 1) AS pepPosition2, \n"
                    + "CASE WHEN p2.sequence IS NULL THEN 1 ELSE (4.0/5.0+(p1.peptide_length/(p1.peptide_length+p2.peptide_length)))/2 END AS score_ratio \n"
                    + " , s.precursor_charge as exp_charge \n"
                    + " , array_agg(pr1.id) as protein1id\n"
                    + " , array_agg(pr2.id) as protein2id\n"
                    + " , p1.id as peptide1id\n"
                    + " , p2.id as peptide2id\n"
                    + " , v.run_name, v.scan_number \n"
                    + " , scorepeptide1matchedconservative AS scoreP1Coverage \n"
                    + " , scorepeptide2matchedconservative AS scoreP2Coverage \n"
                    + " , scoredelta AS deltaScore\n"
                    + " , scorelinksitedelta AS LinkSiteDelta\n"
                    + " , s.precursor_mz AS exp_mz\n"
                    + " , sm.calc_mass\n"
                    + " , sm.precursor_charge as match_charge\n"
                    + " , p1.mass as pep1mass\n"
                    + " , p2.mass as pep2mass\n"
                    + " , sm.search_id\n"
                    + " , sm.spectrum_id\n"
                    + " , s.precursor_charge\n"
                    + " , autovalidated\n"
                    + " , cl.name  AS crosslinker \n"
                    + " , cl.mass  AS clmass \n"
                    + " , sm.rank \n"
                    + " , s.scan_index\n"
                    + " , plf.name as peaklistfile\n"
                    + " , scoremgxrank AS mgxrank\n"
                    + " , scorecleavclpep1fragmatched AS cleavclpep1fragmatched\n"
                    + " , scorecleavclpep2fragmatched AS cleavclpep2fragmatched\n"
                    + " , scores AS subscores\n"
                    + " , scoredelta AS delta\n"
                    + " , s.precursor_intensity\n"
                    + " , s.elution_time_start as retentiontime\n"
                    + " , validated\n"
                    + " \n"
                    + "FROM \n"
                    + "  (SELECT * FROM Spectrum_match WHERE Search_id = " + searchId + (topOnly ? " AND dynamic_rank = 't'":"")
                    + "     AND score>0 "+ (searchLastScore == null ? "": " AND score <= " + searchLastScore) +" ) sm \n"
                    + "  INNER JOIN (\n"
                    + "SELECT ss.name as run_name, s.scan_number, sm.id as spectrum_match_id FROM (select * from spectrum_match where Search_id = " + searchId + (topOnly ? " AND dynamic_rank = 't'":"")+ " AND score>0) sm inner join spectrum s on sm.spectrum_id = s.id INNER JOIN spectrum_source ss on s.source_id = ss.id\n"                            
                    + ") v \n"
                    + "    ON v.spectrum_match_id = sm.id\n"
                    + "    inner join \n"
                    + "   matched_peptide mp1 on sm.id = mp1.match_id and mp1.match_type =1  LEFT OUTER JOIN \n"
                    + "   matched_peptide mp2 on sm.id = mp2.match_id AND mp2.match_type = 2  INNER JOIN \n"
                    + "   peptide p1 on mp1.peptide_id = p1.id LEFT OUTER JOIN  \n"
                    + "   peptide p2 on mp2.peptide_id = p2.id INNER JOIN \n"
                    + "   has_protein hp1 ON mp1.peptide_id = hp1.peptide_id LEFT OUTER JOIN  \n"
                    + "   has_protein hp2 ON mp2.peptide_id = hp2.peptide_id INNER JOIN \n"
                    + "   protein pr1 ON hp1.protein_id = pr1.id LEFT OUTER JOIN \n"
                    + "   protein pr2 ON hp2.protein_id = pr2.id\n"
                    + "   INNER JOIN  \n"
                    + "   spectrum s ON sm.spectrum_id = s.id \n"
                    + " LEFT OUTER JOIN peaklistfile plf on s.peaklist_id = plf.id\n"
                    + " LEFT OUTER JOIN crosslinker cl on mp1.crosslinker_id = cl.id \n"
                    + " WHERE  \n"
                    + " (s.precursor_charge = sm.precursor_charge OR sm.precursor_charge < 6) \n"
                    + " GROUP BY  sm.id, \n"
                    + "p1.sequence, \n"
                    + "p2.sequence, \n"
                    + "p1.peptide_length, \n"
                    + "p2.peptide_length, \n"
                    + "mp1.link_position + 1, \n"
                    + "mp2.link_position + 1, \n"
                    + "pr1.is_decoy, \n"
                    + "pr2.is_decoy, \n"
                    + "sm.precursor_charge, \n"
                    + "sm.score, \n"
                    + "CASE WHEN p2.sequence IS NULL THEN 1 ELSE (4.0/5.0+(p1.peptide_length/(p1.peptide_length+p2.peptide_length)))/2 END\n"
                    + " , s.precursor_charge \n"
                    + " , p1.id "
                    + " , p2.id "
                    + " , v.run_name, v.scan_number \n"
                    + " , scorepeptide1matchedconservative \n"
                    + " , scorepeptide2matchedconservative \n"
                    + " , scoredelta "
                    + " , scorelinksitedelta "
                    + " , s.precursor_mz "
                    + " , sm.calc_mass\n"
                    + " , sm.precursor_charge "
                    + " , p1.mass "
                    + " , p2.mass "
                    + " , sm.search_id\n"
                    + " , sm.spectrum_id\n"
                    + " , s.precursor_charge\n"
                    + " , autovalidated\n"
                    + " , cl.name \n"
                    + " , cl.mass \n"
                    + " , sm.rank \n"
                    + " , s.scan_index\n"
                    + " , plf.name "
                    + " , scoremgxrank "
                    + " , scorecleavclpep1fragmatched "
                    + " , scorecleavclpep2fragmatched "
                    + " , scores "
                    + " , scoredelta "
                    + " , s.precursor_intensity\n"
                    + " , s.elution_time_start \n"
                    + " , mp1.link_site_score \n"
                    + " , mp2.link_site_score \n"
                    + " , validated \n"
                    + " ORDER BY sm.score DESC)  i \n"
                    + (searchfilter == null || searchfilter.isEmpty() ? "" : " WHERE ( " + searchfilter + " )");

            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "read from db");
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, matchQuerry);
            
            Statement stmEA = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
//            ResultSet rsEA = stmEA.executeQuery("EXPLAIN ANALYZE " + matchQuerry);
//            ResultSetMetaData rsEAM = rsEA.getMetaData();
//            int eaCols = rsEAM.getColumnCount();
//            while (rsEA.next()) {
//                for (int c = 1; c<= eaCols; c++) {
//                    System.err.print(rsEA.getString(c)  +  " , ");
//                }
//                System.err.println();
//            }
            
            if (searchIds.length >1)
                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Search {0} of {1}", new Object[]{currentsearch+1, searchIds.length});
            //getDBConnection().setAutoCommit(false);
            Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
            stm.setFetchSize(1000000);
            ResultSet rs;
            this.keepConnection.set(true);
            try {
                flagDBUsage();
                rs = stm.executeQuery(matchQuerry);
                flagDBUsage();
            } finally{
                this.keepConnection.set(false);
            }
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "go through the results");

            int psmIDColumn = rs.findColumn("psmID");
            int pepSeq1Column = rs.findColumn("pepSeq1");
            int pepSeq2Column = rs.findColumn("pepSeq2");
            int peplen1Column = rs.findColumn("peplen1");
            int peplen2Column = rs.findColumn("peplen2");
            int site1Column = rs.findColumn("site1");
            int site2Column = rs.findColumn("site2");
            int isDecoy1Column = rs.findColumn("isDecoy1");
            int isDecoy2Column = rs.findColumn("isDecoy2");
            int calc_chargeColumn = rs.findColumn("calc_charge");
            int scoreColumn = rs.findColumn("score");
            int pepPosition1Column = rs.findColumn("pepPosition1");
            int pepPosition2Column = rs.findColumn("pepPosition2");
            int score_ratioColumn = rs.findColumn("score_ratio");
            int exp_chargeColumn = rs.findColumn("exp_charge");
            int protein1idColumn = rs.findColumn("protein1id");
            int protein2idColumn = rs.findColumn("protein2id");
            int peptide1idColumn = rs.findColumn("peptide1id");
            int peptide2idColumn = rs.findColumn("peptide2id");
            int run_nameColumn = rs.findColumn("run_name");
            int scan_numberColumn = rs.findColumn("scan_number");
            int scoreP1CoverageColumn = rs.findColumn("scoreP1Coverage");
            int scoreP2CoverageColumn = rs.findColumn("scoreP2Coverage");
            int deltaScoreColumn = rs.findColumn("deltaScore");
            int LinkSiteDeltaColumn = rs.findColumn("LinkSiteDelta");
            int exp_mzColumn = rs.findColumn("exp_mz");
            int calc_massColumn = rs.findColumn("calc_mass");
            int pep1massColumn = rs.findColumn("pep1mass");
            int pep2massColumn = rs.findColumn("pep2mass");
            int search_idColumn = rs.findColumn("search_id");
            int spectrum_idColumn = rs.findColumn("spectrum_id");
            int xlColumn = rs.findColumn("crosslinker");
            int avColumn = rs.findColumn("autovalidated");
            int rankColumn = rs.findColumn("rank");
            int scanIndexColumn = rs.findColumn("scan_index");
            int peaklistfileColumn = rs.findColumn("peaklistfile");
            int mgxrankColumn = rs.findColumn("mgxrank");
            int cleavclpep1fragmatchedColumn = rs.findColumn("cleavclpep1fragmatched");
            int cleavclpep2fragmatchedColumn = rs.findColumn("cleavclpep2fragmatched");
            int subscoresColumn = rs.findColumn("subscores");
            int precIntensityColumn = rs.findColumn("precursor_intensity");
            int retentiontimeColumn = rs.findColumn("retentiontime");
            int linkSiteScore1Column = rs.findColumn("link_site_score1");
            int linkSiteScore2Column = rs.findColumn("link_site_score2");
            int clMassColumn = rs.findColumn("clmass");
            int clDeltaColumn = rs.findColumn("delta");
            Pattern modDetect = Pattern.compile(".*[^A-Z].*");
            HashSet<Double> tmmodcount = new HashSet<>();
            HashMap<Double,Double> xlmodmasses = new HashMap<Double,Double>(1);

            int total = 0;
            try {
                while (rs.next()) {
                    total++;

                    long ipsmID = rs.getLong(psmIDColumn);
                    if (skip.contains(ipsmID)) {
                        continue;
                    }

                    String psmID = Long.toString(ipsmID);
                    String pepSeq1 = rs.getString(pepSeq1Column);
                    String pepSeq2 = rs.getString(pepSeq2Column);
                    int peplen1 = rs.getInt(peplen1Column);
                    int peplen2 = rs.getInt(peplen2Column);
                    int site1 = rs.getInt(site1Column);
                    int site2 = pepSeq2 == null ? -1 : rs.getInt(site2Column);
                    boolean isDecoy1 = rs.getBoolean(isDecoy1Column);
                    boolean isDecoy2 = rs.getBoolean(isDecoy2Column);
                    int charge = rs.getInt(calc_chargeColumn);
                    double score = rs.getDouble(scoreColumn);
                    Integer[] pepPosition1 = (Integer[])rs.getArray(pepPosition1Column).getArray();
                    Integer[] pepPosition2 = pepSeq2 == null ? new Integer[]{-1} : (Integer[])rs.getArray(pepPosition2Column).getArray();
                    //            double coverage1  = rs.getDouble(18);
                    //            double coverage2  = rs.getDouble(19);
                    //            double scoreRatio = coverage1/(coverage1+coverage2);
                    double scoreRatio = rs.getDouble(score_ratioColumn);
                    int spectrum_charge = rs.getInt(exp_chargeColumn);
                    Long[] protein1ID = (Long[])rs.getArray(protein1idColumn).getArray();
                    Long[] protein2ID = pepSeq2 == null ? new Long[]{0l} : (Long[])rs.getArray(protein2idColumn).getArray();
                    long pep1ID = rs.getLong(peptide1idColumn);
                    long pep2ID = rs.getLong(peptide2idColumn);
                    String run = rs.getString(run_nameColumn);
                    String scan = rs.getString(scan_numberColumn);
                    boolean autovalidated = rs.getBoolean(avColumn);
                    int rank = rs.getInt(rankColumn);



                    double p1c = rs.getDouble(scoreP1CoverageColumn);
                    double p2c = rs.getDouble(scoreP2CoverageColumn);
                    double pminc = p1c;
                    
                    if (pepSeq2 != null && !pepSeq2.isEmpty()) {
                        pminc = Math.min(p1c,p2c);
                    }
                    
                    double pmz = rs.getDouble(exp_mzColumn);
                    double f = 1;
                    if (pepSeq2 != null && !pepSeq2.isEmpty() && p1c + p2c > 0) {
                        ////                    double max = Math.max(p1c,p2c);
                        ////                    double min = Math.min(p1c,p2c);
                        ////                    f = min/(p1c+p2c);
                        ////                    score = score * f;
                        scoreRatio = (p1c) / (p1c + p2c + 1);
                        //                    if (p1c <3 || p2c <3) 
                        //                        continue;
                    }
                    
                    double peptide1score = score*scoreRatio;
                    double peptide2score = score*(1-scoreRatio);
                    
                    String xl = rs.getString(xlColumn);
                    if (xl == null) {
                        xl = "";
                    }
                    double calc_mass = rs.getDouble(calc_massColumn);
                    double pep1mass = rs.getDouble(pep1massColumn);
                    double pep2mass = rs.getDouble(pep2massColumn);
                    int search_id = rs.getInt(search_idColumn);
                    int scan_id = rs.getInt(spectrum_idColumn);
                    if (pepSeq2 != null && pepSeq2.matches("^X-?[0-9\\.]*$")) {
                        double mass = Double.parseDouble(pepSeq2.substring(1));
                        AminoModification am = new AminoModification(pepSeq2, AminoAcid.A, mass-18.0105647);
                        m_conf.addKnownModification(am);
                        getConfig(searchId).addKnownModification(am);
                        tmmodcount.add(((int)mass*10)/10.0);
                    }

                    int mgxrank = rs.getInt(mgxrankColumn);
                    boolean cleavclpep1fragmatched = rs.getBoolean(cleavclpep1fragmatchedColumn);
                    boolean cleavclpep2fragmatched = rs.getBoolean(cleavclpep2fragmatchedColumn);
                    java.sql.Array subscores = rs.getArray(subscoresColumn);

                    PSM psm = null;
                    if (pepSeq2 != null && !pepSeq2.isEmpty()){
                        for (int p = 0; p< protein1ID.length; p++) {
                            int p1 = pepPosition1[p];
                            long p1id = protein1ID[p];
                            proteinIds.add(p1id);
                            
                            String a1 = Long.toString(p1id);
                            String n1 = Long.toString(p1id);
                            String d1 = Long.toString(p1id);
                            String s1 = Long.toString(p1id);

                            int p2 = pepPosition2[p];
                            long p2id = protein2ID[p];
                            proteinIds.add(p2id);
                            String a2 = Long.toString(p2id);
                            String n2 = Long.toString(p2id);
                            String d2 = Long.toString(p2id);
                            String s2 = Long.toString(p2id);

                            psm = setUpDBPSM(psmID, run, scan, pep1ID, pep2ID, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, p1id, a1, d1, p2id, a2, d2, p1, p2, s1, s2, peptide1score, peptide2score, spectrum_charge, xl, pmz, calc_mass, pep1mass, pep2mass, search_id, scan_id);
                        }
                    } else {
                        String a2 = "";
                        String n2 = "";
                        String d2 = "";
                        int p2 = -1;
                        long p2id = 0;
                        String s2 = null;
                        for (int p = 0; p< protein1ID.length; p++) {
                            int p1 = pepPosition1[p];
                            long p1id = protein1ID[p];
                            proteinIds.add(p1id);
                            String s1 = Long.toString(p1id);
                            String a1 = Long.toString(p1id);
                            String n1 = Long.toString(p1id);
                            String d1 = Long.toString(p1id);
                            if (d2==null || d2.isEmpty()) {
                                if (n2 == null || n2.isEmpty())
                                    d2 = a2;
                                else
                                    d2 = n2;
                            }
                            psm = setUpDBPSM(psmID, run, scan, pep1ID, pep2ID, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, p1id, a1, d1, p2id, a2, d2, p1, p2, s1, s2, peptide1score, peptide2score, spectrum_charge, xl, pmz, calc_mass, pep1mass, pep2mass, search_id, scan_id);
                        }
                    }
                    
                    psm.setRank(rank);
                    
                    Float[] scorevalues = null;
                    if (subscores != null) {
                        scorevalues = (Float[])subscores.getArray();
                        if (cPepCoverage1 >= 0) {
                            if (cPepCoverage2 <0 || Double.isNaN(scorevalues[cPepCoverage2]))
                                psm.addOtherInfo("minPepCoverage",
                                    scorevalues[cPepCoverage1]);
                            else
                                psm.addOtherInfo("minPepCoverage",
                                    Math.min(scorevalues[cPepCoverage1], 
                                            scorevalues[cPepCoverage2])
                            );
                        }
                    }

                    if (cPepStubs >= 0) {
                        int s = scorevalues[cPepStubs].intValue();
                        psm.peptidesWithStubs =  s;
                        psm.addOtherInfo("PeptidesWithStubs", s);
                        if (s >0)  {
                            stubsFound(true);
                        } 

                    }
                    if (cPepDoublets >= 0) {

                        psm.peptidesWithDoublets = scorevalues[cPepDoublets].intValue();
                        psm.addOtherInfo("PeptidesWithDoublets", psm.peptidesWithDoublets);

                    }
                    
                    
                    psm.addOtherInfo("P1Fragments",p1c);
                    psm.addOtherInfo("P2Fragments",p2c);
                    psm.addOtherInfo("MinFragments",pminc);
                    if (pep2mass>0) {
                        psm.addOtherInfo("minPepCoverageAbsolute",
                            Math.min(p1c, p2c));
                    } else {
                        psm.addOtherInfo("minPepCoverageAbsolute", p1c);
                    }
                    
                    psm.addOtherInfo("cleavclpep1fragmatched",cleavclpep1fragmatched ? 1:0);
                    psm.addOtherInfo("cleavclpep2fragmatched",cleavclpep2fragmatched ? 1:0);
                    psm.addOtherInfo("cleavclpepfragmatched",(cleavclpep2fragmatched ? 1:0) + (cleavclpep1fragmatched ? 1:0));
                    psm.addOtherInfo("deltaScore",rs.getDouble(deltaScoreColumn));
                    psm.addOtherInfo("ScoreDivDelta",rs.getDouble(scoreColumn)/rs.getDouble(deltaScoreColumn));
                    psm.addOtherInfo("linkSiteDelta",rs.getDouble(LinkSiteDeltaColumn));
                    psm.addOtherInfo("mgxrank",rs.getDouble(mgxrankColumn));
                    psm.addOtherInfo("PrecursorIntensity",rs.getDouble(precIntensityColumn));
                    psm.addOtherInfo("RetentionTime",rs.getDouble(retentiontimeColumn));
                    Array lss1 = rs.getArray(linkSiteScore1Column);
                    Array lss2 = rs.getArray(linkSiteScore2Column);
                    if (lss1 != null) {
                        psm.addOtherInfo("LinkSiteScore1",(Float[])lss1.getArray());
                    } else {
                        psm.addOtherInfo("LinkSiteScore1", new Float[0]);
                    }
                    if (lss2 != null) {
                        psm.addOtherInfo("LinkSiteScore2",(Float[])lss2.getArray());
                    } else {
                        psm.addOtherInfo("LinkSiteScore2", new Float[0]);
                    }
                    
                    for (Integer p : peaks) {
                        psm.addOtherInfo(scorenames.get(p), scorevalues[p]);
                    }
                    
                    for (Integer p : scoresForwarded) {
                        psm.addOtherInfo(scorenames.get(p), scorevalues[p]);
                    }
                    
                    Double xlModmassPre = rs.getDouble(clMassColumn);
                    Double xlModmass = xlmodmasses.get(xlModmassPre);
                    if (xlModmass == null) {
                        xlmodmasses.put(xlModmassPre,xlModmassPre);
                        xlModmass = xlModmassPre;
                    }
                    psm.setCrosslinkerModMass(rs.getDouble(clMassColumn));
                    psm.setDeltaScore(rs.getDouble(deltaScoreColumn));

                    if  (autovalidated) {
                        psm.setAutoValidated(autovalidated);
                        if (flagAutoValidated) {
                            psm.setPositiveGrouping("AutoValidated");
                        }
                    }

                    String modLoockup = pepSeq1;
                    if (pepSeq2 != null)
                        modLoockup+=pepSeq2;
                    RunConfig psmconfig  =getConfig(psm.getSearchID());
                    Sequence m = new Sequence(modLoockup, psmconfig);
                    for (AminoAcid aa : m) {
                        if (aa instanceof AminoModification) {
                            if (psmconfig.getVariableModifications().contains(aa)) {
                                psm.setHasVarMods(true);
                            }
                            if (psmconfig.getFixedModifications().contains(aa)) {
                                psm.setHasFixedMods(true);
                            }
                        }
                    }
                    if (psm.hasVarMods()) {
                        if (markModifications())
                            psm.addNegativeGrouping("Modified");
                    }

                    if (markRun())
                        psm.getAdditionalFDRGroups().add(psm.getRun());

                    if (markSearchID())
                        psm.getAdditionalFDRGroups().add(""+psm.getSearchID());
                    
                    if (psm.hasVarMods()) {
                        if (markModifications())
                            psm.addNegativeGrouping("Modified");
                    }
                    

                    int scanIndex = rs.getInt(scanIndexColumn);
                    if (!rs.wasNull()) {
                        psm.setFileScanIndex(scanIndex);
                    }
                    String peaklistfile = rs.getString(peaklistfileColumn);
                    if (!rs.wasNull()) {
                        psm.setPeakListName(peaklistfile);
                    }
                    
//                    if (isMultpleSearches) {
//                        psm.addPositiveGrouping("" + searchId);
//                    }
                    if (total % 100 == 0) {
                        flagDBUsage();
                    }
                    
                    if (total % 10000 == 0) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO, 
                                "go through the results ({0}) Search {1} of {2}", 
                                new Object[]{allPSMs.size(), currentsearch+1, searchIds.length});
                    }
                    lastScore.put(searchId, score);
                    skip.add(ipsmID);

                }
                rs.close();
                stm.close();

                if (allProteins.size() > 0) {
                    readProteinInfos(proteinIds);
                }
                        
                m_search_ids.add(searchId);

            } catch (SQLException sex) {
                if (tries < 5) {
                    int nexttry = tries;
                    if (total<10)
                        nexttry++;
                        Logger.getLogger(this.getClass().getName()).log(Level.FINE, "failed to read from the database retrying", sex);
                    readDBSteps(searchIds, filter, topOnly, allProteinIds, skip, nexttry, lastScore);
                } else {
                    Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Repeatedly (" + total + ") failed to read from the database giving up now", sex);
                }
                return;
            }

            // make sure the hashmap/sets are actually using the correct hashes
            SelfAddHashSet<Protein> newProteins = new SelfAddHashSet<Protein>();
            for (Protein p : allProteins) {
                newProteins.add(p);
            }
            allProteins = newProteins;
            // also the psms use hashes for the proteins
            for (PSM psm : allPSMs) {
                for (org.rappsilber.fdr.entities.Peptide p : psm.getPeptides()) {
                    p.renewPositions();
                }
                psm.reTestInternal();
            }

            
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Count results : " + total);

            if (tmmodcount.size()>10) {
                PeptidePair.ISTARGETED = false;
            }
        }
    }

    protected void readProteinInfos(HashSet<Long> proteinIds) throws SQLException {
        Statement stm;
        ResultSet rs;
        HashMap<Long,Protein> id2Protein = new HashMap<>();
        for (Protein p : allProteins) {
            id2Protein.put(p.getId(), p);
        }
        String proteinQuerry = "SELECT id, accession_number, name, description, sequence FROM protein WHERE id IN (" + RArrayUtils.toString(proteinIds, ", ") + ")";
        Logger.getLogger(this.getClass().getName()).log(Level.FINE, "Geting protein information: \n" + proteinQuerry);
        ensureConnection();
        stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
        stm.setFetchSize(100);
        rs = stm.executeQuery(proteinQuerry);
        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "go through the results");
        int idColumn = rs.findColumn("id");
        int accessionColumn = rs.findColumn("accession_number");
        int descriptionColumn = rs.findColumn("description");
        int proteinSequenceColumn = rs.findColumn("sequence");
        int proteinNameColumn = rs.findColumn("name");
        while (rs.next()) {
            long id = rs.getLong(idColumn);
            String sequence =  rs.getString(proteinSequenceColumn).replaceAll("[^A-Z]", "");
            String accession = rs.getString(accessionColumn);
            String name = rs.getString(proteinNameColumn);
            String description = rs.getString(descriptionColumn);
            
            if (accession == null && description !=null) {
                accession = description;
            }
            
            if (description==null || description.isEmpty()) {
                if (name == null || name.isEmpty())
                    description = accession;
                else
                    description = name;
            }
            
            Protein p = id2Protein.get(id);
            p.setSequence(sequence);
            p.setAccession(accession);
            p.setDescription(description);
            p.setName(name);
            
        }
        rs.close();
        stm.close();
    }

    public void readDB(String[] searchIds, String filter, boolean topOnly, ArrayList<Long> skip, int tries) throws SQLException {

        if (!ensureConnection()) {
            return;
        }

        setConfig(new DBRunConfig(getDBConnection()));
        getConfig().readConfig(searchIds);
        boolean isTargted = false;

        for (CrossLinker xl : getConfig().getCrossLinker()) {
            if (xl.getName().toLowerCase().contains("targetmodification")) {
                isTargted = true;
            }
        }

        String dbNameQuerry = "Select id,name from search_sequencedb ss inner join sequence_file sf ON ss.search_id in (" + RArrayUtils.toString(searchIds, ",") + ") and ss.seqdb_id = sf.id";

        Statement dbnSt = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
        ResultSet rsDBn = dbnSt.executeQuery(dbNameQuerry);
        StringBuilder sb = new StringBuilder();
        sequenceDBs = new ArrayList<String>(1);
        HashSet<Integer> dbIds = new HashSet<Integer>();

        while (rsDBn.next()) {
            int id = rsDBn.getInt(1);
            if (!dbIds.contains(id)) {
                sequenceDBs.add(rsDBn.getString(2));
                dbIds.add(id);
            }
        }

        PeptidePair.ISTARGETED = isTargted;
        boolean xi3db = false;
        try {
            Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
            stm.setFetchSize(100);
            ResultSet rs = stm.executeQuery("SELECT * FROM spectrum_source limit 0;");
            rs.close();
            xi3db = true;
        } catch (SQLException ex) {
            xi3db = false;
        }

        if (!m_search_ids.isEmpty()) {
            if ((filter == null && !filterSetting.isEmpty()) || (!filter.contentEquals(filterSetting))) {
                filterSetting = filterSetting + "\n Search IDS:" + RArrayUtils.toString(searchIds, ",") +":" + filter;
            }
        } else if (filter != null && !filter.isEmpty()) {
            filterSetting = filter;
        }
        boolean shownMinPepWarning =false;

        for (int currentsearch = 0 ; currentsearch<searchIds.length;currentsearch++){
            String searchfilter = filter;
            String searchId = searchIds[currentsearch];
            // wee read that one already
            if (m_search_ids.contains(searchId)) {
                continue;
            }
            String sPepCoverage1 = "peptide1 unique matched non lossy coverage";
            String sPepCoverage2 = "peptide2 unique matched non lossy coverage";
            int cPepCoverage1 = -1;
            int cPepCoverage2 = -1;

            ArrayList<Integer> peaks = new ArrayList<>();
            ArrayList<String> scorenames = new ArrayList<>();
            ArrayList<Integer> scoresForwarded = new ArrayList<>();
            
            if (xi3db) {
                // retrive the subscore names
                String scoreNameQuerry = 
                        "SELECT scorenames FROM search WHERE id = " + searchId +";";
                Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
                stm.setFetchSize(500);
                ResultSet rs = stm.executeQuery(scoreNameQuerry);
                if (rs.next()) {
                    java.sql.Array sn = rs.getArray(1);
                    if (sn != null) {
                        for (String s : (String[])sn.getArray()) {
                            additionalInfoNames.add(s);

                            if (s.contentEquals(sPepCoverage1)) {
                                cPepCoverage1 = scorenames.size();
                            } else if (s.contentEquals(sPepCoverage2)) {
                                cPepCoverage2 = scorenames.size();
                            } if (s.startsWith("peak_")) {
                                peaks.add(scorenames.size());
                            }
                            if (subScoresToForward != null && subScoresToForward.matcher(s).matches()) {
                                scoresForwarded.add(scorenames.size());
                            }

                            scorenames.add(s);
                        }
                    }
                }
                cPepCoverage1 = scorenames.indexOf(sPepCoverage1);
                cPepCoverage2 = scorenames.indexOf(sPepCoverage2);
            }

            if (cPepCoverage2 <0 && !shownMinPepWarning) {
                shownMinPepWarning = true;
                Logger.getLogger(this.getClass().getName()).log(Level.WARNING, "Warning - no relative peptide coverage for peptide 2 - bossting on minimum peptide coverage likely not helpfull");
            }
            
            // if any custom scores where filtered - adapt these
            Pattern subscorepattern = Pattern.compile("\\[%([^%]*)%\\]");
            Matcher subscorematcher = subscorepattern.matcher(searchfilter);
            while(subscorematcher.find()) {
                String scorename = subscorematcher.group(1);
                int scoreid=scorenames.indexOf(scorename);
                searchfilter = searchfilter.substring(0, subscorematcher.start()) +"subscores[" + (scoreid+1) + "]" + searchfilter.substring(subscorematcher.end());
                subscorematcher = subscorepattern.matcher(searchfilter);
            }
            
            String matchQuerry;
            matchQuerry
                    = "SELECT * FROM (SELECT sm.id AS psmID, \n"
                    + "p1.sequence AS pepSeq1, \n"
                    + "p2.sequence AS pepSeq2, \n"
                    + "p1.peptide_length as peplen1, \n"
                    + "p2.peptide_length as peplen2, \n"
                    + "mp1.link_position + 1 as site1, \n"
                    + "mp2.link_position + 1 as site2, \n"
                    + "pr1.is_decoy AS isDecoy1, \n"
                    + "pr2.is_decoy AS isDecoy2, \n"
                    + "sm.precursor_charge AS calc_charge, \n"
                    + "sm.score, \n"
                    + "array_agg(pr1.accession_number) AS accession1, \n"
                    + "array_agg(pr1.description) AS  description1, \n"
                    + "array_agg(pr2.accession_number) AS accession2, \n"
                    + "array_agg(pr2.description) AS  description2, \n"
                    + "array_agg(hp1.peptide_position + 1) AS pepPosition1,  \n"
                    + "array_agg(hp2.peptide_position + 1) AS pepPosition2, \n"
                    + "CASE WHEN p2.sequence IS NULL THEN 1 ELSE (4.0/5.0+(p1.peptide_length/(p1.peptide_length+p2.peptide_length)))/2 END AS score_ratio \n"
                    + " , s.precursor_charge as exp_charge \n"
                    + " , array_agg(pr1.id) as protein1id\n"
                    + " , array_agg(pr2.id) as protein2id\n"
                    + " , p1.id as peptide1id\n"
                    + " , p2.id as peptide2id\n"
                    + " , v.run_name, v.scan_number \n"
                    + " , array_agg(pr1.sequence) AS protein1sequence \n"
                    + " , array_agg(pr2.sequence) AS protein2sequence\n"
                    + " , scorepeptide1matchedconservative AS scoreP1Coverage \n"
                    + " , scorepeptide2matchedconservative AS scoreP2Coverage \n"
                    + " , scoredelta AS deltaScore\n"
                    + " , scorelinksitedelta AS LinkSiteDelta\n"
                    + " , s.precursor_mz AS exp_mz\n"
                    + " , sm.calc_mass\n"
                    + " , sm.precursor_charge as match_charge\n"
                    + " , p1.mass as pep1mass\n"
                    + " , p2.mass as pep2mass\n"
                    + " , sm.search_id\n"
                    + " , sm.spectrum_id\n"
                    + " , array_agg(pr1.name) as protein1name\n"
                    + " , array_agg(pr2.name) as protein2name\n"
                    + " , s.precursor_charge\n"
                    + " , autovalidated\n"
                    + " , cl.name  AS crosslinker \n"
                    + " , cl.mass  AS clmass \n"
                    + " , sm.rank \n"
                    + " , s.scan_index\n"
                    + " , plf.name as peaklistfile\n"
                    + " , scoremgxrank AS mgxrank\n"
                    + " , scorecleavclpep1fragmatched AS cleavclpep1fragmatched\n"
                    + " , scorecleavclpep2fragmatched AS cleavclpep2fragmatched\n"
                    + " , scores AS subscores\n"
                    + " , scoredelta AS delta\n"
                    + " , s.precursor_intensity\n"
                    + " , s.elution_time_start as retentiontime\n"
                    + " \n"
                    + "FROM \n"
                    + "  (SELECT * FROM Spectrum_match WHERE Search_id = " + searchId + (topOnly ? " AND dynamic_rank = 't'":"")+ " AND score>0) sm \n"
                    + "  INNER JOIN (\n"
                    + "SELECT ss.name as run_name, s.scan_number, sm.id as spectrum_match_id FROM (select * from spectrum_match where Search_id = " + searchId + (topOnly ? " AND dynamic_rank = 't'":"")+ " AND score>0) sm inner join spectrum s on sm.spectrum_id = s.id INNER JOIN spectrum_source ss on s.source_id = ss.id\n"                            
                    + ") v \n"
                    + "    ON v.spectrum_match_id = sm.id\n"
                    + "    inner join \n"
                    + "   matched_peptide mp1 on sm.id = mp1.match_id and mp1.match_type =1  LEFT OUTER JOIN \n"
                    + "   matched_peptide mp2 on sm.id = mp2.match_id AND mp2.match_type = 2  INNER JOIN \n"
                    + "   peptide p1 on mp1.peptide_id = p1.id LEFT OUTER JOIN  \n"
                    + "   peptide p2 on mp2.peptide_id = p2.id INNER JOIN \n"
                    + "   has_protein hp1 ON mp1.peptide_id = hp1.peptide_id LEFT OUTER JOIN  \n"
                    + "   has_protein hp2 ON mp2.peptide_id = hp2.peptide_id INNER JOIN \n"
                    + "   protein pr1 ON hp1.protein_id = pr1.id LEFT OUTER JOIN \n"
                    + "   protein pr2 ON hp2.protein_id = pr2.id\n"
                    + "   INNER JOIN  \n"
                    + "   spectrum s ON sm.spectrum_id = s.id \n"
                    + " LEFT OUTER JOIN peaklistfile plf on s.peaklist_id = plf.id\n"
                    + " LEFT OUTER JOIN crosslinker cl on mp1.crosslinker_id = cl.id \n"
                    + " WHERE  \n"
                    + " (s.precursor_charge = sm.precursor_charge OR sm.precursor_charge < 6) \n"
                    + " GROUP BY  sm.id, \n"
                    + "p1.sequence, \n"
                    + "p2.sequence, \n"
                    + "p1.peptide_length, \n"
                    + "p2.peptide_length, \n"
                    + "mp1.link_position + 1, \n"
                    + "mp2.link_position + 1, \n"
                    + "pr1.is_decoy, \n"
                    + "pr2.is_decoy, \n"
                    + "sm.precursor_charge, \n"
                    + "sm.score, \n"
                    + "CASE WHEN p2.sequence IS NULL THEN 1 ELSE (4.0/5.0+(p1.peptide_length/(p1.peptide_length+p2.peptide_length)))/2 END\n"
                    + " , s.precursor_charge \n"
                    + " , p1.id "
                    + " , p2.id "
                    + " , v.run_name, v.scan_number \n"
                    + " , scorepeptide1matchedconservative \n"
                    + " , scorepeptide2matchedconservative \n"
                    + " , scoredelta "
                    + " , scorelinksitedelta "
                    + " , s.precursor_mz "
                    + " , sm.calc_mass\n"
                    + " , sm.precursor_charge "
                    + " , p1.mass "
                    + " , p2.mass "
                    + " , sm.search_id\n"
                    + " , sm.spectrum_id\n"
                    + " , s.precursor_charge\n"
                    + " , autovalidated\n"
                    + " , cl.name \n"
                    + " , cl.mass \n"
                    + " , sm.rank \n"
                    + " , s.scan_index\n"
                    + " , plf.name "
                    + " , scoremgxrank "
                    + " , scorecleavclpep1fragmatched "
                    + " , scorecleavclpep2fragmatched "
                    + " , scores "
                    + " , scoredelta "
                    + " , s.precursor_intensity\n"
                    + " , s.elution_time_start \n"
                    + " ORDER BY sm.score DESC)  i \n"
                    + (searchfilter == null || searchfilter.isEmpty() ? "" : " WHERE ( " + searchfilter + " )");

            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "read from db");
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, matchQuerry);
            if (searchIds.length >1)
                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Search {0} of {1}", new Object[]{currentsearch+1, searchIds.length});
            //getDBConnection().setAutoCommit(false);
            Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
            stm.setFetchSize(100);
            ResultSet rs = stm.executeQuery(matchQuerry);
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "go through the results");

            int psmIDColumn = rs.findColumn("psmID");
            int pepSeq1Column = rs.findColumn("pepSeq1");
            int pepSeq2Column = rs.findColumn("pepSeq2");
            int peplen1Column = rs.findColumn("peplen1");
            int peplen2Column = rs.findColumn("peplen2");
            int site1Column = rs.findColumn("site1");
            int site2Column = rs.findColumn("site2");
            int isDecoy1Column = rs.findColumn("isDecoy1");
            int isDecoy2Column = rs.findColumn("isDecoy2");
            int calc_chargeColumn = rs.findColumn("calc_charge");
            int scoreColumn = rs.findColumn("score");
            int accession1Column = rs.findColumn("accession1");
            int description1Column = rs.findColumn("description1");
            int accession2Column = rs.findColumn("accession2");
            int description2Column = rs.findColumn("description2");
            int pepPosition1Column = rs.findColumn("pepPosition1");
            int pepPosition2Column = rs.findColumn("pepPosition2");
            int score_ratioColumn = rs.findColumn("score_ratio");
            int exp_chargeColumn = rs.findColumn("exp_charge");
            int protein1idColumn = rs.findColumn("protein1id");
            int protein2idColumn = rs.findColumn("protein2id");
            int peptide1idColumn = rs.findColumn("peptide1id");
            int peptide2idColumn = rs.findColumn("peptide2id");
            int run_nameColumn = rs.findColumn("run_name");
            int scan_numberColumn = rs.findColumn("scan_number");
            int protein1sequenceColumn = rs.findColumn("protein1sequence");
            int protein2sequenceColumn = rs.findColumn("protein2sequence");
            int scoreP1CoverageColumn = rs.findColumn("scoreP1Coverage");
            int scoreP2CoverageColumn = rs.findColumn("scoreP2Coverage");
            int deltaScoreColumn = rs.findColumn("deltaScore");
            int LinkSiteDeltaColumn = rs.findColumn("LinkSiteDelta");
            int exp_mzColumn = rs.findColumn("exp_mz");
            int calc_massColumn = rs.findColumn("calc_mass");
            int pep1massColumn = rs.findColumn("pep1mass");
            int pep2massColumn = rs.findColumn("pep2mass");
            int search_idColumn = rs.findColumn("search_id");
            int spectrum_idColumn = rs.findColumn("spectrum_id");
            int xlColumn = rs.findColumn("crosslinker");
            int avColumn = rs.findColumn("autovalidated");
            int rankColumn = rs.findColumn("rank");
            int scanIndexColumn = rs.findColumn("scan_index");
            int peaklistfileColumn = rs.findColumn("peaklistfile");
            int mgxrankColumn = rs.findColumn("mgxrank");
            int cleavclpep1fragmatchedColumn = rs.findColumn("cleavclpep1fragmatched");
            int cleavclpep2fragmatchedColumn = rs.findColumn("cleavclpep2fragmatched");
            int subscoresColumn = rs.findColumn("subscores");
            int precIntensityColumn = rs.findColumn("precursor_intensity");
            int retentiontimeColumn = rs.findColumn("retentiontime");
            int clMassColumn = rs.findColumn("clmass");
            int clDeltaColumn = rs.findColumn("delta");
            int prot1nameColumn = rs.findColumn("protein1name");
            int prot2nameColumn = rs.findColumn("protein2name");
            Pattern modDetect = Pattern.compile(".*[^A-Z].*");
            HashSet<Double> tmmodcount = new HashSet<>();
            HashMap<Double,Double> xlmodmasses = new HashMap<Double,Double>(1);
            


            int total = 0;
            try {
                while (rs.next()) {
                    total++;

                    long ipsmID = rs.getLong(psmIDColumn);
                    if (skip.contains(ipsmID)) {
                        continue;
                    }

                    String psmID = Long.toString(ipsmID);
                    String pepSeq1 = rs.getString(pepSeq1Column);
                    String pepSeq2 = rs.getString(pepSeq2Column);
                    int peplen1 = rs.getInt(peplen1Column);
                    int peplen2 = rs.getInt(peplen2Column);
                    int site1 = rs.getInt(site1Column);
                    int site2 = pepSeq2 == null ? -1 : rs.getInt(site2Column);
                    boolean isDecoy1 = rs.getBoolean(isDecoy1Column);
                    boolean isDecoy2 = rs.getBoolean(isDecoy2Column);
                    int charge = rs.getInt(calc_chargeColumn);
                    double score = rs.getDouble(scoreColumn);
                    String[] accession1 = (String[])rs.getArray(accession1Column).getArray();
                    String[] name1 = (String[])rs.getArray(prot1nameColumn).getArray();
                    String[] description1 = (String[])rs.getArray(description1Column).getArray();
                    String[] accession2 = pepSeq2 == null ? new String[]{""} : (String[])rs.getArray(accession2Column).getArray();
                    String[] name2 = pepSeq2 == null ? new String[]{""} : (String[])rs.getArray(prot2nameColumn).getArray();
                    String[] description2 = pepSeq2 == null ? new String[]{""} : (String[])rs.getArray(description2Column).getArray();
                    Integer[] pepPosition1 = (Integer[])rs.getArray(pepPosition1Column).getArray();
                    Integer[] pepPosition2 = pepSeq2 == null ? new Integer[]{-1} : (Integer[])rs.getArray(pepPosition2Column).getArray();
                    //            double coverage1  = rs.getDouble(18);
                    //            double coverage2  = rs.getDouble(19);
                    //            double scoreRatio = coverage1/(coverage1+coverage2);
                    double scoreRatio = rs.getDouble(score_ratioColumn);
                    int spectrum_charge = rs.getInt(exp_chargeColumn);
                    Long[] protein1ID = (Long[])rs.getArray(protein1idColumn).getArray();
                    Long[] protein2ID = pepSeq2 == null ? new Long[]{0l} : (Long[])rs.getArray(protein2idColumn).getArray();
                    long pep1ID = rs.getLong(peptide1idColumn);
                    long pep2ID = rs.getLong(peptide2idColumn);
                    String run = rs.getString(run_nameColumn);
                    String scan = rs.getString(scan_numberColumn);
                    String[] sequence1 = (String[]) rs.getArray(protein1sequenceColumn).getArray();
                    String[] sequence2 = (String[]) rs.getArray(protein2sequenceColumn).getArray();
                    boolean autovalidated = rs.getBoolean(avColumn);
                    int rank = rs.getInt(rankColumn);

                    if (accession1 == null) {
                        accession1 = description1;
                    }

                    if (accession2 == null) {
                        accession2 = description2;
                    }

                    for (int s = 0 ; s < sequence1.length; s++)
                        if (sequence1[s] != null) {
                            sequence1[s] = sequence1[s].replaceAll("[^A-Z]", "");
                        }

                    for (int s = 0 ; s < sequence2.length; s++)
                        if (sequence2[s] != null) {
                            sequence2[s] = sequence2[s].replaceAll("[^A-Z]", "");
                        }

                    double p1c = rs.getDouble(scoreP1CoverageColumn);
                    double p2c = rs.getDouble(scoreP2CoverageColumn);
                    double pminc = p1c;
                    
                    if (pepSeq2 != null && !pepSeq2.isEmpty()) {
                        pminc = Math.min(p1c,p2c);
                    }
                    
                    double pmz = rs.getDouble(exp_mzColumn);
                    double f = 1;
                    if (pepSeq2 != null && !pepSeq2.isEmpty() && p1c + p2c > 0) {
                        ////                    double max = Math.max(p1c,p2c);
                        ////                    double min = Math.min(p1c,p2c);
                        ////                    f = min/(p1c+p2c);
                        ////                    score = score * f;
                        scoreRatio = (p1c) / (p1c + p2c + 1);
                        //                    if (p1c <3 || p2c <3) 
                        //                        continue;
                    }
                    
                    double peptide1score = score*scoreRatio;
                    double peptide2score = score*(1-scoreRatio);
                    

                    String xl = rs.getString(xlColumn);
                    if (xl == null) {
                        xl = "";
                    }
                    double calc_mass = rs.getDouble(calc_massColumn);
                    double pep1mass = rs.getDouble(pep1massColumn);
                    double pep2mass = rs.getDouble(pep2massColumn);
                    int search_id = rs.getInt(search_idColumn);
                    int scan_id = rs.getInt(spectrum_idColumn);
                    if (pepSeq2 != null && pepSeq2.matches("^X-?[0-9\\.]*$")) {
                        double mass = Double.parseDouble(pepSeq2.substring(1));
                        AminoModification am = new AminoModification(pepSeq2, AminoAcid.A, mass-18.0105647);
                        m_conf.addKnownModification(am);
                        getConfig(searchId).addKnownModification(am);
                        tmmodcount.add(((int)mass*10)/10.0);
                    }

                    int mgxrank = rs.getInt(mgxrankColumn);
                    boolean cleavclpep1fragmatched = rs.getBoolean(cleavclpep1fragmatchedColumn);
                    boolean cleavclpep2fragmatched = rs.getBoolean(cleavclpep2fragmatchedColumn);
                    java.sql.Array subscores = rs.getArray(subscoresColumn);

                    PSM psm = null;
                    if (pepSeq2 != null && !pepSeq2.isEmpty()){
                        for (int p = 0; p< accession1.length; p++) {
                            String a1=accession1[p];
                            String n1 = name1[p];
                            String d1 = description1[p];
                            if (d1==null || d1.isEmpty()) {
                                if (n1 == null || n1.isEmpty())
                                    d1 = a1;
                                else
                                    d1 = n1;
                            }
                            String a2 = accession2[p];
                            String n2 = name2[p];
                            String d2 = description2[p];
                            int p2 = pepPosition2[p];
                            long p2id = protein2ID[p];
                            String s2 = sequence2[p];
                            if (d2==null || d2.isEmpty()) {
                                if (n2 == null || n2.isEmpty())
                                    d2 = a2;
                                else
                                    d2 = n2;
                            }
                            String s1 = sequence1[p];
                            int p1 = pepPosition1[p];
                            long p1id = protein1ID[p];
                            psm = setUpDBPSM(psmID, run, scan, pep1ID, pep2ID, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, p1id, a1, d1, p2id, a2, d2, p1, p2, s1, s2, peptide1score, peptide2score, spectrum_charge, xl, pmz, calc_mass, pep1mass, pep2mass, search_id, scan_id);
                        }
                    } else {
                        String a2 = "";
                        String n2 = "";
                        String d2 = "";
                        int p2 = -1;
                        long p2id = 0;
                        String s2 = sequence2[0];
                        for (int p = 0; p< accession1.length; p++) {
                            String a1=accession1[p];
                            String n1 = name1[p];
                            String d1 = description1[p];
                            if (d2==null || d2.isEmpty()) {
                                if (n2 == null || n2.isEmpty())
                                    d2 = a2;
                                else
                                    d2 = n2;
                            }
                            String s1 = sequence1[p];
                            int p1 = pepPosition1[p];
                            long p1id = protein1ID[p];
                            psm = setUpDBPSM(psmID, run, scan, pep1ID, pep2ID, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, p1id, a1, d1, p2id, a2, d2, p1, p2, s1, s2, peptide1score, peptide2score, spectrum_charge, xl, pmz, calc_mass, pep1mass, pep2mass, search_id, scan_id);
                        }
                        
                    }
                    
                    psm.setRank(rank);
                    
                    Float[] scorevalues = (Float[])subscores.getArray();
                    if (cPepCoverage1 > 0) {
                        if (cPepCoverage2 <0 || Double.isNaN(scorevalues[cPepCoverage2]))
                            psm.addOtherInfo("minPepCoverage",
                                scorevalues[cPepCoverage1]);
                        else
                            psm.addOtherInfo("minPepCoverage",
                                Math.min(scorevalues[cPepCoverage1], 
                                        scorevalues[cPepCoverage2])
                        );
                    }

                    
                    psm.addOtherInfo("P1Fragments",p1c);
                    psm.addOtherInfo("P2Fragments",p2c);
                    psm.addOtherInfo("MinFragments",pminc);
                    if (pep2mass>0) {
                        psm.addOtherInfo("minPepCoverageAbsolute",
                            Math.min(p1c, p2c));
                    } else {
                        psm.addOtherInfo("minPepCoverageAbsolute", p1c);
                    }
                    
                    psm.addOtherInfo("cleavclpep1fragmatched",cleavclpep1fragmatched ? 1:0);
                    psm.addOtherInfo("cleavclpep2fragmatched",cleavclpep2fragmatched ? 1:0);
                    psm.addOtherInfo("cleavclpepfragmatched",(cleavclpep2fragmatched ? 1:0) + (cleavclpep1fragmatched ? 1:0));
                    psm.addOtherInfo("deltaScore",rs.getDouble(deltaScoreColumn));
                    psm.addOtherInfo("ScoreDivDelta",rs.getDouble(scoreColumn)/rs.getDouble(deltaScoreColumn));
                    psm.addOtherInfo("linkSiteDelta",rs.getDouble(LinkSiteDeltaColumn));
                    psm.addOtherInfo("mgxrank",rs.getDouble(mgxrankColumn));
                    psm.addOtherInfo("PrecursorIntensity",rs.getDouble(precIntensityColumn));
                    psm.addOtherInfo("RetentionTime",rs.getDouble(retentiontimeColumn));
                    
                    for (Integer p : peaks) {
                        psm.addOtherInfo(scorenames.get(p), scorevalues[p]);
                    }
                    
                    for (Integer p : scoresForwarded) {
                        psm.addOtherInfo(scorenames.get(p), scorevalues[p]);
                    }
                    
                    Double xlModmassPre = rs.getDouble(clMassColumn);
                    Double xlModmass = xlmodmasses.get(xlModmassPre);
                    if (xlModmass == null) {
                        xlmodmasses.put(xlModmassPre,xlModmassPre);
                        xlModmass = xlModmassPre;
                    }
                    psm.setCrosslinkerModMass(rs.getDouble(clMassColumn));
                    psm.setDeltaScore(rs.getDouble(deltaScoreColumn));

                    if  (autovalidated) {
                        psm.setAutoValidated(autovalidated);
                        if (flagAutoValidated) {
                            psm.setPositiveGrouping("AutoValidated");
                        }
                    }

                    String modLoockup = pepSeq1;
                    if (pepSeq2 != null)
                        modLoockup+=pepSeq2;
                    RunConfig psmconfig  =getConfig(psm.getSearchID());
                    Sequence m = new Sequence(modLoockup, psmconfig);
                    for (AminoAcid aa : m) {
                        if (aa instanceof AminoModification) {
                            if (psmconfig.getVariableModifications().contains(aa)) {
                                psm.setHasVarMods(true);
                            }
                            if (psmconfig.getFixedModifications().contains(aa)) {
                                psm.setHasFixedMods(true);
                            }
                        }
                    }
                    if (psm.hasVarMods()) {
                        if (markModifications())
                            psm.addNegativeGrouping("Modified");
                    }

                    if (markRun())
                        psm.getAdditionalFDRGroups().add(psm.getRun());

                    if (markSearchID())
                        psm.getAdditionalFDRGroups().add(""+psm.getSearchID());
                    
                    if (psm.hasVarMods()) {
                        if (markModifications())
                            psm.addNegativeGrouping("Modified");
                    }
                    

                    int scanIndex = rs.getInt(scanIndexColumn);
                    if (!rs.wasNull()) {
                        psm.setFileScanIndex(scanIndex);
                    }
                    String peaklistfile = rs.getString(peaklistfileColumn);
                    if (!rs.wasNull()) {
                        psm.setPeakListName(peaklistfile);
                    }
                    
//                    if (isMultpleSearches) {
//                        psm.addPositiveGrouping("" + searchId);
//                    }
                    if (total % 10000 == 0) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO, 
                                "go through the results ({0}) Search {1} of {2}", 
                                new Object[]{allPSMs.size(), currentsearch+1, searchIds.length});
                    }

                }
                m_search_ids.add(searchId);

            } catch (SQLException sex) {
                if (tries < 5) {
                    readDB(searchIds, searchfilter, topOnly, skip, tries + 1);
                }
                Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Repeatedly (" + total + ") failed to read from the database giving up now", sex);
                return;
            }
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Count results : " + total);
            
            if (tmmodcount.size()>10) {
                PeptidePair.ISTARGETED = false;
            }
        }

    }

    protected PSM setUpDBPSM(String psmID, String run, String scan, long pep1ID, long pep2ID, String pepSeq1, String pepSeq2, int peplen1, int peplen2, int site1, int site2, boolean isDecoy1, boolean isDecoy2, int charge, double score, long protein1ID, String accession1, String description1, long protein2ID, String accession2, String description2, int pepPosition1, int pepPosition2, String sequence1, String sequence2, double peptide1score, double peptide2score, int spectrum_charge, String xl, double pmz, double calc_mass, double pep1mass, double pep2mass, int search_id, int scan_id) {
        PSM psm = addMatch(psmID, run, scan, pep1ID, pep2ID, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, protein1ID, accession1, description1, protein2ID, accession2, description2, pepPosition1, pepPosition2, sequence1, sequence2, peptide1score, peptide2score, spectrum_charge == -1?"Unknow Charge" : null, xl);
        if (spectrum_charge == -1) {
            psm.setNegativeGrouping(" UnknownCharge");
        }
        
        psm.setExperimentalMZ(pmz);
        psm.setCalcMass(calc_mass);
        psm.setExpCharge((byte)spectrum_charge);
        Double pepMass = psm.getPeptide1().getMass();
        if (pepMass == null || pepMass <= 0 || Double.isNaN(pepMass)) {
            if (psm.getPeptide1().getSequence().contentEquals(pepSeq1)) {
                pepMass = pep1mass;
            } else {
                pepMass = pep2mass;
            }

            if (pepMass != null && pepMass != 0 && pepMass != Double.NaN) {

                psm.getPeptide1().setMass(pepMass);
            }
        }
        if (sequence2 != null) {
            pepMass = psm.getPeptide2().getMass();
            if (pepMass == null || pepMass <= 0 || Double.isNaN(pepMass)) {
                if (psm.getPeptide2().getSequence().contentEquals(pepSeq1)) {
                    pepMass = pep1mass;
                } else {
                    pepMass = pep2mass;
                }
                if (pepMass != null && pepMass != 0 && pepMass != Double.NaN) {
                    psm.getPeptide2().setMass(pepMass);
                }
            }
        }

        if (Double.isNaN(psm.getPeptide1().getMass()) || (sequence2 != null && Double.isNaN(psm.getPeptide2().getMass()))) {
            Logger.getLogger(this.getClass().getName()).log(Level.WARNING, "No mass for peptide");
        }

//        if (psm.getCalcMass() - psm.getPeptide1().getMass() - psm.getPeptide2().getMass() <0) {
//            Logger.getLogger(this.getClass().getName()).log(Level.WARNING, "wrong mass");
//        }
        psm.setSearchID(Integer.toString(search_id));
        psm.setScanID(Integer.toString(scan_id));
        return psm;
//                    if (!psm.isLinear()) {
//                        psm.setCrosslinker(rs.getString(xlColumn));
//                    }
    }


    public void getDBSizes() {
        int targetDBsize = 0;
        int decoyDBsize = 0;
        int targetProtDBsize = 0;
        int decoyProtDBsize = 0;
        int targetLinkDBsize = 0;
        int decoyLinkDBsize = 0;
        try {
            for (String id : m_search_ids) {

                Digestion digest = getConfig().getDigestion_method();
                HashSet<hashableXiPeptide> targetPeps = new HashSet<hashableXiPeptide>();
                HashSet<hashableXiPeptide> decoyPeps = new HashSet<hashableXiPeptide>();

                String dbQuerry = "SELECT sequence, is_Decoy FROM protein p inner join (SELECT DISTINCT protein_id FROM "
                        + " spectrum_match sm inner join "
                        + " matched_peptide mp ON sm.search_id = " + id + " AND sm.id = mp.match_id INNER JOIN "
                        + " has_protein hp on mp.peptide_id = hp.peptide_id) i ON p.id = i.protein_id;";

                //getDBConnection().setAutoCommit(false);
                Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
                stm.setFetchSize(100);
                ResultSet rs = stm.executeQuery(dbQuerry);
                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "go through the results");
                int total = 0;
                PeptideTree xltree = new PeptideTree(getConfig().getPrecousorTolerance());
                PeptideTree lineartree = new PeptideTree(getConfig().getPrecousorTolerance());
                while (rs.next()) {

                    Sequence seq = new Sequence(rs.getString(1), getConfig());
                    boolean isdecoy = rs.getBoolean(2);

                    seq.setDecoy(isdecoy);
                    Digestion d = getConfig().getDigestion_method();
                    d.setPeptideLookup(xltree, lineartree);
                    d.setMaxMissCleavages(getConfig().getMaxMissCleavages());

                    seq.digest(getConfig().getDigestion_method(), Double.MAX_VALUE, getConfig().getCrossLinker());
                    ArrayList<rappsilber.ms.sequence.Peptide> peps = digest.digest(seq, getConfig().getCrossLinker());
                    if (isdecoy) {
                        decoyProtDBsize++;
                        seqpos:
                        for (int i = 0; i < seq.length(); i++) {
                            for (CrossLinker cl : getConfig().getCrossLinker()) {
                                if (cl.canCrossLink(seq, i)) {
                                    decoyLinkDBsize++;
                                    continue seqpos;
                                }
                            }
                        }
                    } else {
                        targetProtDBsize++;
                        seqpos:
                        for (int i = 0; i < seq.length(); i++) {
                            for (CrossLinker cl : getConfig().getCrossLinker()) {
                                if (cl.canCrossLink(seq, i)) {
                                    targetLinkDBsize++;
                                    continue seqpos;
                                }
                            }
                        }

                    }
                }
                for (rappsilber.ms.sequence.Peptide peptide : xltree) {
                    if (peptide.isDecoy()) {
                        decoyDBsize++;
                    } else {
                        targetDBsize++;
                    }
                }
            }
            setTargetDBSize(targetDBsize);
            setDecoyDBSize(decoyDBsize);

            setTargetProtDBSize(targetProtDBsize);
            setDecoyProtDBSize(decoyProtDBsize);

            setTargetLinkDBSize(targetLinkDBsize);
            setDecoyLinkDBSize(decoyLinkDBsize);

        } catch (SQLException ex) {
            Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public void clearDBValidation(String validation) throws SQLException {
        if (!ensureConnection()) {
            return;
        }
        Statement stm = getDBConnection().createStatement();

        if (validation == null || validation.isEmpty()) {
            stm.execute("Update spectrum_match SET validated = null WHERE search_id in (" + RArrayUtils.toString(m_search_ids, ",") + ") AND validated is not null ");

        } else {
            String[] aval = validation.split(",");
            for (int i = 0; i < aval.length; i++) {
                aval[i] = aval[i].trim();
            }

            stm.execute("Update spectrum_match SET validated = null WHERE search_id in (" + RArrayUtils.toString(m_search_ids, ",") + ") AND validated in ('" + RArrayUtils.toString(aval, "','") + "')");
        }
        if (!getDBConnection().getAutoCommit()) {
            getDBConnection().commit();
        }

//        for (int sid : m_search_ids) {
//            stm.execute("Update spectrum_match SET validated = null WHERE search_id = " + sid  + " AND validated in ('" + RArrayUtils.toString(aval, "','") + "')");
//            if (!m_db_connection.getAutoCommit())
//                getDBConnection().commit();
//        }
    }

    public void clearDBFDR() throws SQLException {
        if (!ensureConnection()) {
            return;
        }
        Statement stm = getDBConnection().createStatement();
        for (String sid : m_search_ids) {
            stm.execute("Update spectrum_match SET peptide_fdr = null, link_fdr = null, ppi_fdr = null   WHERE search_id = " + sid);
            if (!m_db_connection.getAutoCommit()) {
                getDBConnection().commit();
            }
        }
    }

    public void writeDBPPIConfidence(FDRResult result) throws SQLException {
        //Collection<PeptidePair> peps = getFDRPeptidePairs();
        int count = 0;
        int maxUpdate = 100;
        if (!ensureConnection()) {
            return;
        }
        boolean autocomit = getDBConnection().getAutoCommit();
        getDBConnection().setAutoCommit(false);
        // also write the fdr-values
        for (PeptidePair pp : result.peptidePairFDR.filteredResults()) {
            ProteinGroupPair pgp = pp.getFdrLink().getFdrPPI();
            double confidence = 100 * (1 - pgp.getFDR());
            for (String psmid : pp.getPSMids()) {
                updateSetConfidence.setDouble(1, confidence);
                updateSetConfidence.setLong(2, Long.parseLong(psmid));
                updateSetConfidence.addBatch();
                if ((count = (count + 1) % maxUpdate) == 0) {
                    updateSetConfidence.executeBatch();
                    Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
                }
            }
        }
        if (count > 0) {
            updateSetConfidence.executeBatch();
            Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
        }
        getDBConnection().commit();
        getDBConnection().setAutoCommit(autocomit);
    }

    public void writeDBPepPairConfidence(FDRResult result) throws SQLException {
        //     Collection<PeptidePair> peps = getFDRPeptidePairs();
        int count = 0;
        int maxUpdate = 100;
        if (!ensureConnection()) {
            return;
        }
        boolean autocomit = getDBConnection().getAutoCommit();
        getDBConnection().setAutoCommit(false);
        // also write the fdr-values
        for (PeptidePair pp : result.peptidePairFDR.filteredResults()) {

            double confidence = 100 * (1 - pp.getFDR());
            for (String psmid : pp.getPSMids()) {
                updateSetConfidence.setDouble(1, confidence);
                updateSetConfidence.setLong(2, Long.parseLong(psmid));
                updateSetConfidence.addBatch();
                if ((count = (count + 1) % maxUpdate) == 0) {
                    updateSetConfidence.executeBatch();
                    Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
                }
            }
        }
        if (count > 0) {
            updateSetConfidence.executeBatch();
            Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
        }
        getDBConnection().commit();
        getDBConnection().setAutoCommit(autocomit);
    }

    public void writeDBLinkConfidence(FDRResult result) throws SQLException {
        //Collection<PeptidePair> peps = getFDRPeptidePairs();
        int count = 0;
        int maxUpdate = 100;
        if (!ensureConnection()) {
            return;
        }
        boolean autocomit = getDBConnection().getAutoCommit();
        getDBConnection().setAutoCommit(false);
        // also write the fdr-values
        for (PeptidePair pp : result.peptidePairFDR.filteredResults()) {
            if (!pp.isLinear()) {
                ProteinGroupLink l = pp.getFdrLink();
                double confidence = 100 * (1 - l.getFDR());
                for (String psmid : pp.getPSMids()) {
                    updateSetConfidence.setDouble(1, confidence);
                    updateSetConfidence.setLong(2, Long.parseLong(psmid));
                    updateSetConfidence.addBatch();
                    if ((count = (count + 1) % maxUpdate) == 0) {
                        updateSetConfidence.executeBatch();
                        Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
                    }
                }
            }
        }
        if (count > 0) {
            updateSetConfidence.executeBatch();
            Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
        }
        getDBConnection().commit();
        getDBConnection().setAutoCommit(autocomit);
    }

    public void writeDB(String validate, boolean overwrite, boolean writeFDR, FDRResult result, boolean within, boolean  between) throws SQLException {
        //Collection<PeptidePair> peps = getFDRPeptidePairs();
        if (!ensureConnection()) {
            return;
        }
        int count = 0;
        int maxUpdate = 100;
        boolean autocomit = getDBConnection().getAutoCommit();
        getDBConnection().setAutoCommit(false);
        flagDBUsage();
        if (validate != null) {
            PreparedStatement stVal = null;
            if (overwrite) { // overwrite validation
                // clear all previous validations of this label
                clearDBValidation(validate);
                stVal = updateValidateOverWrite;
            } else {
                stVal = updateValidateNonOverWrite;
            }
            // write out validations - but only if the match was not already validated
            for (PSM psm : result.psmFDR.filteredResults()) {
                Long psmid = Long.parseLong(psm.getPsmID());
//                    for (PeptidePair pp : result.peptidePairFDR) {
//                        for (String psmid : pp.getPSMids()) {
                if (psm.isLinear() || (psm.isBetween() && between) || (psm.isInternal() && within)) {
                    stVal.setString(1, validate);
                    stVal.setLong(2, psmid);
                    stVal.addBatch();
                    if ((count = (count + 1)) % maxUpdate == 0) {
                        stVal.executeBatch();
                        Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
                        count = 0;
                     }
                }
            }

            if (count % maxUpdate > 0) {
                stVal.executeBatch();
                Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
            }
        }
        getDBConnection().commit();
        getDBConnection().setAutoCommit(autocomit);
        flagDBUsage();
    }

//    /**
//     * adds a psm to the list folds up the scores to peptidespairs links
//     * proteinpairs and proteins
//     *
//     * @param psmID
//     * @param pepid1
//     * @param pepid2
//     * @param peplen1
//     * @param peplen2
//     * @param site1
//     * @param site2
//     * @param charge
//     * @param score
//     * @param proteinId1
//     * @param proteinId2
//     * @param pepPosition1
//     * @param pepPosition2
//     * @param scoreRation
//     * @return a peptide pair that is supported by the given match
//     */
//    @Override
//    public PSM addMatch(String psmID, org.rappsilber.fdr.entities.Peptide peptide1, org.rappsilber.fdr.entities.Peptide peptide2, int peplen1, int peplen2, int site1, int site2, int charge, double score, Protein proteinId1, Protein proteinId2, int pepPosition1, int pepPosition2, double peptide1score, double peptide2score, String isSpecialCase, String crosslinker, String run, String scan) {
//        org.rappsilber.fdr.entities.Peptide npepid1;
//        org.rappsilber.fdr.entities.Peptide npepid2;
//        int npeplen1;
//        int npeplen2;
//        byte nsite1;
//        byte nsite2;
//        Protein nproteinId1;
//        Protein nproteinId2;
//        int npepPosition1;
//        int npepPosition2;
//        int protcomp = proteinId1.compareDecoyUnAware(proteinId2);
//        int pepcomp = peptide1.compare(peptide2);
//        int sitecomp = (site1 - site2);
//        double pep1score = peptide1score;
//        double pep2score = peptide2score;
//        double npep1score;
//        double npep2score;
//
////        if (protcomp < 0 || (protcomp == 0 && pepcomp < 0) || (protcomp == 0 && pepcomp == 0 && site1 < site2)) {
//            npepid1 = peptide1;
//            npepid2 = peptide2;
//            npeplen1 = peplen1;
//            npeplen2 = peplen2;
//            nsite1 = (byte)site1;
//            nsite2 = (byte)site2;
//            nproteinId1 = proteinId1;
//            nproteinId2 = proteinId2;
//            npepPosition1 = pepPosition1;
//            npepPosition2 = pepPosition2;
//            npep1score = pep1score;
//            npep2score = pep2score;
//
////        } else {
////            npepid1 = peptide2;
////            npepid2 = peptide1;
////            npeplen1 = peplen2;
////            npeplen2 = peplen1;
////            nsite1 = (byte)site2;
////            nsite2 = (byte)site1;
////            nproteinId1 = proteinId2;
////            nproteinId2 = proteinId1;
////            npepPosition1 = pepPosition2;
////            npepPosition2 = pepPosition1;
////            npep1score = pep2score;
////            npep2score = pep1score;
////        }
//
//        if (!PSMScoreHighBetter) {
//            score = 10 - (10 * score);
//        }
//
//        PSM psm = new PSM(psmID, npepid1, npepid2, nsite1, nsite2, proteinId1.isDecoy(), proteinId2.isDecoy(), (byte)charge, score, npep1score, npep2score);
//        psm.setNegativeGrouping(isSpecialCase);
//        psm.setRun(registerRun(run));
//        psm.setScan(scan);
//        psm.setCrosslinker(crosslinker);
//
//        PSM regpsm DBPSM= getAllPSMs().register(psm);
//
//        return regpsm;
//    }

    public String argList() {
        return "--dbids=A,B,C --connection= " + super.argList() + " --crosslinkonly --Filter --flagModifications";
    }

    public String argDescription() {
        return super.argDescription() + "\n"
                + "--dbids=A,B,C            List of search ids\n"
                + "--connection=            connection string to be used to\n"
                + "                         access the database\n "
                + "--crosslinkonly          only load cross-links\n"
                + "--filter                 a filter to be applied while loading\n"
                + "                         matches from the database\n"
                + "--user=                  db user name (ignore for now)\n"
                + "--password=              db password (ignore for now)\n"
                + "--flagModifications     should modified peptide make their own sub-group\n";

    }

    public String[] parseArgs(String[] argv, FDRSettings settings) {
        ArrayList<String> unknown = new ArrayList<String>();
        String ConnString = null;
        ArrayList<String> filters = new ArrayList<String>();
        String dbUserName = "";
        String dbPassword = "";
        
        argv = super.parseArgs(argv, settings);

        for (String arg : argv) {
            if (arg.toLowerCase().startsWith("--dbids=")) {

                String[] sids = arg.substring(arg.indexOf("=") + 1).trim().split(",");
                searchIDSetting = new String[sids.length];

                for (int i = 0; i < sids.length; i++) {
                    searchIDSetting[i] = sids[i].trim();
                }

            } else if (arg.toLowerCase().startsWith("--filter=")) {
                filters.add(arg.substring(arg.indexOf("=") + 1).trim());

            } else if(arg.toLowerCase().startsWith("--forward=")) {
                String forwardPattern=arg.substring("--forward=".length());
                setSubScoresToForward(forwardPattern);

            } else if (arg.toLowerCase().equals("--crosslinkonly") || arg.equals("-C")) {
                filters.add("p2.accessen_number is not null");

            } else if (arg.toLowerCase().startsWith("--connection=")) {
                ConnString = arg.substring(arg.indexOf("=") + 1).trim();
            } else if (arg.toLowerCase().startsWith("--user=")) {
                dbUserName = arg.substring(arg.indexOf("=") + 1).trim();
            } else if (arg.toLowerCase().startsWith("--password=")) {
                dbPassword = arg.substring(arg.indexOf("=") + 1).trim();
            }else if (arg.toLowerCase().startsWith("--flagmodifications")) {
                setMarkModifications(true);
            }else if (arg.toLowerCase().equals("--lastowner")) {
                lastMzIDowner = true;
            } else if (arg.toLowerCase().startsWith("--validate=")) {
                command_line_auto_validate = arg.substring(arg.indexOf("=") + 1).trim();
            } else {
                unknown.add(arg);
            }
            

        }
        if (ConnString == null) {
            Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, "No db-connection defined");
            printUsage();
            System.exit(0);
        } else {
            Exception connectionException = null;
            try {
                setDBConnection(ConnString, dbUserName, dbPassword);
            } catch (SQLException ex) {
                connectionException = ex;
            }
            if (connectionException != null) {
                DBConnectionConfig dbc = new DBConnectionConfig();
                try {    
                    dbc.readConfig();

                    for (DBConnectionConfig.DBServer s: dbc.getServers()) {

                        if (s.getName().contentEquals(ConnString)) {
                            connectionException = null;
                            ConnString = s.getConnectionString();
                            dbUserName = s.getUser();
                            dbPassword = s.getPassword();

                            try {
                                setDBConnection(ConnString, dbUserName, dbPassword);
                                break;
                            } catch (SQLException ex) {
                                connectionException = ex;
                            }
                        }
                    }
                } catch (FileNotFoundException ex) {
                } catch (ParseException ex) {
                    Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, "Error parsing db config", ex);
                }
                    
            }
            
            if (connectionException != null) {
                Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, null, connectionException);
                printUsage();
                System.exit(0);
            }
        }

        if (searchIDSetting == null) {
            Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, "No search ids defined");
            printUsage();
            System.exit(0);
        }

        filterSetting = RArrayUtils.toString(filters, " AND ");
        String[] ret = new String[unknown.size()];
        ret = unknown.toArray(ret);
        return ret;
    }

    public static void main(String[] argv) throws SQLException, FileNotFoundException {

        DBinFDR ofdr = new DBinFDR();

        FDRSettings settings = new FDRSettingsImpl();
        
        String[] unknowns = ofdr.parseArgs(argv, settings);

        if (unknowns.length > 0) {
            for (String a : unknowns) {
                Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, "Unknown argument : " + a);
            }
            ofdr.printUsage();
            System.exit(1);
        }

        if (ofdr.getCsvOutDirSetting() != null) {
            if (ofdr.getCsvOutBaseSetting() == null) {
                ofdr.setCsvOutBaseSetting("FDR_" + RArrayUtils.toString(ofdr.searchIDSetting, "_"));
            }
        }

        if (ofdr.getCsvOutBaseSetting() != null) {
            if (ofdr.getCsvOutDirSetting() == null) {
                ofdr.setCsvOutDirSetting(".");
            }

            System.out.println("writing results to " + ofdr.getCsvOutDirSetting() + "/" + ofdr.getCsvOutBaseSetting() + "*");
            System.out.flush();
        }

        Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, "Read datafrom db");
        ofdr.readDB();

        final CalculateWriteUpdate cu = new CalculateWriteUpdate() {
            @Override
            public void setStatus(MaximisingStatus state) {
            }

            @Override
            public void setStatusText(String text) {
            }

            @Override
            public void reportError(String text, Exception ex) {
                Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.SEVERE, text, ex);
            }

            @Override
            public void setCurrent(double psm, double peptidepair, double protein, double link, double ppi) {
                Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.INFO, "next round: PSM FDR: " +psm + " PepPairFDR:" + peptidepair + " Protein FDR:"+protein + " LinkFDR:" + link + " PPI FDR:" + ppi);
            }

            @Override
            public void setComplete() {
                Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.INFO, "Calculate Write Finished");
            }

            @Override
            public boolean stopped() {
                return false;
            }
        };
        Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, "Calculate FDR");
        FDRResult result = ofdr.calculateWriteFDR(ofdr.getCsvOutDirSetting(), ofdr.getCsvOutBaseSetting(), ",", settings, cu);

//        Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, "Result-summary:" + ofdr.summaryString());
//        
//        if (ofdr.getCsvOutBaseSetting() != null) {
//            
//            ofdr.writeFiles(ofdr.getCsvOutDirSetting(), ofdr.getCsvOutDirSetting(), ",");
//        }
        //System.out.println(ofdr.summaryString(result));
        if (ofdr.command_line_auto_validate != null) {
            ofdr.writeDB(ofdr.command_line_auto_validate, true, false, result, true, true);
        }
        
        if (ofdr.singleCalculation()) {
            // we should be able to write out an mzIdentML
            String path =ofdr.getCsvOutDirSetting();
            String baseName = ofdr.getCsvOutBaseSetting();
            String out = path + (path.endsWith(File.separator)?"":File.separator) + baseName + ".mzid";
            try {
                MZIdentMLOwner o = new MZIdentMLOwner();
                
                if (ofdr.lastMzIDowner) {
                    o.readLast();
                } else {
                    o = MZIdentMLOwnerGUI.getSetOwner(o);
                }
                Logger.getLogger(MethodHandles.lookup().lookupClass().getName()).log(Level.INFO, "Writing mzIdentML to " + out);
                MZIdentMLExport mze = new MZIdentMLExport(ofdr, result, out, MZIdentMLOwnerGUI.getSetOwner(new MZIdentMLOwner()));
            } catch (Exception ex) {
                Logger.getLogger(MethodHandles.lookup().lookupClass().getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        System.exit(0);

    }

    @Override
    protected ArrayList<String> getLinkOutputLine(ProteinGroupLink l) {
        ArrayList<String> ret = super.getLinkOutputLine(l);
        HashSet<String> searchIDs = new HashSet<>();
        ret.add(getLinkWindow(l, 0, 20));
        ret.add(getLinkWindow(l, 1, 20));
        Collection<PeptidePair> cpeppairs = l.getPeptidePairs();
        HashSet<String> peps = new HashSet<String>();
        for (PeptidePair pp : cpeppairs) {
            for (org.rappsilber.fdr.entities.Peptide p : pp.getPeptides()) {
                if (p.getSequence().matches("(-\\.)?X-?[0-9.]+(\\.-)?")) {
                    peps.add(p.getSequence());
                }
            }
            for (PSM psm : pp.getAllPSMs()) {
                for (PSM upsm : psm.getRepresented()) {
                    searchIDs.add(upsm.getSearchID());
                }
                //searchIDs.add(((DBPSM) psm).getSearchID());
            }
        }

        ret.add(0, RArrayUtils.toString(searchIDs, ";"));

        if (!peps.isEmpty()) {
            ret.add(RArrayUtils.toString(peps, ";"));
        }

        if (l.getProteinGroup2().accessions().matches("'?OpenModification.*")) {
            ret.add(getLinkWindow(l, 0, 20));
            HashMap<Protein, IntArrayList> pos = l.getPosition1();
            Protein p = pos.keySet().iterator().next();
            ret.add(p.getAccession());
            int pp = pos.get(p).get(0);
            ret.add("" + pp);

        } else if (l.getProteinGroup1().accessions().matches("'?OpenModification.*")) {
            ret.add(getLinkWindow(l, 1, 20));
            HashMap<Protein, IntArrayList> pos = l.getPosition2();
            Protein p = pos.keySet().iterator().next();
            ret.add(p.getAccession());
            int pp = pos.get(p).get(0);
            ret.add("" + pp);
        } else {
            ret.add("");
            ret.add("");
            ret.add("");
        }

        return ret;
    }

    protected ArrayList<String> getLinkOutputHeader() {
        ArrayList<String> ret = super.getLinkOutputHeader();
        ret.add(0, "SearchIDs");
        ret.add("LinkWindow1");
        ret.add("LinkWindow2");
        ret.add("TM_OM_Mods");
        ret.add("TM_OM_ModWindow");
        ret.add("TM_OM_RepresentiveProtein");
        ret.add("TM_OM_RepresentiveModSite");
        return ret;
    }

    @Override
    protected ArrayList<String> getPSMHeader() {
        ArrayList<String> ret = super.getPSMHeader();
//        String header = "PSMID" + seperator + "run" + seperator + "scan" 
//                + seperator + "exp charge"  + seperator + "exp m/z" + seperator + "exp mass" + seperator + "exp fractionalmass" 
//                + seperator + "match charge" + seperator + "match mass" + seperator + "match fractionalmass"   
//                
//                + seperator + "Protein1" + seperator + "Description1" + seperator + "Decoy1" + seperator + "Protein2" + seperator + "Description2" + seperator + "Decoy2" + seperator + "Peptide1" + seperator + "Peptide2" + seperator + "PeptidePosition1" + seperator + "PeptidePosition2" + seperator + "PeptideLength1" + seperator + "PeptideLength2" + seperator + "fromSite" + seperator + "ToSite" + seperator + "Charge" + seperator + "Score" + seperator + "isDecoy" + seperator + "isTT" + seperator + "isTD" + seperator + "isDD" + seperator + "fdrGroup" + seperator + "fdr" + seperator + "" + seperator + "PeptidePairFDR" + seperator + "Protein1FDR" + seperator + "Protein2FDR" + seperator + "LinkFDR" + seperator + "PPIFDR" + seperator + seperator + "peptide pair id" + seperator + "link id" + seperator + "ppi id";
        ret.add(0, "SearchID");
        ret.add(3, "exp charge");
        ret.add(4, "exp m/z");
        ret.add(5, "exp mass");
        ret.add(6, "exp fractionalmass");
        ret.add(7, "match charge");
        ret.add(8, "match mass");
        ret.add(9, "match fractionalmass");
        return ret;
    }

    protected ArrayList<String> getPSMOutputLine(PSM pp) {
        ArrayList<String> ret = super.getPSMOutputLine(pp);
        double mz = pp.getExperimentalMZ();
        int charge = pp.getCharge();
        double mass = (mz - 1.00727646677) * (pp.getExpCharge() == -1 ? charge : pp.getExpCharge());
        double fraction = mass - Math.floor(mass);
        double calcfraction = pp.getCalcMass() - Math.floor(pp.getCalcMass());

        ret.add(0, "" + pp.getSearchID());
        ret.add(3, "" + pp.getExpCharge());
        ret.add(4, sixDigits.format(mz));
        ret.add(5, sixDigits.format(mass));
        ret.add(6, sixDigits.format(fraction));
        ret.add(7, "" + charge);
        ret.add(8, sixDigits.format(pp.getCalcMass()));
        ret.add(9, sixDigits.format(calcfraction));

        return ret;
    }

    @Override
    protected ArrayList<String> getPPIOutputHeader() {
        ArrayList<String> ret = super.getPPIOutputHeader();
        ret.add(0, "SearchID");
        return ret;
    }

    @Override
    protected ArrayList<String> getPPIOutputLine(ProteinGroupPair pgp) {
        HashSet<String> searches = new HashSet<>();
        for (PeptidePair pp : pgp.getPeptidePairs()) {
            for (PSM psm : pp.getAllPSMs()) {
                for (PSM upsm : psm.getRepresented()) {
                    searches.add(upsm.getSearchID());
                }
//                searches.add(((DBPSM)psm).getSearchID());
            }
        }
        ArrayList<String> ret = super.getPPIOutputLine(pgp);
        ret.add(0, RArrayUtils.toString(searches, ";"));
        return ret;
    }

    @Override
    protected ArrayList<String> getProteinGroupOutputHeader() {
        ArrayList<String> ret = super.getProteinGroupOutputHeader();
        ret.add(0, "SearchID");
        return ret;
    }

    @Override
    protected ArrayList<String> getProteinGroupOutput(ProteinGroup pg) {
        HashSet<String> searches = new HashSet<>();
        for (PeptidePair pp : pg.getPeptidePairs()) {
            for (PSM psm : pp.getAllPSMs()) {
                for (PSM upsm : psm.getRepresented()) {
                    searches.add(upsm.getSearchID());
                }
//                searches.add(((DBPSM)psm).getSearchID());
            }
        }
        ArrayList<String> ret = super.getProteinGroupOutput(pg);
        ret.add(0, RArrayUtils.toString(searches, ";"));
        return ret;
    }

    protected String getPeptideSequence(org.rappsilber.fdr.entities.Peptide p) {
        if (!m_writePrePostAA) {
            return p.getSequence();
        }

        ProteinGroup pg = p.getProteinGroup();
        HashMap<org.rappsilber.fdr.entities.Protein, HashSet<Integer>> positions = p.getPositions();
        int len = p.length();
        HashSet<String> N = new HashSet<String>(positions.size());
        HashSet<String> C = new HashSet<String>(positions.size());

        for (org.rappsilber.fdr.entities.Protein prot : positions.keySet()) {
            Sequence seq = getSequence(prot);
            if (seq.length() == 0) {
                N.add("-");
                C.add("-");
            } else {
                for (int pos : positions.get(prot)) {
                    if (pos == 1) {
                        N.add("-");
                    } else {
                        N.add(seq.aminoAcidAt(pos - 2).toString());
                    }
                    if (pos + len > seq.length()) {
                        C.add("-");
                    } else {
                        C.add(seq.aminoAcidAt(pos + len - 1).toString());
                    }
                }
            }
        }

        return RArrayUtils.toString(N, "") + "." + p.getSequence() + "." + RArrayUtils.toString(C, "");
    }

    @Override
    protected ArrayList<String> getXlPepeptideOutputLine(PeptidePair pp) {
        ArrayList<String> sret = super.getXlPepeptideOutputLine(pp);
        ArrayList<String> ret = new ArrayList<String>(sret.size() + 5);
        HashSet<String> searches = new HashSet<>();
        for (PSM psm : pp.getAllPSMs()) {
            for (PSM upsm : psm.getRepresented()) {
                searches.add(upsm.getSearchID());
            }
        }
        ret.add(RArrayUtils.toString(searches, ";"));
        ret.addAll(sret);

        ret.add(getLinkWindow(pp, 0, 20));
        ret.add(getLinkWindow(pp, 1, 20));
        if (pp.getPeptide1().getSequence().matches("X[\\-0-9]+(\\.[0-9]+)?")) {
            ret.add(getLinkWindow(pp, 1, 20));
            ret.add(getLinkWindow(pp.getLink(), 1, 20));
        } else if (pp.getPeptide2().getSequence().matches("X[\\-0-9]+(\\.[0-9]+)?")) {
            ret.add(getLinkWindow(pp, 0, 20));
            ret.add(getLinkWindow(pp.getLink(), 0, 20));
        }
        return ret;
    }

    @Override
    protected ArrayList<String> getXLPepsHeader() {
        ArrayList<String> ret = super.getXLPepsHeader();
        ret.add(0, "SearchIDs");
        ret.add("LinkWindow1");
        ret.add("LinkWindow2");
        ret.add("TM_OM_ModWindow");
        ret.add("TM_OM_ModWindow_noVarMod");
        return ret;
    }

    @Override
    protected ArrayList<String> getLinearPepeptideOutputLine(PeptidePair pp) {
        ArrayList<String> ret = super.getLinearPepeptideOutputLine(pp);
        HashSet<String> searches = new HashSet<>();
        for (PSM psm : pp.getAllPSMs()) {
            for (PSM upsm : psm.getRepresented()) {
                searches.add(upsm.getSearchID());
            }
        }
        ret.add(0, RArrayUtils.toString(searches, ";"));
        return ret;
    }

    @Override
    protected ArrayList<String> getLinearPepsHeader() {
        ArrayList<String> ret = super.getLinearPepsHeader();
        ret.add(0, "SearchIDs");
        return ret;
    }

    public int getFastas(Collection<String> searchID, ArrayList<String> dbIDs, ArrayList<String> names) {
        int c = 0;
        ArrayList<String> i = new ArrayList<>();
        ArrayList<String> n = new ArrayList<String>();
        try {
            Statement s = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
            ResultSet rs = s.executeQuery("SELECT  DISTINCT file_name, id FROM search_sequencedb ss INNER JOIN sequence_file sf on ss.seqdb_id = sf.id WHERE search_id in (" + RArrayUtils.toString(m_search_ids, ",") + ");");
            while (rs.next()) {
                i.add(rs.getInt(2)+"");
                n.add(rs.getString(1));
                c++;
            }
        } catch (SQLException ex) {
            Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, null, ex);
            return -1;
        }
        dbIDs.addAll(i);
        names.addAll(n);
        return c;
    }

    public int getFastas(String searchID, ArrayList<String> dbIDs, ArrayList<String> names) {
        ArrayList<String> id = new ArrayList<>(1);
        id.add(searchID);
        return getFastas(id, dbIDs, names);
    }

    
    /**
     * @return the flagAutoValidated
     */
    public boolean isFlagAutoValidated() {
        return flagAutoValidated;
    }

    /**
     * @param flagAutoValidated the flagAutoValidated to set
     */
    public void setFlagAutoValidated(boolean flagAutoValidated) {
        this.flagAutoValidated = flagAutoValidated;
    }

    /**
     * @return the markModifications
     */
    public boolean markModifications() {
        return markModifications;
    }

    /**
     * @param mm should modification be marked
     */
    public void setMarkModifications(boolean mm) {
        this.markModifications = mm;
    }

    @Override
    public void cleanup() {
        super.cleanup(); //To change body of generated methods, choose Tools | Templates.
        m_db_autoclosetimer.cancel();
    }

    /**
     * @return the subScoresToForward
     */
    public Pattern getSubScoresToForward() {
        return subScoresToForward;
    }

    /**
     * @param subScoresToForward the subScoresToForward to set
     */
    public void setSubScoresToForward(Pattern subScoresToForward) {
        this.subScoresToForward = subScoresToForward;
    }
    
    /**
     * @param subScoresToForward the subScoresToForward to set
     */
    public void setSubScoresToForward(String subScoresToForward) {
        this.subScoresToForward = Pattern.compile(subScoresToForward);
    }

    @Override
    public int getxiMajorVersion() {
        return 1;
    }


    @Override
    public Version getXiVersion() {
        if (this.m_xi_versions.size()>0)
            return this.m_xi_versions.values().iterator().next();
        return null;
    }

    private Version read_version(String search_id) {
        Connection con = getDBConnection();
        Version v = null;
        try {
            Statement stm = con.createStatement();
            ResultSet rs = stm.executeQuery("SELECT version from search s left outer join xiversion v on s.xiversion = v.id WHERE s.id = " + search_id);
            if (rs.next()) {
                v = new Version(rs.getString(1));
            }
            rs.close();
            stm.close();
        } catch (Exception e) {
            v = new Version("1.unknown");
        }
        return v;
    }
    
    @Override
    public Version getXiVersion(String search) {
        if (m_xi_versions == null) {
            m_xi_versions = new HashMap<>();
        }
        Version ret = m_xi_versions.get(search);
        if (ret == null) {
            ret = read_version(search);
            m_xi_versions.put(search, ret);
        }
        return ret;
    }    
    
    public HashMap<String, Version> getXiVersions() {
        return this.m_xi_versions;
    }

    public HashMap<String,? extends RunConfig> getConfigs() {
        return this.m_configs;
    }
    

    @Override
    public void add(OfflineFDR other) {
        super.add(other);
        if (other instanceof XiInFDR) {
            this.getSearchIDs().addAll(((XiInFDR)other).getSearchIDs());
            for (Map.Entry<String, ? extends RunConfig> e : ((XiInFDR)other).getConfigs().entrySet()) {
                RunConfig c = e.getValue();
                this.m_configs.put(e.getKey(), c);
            }
            this.getXiVersions().putAll(((XiInFDR)other).getXiVersions());
        }
    }

        
}
