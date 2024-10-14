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
import java.util.UUID;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;
import org.rappsilber.config.LocalProperties;
import org.rappsilber.fdr.dataimport.Xi2Xi1Config;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.fdr.entities.PeptidePair;
import org.rappsilber.fdr.entities.Protein;
import org.rappsilber.fdr.entities.ProteinGroupLink;
import org.rappsilber.fdr.entities.ProteinGroupPair;
import org.rappsilber.fdr.entities.ProteinGroup;
import org.rappsilber.fdr.gui.components.GetSearch;
import org.rappsilber.fdr.gui.components.MZIdentMLOwnerGUI;
import org.rappsilber.fdr.utils.CalculateWriteUpdate;
import org.rappsilber.fdr.utils.MZIdentMLExport;
import org.rappsilber.fdr.utils.MZIdentMLOwner;
import org.rappsilber.fdr.utils.MaximisingStatus;
import rappsilber.ms.sequence.AminoAcid;
import rappsilber.ms.sequence.Peptide;
import rappsilber.ms.sequence.Sequence;
import org.rappsilber.utils.IntArrayList;
import org.rappsilber.utils.RArrayUtils;
import org.rappsilber.utils.SelfAddHashSet;
import org.rappsilber.utils.Unescape;
import org.rappsilber.utils.Version;
import rappsilber.config.DBConnectionConfig;
import rappsilber.config.RunConfig;
import rappsilber.gui.components.db.DatabaseProvider;
import rappsilber.ms.statistics.utils.UpdateableLong;

/**
 *
 * @author lfischer
 */
 public class DB2inFDR extends org.rappsilber.fdr.OfflineFDR implements XiInFDR{
     public static final String summary_marker_long = "%leave_this_here_to_get_the_summary_added_to_the_notes%";
     public static String summary_marker = "%summarys%";
    private HashMap<String, Version> m_xi_versions = new HashMap<>();

    public class Xi2Score {
        int id;
        String name;
        boolean primary_score;
        boolean higher_is_better;
    
        public Xi2Score(int id, String name, boolean primary_score, boolean higher_is_better) {
            this.id = id;
            this.name = name;
            this.primary_score = primary_score;
            this.higher_is_better = higher_is_better;
        }

        @Override
        public String toString() {
            ArrayList<String> qualifiers = new ArrayList<>();
            
            if (primary_score)
                qualifiers.add("*");

            if (higher_is_better)
                qualifiers.add("h>l");
            else
                qualifiers.add("l>h");
            
            
            return name + "(" + RArrayUtils.toString(qualifiers, ",") + ")";
        }
        
        
    }
    
    public class Xi2ScoreList extends ArrayList<Xi2Score> {
        Xi2Score primaryScore;
        HashMap<String, Xi2Score> nameToScore = new HashMap<>();

        public Xi2ScoreList() {
        }

        public Xi2ScoreList(Xi2ScoreList source) {
            this.primaryScore = source.primaryScore;
            this.addAll(source);
        }

        @Override
        public boolean add(Xi2Score e) {
            ensureCapacity(e.id);
            for (int index=size(); index<=e.id; index++) {
                super.add(null);
            }
            super.set(e.id, e);
            nameToScore.put(e.name, e);
            if (e.primary_score) {
                this.primaryScore = e;
            }
            return true;
        }
        
        public Xi2Score get(String name) {
            return nameToScore.get(name);
        }
        
    }
    
    public class ResultMatchData {
        java.sql.Array link_score;
        java.sql.Array link_score_site1;
        java.sql.Array link_score_site2;
        boolean topranking;

        public ResultMatchData(Array link_score, Array link_score_site1, Array link_score_site2, boolean topranking) {
            this.link_score = link_score;
            this.link_score_site1 = link_score_site1;
            this.link_score_site2 = link_score_site2;
            this.topranking = topranking;
        }
        
    }
    private final static class resultMatchColumns{
        final static String querry = "INSERT INTO resultmatch (resultset_id , match_id, search_id, scores, top_ranking, site1, site2, link_score, link_score_site1, link_score_site2, match_group_id) "
                                   + "VALUES                  (? ,            ?,         ?,         ?,      ?,          ?,     ?,      ? ,        ?,                 ?,               ?)";
        final static int resultset_id =1;
        final static int match_id = 2;
        final static int search_id = 3;
        final static int scores = 4;
        final static int top_ranking  =5;
        final static int site1  = 6;
        final static int site2 = 7;
        final static int link_score = 8;
        final static int link_score_site1 = 9;
        final static int link_score_site2 = 10;       
        final static int match_group_id = 11;
    }
    
    
    public static final String subscoreroperty = "xiFDR.FILTER_SUBSCORES";
    private Connection m_db_connection;
    private Connection m_db_connectionAdmin;
    private UpdateableLong    m_db_last_used  =new UpdateableLong(0);
    private AtomicBoolean       keepConnection  = new AtomicBoolean(false);
    private Timer      m_db_autoclosetimer;
    private PreparedStatement updateValidateOverWrite;
    private PreparedStatement updateValidateNonOverWrite;
    private PreparedStatement updateSetConfidence;
    private ArrayList<UUID> m_resultset_ids = new ArrayList<UUID>();
    private String m_connectionString = null;
    private String m_dbuser = null;
    private String m_dbpass = null;
    private String m_connectionStringAdmin = null;
    private String m_dbuserAdmin = null;
    private String m_dbpassAdmin = null;
    private HashMap<String, HashMap<Long, Sequence>> m_proteinSequences = new HashMap<>();
    private Sequence m_noSequence = new Sequence(new AminoAcid[0]);
    //private DBRunConfig m_conf;
    //private HashMap<Integer,DBRunConfig> m_configs;
    Xi2Xi1Config m_config = new Xi2Xi1Config();
    HashMap<String, Xi2Xi1Config> m_configs = new HashMap<>();

    private boolean m_writePrePostAA = false;
    public static DecimalFormat sixDigits = new DecimalFormat("###0.000000");
    
    HashMap<UUID,Xi2ScoreList> subscoreDefs = new HashMap<>();
    
    Xi2ScoreList resultScores = new Xi2ScoreList();

    UUID[] resultSetIDSetting;
    ArrayList<String> m_searchids = new ArrayList<>();
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
    // reuse the last mzidentml owner infos
    boolean lastMzIDowner = false;
    

    public String getLinkWindow(PeptidePair pp, int proteingroup, int window) {
        return "NOT IMPLEMENTED";
    }

    public String getLinkWindow(ProteinGroupLink l, int proteingroup, int window) {
        return "NOT IMPLEMENTED";
    }

    /**
     * @return the m_db_connection
     */
    public Connection getDBConnection() {
        flagDBUsage();
        if (m_db_autoclosetimer == null) {
            m_db_autoclosetimer = new Timer("Timer - auto close db2", true);
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
                        if (timeUnused > 3600000) {
                            if (!keepConnection.get()) {
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

    public void setDBConnectionAdmin(String connection, String user, String password) throws SQLException {
        m_connectionStringAdmin = connection;
        m_dbuserAdmin = user;
        m_dbpassAdmin = password;
        this.m_db_connectionAdmin = DriverManager.getConnection(connection, user, password);
        setupPreparedStatements();
    }

    /**
     * @param m_db_connection the m_db_connection to set
     */
    public void setDBConnection(Connection connection) throws SQLException {
        this.m_db_connection = connection;
        setupPreparedStatements();
    }

    public void setDBConnectionAdmin(Connection connection) throws SQLException {
        this.m_db_connectionAdmin = connection;
        setupPreparedStatements();
    }


    public DatabaseProvider getDatabaseProvider()  {
        return databaseProvider;
    }

    public void setDatabaseProvider(DatabaseProvider getSearch) throws SQLException {
        databaseProvider = getSearch;
        setDBConnection(databaseProvider.getConnection());
        if (getSearch instanceof GetSearch)
            setDBConnectionAdmin(((GetSearch) getSearch).getConnectionAdmin());
    }

    public synchronized void closeConnection()  {
        try {
            if (m_db_connection!=null) {
                m_db_connection.close();
                m_db_connection=null;
            }
        } catch (SQLException ex) {
            Logger.getLogger(DB2inFDR.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void flagDBUsage() {
        synchronized (m_db_last_used) {
            m_db_last_used.value = Calendar.getInstance().getTimeInMillis();
        }
    }

    public synchronized boolean ensureConnection() throws SQLException {
        flagDBUsage();
        boolean isconnected = false;
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
                    Logger.getLogger(DB2inFDR.class.getName()).log(Level.SEVERE, "Currently I don't have a way to reopen the connection");
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

    public synchronized boolean ensureConnectionAdmin() throws SQLException {
        flagDBUsage();
        boolean isconnected = false;
        try {
            Statement st = m_db_connectionAdmin.createStatement();
            ResultSet rs = st.executeQuery("select 1");
            rs.close();
            st.close();
            isconnected = true;
        } catch (Exception sex) {
            //Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, "Database connection is closed/ non-functioning. Will try to reopen", sex);
            try {
                m_db_connectionAdmin.close();
            } catch (Exception e) {
            }

        }
//        }
//        }

        if (!isconnected) {
            Logger.getLogger(this.getClass().getName()).log(Level.WARNING, "Lost connection to database. Will try to reopen");

            try {
                if (m_connectionStringAdmin != null && m_dbuserAdmin != null && m_dbpassAdmin != null) {
                    this.m_db_connectionAdmin = DriverManager.getConnection(m_connectionStringAdmin, m_dbuserAdmin, m_dbpassAdmin);
                } else if (databaseProvider != null && databaseProvider instanceof GetSearch) {
                    setDBConnectionAdmin(((GetSearch) databaseProvider).getConnectionAdmin());
                } else {
                    Logger.getLogger(DB2inFDR.class.getName()).log(Level.SEVERE, "Currently I don't have a way to reopen the connection");
                    return false;
                }

                isconnected = true;
                this.m_db_connectionAdmin.setAutoCommit(true);

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
        updateSetConfidence = this.m_db_connection.prepareStatement("UPDATE spectrum_match rescored = ? SET WHERE id = ?");
    }

    @Override
    public String getSource() {
        if (m_connectionString != null && !m_connectionString.isEmpty()) {
            return "\n\tDB:" + m_connectionString + "\n\tIDs:" + RArrayUtils.toString(m_resultset_ids, ";") + (filterSetting == null || filterSetting.isEmpty() ? "" : "\n\tFilter: " + filterSetting);
        } else {
            return "\n\tIDs:" + RArrayUtils.toString(m_resultset_ids, ",") + (filterSetting == null || filterSetting.isEmpty() ? "" : "\n\tFilter: " + filterSetting);
        }
    }

    public ArrayList<UUID> getResultSetIDs() {
        return m_resultset_ids;
    }

    @Override
    public ArrayList<String> getSearchIDs() {
        return this.m_searchids;
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

    public DB2inFDR() {
        resultScores.add(new Xi2Score(0, "CSM FDR", false, false));
        resultScores.add(new Xi2Score(1, "CSM PEP", false, false));
        resultScores.add(new Xi2Score(2, "CSM Score", false, true));
        resultScores.add(new Xi2Score(3, "PeptidePair FDR", false, false));
        resultScores.add(new Xi2Score(4, "PeptidePair PEP", false, false));
        resultScores.add(new Xi2Score(5, "PeptdePair Score", false, true));
        resultScores.add(new Xi2Score(6, "ResiduePair FDR", false, false));
        resultScores.add(new Xi2Score(7, "ResiduePair PEP", false, false));
        resultScores.add(new Xi2Score(8, "ResiduePair Score", true, true));
        resultScores.add(new Xi2Score(9, "ProteinPair FDR", false, false));
        resultScores.add(new Xi2Score(10, "ProteinPair PEP", false, false));
        resultScores.add(new Xi2Score(11, "ProteinPair Score", false, true));
        resultScores.add(new Xi2Score(12, "Protein1 FDR", false, false));
        resultScores.add(new Xi2Score(13, "Protein1 PEP", false, false));
        resultScores.add(new Xi2Score(14, "Protein1 Score", false, true));
        resultScores.add(new Xi2Score(15, "Protein2 FDR", false, false));
        resultScores.add(new Xi2Score(16, "Protein2 PEP", false, false));
        resultScores.add(new Xi2Score(17, "Protein2 Score", false, true));
    }

    public DB2inFDR(Connection connection) {
        this.m_db_connection = connection;
    }

    public DB2inFDR(Connection connection, int[] peptideLengthGroups) {
        super(peptideLengthGroups);
        this.m_db_connection = connection;
    }

    
    public void read_score_names(UUID resultset_id) throws SQLException{
        Xi2ScoreList scores = new Xi2ScoreList();
        Xi2Score primary_score = null;
        try {
            ensureConnection();
        } catch (SQLException sex) {
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE,"No connection to DB");
            return;
        }
        String querry = 
                "SELECT "
                + "  score_id, "
                + "  name, "
                + "  primary_score, "
                + "  higher_is_better "
                + "FROM ScoreName WHERE resultset_id = '" +
                resultset_id + "' ORDER BY score_id";
        Statement stm = this.getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
        ResultSet rs = stm.executeQuery(querry);
        while (rs.next()) {
            Xi2Score s = new Xi2Score(
                    rs.getInt("score_id"), 
                    rs.getString("name"), 
                    rs.getBoolean("primary_score"), 
                    rs.getBoolean("higher_is_better"));
            if (s.primary_score)
                primary_score = s;
            scores.add(s);
        }
//        if (primary_score != null)
//            scores.add(primary_score);
//        else 
//            scores.add(null);
        this.subscoreDefs.put(resultset_id, scores);
    }
    
    
    public HashMap<UUID,Xi2Xi1Config> readconfig(UUID resultset_id) throws SQLException{
        HashMap<UUID,Xi2Xi1Config> configs = new HashMap<UUID,Xi2Xi1Config>();
        ensureConnection();
        Statement st = this.getDBConnection().createStatement(
                ResultSet.TYPE_FORWARD_ONLY, 
                ResultSet.CONCUR_READ_ONLY);
        ResultSet rs = st.executeQuery("SELECT s.id as id, config "
                + "FROM "
                + " Search s"
                + " INNER JOIN "
                + " ResultSearch rss ON rss.search_id = s.id "
                + "WHERE rss.resultset_id = '" + resultset_id + "'");
        
        while (rs.next()) {
            String config = Unescape.unescape(rs.getString("config"));
            UUID id = (UUID) rs.getObject("id");
            String sid = id.toString();
            configs.put(id, new Xi2Xi1Config(config.substring(1, config.length()-1)));
            if (!this.m_searchids.contains(sid)){
                m_searchids.add(sid);
            }
        }
        return configs;
    }
    
    
    
    public void readDB() throws SQLException {
        readDB(resultSetIDSetting, filterSetting, true);
    }

    public void readDB(UUID resultsetId, String filter, boolean topOnly) throws SQLException {
        readDB(new UUID[]{resultsetId}, filter, topOnly);
    }

    public void readDB(UUID[] resultsetIds, String filter, boolean topOnly) throws SQLException {
        this.readDBSteps(resultsetIds, filter, topOnly, new HashSet<UUID>(), 0, new HashMap<UUID, Double>());
    }

    public void readDBSteps(UUID[] resultset_ids, String filter, boolean topOnly, HashSet<UUID> prev_skip, int tries, HashMap<UUID,Double>  lastScore) throws SQLException {
        
        HashSet<UUID> skip= new HashSet<>();
            
        if (!ensureConnection()) {
            return;
        }
        for (UUID resultset_id : resultset_ids) {
            HashMap<UUID, Xi2Xi1Config> confs = readconfig(resultset_id);
            for (Map.Entry<UUID, Xi2Xi1Config> confe : confs.entrySet()) {
                m_configs.put(confe.getKey().toString(), confe.getValue());
                // TODO need to update ones we have versions in the DB
                m_xi_versions.put(confe.getKey().toString(), new Version("2.0.beta"));
            }
                
            
            for (Xi2Xi1Config conf :  confs.values()) {
                this.m_config.add(conf);
            }
        }
        HashMap<String, HashMap<Long, Protein>> search_proteins = new HashMap<>();
        for (String searchid : m_configs.keySet()) {
            search_proteins.put(searchid, readProteinInfos(searchid));
        }
        
        if (!m_resultset_ids.isEmpty()) {
            if ((filter == null && !filterSetting.isEmpty()) || (!filter.contentEquals(filterSetting))) {
                filterSetting = filterSetting + "\n ResultSet IDS:" + RArrayUtils.toString(resultset_ids, ",") +":" + filter;
            }
        } else if (filter != null && !filter.isEmpty()) {
            filterSetting = filter;
        }
        boolean shownMinPepWarning = false;
        

        for (int currentresultset = 0 ; currentresultset < resultset_ids.length; currentresultset++){
            ensureConnection();
            UUID resultset_id = resultset_ids[currentresultset];
            Statement stmSearch_ids = getDBConnection().createStatement();
            
            if (m_resultset_ids.contains(resultset_id)) {
                continue;
            }
            
            
            
//            while (rsSearch_ids.next()){
//                UUID search_id = rsSearch_ids.getObject(0, UUID.class);
                String searchfilter = filter;
                

                String sPepCoverage1 = "unique_peak_primary_coverage_p1";
                String sPepCoverage2 = "unique_peak_primary_coverage_p2";
                String sPepUniqueMatchedCons1 = "unique_peak_conservative_fragsites_p1";
                String sPepUniqueMatchedCons2 = "unique_peak_conservative_fragsites_p2";
                String sPepStubs = "cc_pep_frag_pp";
                String sDelta = "delta";
                String sPepDoublets = "fragment CCPepDoubletFound";
                String sCleavCLPep1Fragmatched = "cc_pep_frag_p1";
                String sCleavCLPep2Fragmatched = "cc_pep_frag_p2";
                String sAutoValidated = "auto_validation";
                Integer cPepCoverage1 = null;
                Integer cPepCoverage2 = null;
                Integer cDelta = null;
                Integer cPepStubs = null;
                Integer cPepDoublets = null;
                Integer cPrimaryScore = -1;
                Integer cCleavCLPep1Fragmatched = null;
                Integer cCleavCLPep2Fragmatched = null;
                Integer cAutoValidated = null;
                Integer cPepUniqueMatchedCons1 = null;
                Integer cPepUniqueMatchedCons2 = null;

                ArrayList<Integer> peaks = new ArrayList<>();
                HashSet<Integer> scoresForwarded = new HashSet<>();

                read_score_names(resultset_id);
                
                if (subScoresToForward != null ) {
                    for (Xi2Score s : subscoreDefs.get(resultset_id)) {
                            if (subScoresToForward.matcher(s.name).matches()) {
                                scoresForwarded.add(s.id);
                            }

                    }
                }
                if (subscoreDefs.get(resultset_id).get(sPepUniqueMatchedCons1)!=null)
                    cPepUniqueMatchedCons1 = subscoreDefs.get(resultset_id).get(sPepUniqueMatchedCons1).id;
                if (subscoreDefs.get(resultset_id).get(sPepUniqueMatchedCons2)!=null)
                    cPepUniqueMatchedCons2 = subscoreDefs.get(resultset_id).get(sPepUniqueMatchedCons2).id;

                if (subscoreDefs.get(resultset_id).get(sPepCoverage1)!=null) {
                    cPepCoverage1 = subscoreDefs.get(resultset_id).get(sPepCoverage1).id;
                }
                if (subscoreDefs.get(resultset_id).get(sPepCoverage2)!=null) {
                    cPepCoverage2 = subscoreDefs.get(resultset_id).get(sPepCoverage2).id;
                }
                if (subscoreDefs.get(resultset_id).get(sCleavCLPep1Fragmatched) != null)
                    cCleavCLPep1Fragmatched = subscoreDefs.get(resultset_id).get(sCleavCLPep1Fragmatched).id;
                if (subscoreDefs.get(resultset_id).get(sCleavCLPep2Fragmatched) != null)
                    cCleavCLPep2Fragmatched = subscoreDefs.get(resultset_id).get(sCleavCLPep2Fragmatched).id;
                if (subscoreDefs.get(resultset_id).get(sPepStubs) != null)
                    cPepStubs = subscoreDefs.get(resultset_id).get(sPepStubs).id;
                if (subscoreDefs.get(resultset_id).get(sPepDoublets) != null)
                    cPepDoublets = subscoreDefs.get(resultset_id).get(sPepDoublets).id;
                if (subscoreDefs.get(resultset_id).get(sDelta) != null)
                    cDelta = subscoreDefs.get(resultset_id).get(sDelta).id;
                if (subscoreDefs.get(resultset_id).get(sAutoValidated) != null)
                    cAutoValidated = subscoreDefs.get(resultset_id).get(sAutoValidated).id;
                cPrimaryScore = subscoreDefs.get(resultset_id).primaryScore.id;
                

                if (cPepCoverage2 == null && !shownMinPepWarning) {
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
                    for (int i = 0;i<subscoreDefs.get(resultset_id).size();i++) {
                        if (subscoreDefs.get(resultset_id).get(i).name.toLowerCase().equals(scorename.trim().toLowerCase())) {
                            scoreid = i;
                        }
                    }
                    if (scoreid == -1) {
                        throw new SQLException("subscore [" + scorename + "] not found");
                    }

                    scorefilters_used.add(subscorematcher.group(1));

                    searchfilter = searchfilter.substring(0, subscorematcher.start()) 
                            +"subscores[" + (scoreid) + " +  array_lower(subscores, 1)]" + searchfilter.substring(subscorematcher.end());
                    subscorematcher = subscorepattern.matcher(searchfilter);
                    if (!scoresForwarded.contains(scoreid)) {
                        scoresForwarded.add(scoreid);
                    }
                }   
                // make sure the gui shows previously used subscores higher up in the list
                if (!scorefilters_used.isEmpty()) {
                    String prev = LocalProperties.getProperty(DB2inFDR.subscoreroperty, "");
                    if (prev.length() > 0) {
                        for (String s : prev.split(";")) {
                            if (!scorefilters_used.contains(s))
                                scorefilters_used.add(s);
                        }

                    }
                    LocalProperties.setProperty(DB2inFDR.subscoreroperty, RArrayUtils.toString(scorefilters_used,";"));
                }

                Double searchLastScore = lastScore.get(resultset_id);
                String matchQuerry;
                matchQuerry
                        = "SELECT * FROM (SELECT "
                        + "m.id AS psmID, \n"
                        + "mp1.sequence AS pepSeq1, \n"
                        + "mp1.base_sequence AS baseseq1, \n"
                        + "mp1.modification_ids AS pep_mods1, \n"
                        + "mp1.modification_position AS pep_mods_pos1, \n"
                        + "mp2.base_sequence AS baseseq2, \n"
                        + "mp2.sequence AS pepSeq2, \n"
                        + "mp2.modification_ids AS pep_mods2, \n"
                        + "mp2.modification_position AS pep_mods_pos2, \n"
                        + "mp1.length as peplen1, \n"
                        + "mp2.length as peplen2, \n"
//                        + "mp1.modification_ids as modIDs1, \n"
//                        + "mp2.modification_ids as modIDs2, \n"
//                        + "mp1.modification_position as modPos1, \n"
//                        + "mp2.modification_position as modPos2, \n"
//                        + "m.link_score_site1, \n"
//                        + "m.link_score_site2, \n"
                        + "CASE WHEN rm.link_score IS NULL THEN FALSE ELSE TRUE END AS sites_in_resultmatch, \n"
                        + "CASE WHEN rm.link_score IS NULL THEN m.link_score ELSE rm.link_score END as link_score, \n"
                        + "CASE WHEN rm.link_score_site1 IS NULL THEN m.link_score_site1 ELSE rm.link_score_site1 END as link_score_site1, \n"
                        + "CASE WHEN rm.link_score_site2 IS NULL THEN m.link_score_site2 ELSE rm.link_score_site2 END as link_score_site2, \n"
                        + "CASE WHEN rm.site1 IS NULL THEN m.site1 ELSE rm.site1 END as site1, \n"
                        + "CASE WHEN rm.site1 IS NULL THEN m.site2 ELSE rm.site2 END  as site2, \n"
                        + "mp1.is_decoy AS isDecoy1, \n"
                        + "mp2.is_decoy AS isDecoy2, \n"
                        + "m.assumed_prec_charge AS calc_charge, \n"
                        + "m.assumed_prec_mz AS calc_mz, \n"
                        + "m.score as PSMscore, \n"
                        + "array_agg(pp1.start) AS pepPosition1,  \n"
                        + "array_agg(pp2.start) AS pepPosition2, \n"
                        + "CASE WHEN mp2.sequence IS NULL THEN 1 ELSE (4.0/5.0+(mp1.length/(mp1.length+mp2.length)))/2 END AS score_ratio \n"
                        + " , s.precursor_charge as exp_charge \n"
                        + " , array_agg(pp1.protein_id) as protein1id\n"
                        + " , array_agg(pp2.protein_id) as protein2id\n"
                        + " , pep1_id as peptide1id\n"
                        + " , pep2_id as peptide2id\n"
                        + " , r.name as run_name, s.scan_number \n"
                        + " , s.precursor_mz AS exp_mz\n"
                        + " , m.calc_mass\n"
                        + " , m.assumed_prec_charge as match_charge\n"
                        + " , mp1.mass as pep1mass\n"
                        + " , mp2.mass as pep2mass\n"
                        + " , m.search_id\n"
                        + " , ms.spectrum_id\n"
                        + " , s.precursor_charge\n"
                        + " , crosslinker_id\n"
                        + " , s.scan_index\n"
                        + " , plf.name as peaklistfile\n"
                        + " , scores AS subscores\n"
                        + (cPepUniqueMatchedCons1 == null ? ", 0" : " , scores[" + cPepUniqueMatchedCons1 + "  +  array_lower(scores, 1)]") + " as scoreP1Coverage\n"
                        + (cPepUniqueMatchedCons2 == null ? ", 0" : " , scores[" + cPepUniqueMatchedCons2 + "  +  array_lower(scores, 1)]") + " as scoreP2Coverage\n"
                        + (cDelta == null ? ", 0" : " , scores[" + cDelta + "  +  array_lower(scores, 1)]" ) + " as deltaScore\n"
                        + " , scores[" + cPrimaryScore + " +  array_lower(scores, 1)] as score\n"
                        + (cCleavCLPep1Fragmatched == null ? ", 0" : " , scores[" + cCleavCLPep1Fragmatched +"  +  array_lower(scores, 1)]::int") + " as cleavCLPep1Fragmatched \n"
                        + (cCleavCLPep2Fragmatched == null ? ", 0" : " , scores[" + cCleavCLPep2Fragmatched +"  +  array_lower(scores, 1)]::int") + " as cleavCLPep2Fragmatched \n"
                        + (cAutoValidated == null ? ", false" : " , scores[" + cAutoValidated +"  +  array_lower(scores, 1)] <> 0") + " as autovalidated \n"
                        + " , s.precursor_intensity\n"
                        + " , s.retention_time as retentiontime\n"
                        + " , rm.match_group_id"
                        + " , array_lower(scores, 1) as lowest_score_id"
                        + " \n"
                        + "FROM \n"
                        + " resultmatch rm INNER JOIN match m ON rm.match_id = m.id and rm.top_ranking\n" +
                        "  INNER JOIN \n" +
                        " matchedspectrum ms on m.id = ms.match_id INNER JOIN \n" +
                        " spectrum s on ms.spectrum_id = s.id INNER JOIN \n" +
                        " peaklist plf on s.peaklist_id = plf.id  INNER JOIN \n" +
                        " source ss on s.source_id = ss.id INNER JOIN \n" +
                        " run r on s.run_id = r.id \n" +
                        "    INNER JOIN \n" +
                        "   modifiedpeptide mp1 on m.search_id = mp1.search_id AND m.pep1_id = mp1.id \n" +
                        " INNER JOIN \n" +
                        "   peptideposition pp1 on m.search_id = pp1.search_id AND m.pep1_id = pp1.mod_pep_id \n" +
                        " LEFT OUTER JOIN \n" +
                        "   modifiedpeptide mp2 on m.search_id = mp2.search_id AND m.pep2_id = mp2.id \n" +
                        " LEFT OUTER JOIN \n" +
                        "   peptideposition pp2 on m.search_id = pp2.search_id AND m.pep2_id = pp2.mod_pep_id  \n"
                        + " WHERE  \n"
                        + " (s.precursor_charge = m.assumed_prec_charge OR m.assumed_prec_charge < 6) \n"
                        + " AND rm.resultset_id = '" + resultset_id + "'"
                        + " GROUP BY \n"
                        + "m.id, \n"
                        + "mp1.sequence, \n"
                        + "mp1.base_sequence, \n"
                        + "mp1.modification_ids , \n"
                        + "mp1.modification_position , \n"
                        + "mp2.base_sequence , \n"
                        + "mp2.sequence , \n"
                        + "mp2.modification_ids , \n"
                        + "mp2.modification_position , \n"
                        + "mp1.length , \n"
                        + "mp2.length , \n"
                        + "m.link_score_site1, \n"
                        + "m.link_score_site2, \n"
                        + "m.link_score, \n"
                        + "m.site1 + 1 , \n"
                        + "m.site2 + 1 , \n"
                        + "mp1.is_decoy, \n"
                        + "mp2.is_decoy, \n"
                        + "m.assumed_prec_charge , \n"
                        + "m.assumed_prec_mz , \n"
                        + "m.score, \n"
                        + "CASE WHEN mp2.sequence IS NULL THEN 1 ELSE (4.0/5.0+(mp1.length/(mp1.length+mp2.length)))/2 END \n"
                        + " , s.precursor_charge  \n"
                        + " , pep1_id \n"
                        + " , pep2_id \n"
                        + " , r.name, s.scan_number \n"
                        + " , s.precursor_mz \n"
                        + " , m.calc_mass\n"
                        + " , m.assumed_prec_charge \n"
                        + " , mp1.mass \n"
                        + " , mp2.mass \n"
                        + " , m.search_id\n"
                        + " , ms.spectrum_id\n"
                        + " , s.precursor_charge\n"
                        + " , crosslinker_id\n"
                        + " , s.scan_index\n"
                        + " , plf.name \n"
                        + " , scores \n"
                        + (cPepCoverage1 == null ? "" : " , scores[" + cPepCoverage1 + " +  array_lower(scores, 1)]") + " \n"
                        + (cPepCoverage2 == null ? "" : " , scores[" + cPepCoverage2 + "  +  array_lower(scores, 1)]") + " \n"
                        + (cDelta == null ? "" : " , scores[" + cDelta + "  +  array_lower(scores, 1)]" ) + " \n"
                        + " , scores[" + cPrimaryScore + " +  array_lower(scores, 1)] \n"
                        + (cCleavCLPep1Fragmatched == null ? "" : " , scores[" + cCleavCLPep1Fragmatched +" +  array_lower(scores, 1)]::int") + "  \n"
                        + (cCleavCLPep2Fragmatched == null ? "" : " , scores[" + cCleavCLPep2Fragmatched +" +  array_lower(scores, 1)]::int") + " \n"
                        + (cAutoValidated == null ? "" : " , scores[" + cAutoValidated +" +  array_lower(scores, 1)] <> 0") + " \n"
                        + " , s.precursor_intensity\n"
                        + " , s.retention_time \n"
                        + " , rm.link_score\n"
                        + " , rm.link_score_site1\n"
                        + " , rm.link_score_site2\n"
                        + " , rm.site1\n"
                        + " , rm.site2"
                        + " , rm.match_group_id\n"
                        //+ " ORDER BY scores[" + cPrimaryScore + " +  array_lower(scores, 1)] DESC
                        + ")  i \n"
                        + (searchfilter == null || searchfilter.isEmpty() ? "" : " WHERE ( " + searchfilter + " )");

                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "read from db");
                Logger.getLogger(this.getClass().getName()).log(Level.INFO, matchQuerry);
                if (resultset_ids.length >1)
                    Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Resultset {0} of {1}", new Object[]{currentresultset+1, resultset_ids.length});
                //getDBConnection().setAutoCommit(false);
                Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
                stm.setFetchSize(100);
                flagDBUsage();
                keepConnection.set(true);
                ResultSet rs=null;
                try {
                    rs = stm.executeQuery(matchQuerry);
                    flagDBUsage();
                } finally {
                    keepConnection.set(false);
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
                int MatchGrpIDColumn = rs.findColumn("match_group_id");
                //nt LinkSiteDeltaColumn = rs.findColumn("LinkSiteDelta");
                int exp_mzColumn = rs.findColumn("exp_mz");
                int calc_massColumn = rs.findColumn("calc_mass");
                int pep1massColumn = rs.findColumn("pep1mass");
                int pep2massColumn = rs.findColumn("pep2mass");
                int search_idColumn = rs.findColumn("search_id");
                int spectrum_idColumn = rs.findColumn("spectrum_id");
                int xlColumn = rs.findColumn("crosslinker_id");
                int scanIndexColumn = rs.findColumn("scan_index");
                int peaklistfileColumn = rs.findColumn("peaklistfile");
                int subscoresColumn = rs.findColumn("subscores");
                int precIntensityColumn = rs.findColumn("precursor_intensity");
                int retentiontimeColumn = rs.findColumn("retentiontime");
                int cleavclpep1fragmatchedColumn = rs.findColumn("cleavCLPep1Fragmatched");
                int cleavclpep2fragmatchedColumn = rs.findColumn("cleavCLPep2Fragmatched");
                int linkSiteScoreColumn = rs.findColumn("link_score");
                int linkSiteScore1Column = rs.findColumn("link_score_site1");
                int linkSiteScore2Column = rs.findColumn("link_score_site2");
                int sitesInResultmatchColumn = rs.findColumn("sites_in_resultmatch");
                int avColumn = rs.findColumn("autovalidated");
//                int modIDs1Columns = rs.findColumn("modIDs1");
//                int modIDs2Columns = rs.findColumn("modIDs2");
                Pattern modDetect = Pattern.compile(".*[^A-Z].*");
                HashSet<Double> tmmodcount = new HashSet<>();
                HashMap<Double,Double> xlmodmasses = new HashMap<Double,Double>(1);

                int total = 0;
                try {
                    while (rs.next()) {
                        total++;

                        UUID ipsmID = (UUID)rs.getObject(psmIDColumn);
                        if (prev_skip.contains(ipsmID)) {
                            continue;
                        }
                        
                        UUID search_id = (UUID)rs.getObject("search_id");
                        Xi2Xi1Config conf = m_configs.get(search_id.toString());
                        String psmID = ipsmID.toString();
                        String pepSeq1 = rs.getString(pepSeq1Column);
                        String pepSeq2 = rs.getString(pepSeq2Column);
                        int peplen1 = rs.getInt(peplen1Column);
                        int peplen2 = rs.getInt(peplen2Column);
                        int site1 = rs.getInt(site1Column);
                        int site2 = pepSeq2 == null ? -1 : rs.getInt(site2Column);
                        ResultMatchData rmd =  null;
                        if (pepSeq2 != null && site1 >= 0) {
                            java.sql.Array link_score = rs.getArray(linkSiteScoreColumn);
                            java.sql.Array link_score_site1 = rs.getArray(linkSiteScore1Column);
                            java.sql.Array link_score_site2 = rs.getArray(linkSiteScore2Column);
                            rmd = new ResultMatchData(link_score, link_score_site1, link_score_site2, true);
                        }
                        boolean isDecoy1 = rs.getBoolean(isDecoy1Column);
                        boolean isDecoy2 = rs.getBoolean(isDecoy2Column);
                        int charge = rs.getInt(calc_chargeColumn);
                        double score = rs.getDouble(scoreColumn);
                        Integer[] pepPosition1 = (Integer[])rs.getArray(pepPosition1Column).getArray();
                        Integer[] pepPosition2 = pepSeq2 == null ? new Integer[]{-1} : (Integer[])rs.getArray(pepPosition2Column).getArray();
                        double scoreRatio = rs.getDouble(score_ratioColumn);
                        int spectrum_charge = rs.getInt(exp_chargeColumn);
                        Integer[] protein1ID = (Integer[])rs.getArray(protein1idColumn).getArray();
                        Integer[] protein2ID = pepSeq2 == null ? new Integer[]{0} : (Integer[])rs.getArray(protein2idColumn).getArray();
                        long pep1ID = rs.getLong(peptide1idColumn);
                        long pep2ID = rs.getLong(peptide2idColumn);
                        String run = rs.getString(run_nameColumn);
                        String scan = rs.getString(scan_numberColumn);
                        boolean autovalidated = rs.getBoolean(avColumn);
                        UUID match_group_id = (UUID)rs.getObject(MatchGrpIDColumn);
                        //int rank = rs.getInt(rankColumn);



                        double p1c = rs.getDouble(scoreP1CoverageColumn);
                        double p2c = rs.getDouble(scoreP2CoverageColumn);
                        double pminc = p1c;

                        if (pepSeq2 != null && !pepSeq2.isEmpty()) {
                            pminc = Math.min(p1c,p2c);
                        }

                        double pmz = rs.getDouble(exp_mzColumn);
                        double f = 1;
                        if (pepSeq2 != null && !pepSeq2.isEmpty() && p1c + p2c > 0) {
                            scoreRatio = (p1c) / (p1c + p2c + 1);
                        }

                        double peptide1score = score*scoreRatio;
                        double peptide2score = score*(1-scoreRatio);

                        int xl = rs.getInt(xlColumn);
                        
                        double calc_mass = rs.getDouble(calc_massColumn);
                        double pep1mass = rs.getDouble(pep1massColumn);
                        double pep2mass = rs.getDouble(pep2massColumn);
                        UUID scan_id = (UUID)rs.getObject(spectrum_idColumn);

                       // int mgxrank = rs.getInt(mgxrankColumn);
                        boolean cleavclpep1fragmatched = rs.getBoolean(cleavclpep1fragmatchedColumn);
                        boolean cleavclpep2fragmatched = rs.getBoolean(cleavclpep2fragmatchedColumn);
                        Float[] subscores = (Float[]) rs.getArray(subscoresColumn).getArray();

                        PSM psm = null;
                        if (pepSeq2 != null && !pepSeq2.isEmpty()){
                            for (int p = 0; p< protein1ID.length; p++) {
                                int p1 = pepPosition1[p];
                                long p1id = protein1ID[p];
                                HashMap<Long,Protein> sprots = search_proteins.get(search_id.toString());
                                Protein prot1 = sprots.get(p1id);

                                int p2 = pepPosition2[p];
                                long p2id = protein2ID[p];
                                Protein prot2 = sprots.get(p2id);

                                psm = setUpDBPSM(psmID, run, scan, pep1ID, pep2ID, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, p1id, prot1.getAccession(), prot1.getDescription(), p2id, prot2.getAccession(), prot2.getDescription(), p1, p2, prot1.getSequence(), prot2.getSequence(), peptide1score, peptide2score, spectrum_charge, conf.xi2crosslinker.get(xl).name, pmz, calc_mass, pep1mass, pep2mass, search_id.toString(), scan_id);
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

                                HashMap<Long,Protein> sprots = search_proteins.get(search_id.toString());
                                Protein prot1 = sprots.get(p1id);
                                
                                psm = setUpDBPSM(psmID, run, scan, pep1ID, pep2ID, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, p1id, prot1.getAccession(), prot1.getDescription(), p2id, a2, d2, p1, p2, prot1.getSequence(), s2, peptide1score, peptide2score, spectrum_charge, conf.xi2crosslinker.get(xl).name, pmz, calc_mass, pep1mass, pep2mass, search_id.toString(), scan_id);
                            }
                        }

                        //psm.setRank(rank);

                        Float[] scorevalues = null;
                        if (subscores != null) {
                            scorevalues = subscores;
                            if (cPepCoverage1 != null && cPepCoverage1 >= 0) {
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

                        if (cPepStubs != null) {
                            int s = scorevalues[cPepStubs].intValue();
                            psm.peptidesWithStubs =  s;
                            psm.addOtherInfo("PeptidesWithStubs", s);
                            if (s >0)  {
                                stubsFound(true);
                            } 

                        }
                        if (cPepDoublets != null) {

                            psm.addOtherInfo("PeptidesWithDoublets", scorevalues[cPepDoublets].intValue());
                            psm.peptidesWithDoublets =  scorevalues[cPepDoublets].intValue();

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
                        //psm.addOtherInfo("linkSiteDelta",rs.getDouble(LinkSiteDeltaColumn));
                        //psm.addOtherInfo("mgxrank",rs.getDouble(mgxrankColumn));
                        psm.addOtherInfo("PrecursorIntensity",rs.getDouble(precIntensityColumn));
                        psm.addOtherInfo("RetentionTime",rs.getDouble(retentiontimeColumn));
                        //Array lss1 = rs.getArray(linkSiteScore1Column);
                        //Array lss2 = rs.getArray(linkSiteScore2Column);
                        //if (lss1 != null) {
//                            psm.addOtherInfo("LinkSiteScore1",(Float[])lss1.getArray());
//                        } else {
//                            psm.addOtherInfo("LinkSiteScore1", new Float[0]);
//                        }
//                        if (lss2 != null) {
//                            psm.addOtherInfo("LinkSiteScore2",(Float[])lss2.getArray());
//                        } else {
//                            psm.addOtherInfo("LinkSiteScore2", new Float[0]);
//                        }

                        for (Integer p : peaks) {
                            psm.addOtherInfo(subscoreDefs.get(resultset_id).get(p).name, scorevalues[p]);
                        }

                        for (Integer p : scoresForwarded) {
                            psm.addOtherInfo(subscoreDefs.get(resultset_id).get(p).name, scorevalues[p]);
                        }
                        
                        Double xlModmassPre = conf.xi2crosslinker.get(xl).mass;
                        Double xlModmass = xlmodmasses.get(xlModmassPre);
                        if (xlModmass == null) {
                            xlmodmasses.put(xlModmassPre,xlModmassPre);
                            xlModmass = xlModmassPre;
                        }
                        psm.setCrosslinkerModMass(xlModmass);
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
                        
                        if (rmd != null) {
                            psm.nonPublicFlags.put("ResultMatchData", rmd);
                        }
                        psm.nonPublicFlags.put("match_group_id", match_group_id);

    //                    if (isMultpleSearches) {
    //                        psm.addPositiveGrouping("" + searchId);
    //                    }
                        if (total % 100 == 0) {
                            flagDBUsage();
                        }
                        if (total % 10000 == 0) {
                            Logger.getLogger(this.getClass().getName()).log(Level.INFO, 
                                    "go through the results ({0}) Search {1} of {2}", 
                                    new Object[]{allPSMs.size(), currentresultset+1, resultset_ids.length});
                        }
                        lastScore.put(resultset_id, score);
                        skip.add(ipsmID);

                    }
                    rs.close();
                    stm.close();

                    m_resultset_ids.add(resultset_id);

                } catch (SQLException sex) {
                    if (tries < 5) {
                        int nexttry = tries;
                        if (total<10) {
                            nexttry++;
                            Logger.getLogger(this.getClass().getName()).log(Level.FINE, "failed to read from the database retrying", sex);
                        }
                        Logger.getLogger(this.getClass().getName()).log(Level.FINEST, "failed to read from the database retrying", sex);
                        skip.addAll(prev_skip);
                        readDBSteps(resultset_ids, filter, topOnly, skip, nexttry, lastScore);
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
            //}
        }

    }

    protected HashMap<Long, Protein> readProteinInfos(String search_id) throws SQLException {
        Statement stm;
        ResultSet rs;
        ensureConnection();
        HashMap<Long,Protein> id2Protein = new HashMap<>();
        for (Protein p : allProteins) {
            id2Protein.put(p.getId(), p);
        }
        
        HashMap<Long, Protein> ret = new HashMap<>();
        
        String proteinQuerry = "SELECT id, accession, name, description, sequence, is_decoy FROM protein WHERE search_id = '" + search_id + "'";
        Logger.getLogger(this.getClass().getName()).log(Level.FINE, "Geting protein information: \n" + proteinQuerry);
        stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
        stm.setFetchSize(100);
        rs = stm.executeQuery(proteinQuerry);
        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "go through the results");
        int idColumn = rs.findColumn("id");
        int accessionColumn = rs.findColumn("accession");
        int descriptionColumn = rs.findColumn("description");
        int proteinSequenceColumn = rs.findColumn("sequence");
        int proteinNameColumn = rs.findColumn("name");
        int decoyColumn = rs.findColumn("is_decoy");
        while (rs.next()) {
            long id = rs.getLong(idColumn);
            String sequence =  rs.getString(proteinSequenceColumn).replaceAll("[^A-Z]", "");
            String accession = rs.getString(accessionColumn);
            String name = rs.getString(proteinNameColumn);
            String description = rs.getString(descriptionColumn);
            boolean is_decoy = rs.getBoolean(decoyColumn);
            
            if (accession == null && description !=null) {
                accession = description;
            }
            
            if (description==null || description.isEmpty()) {
                if (name == null || name.isEmpty())
                    description = accession;
                else
                    description = name;
            }
            
            if (name == null) {
                name = description;
            }
            Protein p = new Protein(id, accession, name, description, is_decoy, false, false, false);
            ret.put(id, p);
        }
        rs.close();
        stm.close();
        return ret;
    }
    
    
    protected PSM setUpDBPSM(String psmID, String run, String scan, long pep1ID, long pep2ID, String pepSeq1, String pepSeq2, int peplen1, int peplen2, int site1, int site2, boolean isDecoy1, boolean isDecoy2, int charge, double score, long protein1ID, String accession1, String description1, long protein2ID, String accession2, String description2, int pepPosition1, int pepPosition2, String sequence1, String sequence2, double peptide1score, double peptide2score, int spectrum_charge, String xl, double pmz, double calc_mass, double pep1mass, double pep2mass, String search_id, UUID scan_id) {
        PSM psm = (PSM)addMatch(psmID, run, scan, pep1ID, pep2ID, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, protein1ID, accession1, description1, protein2ID, accession2, description2, pepPosition1, pepPosition2, sequence1, sequence2, peptide1score, peptide2score, spectrum_charge == -1?"Unknow Charge" : null, xl);
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
        psm.setSearchID(search_id);
        psm.setScanID(scan_id.toString());
        return psm;
//                    if (!psm.isLinear()) {
//                        psm.setCrosslinker(rs.getString(xlColumn));
//                    }
    }



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
                resultSetIDSetting = new UUID[sids.length];

                for (int i = 0; i < sids.length; i++) {
                    resultSetIDSetting[i] = UUID.fromString(sids[i]);
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
            Logger.getLogger(DB2inFDR.class.getName()).log(Level.SEVERE, "No db-connection defined");
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
                    Logger.getLogger(DB2inFDR.class.getName()).log(Level.SEVERE, "Error parsing db config", ex);
                }
                    
            }
            
            if (connectionException != null) {
                Logger.getLogger(DB2inFDR.class.getName()).log(Level.SEVERE, null, connectionException);
                printUsage();
                System.exit(0);
            }
        }

        if (resultSetIDSetting == null) {
            Logger.getLogger(DB2inFDR.class.getName()).log(Level.SEVERE, "No search ids defined");
            printUsage();
            System.exit(0);
        }

        filterSetting = RArrayUtils.toString(filters, " AND ");
        String[] ret = new String[unknown.size()];
        ret = unknown.toArray(ret);
        return ret;
    }

    public static void main(String[] argv) throws SQLException, FileNotFoundException {

        DB2inFDR ofdr = new DB2inFDR();

        FDRSettings settings = new FDRSettingsImpl();
        
        String[] unknowns = ofdr.parseArgs(argv, settings);

        if (unknowns.length > 0) {
            for (String a : unknowns) {
                Logger.getLogger(DB2inFDR.class.getName()).log(Level.SEVERE, "Unknown argument : " + a);
            }
            ofdr.printUsage();
            System.exit(1);
        }

        if (ofdr.getCsvOutDirSetting() != null) {
            if (ofdr.getCsvOutBaseSetting() == null) {
                ofdr.setCsvOutBaseSetting("FDR_" + RArrayUtils.toString(ofdr.resultSetIDSetting, "_"));
            }
        }

        if (ofdr.getCsvOutBaseSetting() != null) {
            if (ofdr.getCsvOutDirSetting() == null) {
                ofdr.setCsvOutDirSetting(".");
            }

            System.out.println("writing results to " + ofdr.getCsvOutDirSetting() + "/" + ofdr.getCsvOutBaseSetting() + "*");
            System.out.flush();
        }

        Logger.getLogger(DB2inFDR.class.getName()).log(Level.INFO, "Read datafrom db");
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
        Logger.getLogger(DB2inFDR.class.getName()).log(Level.INFO, "Calculate FDR");
        FDRResult result = ofdr.calculateWriteFDR(ofdr.getCsvOutDirSetting(), ofdr.getCsvOutBaseSetting(), ",", settings, cu);

//        Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, "Result-summary:" + ofdr.summaryString());
//        
//        if (ofdr.getCsvOutBaseSetting() != null) {
//            
//            ofdr.writeFiles(ofdr.getCsvOutDirSetting(), ofdr.getCsvOutDirSetting(), ",");
//        }
        //System.out.println(ofdr.summaryString(result));
        if (ofdr.command_line_auto_validate != null) {
            ofdr.writeResult("xiFDR_offline","",(UUID)null, result, true, true, null);
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
        HashSet<String> searchIDs = new HashSet<String>();
        HashSet<UUID> resultsetIDs = new HashSet<UUID>();
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
        HashSet<String> searches = new HashSet<String>();
        for (PeptidePair pp : pgp.getPeptidePairs()) {
            for (PSM psm : pp.getAllPSMs()) {
                for (PSM upsm : psm.getRepresented()) {
                    searches.add(upsm.getSearchID());
                }
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
        HashSet<String> searches = new HashSet<String>();
        for (PeptidePair pp : pg.getPeptidePairs()) {
            for (PSM psm : pp.getAllPSMs()) {
                for (PSM upsm : psm.getRepresented()) {
                    searches.add(upsm.getSearchID());
                }
            }
        }
        ArrayList<String> ret = super.getProteinGroupOutput(pg);
        ret.add(0, RArrayUtils.toString(searches, ";"));
        return ret;
    }

    protected String getPeptideSequence(org.rappsilber.fdr.entities.Peptide p) {
        return p.getSequence();
    }

    @Override
    protected ArrayList<String> getXlPepeptideOutputLine(PeptidePair pp) {
        ArrayList<String> sret = super.getXlPepeptideOutputLine(pp);
        ArrayList<String> ret = new ArrayList<String>(sret.size() + 5);
        HashSet<String> searches = new HashSet<String>();
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
        HashSet<String> searches = new HashSet<String>();
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

    private int getFastas(ArrayList<String> id, ArrayList<String> dbIDs, ArrayList<String> names) {
        int c = 0;
        ArrayList<String> i = new ArrayList<>();
        ArrayList<String> n = new ArrayList<String>();
        try {
            ensureConnectionAdmin();
            Statement s = m_db_connectionAdmin.createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
            String querry =                     "SELECT  DISTINCT \n" +
            " f.file, f.id \n" +
            " FROM \n" +
            "	core_search_sequences css INNER JOIN \n" +
            "	core_sequence cs on css.sequence_id = cs.id INNER JOIN \n" +
            "	core_uploadedfile f on cs.file_id = f.id INNER JOIN \n" +
            "	(SELECT id, task_uuid, resultset_uuid from core_search UNION "
                                        + "SELECT id, task_uuid, resultset_uuid from core_fdr UNION "
                                        + "SELECT id, task_uuid, resultset_uuid from core_postprocessing) s ON css.search_id = s.id  " +
            " WHERE resultset_uuid in  ('" + RArrayUtils.toString(id, "','") + "');";

            ResultSet rs = s.executeQuery(querry);
            while (rs.next()) {
                String name = rs.getString(1);
                name=name.substring(name.lastIndexOf("/")+1);
                i.add(rs.getInt(2)+"");
                n.add(name);
                c++;
            }
        } catch (SQLException ex) {
            Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, null, ex);
            return -1;
        }
        dbIDs.addAll(i);
        names.addAll(n);
        return c;
        //throw new UnsupportedOperationException("not done yet");
    }

    public int getFastas(String searchID, ArrayList<String> dbIDs, ArrayList<String> names) {
        ArrayList<String> id = new ArrayList<String>(1);
        id.add(searchID);
        return getFastas(id, dbIDs, names);
    }

    public int getFastas(ArrayList<String> dbIDs, ArrayList<String> names) {
//        ArrayList<UUID> ids = new ArrayList<UUID>(m_resultset_ids.size());
//        for (UUID i : m_resultset_ids) {
//            ids.add(i);
//        }
//        return getFastas(ids, dbIDs, names);
        throw new UnsupportedOperationException("not done yet");

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
    
    
    public void writeResult(String name, String notes, UUID FDRresultset, FDRResult result, boolean within, boolean  between, Long userID) throws SQLException {
        UUID id;
        boolean autocommitAdmin = false;
        boolean autocomit_result  =false;
        try {
            String summary = this.getSummary(result);
            if (notes.contains(summary_marker)) {
                notes = notes.replace(summary_marker, summary);
            }
            if (notes.contains(summary_marker_long)) {
                notes = notes.replace(summary_marker_long, summary);
            }
            // write new resultset
            if (FDRresultset == null) {
                id = UUID.randomUUID();
                autocommitAdmin = setupAdminResultSet(name, id, notes, userID);
            } else {
                id = FDRresultset;
            }

            int xifdr_rs_type = 0;
            ensureConnection();
            autocomit_result = this.m_db_connection.getAutoCommit();
            this.m_db_connection.setAutoCommit(false);

            // get the xiFDR resultset type
            Statement st = this.m_db_connection.createStatement();
            ResultSet rs = st.executeQuery("SELECT id FROM resultsettype WHERE name = 'xiFDR'");
            if (rs.next()) {
                xifdr_rs_type = rs.getInt(1);
                rs.close();
                st.close();
            } else {
                rs.close();
                rs = st.executeQuery("INSERT INTO resultsettype (name) VALUES ('xiFDR') RETURNING id");
                rs.next();
                xifdr_rs_type = rs.getInt(1);
                rs.close();
                st.close();
            }

            if (ensureConnection()) {
    //            id uuid NOT NULL,
    //            name character varying NOT NULL,
    //            note character varying,
    //            rstype_id integer NOT NULL,
    //            config text NOT NULL,
    //            main_score integer NOT NULL,            
                PreparedStatement pst  = this.m_db_connection.prepareStatement("INSERT INTO resultset (id, name, note, rstype_id, config, main_score)"
                        + " VALUES "
                        + "(?, ?,?,?,?,0)" );
                pst.setObject(1, id);
                pst.setString(2, name);
                pst.setString(3, notes);
                pst.setInt(4, xifdr_rs_type);
                pst.setString(5, summary);
                pst.execute();
                pst.close();
                
                // link to previous resultsets (history)
                pst = this.m_db_connection.prepareStatement("INSERT INTO history (resultset_id , parent_resultset_id) VALUES (? ,?)");
                for (UUID base: this.m_resultset_ids) {
                    pst.setObject(1, id);
                    pst.setObject(2, base);
                    pst.addBatch();
                }
                pst.executeBatch();
                pst.close();

                // link to search_ids
                pst = this.m_db_connection.prepareStatement("INSERT INTO resultsearch (resultset_id , search_id) VALUES (? ,?)");
                for (String search: this.m_configs.keySet()) {
                    pst.setObject(1, id);
                    pst.setObject(2, UUID.fromString(search));
                    pst.addBatch();
                }
                pst.executeBatch();
                pst.close();

                // get the stored addition infos from the first PSM
                Xi2ScoreList outScores = new Xi2ScoreList(resultScores);
                if (PSM.getOtherInfoNames().length > 0) {
                    for (String subscorename : PSM.getOtherInfoNames()) {
                        outScores.add(new Xi2Score(outScores.size(), subscorename, false, true));
                    }
                }


                // write out score names
                pst = this.m_db_connection.prepareStatement("INSERT INTO scorename (resultset_id , score_id, primary_score, higher_is_better, name) VALUES (? ,?, ?, ?, ?)");
                for (Xi2Score score : outScores) {
                    pst.setObject(1, id);
                    pst.setInt(2, score.id);
                    pst.setBoolean(3, score.primary_score);
                    pst.setBoolean(4, score.higher_is_better);
                    pst.setString(5, score.name);
                    pst.addBatch();
                }
                pst.executeBatch();
                pst.close();

                // write resultmatches
                pst = this.m_db_connection.prepareStatement(resultMatchColumns.querry);
                int batch_count=0;
                if (within && between) {
                    for (PSM psm: result.psmFDR.filteredResults()) {
                        if (!psm.isLinear()) {
                            addPSMtoBatch(pst, id, psm, outScores, result);
                            if (++batch_count % 1000 == 0) {
                                pst.executeBatch();
                            }
                        }
                    }
                    for (PSM psm: result.psmFDR.filteredResults()) {
                        if (psm.isLinear()) {
                            addPSMtoBatch(pst, id, psm, outScores, result);
                            if (++batch_count % 1000 == 0) {
                                pst.executeBatch();
                            }
                        }
                    }
                } else if (within) {
                    for (PSM psm: result.psmFDR.filteredResults()) {
                        if (!psm.isBetween()) {
                            addPSMtoBatch(pst, id, psm, outScores, result);
                            if (++batch_count % 1000 == 0) {
                                pst.executeBatch();
                            }
                        }
                    }
                } else if (between) {
                    for (PSM psm: result.psmFDR.filteredResults()) {
                        if (psm.isBetween()) {
                            addPSMtoBatch(pst, id, psm, outScores, result);
                            if (++batch_count % 1000 == 0) {
                                pst.executeBatch();
                            }
                        }
                    }
                }
                if (batch_count % 1000>0)
                    pst.executeBatch();
                pst.close();


            }
            if (FDRresultset == null) {
                this.m_db_connectionAdmin.commit();
                this.m_db_connectionAdmin.setAutoCommit(autocommitAdmin);
            }
            this.m_db_connection.commit();
            this.m_db_connection.setAutoCommit(autocomit_result);

        }catch (SQLException sex) {
            try {
                m_db_connectionAdmin.rollback();
                m_db_connection.rollback();
            } catch(Exception e) {
                
            }
            throw sex;
            
        }
        
        
    }

    protected boolean setupAdminResultSet(String name, UUID id, String notes, Long user_id) throws SQLException {
        boolean autocommitAdmin;
        ensureConnectionAdmin();
        // guess we also should write out to the admin DB
        autocommitAdmin = this.m_db_connectionAdmin.getAutoCommit();
        this.m_db_connectionAdmin.setAutoCommit(false);
        Statement userst = this.m_db_connectionAdmin.createStatement();
        // get the user id for xiFDR
        if (user_id == null) {
            ResultSet userrs = userst.executeQuery("SELECT id FROM core_user WHERE username = 'xiFDR'");
            if (userrs.next()) {
                user_id = userrs.getLong(1);
                userrs.close();
                userst.close();
            } else {
                userrs.close();
                userrs = userst.executeQuery("INSERT INTO core_user (password,is_superuser,username,first_name,last_name,email,is_staff,is_active,date_joined) VALUES ('',false,'xiFDR','xi','FDR','xiFDR@bioanalytik.tu-berlin.de',false,false,now()) RETURNING id");
                userrs.next();
                user_id = userrs.getLong(1);
                userrs.close();
                userst.close();
            }
        } 
       
        // core_fdr
        // status, name, task_uuid, submit_date, notes, config, deleted,
        // owner_id, resultset_uuid
        PreparedStatement pst = this.m_db_connectionAdmin.prepareStatement("INSERT INTO core_fdr "
                + "(status, name, task_uuid, submit_date, notes, config, deleted, owner_id, resultset_uuid) VALUES \n"
                + "(?,      ?,    ?,         now(),       ?,     ?,      false,   ?,        ?) RETURNING id");
        pst.setString(1, "DONE"); // status
        pst.setString(2, name); // name
        pst.setString(3, id.toString()); // task_uuid = resultset_uuid
        pst.setString(4, notes); // notes
        pst.setString(5, ""); // TODO config
        pst.setLong(6, user_id); // owner_id
        pst.setString(7, id.toString()); // resultset_uuid
        ResultSet rsid = pst.executeQuery();
        rsid.next();
        Long fdr_id = rsid.getLong(1);
        rsid.close();
        pst.close();
        // core_fdr_searches
        // fdr_id
        // search_id
        pst = this.m_db_connectionAdmin.prepareStatement("INSERT INTO core_fdr_searches (fdr_id, search_id)"
                + " SELECT " + fdr_id + " as fdr_id , id as search_id from core_search WHERE resultset_uuid = ?");
        for (String searchid : this.m_configs.keySet()) {
            pst.setString(1, searchid);
            pst.execute();
            //pst.addBatch();
        }
        //pst.executeBatch();
        return autocommitAdmin;
    }

    protected void addPSMtoBatch(PreparedStatement pst, UUID id, PSM psm, Xi2ScoreList outScores, FDRResult result) throws SQLException {
        // resultset_id , match_id, search_id, scores
        pst.setObject(resultMatchColumns.resultset_id, id);
        pst.setObject(resultMatchColumns.match_id, UUID.fromString(psm.getPsmID()));
        pst.setObject(resultMatchColumns.search_id, UUID.fromString(psm.getSearchID()));
        // assemble scores
        Double[] scores = new Double[outScores.size()];
        PeptidePair pp = result.peptidePairFDR.filteredGet(psm.getPeptidePair());
        ProteinGroupLink pgl = result.proteinGroupLinkFDR.filteredGet(pp.getLink());
        ProteinGroupPair pgp = pgl == null ? null : result.proteinGroupPairFDR.filteredGet(pgl.getProteinGroupPair());;
        ProteinGroup pg1 = psm.getFdrProteinGroup1();
        ProteinGroup pg2 = psm.getFdrProteinGroup2();
        
        
        scores[0] = psm.getFDR();
        scores[1] = psm.getPEP();
        scores[2] = psm.getScore();
        scores[3] = pp.getFDR();
        scores[4] = pp.getPEP();
        scores[5] = pp.getScore();
        scores[6] = (psm.isLinear() || psm.isNonCovalent() ? Double.NaN : pgl.getFDR());
        scores[7] = (psm.isLinear() || psm.isNonCovalent() ? Double.NaN : pgl.getPEP());
        scores[8] = (psm.isLinear() || psm.isNonCovalent() ? Double.NaN : pgl.getScore());
        scores[9] = (psm.isLinear() || psm.isNonCovalent()  ? Double.NaN : pgp.getFDR());
        scores[10] = (psm.isLinear() || psm.isNonCovalent()  ? Double.NaN : pgp.getPEP());
        scores[11] = (psm.isLinear() || psm.isNonCovalent()  ? Double.NaN : pgp.getScore());
        scores[12] = (pg1 == null ? Double.NaN : pg1.getFDR());
        scores[13] = (pg1 == null ? Double.NaN : pg1.getPEP());
        scores[14] = (pg1 == null ? Double.NaN : pg1.getScore());
        scores[15] = (pg2 == null ? Double.NaN : pg2.getFDR());
        scores[16] = (pg2 == null ? Double.NaN : pg2.getPEP());
        scores[17] = (pg2 == null ? Double.NaN : pg2.getScore());


        for (int i = 18; i< outScores.size(); i++) {
            scores[i] = psm.getDoubleInfo(outScores.get(i).name);
        }
        java.sql.Array scoreArray = pst.getConnection().createArrayOf("float4", scores);
        pst.setObject(resultMatchColumns.scores, scores);
        ResultMatchData rmd = (ResultMatchData) psm.nonPublicFlags.get("ResultMatchData");
        pst.setBoolean(resultMatchColumns.top_ranking, true);
        if (rmd != null) {
            pst.setByte(resultMatchColumns.site1, (byte) psm.getPeptideLinkSite1());
            pst.setByte(resultMatchColumns.site2, (byte) psm.getPeptideLinkSite2());
            pst.setArray(resultMatchColumns.link_score, rmd.link_score);
            pst.setArray(resultMatchColumns.link_score_site1, rmd.link_score_site1);
            pst.setArray(resultMatchColumns.link_score_site2, rmd.link_score_site2);
        }
        pst.setObject(resultMatchColumns.match_group_id, psm.nonPublicFlags.get("match_group_id"));
        pst.addBatch();
    }
    
     @Override
    public int getxiMajorVersion() {
        return 2;
    }

    public RunConfig getConfig() {
        return m_config;
    }

    public RunConfig getConfig(String searchid) {
        return this.m_configs.get(UUID.fromString(searchid));
    }

    public HashMap<String,? extends RunConfig> getConfigs() {
        return this.m_configs;
    }
    
    public Xi2Xi1Config getXi2Config() {
        return this.m_configs.values().iterator().next();
    }

    public Xi2Xi1Config getXi2Config(String search) {
        return this.m_configs.get(search);
    }

    @Override
    public Version getXiVersion() {
        if (this.m_xi_versions.size()>0)
            return this.m_xi_versions.values().iterator().next();
        return null;
    }

    @Override
    public Version getXiVersion(String search) {
        return this.m_xi_versions.get(search);
    }

    public HashMap<String, Version> getXiVersions() {
        return this.m_xi_versions;
    }

        @Override
    public void add(OfflineFDR other) {
        super.add(other);
        if (other instanceof DB2inFDR) {
            this.getSearchIDs().addAll(((XiInFDR)other).getSearchIDs());
            for (Map.Entry<String, Xi2Xi1Config> e : ((DB2inFDR)other).m_configs.entrySet()) {
                this.m_configs.put(e.getKey(), e.getValue());
            }
            this.getXiVersions().putAll(((XiInFDR)other).getXiVersions());
        } else if (other instanceof DBinFDR || other instanceof XiCSVinFDR) {
            this.getSearchIDs().addAll(((XiInFDR)other).getSearchIDs());
            this.getXiVersions().putAll(((XiInFDR)other).getXiVersions());
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, 
                    "Some trouble combining non xi2 database searches. "
                            + "You Could try to read them in the other way around. "
                            + "mzIdenML export might not work correctly.");
        }
    }
}
