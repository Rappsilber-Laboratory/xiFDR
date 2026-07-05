/*
 * Copyright 2015 Lutz Fischer <lfischer at staffmail.ed.ac.uk>.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.rappsilber.fdr;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.rappsilber.data.csv.ColumnAlternatives;
import org.rappsilber.data.csv.CsvParser;
import org.rappsilber.data.csv.condition.CsvCondition;
import org.rappsilber.fdr.dataimport.Xi2Xi1Config;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.fdr.entities.Peptide;
import org.rappsilber.fdr.entities.Protein;
import org.rappsilber.fdr.entities.ProteinGroup;
import org.rappsilber.fdr.gui.FDRGUI;
import org.rappsilber.fdr.gui.components.MZIdentMLOwnerGUI;
import org.rappsilber.fdr.result.FDRResult;
import org.rappsilber.fdr.utils.CalculateWriteUpdate;
import org.rappsilber.fdr.utils.MZIdentMLExport;
import org.rappsilber.fdr.utils.MZIdentMLOwner;
import org.rappsilber.fdr.utils.MaximisingStatus;
import org.rappsilber.utils.UpdatableChar;
import org.rappsilber.utils.Version;
import rappsilber.config.DBRunConfig;
import rappsilber.config.RunConfig;
import rappsilber.config.RunConfigFile;
import rappsilber.ms.sequence.AminoAcid;
import rappsilber.ms.sequence.Sequence;

/**
 *
 * @author lfischer
 */
public class XiCSVinFDR extends CSVinFDR implements XiInFDR{
    
    /** the config file for used to generate the results */
    private RunConfig m_config;
    HashMap<String,Integer> pepids = new HashMap<>();
    HashMap<String,XiSequence> proteins = new HashMap<>();
    HashMap<XiSequence,Integer> proteinIds = new HashMap<>();
    HashMap<String, Integer> scanIDs = new HashMap<>();
    ArrayList<String> searchedFastas = new ArrayList<>();
    HashMap<String,Double> crossLinkerMass = new HashMap<>(1);
    boolean lastMzIDowner = false;

    private boolean markModifications;
    private boolean writeMzID = false;
    private ArrayList<String> m_search_ids = new ArrayList<>();
    private HashMap<String, RunConfig> m_configs = new HashMap<>();
    private HashMap<String, Version> m_xi_versions = new HashMap<>();

    private void markVariableModifiedPSMs() {
        HashSet<Peptide> var_mod_peps=new HashSet<>();
        HashSet<AminoAcid> allvarmods = new HashSet<>();
        for (AminoAcid aa :m_config.getVariableModifications()) 
            allvarmods.add(aa);
        for (AminoAcid aa : m_config.getLinearModifications()) 
            allvarmods.add(aa);
        
        for (Peptide p :allPeptides) {
            Sequence s = new Sequence(p.getSequence(), m_config);
            for (AminoAcid aa : s) {
                if (allvarmods.contains(aa)) {
                    var_mod_peps.add(p);
                    break;
                }
            }
        }
        
        for (PSM psm : allPSMs) {
            if (var_mod_peps.contains(psm.getPeptide1()) || var_mod_peps.contains(psm.getPeptide2()))
                psm.addNegativeGrouping("VarMod");
        }
    }

    private void setWriteMzID(boolean b) {
        this.writeMzID = b;
    }

    private Version parseVersion(String in_file) {
        // find the version string in the file name
        Pattern p = Pattern.compile(".*(?:Xi|xiSEARCH)([0-9](?:\\.[0-9]+(?:\\.[0-9]+)).*)\\.(?-i)(?:csv|tsv|txt)(\\.gz)?$", Pattern.CASE_INSENSITIVE); 
        Matcher m = p.matcher(in_file);
        if (m.matches()) {
            try {
                return new Version(m.group(1));
            } catch (Exception e) {
                return null;
            }
        }
        return null;
    }

    private void addSearchID(String absolutePath) {
        this.m_search_ids.add(absolutePath);
    }

    private void setConfig(String absolutePath, RunConfig conf) {
        if (m_configs == null) {
            m_configs = new HashMap<>();
        }
        for (rappsilber.ms.crosslinker.CrossLinker cl : conf.getCrossLinker()) {
            crossLinkerMass.put(cl.getName(), cl.getCrossLinkedMass());
        }
        m_configs.put(absolutePath, conf);
    }

    @Override
    public Version getXiVersion() {
        if (this.m_xi_versions.size()>0)
            return this.m_xi_versions.values().iterator().next();
        return null;
    }

    @Override
    public Version getXiVersion(String absPath) {
        return m_xi_versions.get(absPath);
    }

    public void setXiVersion(String absolutePath, Version xiVersion) {
        this.m_xi_versions.put(absolutePath, xiVersion);
    }
    
    class XiSequence extends rappsilber.ms.sequence.Sequence {
        
        public XiSequence(String sequence, RunConfig config) {
            super(sequence, config);
        }
        
        public void addSequence(rappsilber.ms.sequence.Sequence s, int pos) {
            if (length()<pos+s.length()) {
                //  pad with Xs
                rappsilber.ms.sequence.AminoAcid[] dummy = new rappsilber.ms.sequence.AminoAcid[pos+s.length()];
                for (int i=length(); i<pos;i++) {
                    dummy[i]=rappsilber.ms.sequence.AminoAcid.X;
                }
                System.arraycopy(m_sequence, 0, dummy, 0, m_sequence.length);
                System.arraycopy(s.m_sequence, 0, dummy, pos, s.m_sequence.length);
            } else {
                for (int i =0; i<s.length(); i++) {
                    rappsilber.ms.sequence.AminoAcid aa = s.m_sequence[i];
                    if (aa instanceof rappsilber.ms.sequence.AminoModification) {
                        aa = ((rappsilber.ms.sequence.AminoModification)aa).BaseAminoAcid;
                    }
                    m_sequence[i+pos]= aa;
                }
            }
        }
        
    }



    public XiCSVinFDR() {
        super();
        HashMap<String,XiSequence> proteins = new HashMap<>();
    }

    public XiCSVinFDR(int[] peptideLengthGroups) {
        super(peptideLengthGroups);
        HashMap<String,XiSequence> proteins = new HashMap<>();
    }
    

    public PSM addMatch(String psmID, String pepSeq1, String pepSeq2, int peplen1, int peplen2, int site1, int site2, boolean isDecoy1, boolean isDecoy2, int charge, double score, String accession1, String description1, String accession2, String description2, int pepPosition1, int pepPosition2, double peptide1score, double peptide2score, String isSpecialCase, String crosslinker, String run, String scan) {

        long pepid1 = m_pepIDs.toIntValue(pepSeq1);
        long pepid2 = m_pepIDs.toIntValue(pepSeq2);
        long protid1 = m_protIDs.toIntValue(accession1);
        long protid2 = m_protIDs.toIntValue(accession2);


        //return addMatch(pepSeq2, pepSeq1, accession1, accession2, protid1, description2, isDecoy1, pepid1, pepPosition1, peplen1, protid2, isDecoy2, pepid2, pepPosition2, peplen2, psmID, site1, site2, charge, score, scoreRatio, isSpecialCase);
        return addMatch(psmID, pepid1, pepid2, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, protid1, accession1, description1, protid2, accession2, description2, pepPosition1, pepPosition2, "","", peptide1score, peptide2score, isSpecialCase, crosslinker,run,scan);
    }    
    
    public String argList() {
        return super.argList() + " --xiconfig=[path to config] --xiversion=[version] --fasta=[path to fasta] --flagModifications --gui --lastowner";
    }
    
    public String argDescription() {
        return super.argDescription() + "\n"
                + "--xiconfig=             what xi config to use to turn find modifications\n"
                + "                        can be used more themn ones and is applied to all\n"
                + "                        input files listed after\n"
                + "--xiversion=            what xiSEARCH version was used\n"
                + "                        can be used more themn ones and is applied to all\n"
                + "                        input files listed after)\n"
                + "--fasta=                fasta file searched\n"
                + "--flagModifications     should modified peptide make their own sub-group\n"
                + "--writemzid             also write out an mzIdentML result file\n"
                + "--lastowner             instead of asking for an mzIdentML document owner\n"
                + "                        reuse the last defined one\n"
                + "--gui                   forward settings to gui\n";
        
    }
    
        
    public String[] parseArgs(String[] argv, FDRSettings settings) {
        ArrayList<String> unknown = new ArrayList<String>();
        argv = super.parseArgs(argv, settings);
        String confpath = null;
        boolean startGUI =  false;
        String xiVersion = null;
        for (String arg : argv) {
            if (arg.toLowerCase().startsWith("--xiconfig=")) {
                try {
                    confpath=arg.substring("--xiconfig=".length());
                    // try to load as xi2 config
                    try {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Try to parse config as xiSEARCH2 config");
                        setConfig(new Xi2Xi1Config(new File(confpath)));
                    } catch (IOException ex) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Try to parse config as xiSEARCH1 config");
                        setConfig(new RunConfigFile(confpath));
                    }
                    
                } catch (IOException ex) {
                    Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.SEVERE, "Can't read config file:", ex);
                    System.exit(-1);
                } catch (ParseException ex) {
                    Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.SEVERE, "Error parsing config file:", ex);
                    System.exit(-1);
                }
            } else if (arg.toLowerCase().startsWith("--xiversion=")) {
                xiVersion = arg.substring(arg.indexOf("=") + 1);
            }else if (arg.toLowerCase().contentEquals("--gui")) {
                startGUI = true;
            }else if (arg.toLowerCase().startsWith("--fasta=")) {
                    String fastapath=arg.substring("--fasta=".length());
                    searchedFastas.add(fastapath);
            }else if (arg.toLowerCase().equals("--lastowner")) {
                    lastMzIDowner = true;
            }else if (arg.toLowerCase().startsWith("--flagmodifications")) {
                setMarkModifications(true);
            }else if (arg.toLowerCase().startsWith("--writemzid")) {
                setWriteMzID(true);
            } else if (arg.startsWith("--mzidfirst=")) {
                mzid_owner.first=arg.substring(arg.indexOf("=")+1);
                setWriteMzID(true);
            } else if (arg.startsWith("--mzidlast=")) {
                mzid_owner.last=arg.substring(arg.indexOf("=")+1);
                setWriteMzID(true);
            } else if (arg.startsWith("--mzidemail=")) {
                mzid_owner.email=arg.substring(arg.indexOf("=")+1);
                setWriteMzID(true);
            } else if (arg.startsWith("--mzidorg=")) {
                mzid_owner.org=arg.substring(arg.indexOf("=")+1);
                setWriteMzID(true);
            } else if (arg.startsWith("--mzidaddr=")) {
                mzid_owner.address=arg.substring(arg.indexOf("=")+1).replace("\\n","\n");
                setWriteMzID(true);
            }  else {
               unknown.add(arg);
               if (confpath != null) {
                   setConfig(arg, getConfig());
               }
               if (xiVersion != null) {
                   setXiVersion(arg, new Version(xiVersion));
               }
            }
            
        }        
        if (startGUI) {
            final FDRGUI fg = FDRGUI.startGUI();
            if (unknown.size() > 0) {
                fg.setInput(unknown.get(0));
            }
            if (searchedFastas.size()>0) {
                fg.setFasta(searchedFastas.get(0));
            }
            if (confpath != null) {
                fg.setXiConfig(confpath);
            }
            fg.setFDRSettings(settings);
            fg.setXiVersion(xiVersion);
            
            String outdir = getCsvOutDirSetting();
            String basename = getCsvOutBaseSetting();
            fg.setCSVOutFile(outdir + ((outdir.endsWith("/") || outdir.endsWith("\\")) ? "":File.separator) + basename +".csv");
//            SwingUtilities.invokeLater(new Runnable() {
//                @Override
//                public void run() {
//                    fg.readAdditionalCSV();
//                }
//            });
//            
            unknown.clear();
            unknown.add("--GUI--");
            
        }
        String[] ret = new String[unknown.size()];
        ret = unknown.toArray(ret);
        
        
        return ret;        
    }
    
    public static void main (String[] argv) throws SQLException, FileNotFoundException {
        
        XiCSVinFDR ofdr = new XiCSVinFDR();
        
        FDRResult res = runCSVFDR(ofdr, argv);
        if (res == null) {
            System.exit(-1);
        }
        //Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.INFO, ofdr.singleCalculation() + 
        //        " " + (ofdr.getConfig() != null) + " " + (ofdr.searchedFastas != null && ofdr.searchedFastas.size()>0) + "\n");
        
        if (ofdr.writeMzID) {
            if (ofdr.singleCalculation() && ofdr.getConfig() != null && (ofdr.searchedFastas != null && ofdr.searchedFastas.size()>0)) {
                // we should be able to write out an mzIdentML
                String path =ofdr.getCsvOutDirSetting();
                String baseName = ofdr.getCsvOutBaseSetting();
                String out = path + (path.endsWith(File.separator)?"":File.separator) + baseName + ".mzid";
                MZIdentMLOwner o = ofdr.getOwner();
                if (o.isEmpty()) {
                    try {


                        if (ofdr.lastMzIDowner) {
                            o.readLast();
                        } else {
                            o = MZIdentMLOwnerGUI.getSetOwner(o);
                        }
                    } catch (Exception ex) {
                        Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
                Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.INFO, "Writing mzIdentML to " + out);
                try {
                    MZIdentMLExport mze = new MZIdentMLExport(ofdr, res, out, o);

                } catch (Exception ex) {
                    Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.SEVERE, "Error exporting mzIdnetML", ex);
                }
            } else 
                System.err.println("could not write out mzIdentML");
        }
            

        System.exit(0);

        
    }

    

    
    /**
     * @return the m_config
     */
    public RunConfig getConfig() {
        return m_config;
    }

    /**
     * @param m_config the m_config to set
     */
    public void setConfig(RunConfig config) {
        this.m_config = config;
        crossLinkerMass.clear();
        for (rappsilber.ms.crosslinker.CrossLinker cl : m_config.getCrossLinker()) {
            crossLinkerMass.put(cl.getName(), cl.getCrossLinkedMass());
        }
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
            return m_config;
        }
        return ret;
    }

    public HashMap<String, Version> getXiVersions() {
        return this.m_xi_versions;
    }

    @Override
    public HashMap<String, ? extends RunConfig> getConfigs() {
        if (m_configs != null) {
            return m_configs;
        }
        HashMap<String, RunConfig> ret = new HashMap<>();
        for (String sid : getSearchIDs()) {
            ret.put(sid,m_config);
        }
        return ret;
    }

    
    @Override
    public ArrayList<String> getSearchIDs() {
        return this.m_search_ids;
    }

    @Override
    public int getFastas(String searchID, ArrayList<String> dbIDs, ArrayList<String> names) {
        int id =0;
        for (String f : searchedFastas) {
            dbIDs.add((id++)+"");
            names.add(f);
        }
        return id;
    }
    
    public void setFastas(ArrayList<String> names) {
        searchedFastas=names;
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

    
    public void matchFastas() throws IOException {
        // parse the fasta - file
        HashMap<String,String> parts2Proteins = new HashMap<>();
        HashMap<String,String> proteinsToSequence = new HashMap<>();
        String protein = null;
        StringBuilder sequnece = new StringBuilder();
        for (String f : this.searchedFastas) {
            int pos = 0;
            try(BufferedReader br = new BufferedReader(new FileReader(f))) {
                for(String line; (line = br.readLine()) != null; ) {
                    if (line.startsWith(">")) {
                        registerFastaProtein(protein, proteinsToSequence, sequnece, parts2Proteins);
                        protein = line.substring(1).trim();
                        sequnece.setLength(0);
                    } else {
                        sequnece.append(line.replaceAll("[^A-Z]", ""));
                    }
                }
            }
            if (protein!=null)
                registerFastaProtein(protein, proteinsToSequence, sequnece, parts2Proteins);
        }
        int countEmptyTarget = 0;
        int countEmptyDecoy = 0;
        for (Protein p : allProteins) {
            if (p.getSequence() == null || p.getSequence().isEmpty()) {
                if (p.isDecoy()) {
                    // try to find a decoy protein
                    String fastaheader =  parts2Proteins.get("rev_" + p.getAccession());
                    if (fastaheader == null)
                        fastaheader =  parts2Proteins.get("ran_" + p.getAccession());
                    if (fastaheader == null)
                        fastaheader =  parts2Proteins.get("decoy:" + p.getAccession());
                    if (fastaheader != null) {
                        // found it - set the sequence
                        p.setSequence(proteinsToSequence.get(fastaheader));
                    } else {
                        // did not find a decoy protein - so will try the target protein
                        fastaheader =  parts2Proteins.get(p.getAccession().toLowerCase());
                        if (fastaheader != null && fastaheader.isEmpty()) {
                            fastaheader = p.getAccession();
                        }
                        if (fastaheader != null && proteinsToSequence.containsKey(fastaheader)) {
                            // only found target protein - assume same size
                            p.setSize(proteinsToSequence.get(fastaheader).length());
                        } else {
                            // nothing found
                            countEmptyDecoy++;
                        }
                    }
                    
                } else {
                    String fastaheader =  parts2Proteins.get(p.getAccession().toLowerCase());
                    if (fastaheader == null || fastaheader.isEmpty()) {
                        fastaheader = p.getAccession();
                    }
                    
                    if (fastaheader != null & proteinsToSequence.containsKey(fastaheader) ) {
                            p.setSequence(proteinsToSequence.get(fastaheader));
                    } else {
                        countEmptyTarget++;
                    }
                }
            }
        }
        if (countEmptyTarget + countEmptyDecoy >0) {
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, (countEmptyTarget + countEmptyDecoy) + " proteins where not assigned a sequence/size (" + countEmptyTarget +" taget," + countEmptyDecoy+ " decoy)");
        }
    }

    protected void registerFastaProtein(String protein, HashMap<String, String> proteinsToSequence, StringBuilder sequnece, HashMap<String, String> parts2Proteins) {
        if (protein != null && !proteinsToSequence.containsKey(protein)) {
            // the last protein needs to be stored
            proteinsToSequence.put(protein, sequnece.toString());
            // split the fasta header in parts
            String[] parts = protein.split("[\\|]");
            HashSet<String> searched = new HashSet<>();
            registerFastaHeaderParts(parts, searched, parts2Proteins, protein);
            parts = protein.split("[\\s\\t]");
            registerFastaHeaderParts(parts, searched, parts2Proteins, protein);
            parts = protein.split("[\\.]");
            registerFastaHeaderParts(parts, searched, parts2Proteins, protein);
            parts = protein.split("[\\|\\s\\t\\.]");
            registerFastaHeaderParts(parts, searched, parts2Proteins, protein);
            parts = protein.split("[\\|\\s\\t\\.:]");
            registerFastaHeaderParts(parts, searched, parts2Proteins, protein);
        }
    }

    public void registerFastaHeaderParts(String[] parts, HashSet<String> searched, HashMap<String, String> parts2Proteins, String protein) {
        for (String p : parts) {
            if (!searched.contains(p)) {
                searched.add(p);
                String pa = parts2Proteins.get(p.toLowerCase());
                if (pa == null) {
                    // we haven't seen that part before so assume it uniquely identifies the protein
                    parts2Proteins.put(p.toLowerCase(), protein);
                } else if (!pa.isEmpty()) {
                    // it was seen before so it does not uniquely identifies any protein
                    // therefore we don't store this link
                    parts2Proteins.put(p.toLowerCase(), "");
                }
            }
        }
    }


    public boolean readCSV(File csv, Version xiVersion) throws FileNotFoundException, IOException, ParseException {
        if (!this.m_search_ids.contains(csv.getAbsolutePath())) {
            this.m_search_ids.add(csv.getAbsolutePath());
        }
        boolean ret = this.readCSV(csv);
        if (xiVersion == null)
            xiVersion = parseVersion(csv.getAbsolutePath());
        this.addSearchID(csv.getAbsolutePath());
        this.setXiVersion(csv.getAbsolutePath(), xiVersion);
        return ret;
    }

    public boolean readCSV(CsvParser csv, CsvCondition filter, Version xiVersion, RunConfig conf) throws FileNotFoundException, IOException, ParseException {
        boolean ret =  this.readCSV(csv, filter);
        String in_file = csv.getInputFile().getAbsolutePath();
        if (xiVersion == null)
            xiVersion = parseVersion(in_file);
        this.addSearchID(csv.getInputFile().getAbsolutePath());
        this.setXiVersion(csv.getInputFile().getAbsolutePath(), xiVersion);
        if (conf != null) {
            this.setConfig(csv.getInputFile().getAbsolutePath(), conf);
        } else if (getConfig(in_file) == null && getConfig() != null) {
            setConfig(in_file, getConfig());
        }
        
        return ret;
    }
    
    @Override
    public boolean readCSV(CsvParser csv, CsvCondition filter) throws FileNotFoundException, IOException, ParseException {
        String infile = csv.getInputFile().getAbsolutePath();
        if (!this.m_search_ids.contains(infile)) {
            this.m_search_ids.add(infile);
        }
        if (!this.m_configs.containsKey(infile) && this.getConfig() != null) {
            this.m_configs.put(infile, this.getConfig());
        }
        boolean ret = super.readCSV(csv, filter); //To change body of generated methods, choose Tools | Templates.
        if (ret)
            matchFastas();
        if (markModifications()) {
            markVariableModifiedPSMs();
        }
        String absPath = csv.getInputFile().getAbsolutePath();
        if (this.getConfig(absPath) == null) {
            this.setConfig(absPath, this.getConfig());
        }
        Version xiVersion = this.getXiVersion(absPath);
        if (xiVersion == null) {
            xiVersion = parseVersion(absPath);
            setXiVersion(absPath, xiVersion);
        }
        return ret;
    }

    @Override
    public int getxiMajorVersion() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    
    @Override
    public void add(OfflineFDR other) {
        super.add(other);
        if (other instanceof XiInFDR) {
            this.getSearchIDs().addAll(((XiInFDR)other).getSearchIDs());
            for (Map.Entry<String, ? extends RunConfig> e : ((XiInFDR)other).getConfigs().entrySet()) {
                this.setConfig(e.getKey(), e.getValue());
            }
            this.getXiVersions().putAll(((XiInFDR)other).getXiVersions());
        }
    }
    
}
