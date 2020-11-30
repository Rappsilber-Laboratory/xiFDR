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
import java.io.IOException;
import java.sql.SQLException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.rappsilber.data.csv.ColumnAlternatives;
import org.rappsilber.data.csv.CsvParser;
import org.rappsilber.data.csv.condition.CsvCondition;
import org.rappsilber.fdr.entities.DBPSM;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.fdr.entities.Peptide;
import org.rappsilber.fdr.entities.Protein;
import org.rappsilber.utils.IntArrayList;
import org.rappsilber.utils.UpdatableChar;
import rappsilber.config.RunConfig;
import rappsilber.config.RunConfigFile;
import rappsilber.ms.sequence.AminoAcid;
import rappsilber.ms.sequence.AminoModification;
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

    private boolean markModifications;




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
    

    /**
     * adds a psm to the list folds up the scores to peptidespairs links
     * proteinpairs and proteins
     *
     * @param psmID
     * @param pepid1
     * @param pepid2
     * @param peplen1
     * @param peplen2
     * @param site1
     * @param site2
     * @param charge
     * @param score
     * @param proteinId1
     * @param proteinId2
     * @param pepPosition1
     * @param pepPosition2
     * @param scoreRation
     * @return a peptide pair that is supported by the given match
     */
    @Override
    public PSM addMatch(String psmID, org.rappsilber.fdr.entities.Peptide peptide1, org.rappsilber.fdr.entities.Peptide peptide2, int peplen1, int peplen2, int site1, int site2, int charge, double score, Protein proteinId1, Protein proteinId2, int pepPosition1, int pepPosition2, double peptide1score, double peptide2score, String isSpecialCase, String crosslinker, String run, String scan) {
        org.rappsilber.fdr.entities.Peptide npepid1;
        org.rappsilber.fdr.entities.Peptide npepid2;
        int npeplen1;
        int npeplen2;
        int nsite1;
        int nsite2;
        Protein nproteinId1;
        Protein nproteinId2;
        int npepPosition1;
        int npepPosition2;
        int protcomp = proteinId1.compareDecoyUnAware(proteinId2);
        int pepcomp = peptide1.compare(peptide2);
        int sitecomp = (site1 - site2);
        double pep1score = peptide1score;
        double pep2score = peptide2score;
        double npep1score;
        double npep2score;

        if (protcomp < 0 || (protcomp == 0 && pepcomp < 0) || (protcomp == 0 && pepcomp == 0 && site1 < site2)) {
            npepid1 = peptide1;
            npepid2 = peptide2;
            npeplen1 = peplen1;
            npeplen2 = peplen2;
            nsite1 = (byte)site1;
            nsite2 = (byte)site2;
            nproteinId1 = proteinId1;
            nproteinId2 = proteinId2;
            npepPosition1 = pepPosition1;
            npepPosition2 = pepPosition2;
            npep1score = pep1score;
            npep2score = pep2score;

        } else {
            npepid1 = peptide2;
            npepid2 = peptide1;
            npeplen1 = peplen2;
            npeplen2 = peplen1;
            nsite1 = (byte)site2;
            nsite2 = (byte)site1;
            nproteinId1 = proteinId2;
            nproteinId2 = proteinId1;
            npepPosition1 = pepPosition2;
            npepPosition2 = pepPosition1;
            npep1score = pep2score;
            npep2score = pep1score;
        }

        if (!PSMScoreHighBetter) {
            score = 10 - (10 * score);
        }


        DBPSM psm = new DBPSM(psmID, npepid1, npepid2, (byte)nsite1, (byte)nsite2, proteinId1.isDecoy(), proteinId2.isDecoy(), (byte)charge, score, npep1score, npep2score);
        psm.setNegativeGrouping(isSpecialCase);
        psm.setRun(registerRun(run));
        psm.setScan(scan);
        psm.setCrosslinker(crosslinker);
        rappsilber.ms.sequence.Sequence ps1= new rappsilber.ms.sequence.Sequence(peptide1.getSequence(), m_config);
        peptide1.mass = ps1.getWeight();
        if (peptide2 != null && peptide2 != Peptide.NOPEPTIDE) {
            rappsilber.ms.sequence.Sequence ps2= new rappsilber.ms.sequence.Sequence(peptide2.getSequence(), m_config);
            peptide2.mass = ps2.getWeight();
            Double xlmass = crossLinkerMass.get(crosslinker);
            psm.setCalcMass(peptide1.mass + peptide2.mass + xlmass);
        } else {
            psm.setCalcMass(peptide1.mass);
        }

        
                    String modLoockup = npepid1.getSequence();
                    if (npepid1 != null)
                        modLoockup+=npepid2.getSequence();
                    Sequence m = new Sequence(modLoockup, m_config);
                    for (AminoAcid aa : m) {
                        if (aa instanceof AminoModification) {
                            if (m_config.getVariableModifications().contains(aa)) {
                                psm.setHasVarMods(true);
                            }
                            if (m_config.getFixedModifications().contains(aa)) {
                                psm.setHasFixedMods(true);
                            }
                        }
                    }
                    if (psm.hasVarMods()) {
//                        psm.setModified(true);
//                        if (markModifications())
//                            psm.addNegativeGrouping("Modified");
                    }
         

        PSM regpsm = getAllPSMs().register(psm);

        return regpsm;
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
        return super.argList() + " --xiconfig=[path to config] --fasta=[path to fasta] --flagModifications";
    }
    
    public String argDescription() {
        return super.argDescription() + "\n"
                + "--xiconfig=             what xi config to use to turn find modifications\n"
                + "--fasta=                fasta file searched"
                + "--flagModifications     should modified peptide make their own sub-group\n";
        
    }
    
        
    public String[] parseArgs(String[] argv, FDRSettings settings) {
        ArrayList<String> unknown = new ArrayList<String>();
        argv = super.parseArgs(argv, settings);
        
        for (String arg : argv) {
            if (arg.toLowerCase().startsWith("--xiconfig=")) {
                try {
                    String confpath=arg.substring("--xiconfig=".length());
                    m_config = new RunConfigFile(confpath);
                } catch (IOException ex) {
                    Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.SEVERE, "Can't read config file:", ex);
                    System.exit(-1);
                } catch (ParseException ex) {
                    Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.SEVERE, "Error parsing config file:", ex);
                    System.exit(-1);
                }
                
            }else if (arg.toLowerCase().startsWith("--fasta=")) {
                    String confpath=arg.substring("--fasta=".length());
                    searchedFastas.add(confpath);
                    
            }else if (arg.toLowerCase().startsWith("--flagmodifications")) {
                setMarkModifications(true);
            }  else {
               unknown.add(arg);
            }
            
        }        
        String[] ret = new String[unknown.size()];
        ret = unknown.toArray(ret);
        return ret;        
    }
    
    public static void main (String[] argv) throws SQLException, FileNotFoundException {
        
        XiCSVinFDR ofdr = new XiCSVinFDR();
        
        FDRSettings settings = new FDRSettingsImpl();
        String[] files = ofdr.parseArgs(argv, settings);
        
        readFromCommandline(ofdr, files);

        if (ofdr.cofilterArgs != null) {
            ofdr.cofilter = new XiCSVinFDR();
            String[] cofiles = ofdr.cofilter.parseArgs(argv, settings);
            readFromCommandline((XiCSVinFDR)ofdr.cofilter, cofiles);
        }
        
        Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.INFO, "Calculate FDR");
        ofdr.calculateWriteFDR(ofdr.getCsvOutDirSetting(), ofdr.getCsvOutBaseSetting(), ",", settings);



        System.exit(0);

        
    }

    protected static void readFromCommandline(XiCSVinFDR ofdr, String[] files) {
        // assume that everything that was not matched to an argument is a file
        if (files.length == 0) {
            ofdr.printUsage();
            System.exit(1);
        }
        
        
        if (ofdr.getCsvOutDirSetting() != null) {
            if (ofdr.getCsvOutBaseSetting() == null) {
                ofdr.setCsvOutBaseSetting("FDR");
            }
        }
        
        if (ofdr.getCsvOutBaseSetting() != null) {
            if (ofdr.getCsvOutDirSetting() == null)
                ofdr.setCsvOutDirSetting(".");
            
            System.out.println("writing results to " + ofdr.getCsvOutDirSetting() + "/" + ofdr.getCsvOutBaseSetting() + "*");
            System.out.flush();
        }

        CsvParser csv = new CsvParser();
        if (ofdr.getCommandLineColumnMapping() != null) {
            for (String[] map : ofdr.getCommandLineColumnMapping()) {
                for (int i = 1; i<map.length;i++) {
                    csv.setAlternative(map[0], map[i]);
                }
            }
        } else {
            ColumnAlternatives.setupAlternatives(csv,CSVinFDR.DEFAULT_COLUMN_MAPPING);
        }
        
        
        // read in all files
        for (String f : files) {
            Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.INFO, "seeting up csv input");
            
            UpdatableChar delimChar = new UpdatableChar(',');
            UpdatableChar quoteChar = new UpdatableChar('"');
            if (ofdr.getDelimiter() != null) {
                delimChar.value = ofdr.getDelimiter();
            }
            if (ofdr.getQuote() != null) {
                quoteChar.value = ofdr.getQuote();
            }
            if (ofdr.getDelimiter() == null || ofdr.getQuote() == null) {
                try {
                    csv.guessDelimQuote(new File(f), 50, delimChar, quoteChar);
                } catch (IOException ex) {
                    Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.SEVERE, "error while quessing csv-definitions", ex);
                }
            }
            try {
                csv.openFile(new File(f), true);
            } catch (IOException ex) {
                Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.SEVERE, "Could not read the file", ex);
                System.exit(-1);
            }
            
            Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.INFO, "Read datafrom CSV");
            try {
                ofdr.readCSV(csv,(CsvCondition)null);
            } catch (IOException ex) {
                Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.SEVERE, "Error while reading file: " + f, ex);
                System.exit(-1);
            } catch (ParseException ex) {
                Logger.getLogger(XiCSVinFDR.class.getName()).log(Level.SEVERE, "Error parsing file: " + f, ex);
                System.exit(-1);
            }
        }
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
    public void setConfig(RunConfig m_config) {
        this.m_config = m_config;
        crossLinkerMass.clear();
        for (rappsilber.ms.crosslinker.CrossLinker cl : m_config.getCrossLinker()) {
            crossLinkerMass.put(cl.getName(), cl.getCrossLinkedMass());
        }
    }


    @Override
    public RunConfig getConfig(int searchid) {
        return m_config;
    }

    @Override
    public IntArrayList getSearchIDs() {
        return new IntArrayList(new int[]{0});
    }

    @Override
    public int getFastas(int searchID, ArrayList<Integer> dbIDs, ArrayList<String> names) {
        int id =0;
        for (String f : searchedFastas) {
            dbIDs.add(id++);
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
                        fastaheader =  parts2Proteins.get(p.getAccession());
                        if (fastaheader != null) {
                            // only found target protein - assume same size
                            p.setSize(proteinsToSequence.get(fastaheader).length());
                        } else {
                            // nothing found
                            countEmptyDecoy++;
                        }
                    }
                    
                } else {
                    String fastaheader =  parts2Proteins.get(p.getAccession());
                    if (fastaheader != null) {
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
            String[] parts = protein.split("[\\|\\s\\t\\.]");
            for (String p : parts) {
                String pa = parts2Proteins.get(p.toLowerCase());
                if (pa == null) {
                    // we haven't seen that part before so assume it uniquely identifies the protein
                    parts2Proteins.put(p.toLowerCase(), protein);
                } else if (!pa.isEmpty()) {
                    // it was seen before so it does not uniquely identifies any protein
                    // therefore we don't store this link
                    parts2Proteins.put(p, "");
                }
            }
        }
    }

    @Override
    public void readCSV(CsvParser csv, CsvCondition filter) throws FileNotFoundException, IOException, ParseException {
        super.readCSV(csv, filter); //To change body of generated methods, choose Tools | Templates.
        matchFastas();
    }
    
    
    
}
