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

import org.rappsilber.fdr.result.FDRResult;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.fdr.entities.PeptidePair;
import org.rappsilber.fdr.utils.StreamReplaceWriter;
import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.FragmentationTable;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.jmzidml.model.mzidml.Param;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidenceRef;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionHypothesis;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionList;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectraData;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.model.mzidml.SubstitutionModification;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.jmzidml.model.utils.MzIdentMLVersion;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

public class MZIdentXLFDR extends OfflineFDR {

    private MzIdentMLUnmarshaller unmarshaller;
    //private URL xmlFileURL = JmzIdentMLParser.class.getClassLoader().getResource("CPTAC_Progenesis_Identifications.mzid");
    //private URL xmlFileURL = JmzIdentMLParser.class.getClassLoader().getResource("55merge_mascot_full.mzid");
    private String xmlFile;
    private List<String> peptideList = new ArrayList();
    private List<String> dbSequenceList = new ArrayList();
    private List<String> matrix = new ArrayList();
    private List<String> matrix_sii = new ArrayList();
    private List<String> matrix_pe = new ArrayList();
    private HashMap<String, PeptideEvidence> peptideEvidenceIdHashMap = new HashMap<String, PeptideEvidence>();
    private HashMap<String, SpectraData> spectraDataIdHashMap = new HashMap<String, SpectraData>();
    private HashMap<String, SpectrumIdentificationItem> siiIdHashMap = new HashMap<String, SpectrumIdentificationItem>();
    private HashMap<SpectrumIdentificationItem, Double> psmScores = new HashMap<SpectrumIdentificationItem, Double>();
    private HashMap<String, SpectrumIdentificationResult> siiIdToSirHashMap = new HashMap<String, SpectrumIdentificationResult>();
    private HashMap<Integer, String> columnToScoreMap = new HashMap<Integer, String>();
    private HashMap<Integer, String> columnToProtScoreMap = new HashMap<Integer, String>();
    private HashMap<String, DBSequence> dbSequenceIdHashMap = new HashMap<String, DBSequence>();
    private HashMap<String, ProteinDetectionHypothesis> pdhIdHashMap = new HashMap<String, ProteinDetectionHypothesis>();
    private HashMap<String, List<ProteinDetectionHypothesis>> peptide_pdh_HashMap = new HashMap<String, List<ProteinDetectionHypothesis>>();
    private HashMap<String, ArrayList<SpectrumIdentificationItem>> crosslinkedPSM = new HashMap<String, ArrayList<SpectrumIdentificationItem>>();
    private HashMap<String, ArrayList<SpectrumIdentificationItem>> PSMidToCrosslinkedPSM = new HashMap<String, ArrayList<SpectrumIdentificationItem>>();
    private HashMap<String, uk.ac.ebi.jmzidml.model.mzidml.Peptide> peptideIdHashMap = new HashMap<String, uk.ac.ebi.jmzidml.model.mzidml.Peptide>();
    private ArrayList<SpectrumIdentificationItem> linearPSM = new ArrayList<SpectrumIdentificationItem>();
//    private List<ProteinDetectionHypothesis> proteinDetectionHypothesisList = new ArrayList<ProteinDetectionHypothesis>();
//    private List<ProteinAmbiguityGroup> proteinAmbiguityGroupList = new ArrayList<ProteinAmbiguityGroup>();
    private ArrayList<SpectrumIdentificationResult> sirList = new ArrayList();
    
    ProteinDetectionList proteinDetectionList = new ProteinDetectionList();
    private String sep = ",";
    private String pagHeader = "PAG ID" + sep + "PAG score" + sep + "protein accession" + sep + "Pass Threshold (Protein)" + sep + "description" + sep + "group membership" + sep;
    private String spectrumHeader = "Raw data location" + sep + "Spectrum ID" + sep + "Spectrum Title" + sep + "Retention Time (s)" + sep;
    private String psmHeader = "PSM_ID" + sep + "rank" + sep + "Pass Threshold" + sep + "Calc m/z" + sep + "Exp m/z" + sep + "Charge" + sep + "Sequence" + sep + "Modifications";
    private String pScoreHeader = "";     //Protein score header will be set only after reading the file
    private String scoreHeader = "";     //This will be set only after reading the file
    private String endPsmHeader = sep + "proteinacc_start_stop_pre_post_;" + sep + "Is decoy";
    private String representativeProteinAcc = "MS:1001591";     //Used to identify the representative of each group - only used for the one line export of PAGS

    
    //    protected String crosslinkedModAcc = "MS:8888888";     
    /**
     * cvTerm used to identify modifications, that span several peptides
     * This is the cvTerm for the Modification, that holds the mass of the cross-linker
     */
    private String crosslinkedDonorModAcc = "MS:1002509";     
    /**
     * cvTerm used to identify modifications, that span several peptides
     * This is the cvTerm for the Modification, that holds a zero mass to denote the second (third, forth ...) site a cross-linker is attached to
     */
    private String crosslinkedReceptorModAcc = "MS:1002510";     
//    protected String crosslinkedSIIAcc = "MS:9999999";     
    /** 
     * cvTerm used to identify members of cross-linked PSMs 
     */
    private String crosslinkedSIIAcc = "MS:1002511";     
//    /** 
//     * cvTerm used to identify members of cross-linked PSMs 
//     * This one denotes the "beta" peptide
//     */
//    private String crosslinkedReceptorSIIAcc = "MS:9999992";     
    /** Identifies a mzIdentML-file, that represents a cross-link search */
    private String crosslinkedSearchAcc = "MS:9999XXX";     

    /** XiFDR mzIdenML - id*/
    private String analysesSoftwareFDR =
                    "      <AnalysisSoftware version=\"%XIFDRVERSION%\" name=\"XiFDR\" id=\"fdr_software\">\n" +
                    "      </AnalysisSoftware>\n";

    
    
    /**Used to identify modifications, that span several peptides*/
    protected String PSMScore = "score";     
    
    private boolean scoreFound = false;
    
    private Boolean isVerbose = true;
    
    private boolean deletePassThreshold = true;
//    
//    public static void main(String[] args) {
//        MZIdentXLFDR mzidToCsv = new MZIdentXLFDR();
//        
//        
//        //TODO - Undecided which if any command line arguments to include - minimally need to know whether to export Peptides or PAGs
//
//        if (args != null && args.length == 3) {
//            mzidToCsv.unmarshaller = new MzIdentMLUnmarshaller(new File(args[0]));
//            mzidToCsv.init(args[1],args[2]);
//
//        } else {
//
//            mzidToCsv.unmarshaller = new MzIdentMLUnmarshaller(new File(mzidToCsv.xmlFile));
//            mzidToCsv.init("out.csv","exportPSMs");
//            //mzidToCsv.init("out.csv","exportProteinGroups");
//
//            //System.out.println("Error - correct usage MzIdentMLToCSV inputFile.mzid outputFile.csv [exportProteinGroups|exportPSMs|exportProteinsOnly]");
//            //System.exit(1);
//        }
//        
//        
//    }
//    

    public void readMzIdentML(File f, boolean passThreshHoldOnly) throws FileNotFoundException, IOException, ParseException {
        MzIdentMLUnmarshaller mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(f);
        //MZIdentXLFDR mzidToCsv = new MZIdentXLFDR();
        this.unmarshaller = mzIdentMLUnmarshaller;
        this.init(passThreshHoldOnly);
    }


    private void init(boolean passThreshHoldOnly) {

        try {

//            BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));            

            //Read all the objects we will need into hashes that are not automatically resolved by object reference

            if (isVerbose) {
                System.out.print("About to iterate over PepEvid...");
            }
            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"About to iterate over PepEvid...");

            Iterator<PeptideEvidence> iterPeptideEvidence = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.PeptideEvidence);
            while (iterPeptideEvidence.hasNext()) {
                PeptideEvidence peptideEvidence = iterPeptideEvidence.next();
                peptideEvidenceIdHashMap.put(peptideEvidence.getId(), peptideEvidence);
            }

            if (isVerbose) {
                System.out.println("...done");
                System.out.print("About to iterate over uk.ac.ebi.jmzidml.model.mzidml.Peptide");
            }
            
            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"About to iterate over uk.ac.ebi.jmzidml.model.mzidml.Peptide");
            
            Iterator<uk.ac.ebi.jmzidml.model.mzidml.Peptide> iterPeptide = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
            while (iterPeptide.hasNext()) {
                uk.ac.ebi.jmzidml.model.mzidml.Peptide peptide = iterPeptide.next();
                peptideIdHashMap.put(peptide.getId(), peptide);
            }

            if (isVerbose) {
                System.out.println("...done");
                System.out.print("About to iterate over Spectra Data");
            }
            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"About to iterate over Spectra Data");
            
            Iterator<SpectraData> iterSpectraData = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectraData);
            while (iterSpectraData.hasNext()) {
                SpectraData spectraData = iterSpectraData.next();
                spectraDataIdHashMap.put(spectraData.getId(), spectraData);
            }


            if (isVerbose) {
                System.out.println("...done");
                System.out.print("About to iterate over DBsequence");
            }
            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"About to iterate over DBsequence");
            
            Iterator<DBSequence> iterDBSequence = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
            while (iterDBSequence.hasNext()) {
                DBSequence dbSequence = iterDBSequence.next();
                dbSequenceIdHashMap.put(dbSequence.getId(), dbSequence);
            }

            if (isVerbose) {
                System.out.println("...done");
                System.out.print("About to iterate over PDH");
            }
            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"About to iterate over PDH");
            
            Iterator<ProteinDetectionHypothesis> iterPDH = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.ProteinDetectionHypothesis);
            Integer pCounter = 0;
            while (iterPDH.hasNext()) {
                ProteinDetectionHypothesis pdh = iterPDH.next();
                pdhIdHashMap.put(pdh.getId(), pdh);

                for (CvParam cvParam : pdh.getCvParam()) {
                    if (cvParam.getAccession().equals("MS:1001591") || cvParam.getAccession().equals("MS:1001592") || cvParam.getAccession().equals("MS:1001593")
                            || cvParam.getAccession().equals("MS:1001594") || cvParam.getAccession().equals("MS:1001595") || cvParam.getAccession().equals("MS:1001596")
                            || cvParam.getAccession().equals("MS:1001597")
                            || cvParam.getAccession().equals("MS:1001598")
                            || cvParam.getAccession().equals("MS:1001599")) {       
                        //do nothing - these are specifically handled
                        //ToDO this code could be improved using an array of values...
                    } else if (cvParam.getValue() != null) {
                        if (!columnToProtScoreMap.containsValue(cvParam.getName())) {
                            columnToProtScoreMap.put(pCounter, cvParam.getName());
                            pCounter++;
                        }
                    }
                }

                for (UserParam userParam : pdh.getUserParam()) {
                    if (!columnToProtScoreMap.containsValue(userParam.getName())) {
                        columnToProtScoreMap.put(pCounter, userParam.getName());
                        pCounter++;
                    }

                }


            }

            for (int i = 0; i < pCounter; i++) {
                pScoreHeader += columnToProtScoreMap.get(i) + sep;
            }

            //Now let's see what scores we have in the file
            //TODO - I'm not sure this is the fastest way to parse the files; these are unmarshalled again below - inefficient?
            //Iterator<SpectrumIdentificationItem> iterSII = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationItem);

            Integer counter = 0;

            if (isVerbose) {
                System.out.println("...done");
                System.out.print("About to iterate over SIR");
            }
            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"About to iterate over SIR");
            
            Iterator<SpectrumIdentificationResult> iterSIR = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);

            while (iterSIR.hasNext()) {
                SpectrumIdentificationResult sir = iterSIR.next();
                sirList.add(sir);

                List<SpectrumIdentificationItem> listSII = sir.getSpectrumIdentificationItem();

                for (SpectrumIdentificationItem sii : listSII) {
                    siiIdHashMap.put(sii.getId(), sii);
                    siiIdToSirHashMap.put(sii.getId(), sir);
                    boolean isCrosslinked = false;
                    for (CvParam cvParam : sii.getCvParam()) {
                        String cvpv = cvParam.getValue();
                        Double score = null;
                        if (cvpv != null) {
                            String accession = cvParam.getAccession();
                            if (accession.contentEquals(getCrosslinkedSIIAcc())) {
                                isCrosslinked = true;
                                ArrayList<SpectrumIdentificationItem> xlSIIs =  crosslinkedPSM.get(cvpv);
                                if (xlSIIs == null) {
                                    xlSIIs = new ArrayList<SpectrumIdentificationItem>();
                                    crosslinkedPSM.put(cvpv, xlSIIs);
                                    PSMidToCrosslinkedPSM.put(sii.getId(), xlSIIs);
                                }
                                xlSIIs.add(sii);
                            }
                            if (!columnToScoreMap.containsValue(cvParam.getName())) {
                                columnToScoreMap.put(counter, cvParam.getName());
                                counter++;
                            }
                            if (cvParam.getAccession().contentEquals(getPSMScore()) || 
                                    cvParam.getName().contentEquals(getPSMScore()) || 
                                    (score == null &&  cvParam.getName().toLowerCase().endsWith(":" + getPSMScore().toLowerCase()))) {
                                score = Double.parseDouble(cvParam.getValue());
                                psmScores.put(sii, score);
                                scoreFound =true;
                            }
                        }
                    }
                    if (!isCrosslinked) {
                        if (sii.isPassThreshold() || !passThreshHoldOnly)
                            linearPSM.add(sii);
                    }
                }
            }

            

            for (int i = 0; i < counter; i++) {
                scoreHeader += sep + columnToScoreMap.get(i);
            }

            if (isVerbose) {
                System.out.println("...done");
                if (scoreFound) {
                    System.out.println("Scores found");
                } else {
                    System.err.println("!!!!!!!!!!!!!! NO SCORES FOUND !!!!!!!!!!!!!!!!!");
                }
                System.out.print("register linear matches");
            }
            if (scoreFound) {
                Logger.getLogger(this.getClass().getName()).log(Level.INFO,"Scores found");
            } else {
                Logger.getLogger(this.getClass().getName()).log(Level.INFO,"!!!!!!!!!!!!!! NO SCORES FOUND !!!!!!!!!!!!!!!!!");
            }
            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"register linear matches");
            
            boolean linearDecoys = false;
            
            for (SpectrumIdentificationItem sii : linearPSM) {
                uk.ac.ebi.jmzidml.model.mzidml.Peptide pep = peptideIdHashMap.get(sii.getPeptideRef());
                List<PeptideEvidenceRef> pepevs =sii.getPeptideEvidenceRef();
                
                for (PeptideEvidenceRef pepevref : pepevs) {
                    PeptideEvidence pepev = peptideEvidenceIdHashMap.get(pepevref.getPeptideEvidenceRef());
                    String seqref =  pepev.getDBSequenceRef();
                    
                    int pepstart = pepev.getStart();
                    //int pepend = pepev.getEnd();
                    String pepSeq = pep.getPeptideSequence();
                    // if we have a loop link we need to get the link positions
                    ArrayList<Integer> linkPositions = new ArrayList<Integer>(2);
                    for (Modification m : pep.getModification()) {
                        pepSeq += " - " + m.getLocation() + " " + modToString(m);
                        // check wether we have a loop link
                        for (CvParam cvp : m.getCvParam()) {
                            String AccString =  cvp.getAccession();
                            if (AccString.contentEquals(getCrosslinkedDonorModAcc()) || AccString.contentEquals(getCrosslinkedReceptorModAcc())) {
                                // ok we have a loop link
                                linkPositions.add(m.getLocation());
                            }
                        }
                    }
                    
                    DBSequence seq = dbSequenceIdHashMap.get(seqref);
                    String acc = seq.getAccession();
                    String desc = seq.getName();
                    Double score = psmScores.get(sii);
                    if (score == null)
                        score = 0.0;
                    //boolean isDecoy = seq.getSearchDatabase().
                            
                         
                    if ( pepev.isIsDecoy())
                        linearDecoys = true;
                    
                    
                    if (linkPositions.size() > 0) {
                        if (linkPositions.size() == 2) {
                            addMatch(sii.getId(), pepSeq, null, pep.getPeptideSequence().length(), 0, linkPositions.get(0), linkPositions.get(1), pepev.isIsDecoy(), false, sii.getChargeState(),score, acc, desc, null, null, pepstart, pepstart, score,0, null);
                        } else {
                            System.err.println(sii.getId() + ": Currently only loop links with exactly two links within the peptide are supported - will add this match as linear (non-cross-linked) match");
                            Logger.getLogger(this.getClass().getName()).log(Level.WARNING, "{0}Currently only loop links with exactly two links within the peptide are supported - will add this match as linear (non-cross-linked) match", sii.getId());
                            addMatch(sii.getId(), pepSeq, null, pep.getPeptideSequence().length(), 0, -1, -1, pepev.isIsDecoy(), false, sii.getChargeState(),score, acc, desc, null, null, pepstart, pepstart, score,0, null);
                        }
                    } else {
                        addMatch(sii.getId(), pepSeq, null, pep.getPeptideSequence().length(), 0, -1, -1, pepev.isIsDecoy(), false, sii.getChargeState(),score, acc, desc, null, null, pepstart, pepstart, score,0, null);
                    }
                    
                }
            }
            
            if (isVerbose) {
                System.out.println("...done");
                if (linearDecoys) {
                    System.out.println("linear decoys found");
                } else {
                    System.err.println("!!!!!!!!!! no linear decoys found!!!!!!!!!!!!!!");
                }
                System.out.print("register cross-linked matches");
            }
            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"register cross-linked matches");
            
            boolean crosslinkedDecoys = false;
            int countCrosslinkedPSM = 0;
            for (ArrayList<SpectrumIdentificationItem> xlsiis : crosslinkedPSM.values()) {
                SpectrumIdentificationItem sii1 = xlsiis.get(0);
                SpectrumIdentificationItem sii2 = xlsiis.get(1);
                uk.ac.ebi.jmzidml.model.mzidml.Peptide pep1 = peptideIdHashMap.get(sii1.getPeptideRef());
                List<PeptideEvidenceRef> pepevs1 =sii1.getPeptideEvidenceRef();
                int peplinksite1 =-1;
                        
                String pepSeq1 = pep1.getPeptideSequence();
                for (Modification m : pep1.getModification()) {
                    pepSeq1 += " - " + m.getLocation() + " " + modToString(m);
                    for (CvParam cvp : m.getCvParam()) {
                        if (cvp.getAccession().contentEquals(getCrosslinkedReceptorModAcc()) || cvp.getAccession().contentEquals(getCrosslinkedDonorModAcc())) {
                            peplinksite1 = m.getLocation();
                            break;
                        }
                    }
                }
                

                uk.ac.ebi.jmzidml.model.mzidml.Peptide pep2 = peptideIdHashMap.get(sii2.getPeptideRef());
                List<PeptideEvidenceRef> pepevs2 =sii2.getPeptideEvidenceRef();
                int peplinksite2 =-1;
                
                String pepSeq2 = pep2.getPeptideSequence();
                for (Modification m : pep2.getModification()) {
                    pepSeq2 += " - " + m.getLocation() + " " + modToString(m);
                    for (CvParam cvp : m.getCvParam()) {
                        if (cvp.getAccession().contentEquals(getCrosslinkedReceptorModAcc()) || cvp.getAccession().contentEquals(getCrosslinkedDonorModAcc())) {
                            peplinksite2 = m.getLocation();
                            break;
                        }
                    }
                }
                
                
                
                for (PeptideEvidenceRef pepevref1 : pepevs1) {
                    for (PeptideEvidenceRef pepevref2 : pepevs2) {
                        PeptideEvidence pepev1 = peptideEvidenceIdHashMap.get(pepevref1.getPeptideEvidenceRef());
                        String seqref1 =  pepev1.getDBSequenceRef();
                        PeptideEvidence pepev2 = peptideEvidenceIdHashMap.get(pepevref2.getPeptideEvidenceRef());
                        String seqref2 =  pepev2.getDBSequenceRef();

                        int pepstart1 = pepev1.getStart();
//                        int pepend1 = pepev1.getEnd();

                        int pepstart2 = pepev2.getStart();
//                        int pepend2 = pepev2.getEnd();
                        
                        DBSequence seq1 = dbSequenceIdHashMap.get(seqref1);
                        String acc1 = seq1.getAccession();
                        String desc1 = seq1.getName();
                        Double score1 = psmScores.get(sii1);
                        DBSequence seq2 = dbSequenceIdHashMap.get(seqref2);
                        String acc2 = seq2.getAccession();
                        String desc2 = seq2.getName();
                        Double score2 = psmScores.get(sii2);
                        //boolean isDecoy = seq.getSearchDatabase().
                        double score  = 0;
                        double peptide1score = 0;
                        double peptide2score = 0;
                        double scoreRatio = 0;
                        if (score1 == score2)  {
                            score  = score1;
                            scoreRatio = (4/5+pepSeq1.length()/pepSeq2.length())/2;
                            peptide1score = score *scoreRatio;
                            peptide1score = score * (1-scoreRatio);
                        } else  {
                            // @TODO need to change that
                            score  = score1+score2;
                            peptide1score = score1;
                            peptide2score = score2;
                        }

                        if (pepev1.isIsDecoy() || pepev2.isIsDecoy())
                            crosslinkedDecoys = true;
                        
                        countCrosslinkedPSM++;
                        addMatch(sii1.getId(), pepSeq1, pepSeq2, pep1.getPeptideSequence().length(), pep2.getPeptideSequence().length(), peplinksite1, peplinksite2, pepev1.isIsDecoy(), pepev2.isIsDecoy(), sii1.getChargeState(),score, acc1, desc1, acc2, desc2, pepstart1, pepstart2, peptide1score, peptide2score, null);
                    }
                    
                }
            }
            if (isVerbose) {
                System.out.println("...done");
                
            }
            
            if (crosslinkedDecoys == false) {
                System.err.println("!!!!!!!!!!! NO crosslinked decoys found !!!!!!!!!!!!!!!!!");
                Logger.getLogger(this.getClass().getName()).log(Level.WARNING,"!!!!!!!!!!! NO crosslinked decoys found !!!!!!!!!!!!!!!!!");
            }
            if (countCrosslinkedPSM == 0) {
                System.err.println("!!!!!!!!!!! NO crosslinked PSM found !!!!!!!!!!!!!!!!!");
                Logger.getLogger(this.getClass().getName()).log(Level.WARNING,"!!!!!!!!!!! NO crosslinked PSM found !!!!!!!!!!!!!!!!!");
            }
                
                
            
        } catch (Exception e) {
            System.err.println(e.getMessage());
            e.printStackTrace();
        }

    }






    /*
     * Method to convert an mzid Mod element into a string of type ModName@location
     */
    private String modToString(Modification mod) {

        String modString = "";

        if (mod.getCvParam() != null) {

            for (CvParam cvParam : mod.getCvParam()) {

                modString += cvParam.getName();
            }
        } else {

            if (mod.getMonoisotopicMassDelta() != null) {
                modString += mod.getMonoisotopicMassDelta();
            } else if (mod.getAvgMassDelta() != null) {
                modString += mod.getAvgMassDelta();
            }
        }

        if (mod.getLocation() != null) {
            modString += ":" + mod.getLocation();
        }
        return modString;

    }

    /*
     * Method to create and return a string representation of a substitution modification
     */
    private String subModToString(SubstitutionModification subMod) {
        return subMod.getOriginalResidue() + "=>" + subMod.getReplacementResidue() + "@" + subMod.getLocation();

    }
    
    
    

//    /**
//     * Used to identify members of cross-linked PSMs
//     * @return the crosslinkedSIIAcc
//     */
//    public String getCrosslinkedSIIAcc() {
//        return crosslinkedSIIAcc;
//    }

    /**
//     * Used to identify members of cross-linked PSMs
//     * @param crosslinkedSIIAcc the crosslinkedSIIAcc to set
//     */
//    public void setCrosslinkedSIIAcc(String crosslinkedSIIAcc) {
//        this.crosslinkedSIIAcc = crosslinkedSIIAcc;
//    }

//    /**
//     * CvTerm that is used to identify modifications, that span several peptides
//     * @return the crosslinkedModAcc
//     */
//    public String getCrosslinkedModAcc() {
//        return crosslinkedModAcc;
//    }
//
//    /**
//     * CvTerm that is used to identify modifications, that span several peptides
//     * @param crosslinkedModAcc the crosslinkedModAcc to set
//     */
//    public void setCrosslinkedModAcc(String crosslinkedModAcc) {
//        this.crosslinkedModAcc = crosslinkedModAcc;
//    }

    /**
     * Used to identify modifications, that span several peptides
     * @return the PSMScore
     */
    public String getPSMScore() {
        return PSMScore;
    }

    /**
     * Used to identify modifications, that span several peptides
     * @param PSMScore the PSMScore to set
     */
    public void setPSMScore(String PSMScore) {
        this.PSMScore = PSMScore;
    }
    
    
    protected void assignFDR(SpectrumIdentificationItem sii, String psmID) {
        
    }
    
    public void writeMZIdentML(String  mzidFileName, FDRResult result) {
        if (deletePassThreshold) {
            for (SpectrumIdentificationItem sii : siiIdHashMap.values()) {
                sii.setPassThreshold(false);
            }
        }
        
//        for (PSM psm : getFDRLinearPSMs()) {
//            String id = psm.getPsmID();
//            SpectrumIdentificationItem sii = siiIdHashMap.get(id);
//            sii.setPassThreshold(true);
//            PeptidePair pp = psm.getFdrPeptidePair();
//            if (pp != null) {
//                CvParam cvParamfdrscore = new CvParam();
//                Cv cv = new Cv();
//                cv.setId("PSI-MS");
//                cv.setUri("http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");
//                cv.setVersion("2.25.0");
//                cv.setFullName("PSI-MS");
//
//                cvParamfdrscore.setAccession("MS:1001364");
//                cvParamfdrscore.setName("pep:global FDR");
//                cvParamfdrscore.setValue(Double.toString(pp.getFDR()));
//                cvParamfdrscore.setCv(cv);
//                sii.getCvParam().add(cvParamfdrscore);
//            }
//        }
        
        
        for (PSM psm : result.getPsmFDR()) {
            PeptidePair pp = psm.getFdrPeptidePair();
            if (pp != null) {
                if (pp.isLinear()) {
                    SpectrumIdentificationItem sii = siiIdHashMap.get(psm.getPsmID());
                    sii.setPassThreshold(true);
                    CvParam cvParamfdrscore = new CvParam();
                    Cv cv = new Cv();
                    cv.setId("PSI-MS");
                    cv.setUri("http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");
                    cv.setVersion("2.25.0");
                    cv.setFullName("PSI-MS");

                    cvParamfdrscore.setAccession("MS:1001364");
                    cvParamfdrscore.setName("pep:global FDR");
                    cvParamfdrscore.setValue(Double.toString(pp.getFDR()));
                    cvParamfdrscore.setCv(cv);
                    sii.getCvParam().add(cvParamfdrscore);
                } else {
                    ArrayList<SpectrumIdentificationItem> siiList = PSMidToCrosslinkedPSM.get(psm.getPsmID());
                    for (SpectrumIdentificationItem sii : siiList) {
                        CvParam cvParamfdrscore = new CvParam();
                        Cv cv = new Cv();
                        cv.setId("PSI-MS");
                        cv.setUri("http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");
                        cv.setVersion("2.25.0");
                        cv.setFullName("PSI-MS");

                        cvParamfdrscore.setAccession("MS:1001364");
                        cvParamfdrscore.setName("pep:global FDR");
                        cvParamfdrscore.setValue(Double.toString(pp.getFDR()));
                        cvParamfdrscore.setCv(cv);
                        sii.getCvParam().add(cvParamfdrscore);
                        sii.setPassThreshold(true);
                    }
                }
            }
        }
        
        
        try {
            String outFile = mzidFileName;
            if (!outFile.endsWith(".mzid")) {
                outFile = outFile + ".mzid";
            }
            FileOutputStream fwriter = new FileOutputStream(outFile);
            OutputStreamWriter writer = new OutputStreamWriter(fwriter, "UTF-8");
            
            MzIdentMLMarshaller marshaller;
            marshaller = new MzIdentMLMarshaller(MzIdentMLVersion.Version_1_3);
    
            writer.write(marshaller.createXmlHeader() + "\n");

            
            AnalysisSoftwareList analysisSoftwareList;
            AuditCollection auditCollection;
            Provider provider;
            AnalysisProtocolCollection analysisProtocolCollection;
            CvList cvList;
            AnalysisCollection analysisCollection;
            Inputs inputs;
            String searchDatabase_Ref;
            cvList = unmarshaller.unmarshal(MzIdentMLElement.CvList);
            //analysisSoftwareList = mzIdentML.getAnalysisSoftwareList();
            analysisSoftwareList = unmarshaller.unmarshal(MzIdentMLElement.AnalysisSoftwareList);
            //auditCollection = mzIdentML.getAuditCollection();
            auditCollection = unmarshaller.unmarshal(MzIdentMLElement.AuditCollection);
            //provider = mzIdentML.getProvider();
            provider = unmarshaller.unmarshal(MzIdentMLElement.Provider);
            // analysisProtocolCollection = mzIdentML.getAnalysisProtocolCollection();
            analysisProtocolCollection = unmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);
            //analysisCollection = mzIdentML.getAnalysisCollection();
            analysisCollection = unmarshaller.unmarshal(MzIdentMLElement.AnalysisCollection);
            //inputs = mzIdentML.getDataCollection().getInputs();
            inputs = unmarshaller.unmarshal(MzIdentMLElement.Inputs);
            searchDatabase_Ref = inputs.getSearchDatabase().get(0).getId();
            
            
            
            writer.write(marshaller.createMzIdentMLStartTag("12345") + "\n");
            
            
            if (cvList != null) {
                
                marshaller.marshal(cvList, writer);
            }
            writer.write("\n");

            AnalysisSoftware asXi = new AnalysisSoftware();
            Param xiNameParam = new Param();
            UserParam xiUP = new UserParam();
            xiUP.setName("XiFDR");
            xiNameParam.setParam(xiUP);
            asXi.setSoftwareName(xiNameParam);
            asXi.setVersion(OfflineFDR.getXiFDRVersion().toString());
            
            if (analysisSoftwareList != null) {
                analysisSoftwareList.getAnalysisSoftware().add(asXi);
                
                marshaller.marshal(analysisSoftwareList, writer);
            } else {
                AnalysisSoftwareList asl = new AnalysisSoftwareList();
                asl.getAnalysisSoftware().add(asXi);
                marshaller.marshal(analysisSoftwareList, writer);
            }
            writer.write("\n");
            
            if (provider != null) {
                marshaller.marshal(provider, writer);
            }
            writer.write("\n");
            
            if (auditCollection != null) {
                marshaller.marshal(auditCollection, writer);
            }
            writer.write("\n");
            
            SequenceCollection sequenceCollection = unmarshaller.unmarshal(MzIdentMLElement.SequenceCollection);
            
            if (sequenceCollection != null) {
                marshaller.marshal(sequenceCollection, writer);
            }
            writer.write("\n");
            
            
            
            if (analysisCollection != null) {
                marshaller.marshal(analysisCollection, writer);
            }
            writer.write("\n");
            
            
            if (analysisProtocolCollection != null) {
                marshaller.marshal(analysisProtocolCollection, writer);
            }
            writer.write("\n");
            
            
            writer.write(marshaller.createDataCollectionStartTag() + "\n");
            
            writer.write("\n");
            
            if (inputs != null) {
                marshaller.marshal(inputs, writer);
            }
            writer.write("\n");
            
            writer.write(marshaller.createAnalysisDataStartTag() + "\n");
            
            
            String spectrumIdentificationListRef = "";
            if (analysisCollection.getSpectrumIdentification().size() > 0) {
                spectrumIdentificationListRef = analysisCollection.getSpectrumIdentification().get(0).getSpectrumIdentificationListRef();
            }
            SpectrumIdentificationList siList;
            siList = new SpectrumIdentificationList();
            
            siList.setId(spectrumIdentificationListRef);
            
            Iterator<FragmentationTable> iterFragmentationTable = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.FragmentationTable);
            while (iterFragmentationTable.hasNext()) {
                
                FragmentationTable fr = iterFragmentationTable.next();
                siList.setFragmentationTable(fr);
            }
            
            
            HashMap<String, String> sii_stringMap = new HashMap();
            
            
            Iterator<SpectrumIdentificationResult> iterSpectrumIdentificationResult = sirList.iterator();
            while (iterSpectrumIdentificationResult.hasNext()) {
                
                SpectrumIdentificationResult sr = iterSpectrumIdentificationResult.next();
                

                siList.getSpectrumIdentificationResult().add(sr);
            }

//            if (siListList != null) {
//                
//            }
            
            marshaller.marshal(siList, writer);
            writer.write("\n");
            
            
            writer.write(marshaller.createProteinDetectionListStartTag("PDL_1", null) + "\n");
            
            
            
            writer.write(marshaller.createProteinDetectionListClosingTag() + "\n");
            writer.write(marshaller.createAnalysisDataClosingTag() + "\n");
            writer.write(marshaller.createDataCollectionClosingTag() + "\n");
            
            writer.write(marshaller.createMzIdentMLClosingTag());
            
            writer.close();
            
            System.out.println("Output written to " + outFile);
            
        } catch (IOException e) {
            e.printStackTrace();
        }
    }    

    @Override
    public String getSource() {
        return xmlFile;
    }

    /**
     * cvTerm used to identify modifications, that span several peptides
     * This is the cvTerm for the Modification, that holds the mass of the cross-linker
     * @return the crosslinkedDonorModAcc
     */
    public String getCrosslinkedDonorModAcc() {
        return crosslinkedDonorModAcc;
    }

    /**
     * cvTerm used to identify modifications, that span several peptides
     * This is the cvTerm for the Modification, that holds the mass of the cross-linker
     * @param crosslinkedDonorModAcc the crosslinkedDonorModAcc to set
     */
    public void setCrosslinkedDonorModAcc(String crosslinkedDonorModAcc) {
        this.crosslinkedDonorModAcc = crosslinkedDonorModAcc;
    }

    /**
     * cvTerm used to identify modifications, that span several peptides
     * This is the cvTerm for the Modification, that holds a zero mass to denote the second (third, forth ...) site a cross-linker is attached to
     * @return the crosslinkedReceptorModAcc
     */
    public String getCrosslinkedReceptorModAcc() {
        return crosslinkedReceptorModAcc;
    }

    /**
     * cvTerm used to identify modifications, that span several peptides
     * This is the cvTerm for the Modification, that holds a zero mass to denote the second (third, forth ...) site a cross-linker is attached to
     * @param crosslinkedReceptorModAcc the crosslinkedReceptorModAcc to set
     */
    public void setCrosslinkedReceptorModAcc(String crosslinkedReceptorModAcc) {
        this.crosslinkedReceptorModAcc = crosslinkedReceptorModAcc;
    }

    /**
     * cvTerm used to identify members of cross-linked PSMs
     * @return the crosslinkedDonorSIIAcc
     */
    public String getCrosslinkedSIIAcc() {
        return crosslinkedSIIAcc;
    }

    /**
     * cvTerm used to identify members of cross-linked PSMs
     * @param crosslinkedDonorSIIAcc the crosslinkedDonorSIIAcc to set
     */
    public void setCrosslinkedSIIAcc(String crosslinkedDonorSIIAcc) {
        this.crosslinkedSIIAcc = crosslinkedDonorSIIAcc;
    }


    /**
     * Identifies a mzIdentML-file, that represents a cross-link search
     * @return the crosslinkedSearchAcc
     */
    public String getCrosslinkedSearchAcc() {
        return crosslinkedSearchAcc;
    }

    /**
     * Identifies a mzIdentML-file, that represents a cross-link search
     * @param crosslinkedSearchAcc the crosslinkedSearchAcc to set
     */
    public void setCrosslinkedSearchAcc(String crosslinkedSearchAcc) {
        this.crosslinkedSearchAcc = crosslinkedSearchAcc;
    }
}