/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr.utils;

//import uk.ac.liv.mzidlib.util.CVUtils;
//import uk.ac.liv.mzidlib.util.Utils;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.lang.reflect.Method;
import java.text.ParseException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;
import org.rappsilber.fdr.XiInFDR;

import uk.ac.ebi.jmzidml.model.mzidml.AbstractContact;
import uk.ac.ebi.jmzidml.model.mzidml.Affiliation;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSampleCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.ContactRole;
import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.Enzyme;
import uk.ac.ebi.jmzidml.model.mzidml.Enzymes;
import uk.ac.ebi.jmzidml.model.mzidml.FileFormat;
import uk.ac.ebi.jmzidml.model.mzidml.InputSpectra;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.IonType;
import uk.ac.ebi.jmzidml.model.mzidml.Measure;
import uk.ac.ebi.jmzidml.model.mzidml.ModificationParams;
import uk.ac.ebi.jmzidml.model.mzidml.Organization;
import uk.ac.ebi.jmzidml.model.mzidml.Param;
import uk.ac.ebi.jmzidml.model.mzidml.ParamList;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidenceRef;
import uk.ac.ebi.jmzidml.model.mzidml.Person;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;
import uk.ac.ebi.jmzidml.model.mzidml.Role;
import uk.ac.ebi.jmzidml.model.mzidml.SearchDatabase;
import uk.ac.ebi.jmzidml.model.mzidml.SearchDatabaseRef;
import uk.ac.ebi.jmzidml.model.mzidml.SearchModification;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectraData;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIDFormat;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentification;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationProtocol;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.model.mzidml.Tolerance;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.model.utils.MzIdentMLVersion;

//import uk.ac.liv.mzidconverters.ReadUnimod;
//import uk.ac.liv.unimod.ModT;

import org.rappsilber.fdr.OfflineFDR;
import org.rappsilber.fdr.dataimport.Xi2Xi1Config;
//import org.rappsilber.fdr.entities.DBPSM;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.fdr.entities.PeptidePair;
import org.rappsilber.fdr.entities.Protein;
import org.rappsilber.fdr.entities.ProteinGroupLink;
import org.rappsilber.fdr.entities.ProteinGroupPair;
import org.rappsilber.fdr.entities.ProteinGroup;
import org.rappsilber.fdr.result.FDRResult;
import org.rappsilber.xlmod.XLMOD;
import org.rappsilber.xlmod.XLModEntry;
import org.rappsilber.xlmod.XLModQuery;
import rappsilber.config.RunConfig;
import rappsilber.ms.crosslinker.NonCovalentBound;
import rappsilber.ms.sequence.AminoAcid;
import rappsilber.ms.sequence.AminoModification;
import rappsilber.ms.sequence.ions.AIon;
import rappsilber.utils.Util;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisData;
import uk.ac.ebi.jmzidml.model.mzidml.DataCollection;
import uk.ac.ebi.jmzidml.model.mzidml.MzIdentML;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideHypothesis;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinAmbiguityGroup;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionHypothesis;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionList;
import uk.ac.ebi.jmzidml.model.mzidml.SpecificityRules;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItemRef;
import uk.ac.liv.mzidlib.converters.ReadUnimod;

/**
 *
 * @author lfischer
 */
public class MZIdentMLExport {

    String[][] psiEnzymes = new String[][]{
        new String[]{"MS:1001251", "Trypsin"},
        new String[]{"MS:1001313", "Trypsin/P"},
        new String[]{"MS:1001303", "Arg-C"},
        new String[]{"MS:1001304", "Asp-N"},
        new String[]{"MS:1001305", "Asp-N_ambic"},
        new String[]{"MS:1001306", "Chymotrypsin"},
        new String[]{"MS:1001307", "CNBr"},
        new String[]{"MS:1001308", "Formic_acid"},
        new String[]{"MS:1001309", "Lys-C"},
        new String[]{"MS:1001310", "Lys-C/P"},
        new String[]{"MS:1001311", "PepsinA"},
        new String[]{"MS:1001312", "TrypChymo"},
        new String[]{"MS:1001314", "V8-DE"},
        new String[]{"MS:1001315", "V8-E"},
        new String[]{"MS:1001915", "leukocyte elastase"},
        new String[]{"MS:1001916", "proline endopeptidase"},
        new String[]{"MS:1001917", "glutamyl endopeptidase"},
        new String[]{"MS:1001918", "2-iodobenzoate"},
        new String[]{"MS:1001955", "no cleavage"},
        new String[]{"MS:1001956", "unspecific cleavage"}};
    


    
    HashMap<String,String[]> name2Enzyme=new HashMap<String, String[]>();
    {
        for (String[] e : psiEnzymes) {
            name2Enzyme.put(e[0],e);
            name2Enzyme.put(e[1],e);
            name2Enzyme.put(e[1].toLowerCase(),e);
            name2Enzyme.put(e[1].toLowerCase().replace(" ", "").replace("_", "").replace("-", ""),e);
        }
    }

    OfflineFDR fdr;
    FDRResult fdrResult;
    /**
     * Used to identify members of cross-linked PSMs
     */
    private static String crosslinkedSIIAcc = "MS:1002511";
    private static String crosslinkedSIIName = "cross-link spectrum identification item";
    /**
     * Used to identify members of non-covalent PSMs
     */
    private static String noncovalentSIIAcc = "MS:1003331";
    private static String noncovalentSIIName = "noncovalently associated peptides spectrum identification item";
    /**
     * Used to identify modifications, that span several peptides.
     * This is denotes the donor
     */
    private static String crosslinkedDonorModAcc = "MS:1002509";
    private static String crosslinkedDonorModName = "cross-link donor";
    /**
     * Used to identify modifications, that span several peptides.
     * This is denotes the receiver
     */
    private static String crosslinkedAcceptorModAcc = "MS:1002510";
    private static String crosslinkedAcceptorModName = "cross-link acceptor";
    // crosslinker stub    
    private static String crosslinkerStubAcc= "MS:1003390";
    private static String crosslinkerStubName = "crosslinker cleavage characteristics";
    // name to xi crosslinker lookup
    private HashMap<String,rappsilber.ms.crosslinker.CrossLinker> name2xiXL = new HashMap<>();
    //These are the main structures to be output by the main writing method
    private SequenceCollection sequenceCollection;
    private HashMap<String,HashMap<String,SpectrumIdentificationList>> siList = new HashMap<>();
//    private HashMap<String,SpectrumIdentificationProtocol> siProtocol;
//    CvList cvList;
    private AnalysisSoftwareList analysisSoftwareList;
//    Provider provider;
    private AuditCollection auditCollection;
    private AnalysisSampleCollection analysisSampleCollection;
    private AnalysisCollection analysisCollection;
    private AnalysisProtocolCollection analysisProtocolCollection;
    private ProteinDetectionList proteinDetectionList;
    private HashMap<org.rappsilber.fdr.entities.Protein, DBSequence> protein2DBSequence = new HashMap<Protein, DBSequence>();

    private Inputs inputs;
    //Some IDs to be used throughout;
    private static String siiListID = "SII_LIST_1";
    private static String spectraDataID = "SID_1";
    private static String psiCvID = "PSI-MS";
    private static String siProtocolID = "SearchProtocol_1";
//    private static String searchDBID = "SearchDB_1";
//    private static String pepEvidenceListID = "PepEvidList_1";
    private static String specIdentID = "SpecIdent_";
    private static String unimodID = "UNIMOD";
    private static String xlmodID = "XLMOD";
    private static String unitCvID = "UO";
    private static String measureMzID = "Measure_MZ";
    private static String measureIntID = "Measure_Int";
    private static String measureErrorID = "Measure_Error";
    private static String sourceFileID = "SourceFile_1";
//    private String decoyRegex = null;    //Added by ARJ for setting is decoy
    //Some objects we will need globally
    private CvList cvList;
    private Cv xlmodCV;
    private Cv unimodCV;
    private Cv psiCV;
    private Cv unitCV;
    private Provider provider;
    private HashMap<String,SpectrumIdentificationProtocol> siProtocol = new HashMap<>();
//    private SearchDatabase searchDB;
//    private SpectraData spectraData;
    //private List<SpectraData> spectraDataList;
    private Person docOwner;
    private AnalysisSoftware analysisSoftware;
    private Map<String, String> cvMap = null;
    private HashMap<org.rappsilber.fdr.entities.Peptide, uk.ac.ebi.jmzidml.model.mzidml.Peptide> uniquePeps;
    //HashMap<org.rappsilber.fdr.entities.PeptidePair, HashMap<org.rappsilber.fdr.entities.Peptide, uk.ac.ebi.jmzidml.model.mzidml.Peptide>> uniquePeptidePairs;
    private HashMap<org.rappsilber.fdr.entities.PeptidePair,uk.ac.ebi.jmzidml.model.mzidml.Peptide[]> uniquePeptidePairs;
    private HashMap<org.rappsilber.fdr.entities.PeptidePair,PeptidePair> firstFoundPeptidePairs;

    //int pepCounter = 0;
    int pepEvidCounter = 0;
    //List<SpectrumIdentificationResult> specIdentResults;
    ReadUnimod unimodDoc;
    double unimodMassError = 0.01;              //TODO This parameter is hard-coded (ARJ changed from 0.001 to 0.01 - Aug2012; perhaps should be set dynamically from search params)
    //defaults:
    private boolean isMs2SpectrumIdStartingAtZero = false;
    private String databaseFileFormatID;
    private String databaseFileFormatName;
    private String massSpecFileFormatID;
    private String massSpecFileFormatName;
    boolean outputFragmentation = true;
    
    private String forceExtension;
    
    XLMOD xlmod;
    
    
    private MZIdentMLOwner owner = new MZIdentMLOwner("fisrt", "last", "email", "org", "address");


    private PSMToMzIdentScanId mgfScanID = new PSMToMzIdentScanId() {
        @Override
        public String getID(PSM psm) {
            return "index=" + psm.getFileScanIndex();
        }
    };

    private MzMLScanTranslation mzMLscanIDTranslation = new MzMLScanTranslation("controllerType=0 controllerNumber=1 scan=%s%");

    private PSMToMzIdentScanId scanIDTranslation = new PSMToMzIdentScanId() {
        @Override
        public String getID(PSM psm) {
            return "index=" + psm.getScan();
        }
    };



    private PSMToMzIdentScanId thermoRawScaNumber = new PSMToMzIdentScanId() {
        @Override
        public String getID(PSM psm) {
            return "index=" + (Long.parseLong(psm.getScan()) - 1);
        };
    };
    
    
    /**
     *
     * @param inputfile	: the X!Tandem output file (BIOML xml format)
     * @param outputfile : the name of the new mzIdentML file to create
     * @param databaseFileFormatID : optional. If set to null, then we try to
     * find the right code based on the file extension found in the xtandem
     * parameter "list path, sequence source #1" found in the xtandem file,
     * falling back to "MS:1001348","FASTA format" if it cannot infer the format
     * based on the extension .
     * @param massSpecFileFormatID : optional. If set to null, then we try to
     * find the right code based on the file extension found in the xtandem
     * parameter "spectrum, path" found in the xtandem file.
     * @param isMs2SpectrumIdStartingAtZero : set this to true if the spectra
     * file originally submitted to X!Tandem had spectrum id numbering staring
     * at 0 (e.g. is the case with mzML format). This is important because then
     * we can calculate the correct number for the new mzIdentML file which
     * always has to start from 0 (this conforms to the specifications of the
     * controlled vocabulary item MS:1000774 where spectrumID should start from
     * 0). CV available at: /resources/CV_psi-ms.obo.txt
     * @param decoyRegularExpression : optional. if a the referenced protein
     * accession from PeptideEvidence contains this string value, the attribute
     * isDecoy will be set to true
     * @param outputFragmentation : optional. If set to null, by default
     * fragment ions are output (produces much larger mzid files). If set to
     * false, fragment ions are not exported.
     *
     *
     * @throws Exception
     */
    //CV also available at http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo
    public MZIdentMLExport(OfflineFDR fdr, FDRResult fdrResult, String outputfile, String databaseFileFormatID, String massSpecFileFormatID, MZIdentMLOwner owner) throws Exception {
        //XTandemFile xfile = new XTandemFile(inputfile);
        this.owner = owner;
        initializeVariables(fdr, fdrResult, databaseFileFormatID, massSpecFileFormatID);
        convertFile(fdr, fdrResult, outputfile);
    }

    
    /**
     *
     * @param inputfile	: the X!Tandem output file (BIOML xml format)
     * @param outputfile : the name of the new mzIdentML file to create
     * @param databaseFileFormatID : optional. If set to null, then we try to
     * find the right code based on the file extension found in the xtandem
     * parameter "list path, sequence source #1" found in the xtandem file,
     * falling back to "MS:1001348","FASTA format" if it cannot infer the format
     * based on the extension .
     * @param massSpecFileFormatID : optional. If set to null, then we try to
     * find the right code based on the file extension found in the xtandem
     * parameter "spectrum, path" found in the xtandem file.
     * @param isMs2SpectrumIdStartingAtZero : set this to true if the spectra
     * file originally submitted to X!Tandem had spectrum id numbering staring
     * at 0 (e.g. is the case with mzML format). This is important because then
     * we can calculate the correct number for the new mzIdentML file which
     * always has to start from 0 (this conforms to the specifications of the
     * controlled vocabulary item MS:1000774 where spectrumID should start from
     * 0). CV available at: /resources/CV_psi-ms.obo.txt
     * @param decoyRegularExpression : optional. if a the referenced protein
     * accession from PeptideEvidence contains this string value, the attribute
     * isDecoy will be set to true
     * @param outputFragmentation : optional. If set to null, by default
     * fragment ions are output (produces much larger mzid files). If set to
     * false, fragment ions are not exported.
     *
     *
     * @throws Exception
     */
    //CV also available at http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo
    public MZIdentMLExport(OfflineFDR fdr, FDRResult fdrResult, String databaseFileFormatID, String massSpecFileFormatID, MZIdentMLOwner owner) throws Exception {
        //XTandemFile xfile = new XTandemFile(inputfile);
        this.owner = owner;
        initializeVariables(fdr, fdrResult, databaseFileFormatID, massSpecFileFormatID);
    }    
    
    public MZIdentMLExport(OfflineFDR fdr, FDRResult fdrResult, String outputfile, MZIdentMLOwner owner) throws Exception {
        this(fdr, fdrResult, outputfile, "MS:1001348", "MS:1001062", owner);
    }

    public MZIdentMLExport(OfflineFDR fdr, FDRResult fdrResult, MZIdentMLOwner owner) throws Exception {
        this(fdr, fdrResult, "MS:1001348", "MS:1001062", owner);
    }
    
    public void readXLMOD() throws IOException, ParseException {
        xlmod = new XLMOD();
        xlmod.read();

    }
    
    
    
    
    /**
     * This method initializes some of the global variables based on the values
     * given in the parameters.
     *
     * @param databaseFileFormatID : optional. If set to null, then we try to
     * find the right code based on the file extension found in the xtandem
     * parameter "list path, sequence source #1" found in the xtandem file,
     * falling back to "MS:1001348","FASTA format" if it cannot infer the format
     * based on the extension .
     * @param massSpecFileFormatID : optional. If set to null, then we try to
     * find the right code based on the file extension found in the xtandem
     * parameter "spectrum, path" found in the xtandem file.
     * @param isMs2SpectrumIdStartingAtZero
     * @param decoyRegularExpression
     * @param outputFragmentation
     * @throws IOException
     */
    private void initializeVariables(OfflineFDR fdr, FDRResult result, String databaseFileFormatID, String massSpecFileFormatID) throws IOException {
        this.fdr = fdr;
        this.fdrResult  = result;
//        if (decoyRegularExpression != null && !decoyRegularExpression.equals("")) {
//            this.decoyRegex = decoyRegularExpression;
//        }

        this.outputFragmentation = outputFragmentation;

        this.databaseFileFormatID = databaseFileFormatID;
        this.databaseFileFormatName = getCVName(databaseFileFormatID);
        this.massSpecFileFormatID = massSpecFileFormatID;
        this.massSpecFileFormatName = getCVName(massSpecFileFormatID);

//        if (isMs2SpectrumIdStartingAtZero == null) {
            //if file format is mzML (MS:1000584), then spectrum starts at 0, otherwise it starts at 1
            if (this.massSpecFileFormatID.equalsIgnoreCase("MS:1000584")) {
                this.isMs2SpectrumIdStartingAtZero = true;
            } else {
                this.isMs2SpectrumIdStartingAtZero = false;
            }
//        } else {
//            this.isMs2SpectrumIdStartingAtZero = isMs2SpectrumIdStartingAtZero;
//        }
    }

    /**
     * Gets the Controlled Vocabulary (CV) item name for the given item ID
     *
     * @param cvItemID
     * @return
     * @throws IOException
     */
    private String getCVName(String cvItemID) throws IOException {
        //If CV map is not yet initialized, do it:
        if (this.cvMap == null) {
            this.cvMap = uk.ac.liv.mzidlib.util.Utils.getInitializedCVMap();
        }

        //validate:
        if (this.cvMap.get(cvItemID) == null) {
            throw new RuntimeException("Given item not found in Controlled Vocabulary : " + cvItemID);
        } else {
            return this.cvMap.get(cvItemID);
        }
    }

    public void convertFile(String outputfile) throws Exception {
        convertFile(fdr, fdrResult, outputfile);
    }

    private void convertFile(OfflineFDR fdr, FDRResult fdrResult, String outputfile) throws Exception {
        try {
            unimodDoc = new ReadUnimod();
            readInFDR(fdr, fdrResult, 20);
            writeMzidFile(outputfile);
        } catch (Exception e) {
            throw e;
        }
    }

    public void readInFDR(OfflineFDR fdr, FDRResult result, double massError) {
        xlmod = new XLMOD();
        try {
            xlmod.read();
        } catch (IOException ex) {
            Logger.getLogger(MZIdentMLExport.class.getName()).log(Level.SEVERE, null, ex);
        } catch (ParseException ex) {
            Logger.getLogger(MZIdentMLExport.class.getName()).log(Level.SEVERE, null, ex);
        }
        rappsilber.config.RunConfig conf = ((XiInFDR)fdr).getConfig();

        // Setup the mzid objects
        handleCVs();
        ArrayList<CvParam> stubs = parseStubs(conf);
        

        HashMap<String, DBSequence> foundProts = new HashMap<String, DBSequence>();
        HashMap<String, PeptideEvidence> pepEvidLookup = new HashMap<String, PeptideEvidence>();
        //peptideLookup= new HashMap<String, uk.ac.ebi.jmzidml.model.mzidml.Peptide>();   //lookup to get a peptide by peptideseq_varmods_fixedmods_start_stop
        uniquePeps = new HashMap<org.rappsilber.fdr.entities.Peptide, uk.ac.ebi.jmzidml.model.mzidml.Peptide>();      //lookup to get a peptide by peptideseq_varmods_fixedmods (i.e. uniqueness check)
        uniquePeptidePairs = new HashMap<PeptidePair, Peptide[]>();
        firstFoundPeptidePairs = new HashMap<PeptidePair, PeptidePair>();
        sequenceCollection = new SequenceCollection();

        proteinDetectionList = new ProteinDetectionList();
        HashMap<ProteinGroup, ProteinAmbiguityGroup> pg2Pag = new HashMap<ProteinGroup, ProteinAmbiguityGroup>();
        List<ProteinAmbiguityGroup> pagList = proteinDetectionList.getProteinAmbiguityGroup();
        int pagID =0;
        
        if (OfflineFDR.getXiFDRVersion().toString() != null)
            handleAnalysisSoftware(OfflineFDR.getXiFDRVersion().toString(), (XiInFDR)fdr);
        handleAuditCollection(owner.first, owner.last, owner.email, owner.address, owner.org);
        handleProvider();                //Performed after auditcollection, since contact is needed for provider

        boolean fragmentIsMono = handleAnalysisProtocolCollection(fdr, result, stubs);
        inputs = new Inputs();
//        handleInputs(fdr, "", "", result.proteinGroupFDR.getInputCount(), "");
//        handleAnalysisCollection("");


        // Create fragmentation table
        /*
         <Measure id="m_mz">
         <cvParam cvRef="PSI-MS" accession="MS:1001225" name="product ion m/z"/>
         </Measure>
         <Measure id="m_intensity">
         <cvParam cvRef="PSI-MS" accession="MS:1001226" name="product ion intensity"/>
         </Measure>
         <Measure id="m_error">
         <cvParam cvRef="PSI-MS" accession="MS:1001227" name="product ion m/z error" unitAccession="MS:1000040" unitName="m/z" unitCvRef="PSI-MS"/>
         </Measure>
         */

        Measure mzMeasure = null;
        Measure intMeasure = null;
        Measure errorMeasure = null;



        //List<SpectrumIdentificationResult> specIdentResults = siList.getSpectrumIdentificationResult();
        List<uk.ac.ebi.jmzidml.model.mzidml.Peptide> peptideList = sequenceCollection.getPeptide();

        
        HashMap<Integer,CvParam> xlMass2DonorModParam = new HashMap<Integer, CvParam> (conf.getCrossLinker().size());
        HashMap<String,CvParam> xlMass2DonorModParamByName = new HashMap<String, CvParam> (conf.getCrossLinker().size());
        
        for (rappsilber.ms.crosslinker.CrossLinker cl : conf.getCrossLinker()) {
            if (cl.getName().trim().replaceAll("[- ]*", "").toLowerCase().contains("noncovalent"))
                continue;
            double xlMass = cl.getCrossLinkedMass();
            CvParam XLDonorModParam = new CvParam();
            
            xlspecificityloop:for (rappsilber.ms.sequence.AminoAcid aa : rappsilber.ms.sequence.AminoAcid.getRegisteredAminoAcids()) {
                for (rappsilber.ms.sequence.AminoAcid aa2 : rappsilber.ms.sequence.AminoAcid.getRegisteredAminoAcids()) {
                    rappsilber.ms.sequence.Sequence s = new rappsilber.ms.sequence.Sequence(new rappsilber.ms.sequence.AminoAcid[]{aa,aa,aa,aa,aa,aa,aa,aa,aa});
                    rappsilber.ms.sequence.Peptide p = new rappsilber.ms.sequence.Peptide(s, 2, 5);
                    rappsilber.ms.sequence.Sequence s2 = new rappsilber.ms.sequence.Sequence(new rappsilber.ms.sequence.AminoAcid[]{aa2,aa2,aa2,aa2,aa2,aa2,aa2,aa2,aa2});
                    rappsilber.ms.sequence.Peptide p2 = new rappsilber.ms.sequence.Peptide(s2, 2, 5);
                    if (cl.canCrossLink(p,p2)) {
                        // querry xlmod
                        boolean[] term = new boolean[]{false,false};
                        XLModEntry e = xlmod.guessModification(new XLModQuery(xlMass, cl.getName(), new String[] {aa.toString(),aa2.toString()}, term,term,term,term, true));
                        if (e == null) {
                            e = xlmod.guessModification(new XLModQuery(xlMass, new String[] {aa.toString(),aa2.toString()}, term,term,term,term, false));
                        }

                        if (e != null) {
                            XLDonorModParam.setAccession(e.getId());
                            XLDonorModParam.setCv(xlmodCV);
                            XLDonorModParam.setName(e.getName());
                        } else {
                            uk.ac.liv.unimod.ModT xlMod = unimodDoc.getModByMass(cl.getCrossLinkedMass(), unimodMassError, fragmentIsMono, aa.SequenceID.charAt(0));
                            //Set the found details to modParam. If no unimod record was found, set the modification as "unknown" 
                            if (xlMod != null) {
                                XLDonorModParam.setAccession("UNIMOD:" + xlMod.getRecordId());
                                XLDonorModParam.setCv(unimodCV);
                                XLDonorModParam.setName(xlMod.getTitle());
                            } else {
                                //modification with mass not recognized:
                                XLDonorModParam.setName("unknown modification");
                                XLDonorModParam.setCv(psiCV);
                                XLDonorModParam.setAccession("MS:1001460");
                            }
                        }

                        break xlspecificityloop;
                    }
                }
            }
            xlMass2DonorModParam.put((int)(massNormalise(xlMass)*1000), XLDonorModParam);
            xlMass2DonorModParamByName.put(cl.getName(), XLDonorModParam);
        }
        CvParam XLReceiverModParam = new CvParam();
        XLReceiverModParam.setName("unknown modification");
        XLReceiverModParam.setCv(psiCV);
        XLReceiverModParam.setAccession("MS:1001460");
        
        
        int sirCounter = 1; //Counter used to create unique ID for SpectrumIdentificationResult
        
        HashMap<String, ArrayList<PSM>> allSpectra = new HashMap<String, ArrayList<PSM>>();

        
        
        // join PSMs by scans
        for (PSM ipsm : result.getInput()) {
            for (PSM p : ipsm.getRepresented()) {
                PSM psm = (PSM)p;
                String spectrumID = psm.getScan() + " - " + psm.getRun();
                ArrayList<PSM> scanpsms = allSpectra.get(spectrumID);
                if (scanpsms == null) {
                    scanpsms = new ArrayList<PSM>();
                    allSpectra.put(spectrumID, scanpsms);
                }
                scanpsms.add(psm);
            }
        }
        int xlModId = 0;

        
        
        int proteingrouppairID = 1;
        HashMap<ProteinAmbiguityGroup,HashMap<ProteinAmbiguityGroup,Integer>> protPair2ID = new HashMap<ProteinAmbiguityGroup,HashMap<ProteinAmbiguityGroup,Integer>>();
        HashMap<ProteinGroupPair,CvParam> protein_pair_globa_fdr_term1 = new HashMap<>();
        HashMap<ProteinGroupPair,CvParam> protein_pair_globa_fdr_term2 = new HashMap<>();
        int factor = (int) Math.pow(10,Math.round(Math.log10(fdrResult.getProteinGroupLinkFDR().size()+1)+0.5));
        
        for (ArrayList<PSM> psms : allSpectra.values()) {
            // Get the next spectrum.
//            Spectrum spectrum = iter.next();
            PSM f = psms.get(0);
            int spectrumNumber = Integer.parseInt(f.getScan());//note: spectrum number seems to be a sequential index. For the spectrum number as found in xtandem file use spectrumId
            String run = f.getRun();

            /**
             * ***********************************************
             *  *** Setup SpectrumIdentificationResult ****
             * *********************************************
             */
            SpectrumIdentificationResult specIdentRes = new SpectrumIdentificationResult();
            


            //int pepEvidCounter = 1;
            int siiCounter = 1; //Counter used to create unique ID for SpectrumIdentificationItem

            SpectrumIdentificationItem sii = null;
            //uk.ac.ebi.jmzidml.model.mzidml.Peptide mzidPep = null;

            List<IonType> ionTypeList = null;

            HashMap<String, SpectrumIdentificationItem> sIIMap = new HashMap<String, SpectrumIdentificationItem>();
            HashSet<org.rappsilber.fdr.entities.Peptide> allPeps = new HashSet<org.rappsilber.fdr.entities.Peptide>();
            HashMap<String,SpectraData> runData = new HashMap<String,SpectraData>();
            boolean firstSpecPSM = true;
            for (PSM psm : psms) {
                boolean passed = result.getPsmFDR().filteredContains(psm) || (psm.getPartOfUniquePSM() != null && psm.getPartOfUniquePSM() == psm && result.getPsmFDR().filteredContains(psm.getPartOfUniquePSM()));
                xlModId++;
                org.rappsilber.fdr.entities.PeptidePair peppair = psm.getPeptidePair();
                
                SpectraData specData =  getSpectraData(psm.getSearchID(), psm.getPeakListName());

                double calcMass = psm.getCalcMass();
                double psmXLMass = psm.getCalcMass();
                for (org.rappsilber.fdr.entities.Peptide pep : peppair.getPeptides()) {
                    psmXLMass-=pep.mass;
                }
                String xl_name = psm.getCrosslinker();
                rappsilber.ms.crosslinker.CrossLinker selected_xl = name2xiXL.get(xl_name);
                if (selected_xl == null) {
                    if (conf.getCrossLinker().size() == 1) {
                        selected_xl = conf.getCrossLinker().get(0);
                    } else {
                        double massDiff = Double.MAX_VALUE;

                        for (rappsilber.ms.crosslinker.CrossLinker xl : conf.getCrossLinker()) {
                            if (xl.getName().contentEquals(xl_name)) {
                                selected_xl = xl;
                                break;
                            } else {
                                double td = Math.abs(xl.getCrossLinkedMass() - psmXLMass);
                                if (td < massDiff){
                                    massDiff = td;
                                    selected_xl = xl;
                                }
                            }
                        }
                    }
                    name2xiXL.put(xl_name, selected_xl);

                }
                if (selected_xl != null)
                    psmXLMass = massNormalise(selected_xl.getCrossLinkedMass());
                //In mzIdentML we have 1 SepctrumIdentificationItem (SII) linked to 1 Peptide item
                //via the peptide_ref attribute. Each Peptide item is a unique combination of:
                // peptidesequence + modifications + substitutionModifications.
                //String uniquePepKey = getPeptideKey(domain, iXTandemFile);
                int pi=-1;
                PeptidePair pp = psm.getPeptidePair();
                ArrayList<org.rappsilber.fdr.entities.Peptide> pp_peptides = pp.getPeptides();
                
                // did we already get a entry for the peptide-pair
                Peptide[] psmpeps = uniquePeptidePairs.get(peppair);
                // and if so we also want to know the orginial peptide pair
                PeptidePair fpp  = firstFoundPeptidePairs.get(peppair);
                
                String peppairkey;
                boolean pepsNew = false;
                if (psmpeps == null) {
                    pepsNew = true;
                    // first time that we found the peptide pair
                    // so we need to generate new mzIdenML-peptides
                    
                    if (pp.isLinear()) {
                        // linear match -> only one peptide
                        psmpeps = new Peptide[1];
                        peppairkey = pp.getPeptide1().getId() + "_" + pp.getPeptide1().getSequence();
                        
                    } else {
                        // crosslink -> two peptides
                        psmpeps = new Peptide[2];
                        peppairkey = pp.getPeptide1().getId() +"_"+pp.getPeptide1().getSequence() + "_" + pp.getPeptide2().getId()  +"_"+pp.getPeptide2().getSequence()+
                                 "_" + pp.getPeptideLinkSite1() + "_" + pp.getPeptideLinkSite2() + "_" + psm.getCrosslinker();
                    }
                    //psmpeps = new Peptide[pp_peptides.size()];
                    
                    for (pi = pp_peptides.size() - 1; pi>=0; pi --) {
                        // get the peptide
                        org.rappsilber.fdr.entities.Peptide pep = pp_peptides.get(pi);

                        int link = pp.getPeptideLinkSite(pi)-1;

                        // new mz identml peptide
                        uk.ac.ebi.jmzidml.model.mzidml.Peptide mzidPep = new uk.ac.ebi.jmzidml.model.mzidml.Peptide();
                        mzidPep.setPeptideSequence(pep.getUppercaseSequence());
                        // unique key for the peptide
                        // it's the peptide-pair  key + wether it is the first or the second peptide
                        String uniquePepKey = peppairkey + "_p" + pi;
                        mzidPep.setId(uniquePepKey);

                        //Parse the modifications and add them to the mzidPep object:
                        parseModificationsAndSubstitutions(mzidPep, pep, fdr, fragmentIsMono);

                        
                        if (!uniquePeps.containsKey(uniquePepKey)) {
                            peptideList.add(mzidPep);
                            uniquePeps.put(pep, mzidPep);
                        }
                        
                        
                        // the first peptide get the crosslinker as modification added
                        if (pi == 0 ) {
                            if (!psm.isLinear()) {
                                if (!psm.getCrosslinker().replaceAll("[-\\s]","").toLowerCase().contains("noncovalent")) {
                                    CvParam  XLDonorModParam =  xlMass2DonorModParam.get((int)(psmXLMass*1000));
                                    if (XLDonorModParam == null)
                                        XLDonorModParam =  xlMass2DonorModParamByName.get(xl_name);
                                    uk.ac.ebi.jmzidml.model.mzidml.Modification mod = getCrosslinkerDonorModification(XLDonorModParam, link, psmXLMass, fragmentIsMono);
                                    CvParam xlModParam = new CvParam();
                                    xlModParam.setAccession(crosslinkedDonorModAcc);
                                    xlModParam.setName(crosslinkedDonorModName);
                                    xlModParam.setCv(psiCV);
                                    xlModParam.setValue(Integer.toString(xlModId));
                                    mod.getCvParam().add(xlModParam);
                                    mod.getCvParam().addAll(stubs);
                                    mzidPep.getModification().add(mod);
                                }
                            } else if (psm.isLoop()) {
                                CvParam  XLDonorModParam =  xlMass2DonorModParam.get((int)(psmXLMass*1000));
                                if (XLDonorModParam == null)
                                    XLDonorModParam =  xlMass2DonorModParamByName.get(xl_name);
                                uk.ac.ebi.jmzidml.model.mzidml.Modification mod = getCrosslinkerDonorModification(XLDonorModParam, link, psmXLMass, fragmentIsMono);
                                CvParam xlModParam = new CvParam();
                                xlModParam.setAccession(crosslinkedDonorModAcc);
                                xlModParam.setName(crosslinkedDonorModName);
                                xlModParam.setCv(psiCV);
                                xlModParam.setValue(Integer.toString(xlModId));
                                mod.getCvParam().add(xlModParam);
                                mod.getCvParam().addAll(stubs);
                                mzidPep.getModification().add(mod);
                                link = pp.getPeptideLinkSite(pi+1);
                                mod = getCrosslinkerReceptorModification(link, 0, fragmentIsMono);
                                xlModParam = new CvParam();
                                xlModParam.setAccession(getCrosslinkedAcceptorModAcc());
                                xlModParam.setName(getCrosslinkedAcceptorModName());
                                xlModParam.setCv(psiCV);
                                xlModParam.setValue(Integer.toString(xlModId));
                                mod.getCvParam().add(xlModParam);
                                mod.getCvParam().addAll(stubs);
                                mod.getCvParam().addAll(stubs);
                                mzidPep.getModification().add(mod);
                            }
                        } else {
                            if (!psm.getCrosslinker().replaceAll("[-\\s]","").toLowerCase().contains("noncovalent")) {
                                // the second peptide becomes a zero-mass modification added, that denote the linkage site on this peptide
                                uk.ac.ebi.jmzidml.model.mzidml.Modification mod = getCrosslinkerReceptorModification(link, 0, fragmentIsMono);
                                CvParam xlModParam = new CvParam();
                                xlModParam.setAccession(getCrosslinkedAcceptorModAcc());
                                xlModParam.setName(getCrosslinkedAcceptorModName());
                                xlModParam.setCv(psiCV);
                                xlModParam.setValue(Integer.toString(xlModId));
                                mod.getCvParam().add(xlModParam);
                                mod.getCvParam().addAll(stubs);
                                mzidPep.getModification().add(mod);
                            }
                        }
                        psmpeps[pi] = mzidPep;
                    }
                    fpp = pp;
                    uniquePeptidePairs.put(pp, psmpeps);
                    firstFoundPeptidePairs.put(fpp, fpp);
                    
                    
                    
                } else {
                    pp_peptides = fpp.getPeptides();
                }
                    
                
                StringBuffer pepkey=new StringBuffer(psmpeps[0].getId());
                for (int pio =  1; pio  < psmpeps.length; pio ++) {    
                    pepkey.append("_").append(psmpeps[pio].getId());
                }
                for (pi = psmpeps.length - 1; pi>=0; pi --) {
                    //pi++;
                    Peptide mzidPep = psmpeps[pi];
                    //If it is a new global peptide, then initialize a new mzIdentML Peptide object: 
                    

                    String sIIKey = psm.getPsmID() + "_" + pi;
                    org.rappsilber.fdr.entities.Peptide pep = pp_peptides.get(pi);

                    /*
                    ****************************************************
                    ****** Create new SpectrumIdentificationItem *******
                    ****************************************************
                    */
                    sii = new SpectrumIdentificationItem();

                    sii.setPassThreshold(passed);

                    //add sii to sir:
                    specIdentRes.getSpectrumIdentificationItem().add(sii);
                    sIIMap.put(sIIKey, sii);


                    sii.setPeptide(mzidPep);
                    if (psmpeps.length>1) {
                        if (!psm.getCrosslinker().replaceAll("[-\\s]","").toLowerCase().contains("noncovalent")) {
                            CvParam xlModParam = new CvParam();
                            xlModParam.setAccession(crosslinkedSIIAcc);
                            xlModParam.setValue(Integer.toString(xlModId));
                            xlModParam.setName(crosslinkedSIIName);
                            xlModParam.setCv(psiCV);
                            sii.getCvParam().add(xlModParam);
                        } else {
                            CvParam xlModParam = new CvParam();
                            xlModParam.setAccession(noncovalentSIIAcc);
                            xlModParam.setValue(Integer.toString(xlModId));
                            xlModParam.setName(noncovalentSIIName);
                            xlModParam.setCv(psiCV);
                            sii.getCvParam().add(xlModParam);
                        }
                    }

                    //parseScoresAndOtherSIIAttributes(sii, domain, spectrum);
                    parseScoresAndOtherSIIAttributes(sii, pep.getProteinGroup(), psm);

                    sii.setId("SII_" + sirCounter + "_" + siiCounter);
                    siiCounter++;

                    ProteinGroup pg = pep.getProteinGroup();
                    ProteinAmbiguityGroup pag = pg2Pag.get(pg);
                    List<ProteinDetectionHypothesis> pdhList;

                    boolean groupIsNew = false;
                    if (pag == null) {
                        pag = new ProteinAmbiguityGroup();
                        pag.setId("PAG_"+pagID++);
                        pagList.add(pag);
                        pg2Pag.put(pg, pag);
                        groupIsNew = true;
                        CvParam cvp = makeCvParam("MS:1002415", "protein group passes threshold", psiCV,""+fdrResult.getProteinGroupFDR().filteredContains(pg));
                        pag.getCvParam().add(cvp);
                     
                    }
                    pdhList = pag.getProteinDetectionHypothesis();                        
                    
                   //Parse protein details into DBSequence objects:
                    HashMap<org.rappsilber.fdr.entities.Protein, HashSet<Integer>> protPositions = pep.getPositions();
                    for (org.rappsilber.fdr.entities.Protein prot : protPositions.keySet())  {
                        String protKey =  (prot.isDecoy() ? "decoy" : "") + prot.getAccession();
                        String dbSeqId = "dbseq_" + protKey;
                        DBSequence dbSeq = foundProts.get(protKey);
                        if (dbSeq == null) {
                            dbSeq = new DBSequence();
                            foundProts.put(protKey, dbSeq);
                            dbSeq.setAccession(prot.getAccession());
                            if (prot.isDecoy()) {
                                dbSeq.setName("decoy");
                            } else {
                                dbSeq.setName(prot.getName() == null || prot.getName().trim().length() == 0 ? prot.getAccession() : prot.getName());
                            }
                            if (prot.getDescription() != null  && prot.getDescription().trim().length()>0)
                                dbSeq.getCvParam().add(makeCvParam("MS:1001088", "protein description", psiCV,prot.getDescription()));
                            if (prot.getSize() >0)
                                dbSeq.setLength(prot.getSize());
                            if (prot.getSequence() != null && !prot.getSequence().isEmpty())
                                dbSeq.setSeq(prot.getSequence());
                            //dbSeq.setSeq(protSeq);
//                            dbSeq.setLength(prot.length());
                            dbSeq.setId(dbSeqId);
                            dbSeq.setSearchDatabase(getDatabases(psm.getSearchID()).get(0));
                            //dbSeq.setSearchDatabase(searchDB);
                            sequenceCollection.getDBSequence().add(dbSeq);                        
                            foundProts.put(protKey, dbSeq);
                            protein2DBSequence.put(prot, dbSeq);
                        }
                        List<PeptideHypothesis> pepHList = null;
                        // if this is the first time we encounter these peptides we have to register the support for ProteinAmbiguityGroups
                        ProteinDetectionHypothesis pdh = null;
                        if (groupIsNew) {
                            // if this is a new protein group then add all the proteins to it;
                            pdh = new ProteinDetectionHypothesis();
                            pdh.setId(pag.getId()+"_PDH_"+pdhList.size());
                            pdhList.add(pdh);
                            pdh.setDBSequence(dbSeq);
                            pepHList = pdh.getPeptideHypothesis();
                            if (pdhList.size() == 1) {
                                pdh.getCvParam().add(makeCvParam("MS:1002403", "group representative", psiCV));
                            }
                            pdh.getCvParam().add(makeCvParam("MS:1001593", "group member with undefined relationship OR ortholog protein", psiCV));
                            
                        } else {
                            for (ProteinDetectionHypothesis tpdh : pdhList) {
                                if (tpdh.getDBSequence() == dbSeq) {
                                    pepHList = tpdh.getPeptideHypothesis();
                                    pdh=tpdh;
                                    break;
                                }
                            }
                            if (pepHList==null) {
                                Logger.getLogger(this.getClass().getName()).log(Level.SEVERE,"PeptideHypothesis is missing");
                            }
                        }
                        
                        for (int pos : protPositions.get(prot)) {
                            String pepEvidKey = "pepevid_pep_" + pep.getId() +"_" + pep.getSequence() + "_prot" + protKey + "_" + pos + "_pepkey_" + pi+ "_"+pepkey.toString();
                            PeptideEvidence pepEvid = pepEvidLookup.get(pepEvidKey);
                            if (pepEvid == null) {
                                pepEvid = new PeptideEvidence();
                                pepEvid.setDBSequence(dbSeq);
                                pepEvid.setStart(pos);
                                pepEvid.setEnd(pos+pep.length()-1);
                                pepEvid.setIsDecoy(pep.isDecoy());
                                pepEvid.setPeptide(mzidPep);
                                pepEvid.setId(pepEvidKey);
                                pepEvidLookup.put(pepEvidKey, pepEvid);
                                //link it to our sequenceCollection list:
                                sequenceCollection.getPeptideEvidence().add(pepEvid);
                                PeptideHypothesis pepH = new PeptideHypothesis();
                                pepH.setPeptideEvidence(pepEvid);
                                pepHList.add(pepH);
                                SpectrumIdentificationItemRef siir = new SpectrumIdentificationItemRef();
                                siir.setSpectrumIdentificationItem(sii);
                                pepH.getSpectrumIdentificationItemRef().add(siir);
                            } else {
                                // find the peptide hypothesis for this peptide
                                for (PeptideHypothesis pepH : pepHList) {
                                    if (pepH.getPeptideEvidence().equals(pepEvid) ) {
                                        SpectrumIdentificationItemRef siir = new SpectrumIdentificationItemRef();
                                        siir.setSpectrumIdentificationItem(sii);
                                        pepH.getSpectrumIdentificationItemRef().add(siir);
                                        break;
                                    }
                                }
                            }
                            PeptideEvidenceRef peptideEvidenceRef = new PeptideEvidenceRef();
                            peptideEvidenceRef.setPeptideEvidence(pepEvid);
                            sii.getPeptideEvidenceRef().add(peptideEvidenceRef);

                        }
                    }
                }

                if (firstSpecPSM) {
                    specIdentRes.setSpectrumID(this.scanIDTranslation.getID(f));
                    specIdentRes.setId("SIR_" + sirCounter);
                    specIdentRes.setSpectraData(specData);

                    specIdentRes.getCvParam().add(makeCvParam("MS:1000797", "peak list scans", psiCV,psm.getScan()+""));
                    firstSpecPSM = false;
                }
                

                {
                    
                    ProteinAmbiguityGroup pag1 = pg2Pag.get(psm.getProteinGroup1());
                    if (pag1==null) {
                        try {
                            //  I don't know why!?!
                            Thread.sleep(10);
                        } catch (InterruptedException ex) {
                            Logger.getLogger(MZIdentMLExport.class.getName()).log(Level.SEVERE, null, ex);
                        }
                        pag1 = pg2Pag.get(psm.getProteinGroup1());
                    }
                    ProteinAmbiguityGroup pag2 = null;
                    if (!psm.isLinear()) {
                        pag2=pg2Pag.get(psm.getProteinGroup2());
                        if (pag2==null) {
                            try {
                                //  I don't know why!?!
                                Thread.sleep(10);
                            } catch (InterruptedException ex) {
                                Logger.getLogger(MZIdentMLExport.class.getName()).log(Level.SEVERE, null, ex);
                            }
                            pag2 = pg2Pag.get(psm.getProteinGroup2());
                        }
                    }
                    HashMap<ProteinAmbiguityGroup,Integer> prot22id = protPair2ID.get(pag1);
                    if (prot22id == null) {
                        prot22id = new HashMap<ProteinAmbiguityGroup, Integer>();
                        protPair2ID.put(pag1, prot22id);
                    }
                    
                    Integer currentppiID = prot22id.get(pag2);
                    
                    if (currentppiID == null) {
                        currentppiID = proteingrouppairID;
                        prot22id.put(pag2, currentppiID);
                        proteingrouppairID++;
                    }


                    if (!pag1.equals(pag2)) {
                        prot22id = protPair2ID.get(pag2);
                        if (prot22id == null) {
                            prot22id = new HashMap<ProteinAmbiguityGroup, Integer>();
                            protPair2ID.put(pag2, prot22id);
                        }
                        prot22id.put(pag1, currentppiID);
                    }
                    CvParam cvp1 = null;
                    CvParam cvp2 = null;
                    if (pag2 != null) {
                        ProteinGroupPair pgp = psm.getLinks().iterator().next().getProteinGroupPair();
                        ProteinGroupPair pgp_final =  result.getProteinGroupPairFDR().filteredGet(pgp);
                        double pgp_fdr = (pgp_final == null?pgp.getFDR(): pgp_final.getFDR());
                        if (passed) {
                            Iterator<ProteinGroupLink> linkiter =  psm.getLinks().iterator();                                    
                            CvParam prev_cvp1 = protein_pair_globa_fdr_term1.get(pgp);
                            if (prev_cvp1 != null) {
                                prev_cvp1.setValue((currentppiID*factor)+".a:null:"+pgp_fdr+":"+passed);
                                protein_pair_globa_fdr_term2.get(pgp).setValue((currentppiID*factor)+".b:null:"+pgp_fdr+":"+passed);
                            } else {
                                cvp1 = makeCvParam("MS:1002676", "protein-pair-level global FDR", psiCV,(currentppiID*factor)+".a:null:"+pgp_fdr+":"+passed);
                                cvp2 = makeCvParam("MS:1002676", "protein-pair-level global FDR", psiCV,(currentppiID*factor)+".b:null:"+pgp_fdr+":"+passed);
                                protein_pair_globa_fdr_term1.put(pgp,cvp1);
                                protein_pair_globa_fdr_term2.put(pgp,cvp2);
                            }
                        } else {
                            if (protein_pair_globa_fdr_term1.get(pgp)== null) {
                                if (pgp.getFDR() != -1) {
                                    cvp1 = makeCvParam("MS:1002676", "protein-pair-level global FDR", psiCV,(currentppiID*factor)+".a:null:"+pgp_fdr+":"+passed);
                                    cvp2 = makeCvParam("MS:1002676", "protein-pair-level global FDR", psiCV,(currentppiID*factor)+".b:null:"+pgp_fdr+":"+passed);
                                    protein_pair_globa_fdr_term1.put(pgp,cvp1);
                                    protein_pair_globa_fdr_term2.put(pgp,cvp2);
                                }
                            }
                        }
                    }
                    if (cvp1 != null) {
                        for (ProteinDetectionHypothesis pdh : pag1.getProteinDetectionHypothesis()) {
                            pdh.getCvParam().add(cvp1);
                        }
                    }
                    if (cvp2 != null) {
                        for (ProteinDetectionHypothesis pdh : pag2.getProteinDetectionHypothesis()) {
                            pdh.getCvParam().add(cvp2);
                        }
                    }
                }
            }
            SpectrumIdentificationList sil =  getSpectrumIdentificationList(psms.get(0));
            
            sil.getSpectrumIdentificationResult().add(specIdentRes);
            
            
            sirCounter++;

            
            
            
        }

        
        int residuepairID = 1;
        for (ProteinGroupLink pgl : fdrResult.getProteinGroupLinkFDR()) {
            boolean passed = false;
            if (fdrResult.getProteinGroupLinkFDR().filteredContains(pgl)) {
                passed = true;
            }

            ProteinAmbiguityGroup pag1 = pg2Pag.get(pgl.getProteinGroup1());
            ProteinAmbiguityGroup pag2 = pg2Pag.get(pgl.getProteinGroup2());
            Integer ppid = null;
            HashMap<ProteinAmbiguityGroup,Integer> prot22id  = protPair2ID.get(pag1);
            if (prot22id == null) {
                prot22id = protPair2ID.get(pag2);
                ppid = prot22id.get(pag1);
            } else {
                ppid = prot22id.get(pag2);
            }
            
            int linkid = (int) (ppid*factor+residuepairID);
            
            
            for (org.rappsilber.fdr.entities.Protein prot : pgl.getPosition1().keySet()) {
                for (ProteinDetectionHypothesis pdh : pag1.getProteinDetectionHypothesis()) {
                    if (protein2DBSequence.get(prot) == pdh.getDBSequence()) {
                        for (int pos :  pgl.getPosition1().get(prot)) {
                            CvParam cvp = makeCvParam("MS:1002677", "residue-pair-level global FDR", psiCV,linkid+".a:"+pos+":"+pgl.getFDR()+":"+passed);
                            pdh.getCvParam().add(cvp);
                        }
                        break;
                    }
                }
            }
            for (org.rappsilber.fdr.entities.Protein prot : pgl.getPosition2().keySet()) {
                for (ProteinDetectionHypothesis pdh : pag2.getProteinDetectionHypothesis()) {
                    if (protein2DBSequence.get(prot) == pdh.getDBSequence()) {
                        for (int pos :  pgl.getPosition2().get(prot)) {
                            CvParam cvp = makeCvParam("MS:1002677", "residue-pair-level global FDR", psiCV,linkid+".b:"+pos+":"+pgl.getFDR()+":"+passed);
                            pdh.getCvParam().add(cvp);
                        }
                        break;
                    }
                }
            }
            residuepairID++;
        }        
        
    }

    private ArrayList<CvParam> parseStubs(RunConfig conf) {
        // look up stubs
        ArrayList<CvParam> stubs = new ArrayList();
        int ast = 0;
        CvParam stubA = null;
        CvParam stubS = null;
        CvParam stubT = null;
        for (rappsilber.ms.sequence.ions.CrossLinkedFragmentProducer cfp: conf.getCrossLinkedFragmentProducers()) {
            if (cfp instanceof rappsilber.ms.sequence.ions.loss.CleavableCrossLinkerPeptide) {
                rappsilber.ms.sequence.ions.loss.CleavableCrossLinkerPeptide p = (rappsilber.ms.sequence.ions.loss.CleavableCrossLinkerPeptide) cfp;
                String name = p.getName().substring(0,1);
                Double mass = p.getDeltamass();
                CvParam cvp = new CvParam();
                cvp.setAccession(crosslinkerStubAcc);
                cvp.setCv(psiCV);
                cvp.setName(crosslinkerStubName);
                cvp.setValue(name +":" + mass + ":");
                if (name.toLowerCase().contentEquals("a")) {
                    ast++;
                    stubA = cvp;
                } else if (name.toLowerCase().contentEquals("s")) {
                    ast++;
                    stubS = cvp;
                } else if  (name.toLowerCase().contentEquals("t")) {
                    ast++;
                    stubT = cvp;
                }
                stubs.add(cvp);
            }
        }
        // make some eductated guesses for stub pairing
        if (ast == 3) {
            char aName =  stubA.getValue().charAt(0);
            char sName =  stubS.getValue().charAt(0);
            char tName =  stubT.getValue().charAt(0);
            stubA.setValue(stubA.getValue()+sName + tName);
            stubT.setValue(stubT.getValue()+aName);
            stubS.setValue(stubS.getValue()+aName);
        } else if (stubs.size() == 2) {
            String s0 = stubs.get(0).getValue();
            String s1 = stubs.get(1).getValue();
            String n0 = s0.substring(0,1);
            String n1 = s1.substring(0,1);
            stubs.get(0).setValue(s0+n1);
            stubs.get(1).setValue(s1+n0);
        }
        return stubs;
    }
    

    /**
     * DOCUMENT ME!
     *
     * @param sii
     * @param domain
     * @param spectrum
     */
    private void parseScoresAndOtherSIIAttributes(SpectrumIdentificationItem sii, ProteinGroup pg, PSM psm) {
        double pgscore = pg.getScore();

        int precursorCharge = psm.getCharge();

        //Get precursorMh reported by X!Tandem and convert it back
        //to m/z using:
        //     m/z=(mh + z*1.007276 - 1*1.007276)/z
        //where mh is M+H    and H=1.007276466812 (proton mass rounded from 1.007276466812)
        //  X!Tandem is using H=1.007276 (proton mass rounded from 1.007276466812) so we will use that
//        double expMZ = (spectrum.getPrecursorMh() + precursorCharge * protonMass - protonMass) / precursorCharge;
        //round at 6 decimals:
//        expMZ = Utils.round(expMZ, 6);
        sii.setExperimentalMassToCharge(psm.getExperimentalMZ());
//
//        double calcMZ = (domain.getDomainMh() + precursorCharge * protonMass - protonMass) / precursorCharge;
//        calcMZ = Utils.round(calcMZ, 6);
        if (psm.getCalcMass()>0) {
            sii.setCalculatedMassToCharge(psm.getCalcMass()/psm.getCharge()+Util.PROTON_MASS);
        }
        sii.setChargeState(precursorCharge);
        sii.setRank(psm.getRank());

        List<CvParam> cvParamList = sii.getCvParam();
        cvParamList.add(makeCvParam("MS:1002545", "xi:score", psiCV, "" + psm.getScore()));
    }


    /**
     * Parses the modifications and substitutions storing them in the respective
     * mzidPep modifications and substitutions list.
     *
     * @param mzidPep
     * @param domain
     * @param iXTandemFile
     * @param fragmentIsMono : whether fragment masses are based on a single
     * isotope (12C) mass or on the weighted average of all possible
     * 13C-containing molecular masses. Most modern mass spectrometers have
     * sufficient mass resolution to measure the 12C mass separately and thus
     * use this to report the fragment masses.
     *
     */
    private void parseModificationsAndSubstitutions(Peptide mzidPep, org.rappsilber.fdr.entities.Peptide fdrPep, OfflineFDR fdr, boolean fragmentIsMono) {

        //Parse the modifications
        rappsilber.config.RunConfig conf = ((XiInFDR)fdr).getConfig();
        rappsilber.ms.sequence.Sequence seq = null;
        if (conf instanceof Xi2Xi1Config && ((Xi2Xi1Config) conf).isModX) {
            seq = new rappsilber.ms.sequence.Sequence(fdrPep.getSequence(), conf, true);
        } else {
            seq = new rappsilber.ms.sequence.Sequence(fdrPep.getSequence(), conf);
        }
        rappsilber.ms.sequence.AminoAcid[] aaseq = seq.toArray();

        for (int a = 0; a < aaseq.length; a++) {
            if (aaseq[a] instanceof rappsilber.ms.sequence.AminoModification) {
                if (!aaseq[a].SequenceID.contains("X")) {
                    rappsilber.ms.sequence.AminoModification am = (rappsilber.ms.sequence.AminoModification) aaseq[a];
                    uk.ac.ebi.jmzidml.model.mzidml.Modification mzidMod = translateToMzidModification(am, a, seq, fragmentIsMono);
                    mzidPep.getModification().add(mzidMod);
                }
            }
        }


    }

//    /**
//     * Parses and translates the X!tandem reported modification to the
//     * uk.ac.ebi.jmzidml.model.mzidml.SubstitutionModification object format.
//     *
//     * @param reportedMod: X!tandem parser modification object ( for with
//     * reportedMod.isSubstitution() is true)
//     * @param domain : the domain X!Tandem parser object in which the
//     * modification was reported
//     * @return
//     */
//    private SubstitutionModification translateToMzidSubstitution(Modification reportedMod, Domain domain, boolean fragmentIsMono) {
//        SubstitutionModification mzidSubs = new SubstitutionModification();
//
//        int loc = Integer.parseInt(reportedMod.getLocation());
//
//        int pepLoc = loc - domain.getDomainStart(); //location in Tandem is given as location within the whole protein
//        mzidSubs.setLocation(pepLoc + 1);
//
//        //Note: in X!Tandem, like in mzIdentML: sequence reported is the original:
//        char reportedAminoAcid = domain.getDomainSequence().charAt(pepLoc);
//        mzidSubs.setOriginalResidue(reportedAminoAcid + "");
//
//        mzidSubs.setReplacementResidue(reportedMod.getSubstitutedAminoAcid());
//
//
//        //TODO - the current storing of the modification mass in either
//        //setMonoisotopicMassDelta or setAvgMassDelta seems wrong as the massDelta is 
//        //a delta on the parent ion mass which is always monoisotopic in xtandem. On the other hand, 
//        //a modification has to have one or more matching fragments as well which all are 
//        //shifted by x, x being the modification mass.
//        //So check with X!Tandem developer whether the modification (and substitution) masses 
//        //follow the fragment mass type or the parent mass type.
//
//        double mass = reportedMod.getMass();
//
//        if (fragmentIsMono) {
//            mzidSubs.setMonoisotopicMassDelta(mass);
//        } else {
//            mzidSubs.setAvgMassDelta(mass);
//        }
//
//        return mzidSubs;
//
//    }

    /**
     * Parses and translates the X!tandem reported modification to the
     * uk.ac.ebi.jmzidml.model.mzidml.Modification object format.
     *
     * @param reportedMod : X!tandem parser modification object
     * @param domain : the domain X!Tandem parser object in which the
     * modification was reported
     * @param fragmentIsMono :
     * @return
     */
    private uk.ac.ebi.jmzidml.model.mzidml.Modification translateToMzidModification(rappsilber.ms.sequence.AminoModification reportedMod, int location,
            rappsilber.ms.sequence.Sequence seq, boolean fragmentIsMono) {
        uk.ac.ebi.jmzidml.model.mzidml.Modification mzidMod = new uk.ac.ebi.jmzidml.model.mzidml.Modification();
        double mass = massNormalise(reportedMod.weightDiff);
        int loc = location;

        CvParam modParam = null;

        if (fragmentIsMono) {
            mzidMod.setMonoisotopicMassDelta(mass);
        } else {
            mzidMod.setAvgMassDelta(mass);
        }

        int pepLoc = location; //location in Tandem is given as location within the whole protein
        char reportedAminoAcid = reportedMod.BaseAminoAcid.SequenceID.charAt(0);
        mzidMod.setLocation(pepLoc + 1);        //mzid starts counting from 1, except for NTerm/CTerm mods which are 0 
        List<String> residueList = mzidMod.getResidues();
        residueList.add("" + reportedAminoAcid);

        boolean[] term =new boolean[]{false};
        XLModEntry e = xlmod.guessModificationCached(new XLModQuery(mass, reportedMod.SequenceID.substring(1) , new String[]{reportedMod.BaseAminoAcid.SequenceID}, term, term, term, term, false));
        if  (e == null) {
            e = xlmod.guessModificationCached(new XLModQuery(mass, new String[]{reportedMod.BaseAminoAcid.SequenceID}, term, term, term, term, false));
        }

        if  (e == null) {
            //check if we can find a known modification type that matches what is reported in the given mass and aminoAcid, within
            //a given mass error tolerance: 
            uk.ac.liv.unimod.ModT unimod = unimodDoc.getModByMass(mass, unimodMassError, fragmentIsMono, reportedAminoAcid);

            //Part below is needed because it is not clear in XTandem output what are N or C-term mods, 
            //these have the same location as mods on the first aa or last aa in the peptide.
            //So if no unimod was found above, perhaps this is because it is really a N or C-term modification. 
            //Below we try to see if it fits a known N or C-term modification:
            if (unimod == null && pepLoc == 0) {
                //See if this is a possible N-terminal mod
                unimod = unimodDoc.getModByMass(mass, unimodMassError, fragmentIsMono, '[');
                mzidMod.setLocation(0);
            }
            if (unimod == null && loc == seq.length() - 1) {
                //See if this is a possible C-terminal mod
                unimod = unimodDoc.getModByMass(mass, unimodMassError, fragmentIsMono, ']');
                mzidMod.setLocation(0); //also 0?
            }

            //Set the found details to modParam. If no unimod record was found, set the modification as "unknown" 
            if (unimod != null) {
                modParam = makeCvParam("UNIMOD:" + unimod.getRecordId(), unimod.getTitle(), unimodCV);
            } else {
                //modification with mass not recognized:
                modParam = makeCvParam("MS:1001460", "unknown modification", psiCV);
            }
        } else {
            modParam = makeCvParam(e.getId(), e.getName(), xlmodCV);
        }
    
        if (modParam != null)
            mzidMod.getCvParam().add(modParam);
        return mzidMod;
    }

    
   /**
     * Parses and translates the X!tandem reported modification to the
     * uk.ac.ebi.jmzidml.model.mzidml.Modification object format.
     *
     * @param reportedMod : X!tandem parser modification object
     * @param domain : the domain X!Tandem parser object in which the
     * modification was reported
     * @param fragmentIsMono :
     * @return
     */
    private uk.ac.ebi.jmzidml.model.mzidml.Modification getCrosslinkerDonorModification(CvParam xlUnimod, int location, double xlMass, boolean fragmentIsMono) {
        uk.ac.ebi.jmzidml.model.mzidml.Modification mzidMod = new uk.ac.ebi.jmzidml.model.mzidml.Modification();
        double mass = xlMass;
        int loc = location;


        if (fragmentIsMono) {
            mzidMod.setMonoisotopicMassDelta(mass);
        } else {
            mzidMod.setAvgMassDelta(mass);
        }

        int pepLoc = location; //location in Tandem is given as location within the whole protein
        mzidMod.setLocation(pepLoc + 1);        //mzid starts counting from 1, except for NTerm/CTerm mods which are 0 

        mzidMod.getCvParam().add(xlUnimod);
        return mzidMod;
    }    

   /**
     * Parses and translates the X!tandem reported modification to the
     * uk.ac.ebi.jmzidml.model.mzidml.Modification object format.
     *
     * @param reportedMod : X!tandem parser modification object
     * @param domain : the domain X!Tandem parser object in which the
     * modification was reported
     * @param fragmentIsMono :
     * @return
     */
    private uk.ac.ebi.jmzidml.model.mzidml.Modification getCrosslinkerReceptorModification(int location, double xlMass, boolean fragmentIsMono) {
        uk.ac.ebi.jmzidml.model.mzidml.Modification mzidMod = new uk.ac.ebi.jmzidml.model.mzidml.Modification();
        double mass = xlMass;
        int loc = location;


        if (fragmentIsMono) {
            mzidMod.setMonoisotopicMassDelta(mass);
        } else {
            mzidMod.setAvgMassDelta(mass);
        }

        int pepLoc = location; 
        mzidMod.setLocation(pepLoc + 1);        //mzid starts counting from 1, except for NTerm/CTerm mods which are 0 

//        mzidMod.getCvParam().add(xlUnimod);
        return mzidMod;
    }    
    
//    /**
//     * This method will check if the given X!tandem domain object corresponds to
//     * a new SpectrumIdentificationItem (SII) in mzIdentML
//     *
//     * @param sIIKey : the siiKey
//     * @param siiMap : hashmap containing the sii found until now with key =
//     * peptidesequence + modifications + substitutionModifications
//     *
//     * @return returns true in case this is really a new
//     * SepctrumIdentificationItem, false otherwise (i.e. it is just a new
//     * PeptideEvidence for an existing SII)
//     */
//    private boolean isNewSII(String sIIKey, HashMap<String, SpectrumIdentificationItem> sIIMap) {
//        if (sIIMap.get(sIIKey) != null) {
//            return false;
//        } else {
//            return true;
//        }
//    }

//    /**
//     * In mzIdentML, the same Peptide object which is a combination of
//     * peptidesequence + modifications + substitutionModifications can be
//     * matched to multiple SeptrumIdentificationItem objects if these are in
//     * different SepctrumIdentificationResult themselves.
//     *
//     * @param domain
//     * @param iXTandemFile
//     * @return key = peptidesequence + modifications + substitutionModifications
//     */
//    private String getPeptideKey(Domain domain, XTandemFile iXTandemFile) {
//        //is really the same as in getSIIKey, but the context of both maps is different (peptide map is global and siimap is local within 
//        //a SepctrumIdentificationResult ):
//        return getSIIKey(domain, iXTandemFile);
//    }

//    /**
//     * This method returns the unique identifier that maps a X!tandem domain
//     * object to a mzIdentML SpectrumIdentificationItem(SII) which is key =
//     * peptidesequence + modifications + substitutionModifications
//     *
//     * @param domain
//     * @param iXTandemFile
//     * @return key = peptidesequence + modifications + substitutionModifications
//     */
//    private String getSIIKey(Domain domain, XTandemFile iXTandemFile) {
//        ArrayList<de.proteinms.xtandemparser.interfaces.Modification> fixModList = iXTandemFile.getModificationMap().getFixedModifications(domain.getDomainKey());
//        ArrayList<de.proteinms.xtandemparser.interfaces.Modification> varModList = iXTandemFile.getModificationMap().getVariableModifications(domain.getDomainKey());
//
//        String fixMods = "";
//
//        for (de.proteinms.xtandemparser.interfaces.Modification fixMod : fixModList) {
//            String name = fixMod.getName();
//            if (fixMod.isSubstitution()) {
//                name += "_subs_" + fixMod.getSubstitutedAminoAcid();
//            }
//            int loc = Integer.parseInt(fixMod.getLocation());
//            int pepLoc = loc - domain.getDomainStart() + 1;
//            fixMods += name + "$" + pepLoc + ";";
//        }
//
//        String varMods = "";
//        for (de.proteinms.xtandemparser.interfaces.Modification varMod : varModList) {
//            String name = varMod.getName();
//            if (varMod.isSubstitution()) {
//                name += "_subs_" + varMod.getSubstitutedAminoAcid();
//            }
//
//            int loc = Integer.parseInt(varMod.getLocation());
//            int pepLoc = loc - domain.getDomainStart() + 1;
//            varMods += name + "$" + pepLoc + ";";
//        }
//        String sIIKey = domain.getDomainSequence() + "_" + varMods + "_" + fixMods + "_";
//        return sIIKey;
//
//    }

    public void handleCVs() {


        //<cv id="PSI-MS" fullName="PSI-MS" URI="http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo" version="2.25.0"/>
        //<cv id="UNIMOD" fullName="UNIMOD" URI="http://www.unimod.org/obo/unimod.obo" />
        //<cv id="UO" fullName="UNIT-ONTOLOGY" URI="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo"></cv>

        cvList = new CvList();
        List<Cv> localCvList = cvList.getCv();
        psiCV = new Cv();

        //https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo
        psiCV.setUri("https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo");
        psiCV.setId(psiCvID);
        psiCV.setVersion("4.0.0");https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo
        psiCV.setFullName("PSI-MS");

        unimodCV = new Cv();
        unimodCV.setUri("http://www.unimod.org/obo/unimod.obo");
        unimodCV.setId(unimodID);
        unimodCV.setFullName("UNIMOD");

        xlmodCV = new Cv();
        xlmodCV.setUri("https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/cv/XLMOD.obo");
        xlmodCV.setId(xlmodID);
        xlmodCV.setFullName("XLMOD");
        
        unitCV = new Cv();
        unitCV.setUri("https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo");
        //unitCV.setUri("http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo");
        unitCV.setId(unitCvID);
        unitCV.setFullName("UNIT-ONTOLOGY");

        localCvList.add(psiCV);
        localCvList.add(xlmodCV);
        localCvList.add(unimodCV);
        localCvList.add(unitCV);
    }

    /**
     * Helper method to create and return a CvParam from accession, name and CV
     *
     * @return CvParam
     */
    public CvParam makeCvParam(String accession, String name, Cv cv) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        return cvParam;
    }

    /**
     * Helper method to create and return a CvParam from accession, name and CV
     *
     * @return CvParam
     */
    public CvParam makeCvParam(String accession, String name, Cv cv, String value) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        cvParam.setValue(value);
        return cvParam;
    }

    /**
     * Helper method to create and return a CvParam from accession, name, CV,
     * unitAccession and unitName (unitCV is automatically provided)
     *
     * @return CvParam
     */
    public CvParam makeCvParam(String accession, String name, Cv cv, String unitAccession, String unitName) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        cvParam.setUnitAccession(unitAccession);
        cvParam.setUnitCv(unitCV);
        cvParam.setUnitName(unitName);
        return cvParam;
    }

    /**
     * Helper method to create and return a CvParam from accession, name, CV,
     * unitAccession, unitName and unitCV
     *
     * @return CvParam
     */
    public CvParam makeCvParam(String accession, String name, Cv cv, String unitAccession, String unitName, Cv alternateUnitCV) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        cvParam.setUnitAccession(unitAccession);
        cvParam.setUnitCv(alternateUnitCV);
        cvParam.setUnitName(unitName);
        return cvParam;
    }

    /**
     *
     * Aim is to write out set up the analysisSoftwareList following this
     * structure:
     * <AnalysisSoftware id="ID_software" name="xtandem" version="2008.12.1.1" >
     * <SoftwareName>
     * <cvParam accession="MS:1001476" name="xtandem" cvRef="PSI-MS" />
     * </SoftwareName>
     *
     */
    public void handleAnalysisSoftware(String version, XiInFDR fdr) {
        analysisSoftwareList = new AnalysisSoftwareList();
        List<AnalysisSoftware> analysisSoftwares = analysisSoftwareList.getAnalysisSoftware();
        analysisSoftware = new AnalysisSoftware();
        analysisSoftware.setName("xi");

        analysisSoftware.setId("xi_id");
        Param p = new Param();
        // TODO need to ask for version/software - or get from cmd param
        p.setParam(makeCvParam("MS:1002544","xi",psiCV));
        analysisSoftware.setSoftwareName(p);
        if (fdr.getXiVersion() != null)
            analysisSoftware.setVersion(fdr.getXiVersion().toString());

        analysisSoftwares.add(analysisSoftware);

        analysisSoftware = new AnalysisSoftware();
        analysisSoftware.setName("XiFDR");
        p = new Param();
        p.setParam(makeCvParam("MS:1002543","xiFDR",psiCV));
        analysisSoftware.setSoftwareName(p);

        analysisSoftware.setId("xiFDR_id");
        analysisSoftware.setVersion(version);

        analysisSoftwares.add(analysisSoftware);

    }

    /**
     * Setup Provider element as follows
     * <Provider id="PROVIDER">
     * <ContactRole Contact_ref="PERSON_DOC_OWNER">
     * <role>
     * <cvParam accession="MS:1001271" name="researcher" cvRef="PSI-MS"/>
     * </role>
     * </ContactRole>
     * </Provider>
     *
     */
    public void handleProvider() {
        provider = new Provider();
        provider.setId("PROVIDER");

        ContactRole contactRole = new ContactRole();
        contactRole.setContact(docOwner);


        Role role = new Role();
        role.setCvParam(makeCvParam("MS:1001271", "researcher", psiCV));
        contactRole.setRole(role);

        provider.setContactRole(contactRole);

    }

    /**
     * TO DO Capture name and email of the user
     * <AuditCollection>
     * <Person id="PERSON_DOC_OWNER" firstName="Andy" lastName="Jones"
     * email="someone@someuniversity.com">
     * <affiliations Organization_ref="ORG_DOC_OWNER"/>
     * </Person>
     * <Organization id="ORG_DOC_OWNER" address="Some address" name="Some place"
     * />
     * </AuditCollection>
     *
     *
     */
    public void handleAuditCollection(String firstName, String lastName, String email, String address, String affiliationName) {
        auditCollection = new AuditCollection();
        //List<Contact> contactList = auditCollection.getContactGroup();
        List<AbstractContact> contactList = auditCollection.getPersonOrOrganization();
        docOwner = new Person();
        docOwner.setId("PERSON_DOC_OWNER");
        docOwner.setFirstName(firstName);
        docOwner.setLastName(lastName);
        docOwner.getCvParam().add(makeCvParam("MS:1000587", "contact address", psiCV, address));
        docOwner.getCvParam().add(makeCvParam("MS:1000589", "contact email", psiCV, email));

        //docOwner.setEmail(email);

        Organization org = new Organization();
        org.setId("ORG_DOC_OWNER");
        org.setName(affiliationName);
        org.getCvParam().add(makeCvParam("MS:1000586", "contact name", psiCV, affiliationName));



        List<Affiliation> affList = docOwner.getAffiliation();
        Affiliation aff = new Affiliation();
        aff.setOrganization(org);
        affList.add(aff);
        contactList.add(docOwner);
        contactList.add(org);

    }



    
    public SpectrumIdentificationList getSpectrumIdentificationList(PSM psm) {
        return getSpectrumIdentificationList(psm.getSearchID(),psm.getPeakListName());
    }
    
    public SpectrumIdentificationList getSpectrumIdentificationList(String searchid,String run) {
        HashMap<String,SpectrumIdentificationList> run2list = siList.get(searchid);
        if (run2list == null) {
            run2list = new HashMap<String, SpectrumIdentificationList>();
            siList.put(searchid, run2list);
        }
        SpectrumIdentificationList sil =  run2list.get(run);
        if (sil ==null) {
            sil = new SpectrumIdentificationList();
            run2list.put(run, sil);
            sil.setId(siiListID+"_"+siList.size()+"_"+searchid+"_"+run);
            //sil.getCvParam().add(makeCvParam("MS:1002439", "final PSM list", psiCV));
            getSpectrumIdentification(searchid, run);
        }
        return sil;
    }
            
    
    HashMap<String,HashMap<String,SpectrumIdentification>> spectrumIdentifications =  new HashMap<>();
    public SpectrumIdentification getSpectrumIdentification(String search_id, String run_name)  {
        if (analysisCollection ==null)
            analysisCollection = new AnalysisCollection();
        HashMap<String,SpectrumIdentification> siruns = spectrumIdentifications.get(search_id);
        if (siruns == null) {
            siruns = new HashMap<String, SpectrumIdentification>();
            spectrumIdentifications.put(search_id, siruns);
        }
        SpectrumIdentification specIdent = siruns.get(run_name);
        if (specIdent == null) {
            List<SpectrumIdentification> specIdentList = analysisCollection.getSpectrumIdentification();
            specIdent = new SpectrumIdentification();
            specIdent.setId(specIdentID+"_"+search_id+"_"+run_name);
            specIdentList.add(specIdent);
            siruns.put(run_name, specIdent);
           // analysisCollection.getSpectrumIdentification().add(specIdent);
            
            specIdent.setSpectrumIdentificationList(getSpectrumIdentificationList(search_id, run_name));
            specIdent.setSpectrumIdentificationProtocol(siProtocol.get(search_id));
            List<SearchDatabaseRef> searchDBRefList = specIdent.getSearchDatabaseRef();
            for (SearchDatabase sb : getDatabases(search_id)) {
                SearchDatabaseRef searchDBRef = new SearchDatabaseRef();
                searchDBRef.setSearchDatabase(sb);
                searchDBRefList.add(searchDBRef);
            }

            List<InputSpectra> inputSpecList = specIdent.getInputSpectra();
            InputSpectra inputSpec = new InputSpectra();
            inputSpec.setSpectraData(getSpectraData(search_id, run_name));
            inputSpecList.add(inputSpec);

        }
        return specIdent;
    }
    
    HashMap<String,HashMap<String,SpectraData>> spectraDataList = new HashMap<String, HashMap<String, SpectraData>>();
    public SpectraData getSpectraData(String searchid, String run) {
        HashMap<String,SpectraData> run2SpecData = spectraDataList.get(searchid);
        if (run2SpecData == null) {
            run2SpecData = new HashMap<String, SpectraData>();
            spectraDataList.put(searchid, run2SpecData);
        }
        SpectraData sd = run2SpecData.get(run);
        if (sd == null) {
            sd = new SpectraData();
            SpectrumIDFormat sif = new SpectrumIDFormat();
            if (forceExtension != null && forceExtension.toLowerCase().contentEquals("mzml")) {
                sif.setCvParam(makeCvParam("MS:1001530",  "mzML unique identifier",psiCV));
            } else {
                sif.setCvParam(makeCvParam("MS:1000774", "multiple peak list nativeID format", psiCV));
            }
            sd.setSpectrumIDFormat(sif);
            String runOut = run;
            String extension = "";
            // get the file extension of the run-name
            if (runOut.contains(".")) {
                extension = runOut.substring(runOut.lastIndexOf(".")+1);
            }
            
            if (forceExtension != null) {
                if (runOut.contains("."))
                    runOut = runOut.substring(0,runOut.lastIndexOf("."))+"."+forceExtension;
                else 
                    runOut = runOut + "." + forceExtension;
                extension  = forceExtension;
            }
            sd.setLocation(runOut);
            
            sd.setId("SD_"+searchid + "_"+run);

            /*[Term]
            id: MS:1001062
            name: Mascot MGF format
            def: "Mascot MGF file format." [PSI:MS]
            is_a: MS:1000560 ! mass spectrometer file format
            */
            
            FileFormat ff = new FileFormat();
            if (extension.toLowerCase().contentEquals("mzml")) {
                ff.setCvParam(makeCvParam("MS:1000584",  "mzML format",psiCV));
                scanIDTranslation= this.mzMLscanIDTranslation;
            } else if (extension.toLowerCase().contentEquals("raw")) {
                ff.setCvParam(makeCvParam("MS:1000563","Thermo RAW format", psiCV));
                scanIDTranslation= this.thermoRawScaNumber;
            } else if (extension.toLowerCase().contentEquals("apl")) {
                ff.setCvParam(makeCvParam("MS:1002996","Andromeda:apl file format", psiCV));
                scanIDTranslation= this.thermoRawScaNumber;
            } else {
                ff.setCvParam(makeCvParam("MS:1001062","Mascot MGF format", psiCV));
                scanIDTranslation= this.mgfScanID;
            }
                    
            sd.setFileFormat(ff);
            inputs.getSpectraData().add(sd);
            run2SpecData.put(run, sd);
        }
        return sd;
    }    
    
    
    ArrayList<Enzyme> enzymes = new ArrayList<Enzyme>();
    public ArrayList<Enzyme> getEnzyme(String enzymeName, final int missedCleavage) {
        
        final ArrayList<Enzyme> ret = new ArrayList<Enzyme>(1);
        final String enzymeNameI = enzymeName.replaceAll("(?i)NAME=","");
        String name=enzymeNameI;

        
        String[] e = name2Enzyme.get(name);
        if (e== null) {
            name=name.toLowerCase();
            e=name2Enzyme.get(name);
        }
        if (e==null) {
            name=name.replaceAll("[_\\s-]", "");
            e=name2Enzyme.get(name);
        }
        if (e==null) {
            name=name.replaceAll("\\\\", "/");
            e=name2Enzyme.get(name);
        }

        if (e == null) {
            if (name.contains("+")) {
                for (String s : name.split("\\+")) {
                    ArrayList<Enzyme> r = getEnzyme(s, missedCleavage);
                    ret.addAll(r);
                }
            } else {
                final JDialog dialog = new JDialog();
                dialog.setModal(true);
                dialog.setTitle("Select Enzyme");
                dialog.getContentPane().setLayout(new BoxLayout(dialog.getContentPane(), BoxLayout.Y_AXIS));
                String[] en = new String[psiEnzymes.length];
                int i =0;
                for (String[] s : psiEnzymes) {
                    en[i++] = s[1];
                }
                JLabel l1 = new JLabel("Found a enzyme I don't know:");
                JLabel l2 = new JLabel(enzymeName);
                JLabel l3 = new JLabel("Please select one from the known enzymes:");
                final JList cbEnzymes = new JList(en);
                cbEnzymes.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
                JScrollPane spEzymes = new JScrollPane(cbEnzymes);
                JPanel okcancel = new JPanel();
                JButton ok = new JButton("OK");
                ok.addActionListener(new ActionListener() {

                    public void actionPerformed(ActionEvent e) {
                        String es;
                        for (Object o : cbEnzymes.getSelectedValuesList()) {
                            es = o.toString();
                            ret.addAll(getEnzyme(es, missedCleavage));                        
                            
                        }
                        if (cbEnzymes.getSelectedValuesList().size() == 1) {
                            for (String[] s:psiEnzymes) {
                                if (s[0] == ret.get(0).getEnzymeName().getCvParam().get(0).getAccession()) {
                                    name2Enzyme.put(enzymeNameI, s);
                                }
                            }
                        }
                        dialog.setVisible(false);
                    }
                });
                dialog.getContentPane().add(l1);
                dialog.getContentPane().add(l2);
                dialog.getContentPane().add(l3);
                dialog.getContentPane().add(spEzymes);
                dialog.getContentPane().add(ok);
                dialog.pack();
                dialog.setVisible(true);
       

            }
        } else {
            Enzyme enzyme = new Enzyme();
            enzyme.setId(e[1]+"_"+enzymes.size());
            enzyme.setCTermGain("OH");
            enzyme.setNTermGain("H");   
            enzyme.setMissedCleavages(missedCleavage);
            enzyme.setSemiSpecific(false);
            ParamList paramList = enzyme.getEnzymeName();
            if (paramList == null) {
                paramList = new ParamList();
                enzyme.setEnzymeName(paramList);
            }
            List<CvParam> cvParamList = paramList.getCvParam();
            cvParamList.add(makeCvParam(e[0], e[1], psiCV));    
            enzymes.add(enzyme);
            ret.add(enzyme);
        }
        return ret;
    }
    
    /**
     * <AnalysisProtocolCollection>
     * <SpectrumIdentificationProtocol id="SearchProtocol"
     * AnalysisSoftware_ref="ID_software">
     * <SearchType>
     * <cvParam accession="MS:1001083" name="ms-ms search" cvRef="PSI-MS"/>
     * </SearchType>
     * <AdditionalSearchParams>
     * <cvParam accession="MS:1001211" name="parent mass type mono"
     * cvRef="PSI-MS"/>
     * <cvParam accession="MS:1001256" name="fragment mass type mono"
     * cvRef="PSI-MS"/>
     * </AdditionalSearchParams>
     * <ModificationParams>
     * <SearchModification fixedMod="true">
     * <ModParam massDelta="57.021464" residues="C">
     * <cvParam accession="UNIMOD:4" name="Carbamidomethyl" cvRef="UNIMOD" />
     * </ModParam>
     * </SearchModification>
     * <SearchModification fixedMod="false">
     * <ModParam massDelta="15.994919" residues="M">
     * <cvParam accession="UNIMOD:35" name="Oxidation" cvRef="UNIMOD" />
     * </ModParam>
     * </SearchModification>
     * </ModificationParams>
     * <Enzymes independent="0">
     * <Enzyme id="ENZ_1" CTermGain="OH" NTermGain="H" missedCleavages="1"
     * semiSpecific="0">
     * <EnzymeName>
     * <cvParam accession="MS:1001251" name="Trypsin" cvRef="PSI-MS" />
     * </EnzymeName>
     * </Enzyme>
     * </Enzymes>
     * <MassTable id="0" msLevel="2">
     * </MassTable>
     * <FragmentTolerance>
     * <cvParam accession="MS:1001412" name="search tolerance plus value"
     * value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton"
     * unitCvRef="UO" />
     * <cvParam accession="MS:1001413" name="search tolerance minus value"
     * value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton"
     * unitCvRef="UO" />
     * </FragmentTolerance>
     * <ParentTolerance>
     * <cvParam accession="MS:1001412" name="search tolerance plus value"
     * value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton"
     * unitCvRef="UO" />
     * <cvParam accession="MS:1001413" name="search tolerance minus value"
     * value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton"
     * unitCvRef="UO" />
     * </ParentTolerance>
     * <Threshold>
     * <cvParam accession="MS:1001494" name="no threshold" cvRef="PSI-MS" />
     * </Threshold>
     * </SpectrumIdentificationProtocol>
     * </AnalysisProtocolCollection>
     *
     *
     */
    public boolean handleAnalysisProtocolCollection(OfflineFDR fdr, FDRResult result, ArrayList<CvParam> stubs) {//InputParams inputParams){

        //boolean (parentIsMono, boolean fragmentIsMono, SearchModification[] searchMods, String enzymeName, double parTolPlus, double parTolMinus, double fragTolPlus, double fragTolMinus);
        if (analysisProtocolCollection == null) {
            analysisProtocolCollection = new AnalysisProtocolCollection();
        }
        List<SpectrumIdentificationProtocol> sipList = analysisProtocolCollection.getSpectrumIdentificationProtocol();
        for (String sID : ((XiInFDR)fdr).getSearchIDs()) {
            RunConfig conf =  ((XiInFDR)fdr).getConfig(sID);
            SpectrumIdentificationProtocol siProt = new SpectrumIdentificationProtocol();
            siProt.setId(siProtocolID+"_"+sID);
            siProt.setAnalysisSoftware(analysisSoftware);
            
            //<cvParam accession="MS:1001083" name="ms-ms search" cvRef="PSI-MS"/>
            //siProtocol.setSearchType(makeCvParam("MS:1001083","ms-ms search",psiCV));
            Param tempParam = new Param();
            tempParam.setParam(makeCvParam("MS:1001083", "ms-ms search", psiCV));
            siProt.setSearchType(tempParam);

            ParamList paramList = siProt.getAdditionalSearchParams();
            if (paramList == null) {
                paramList = new ParamList();
                siProt.setAdditionalSearchParams(paramList);
            }

            List<CvParam> cvParamList = paramList.getCvParam();

            boolean parentIsMono = true;  //does not appear to be a way in Tandem of specifying parent mass is average
            if (parentIsMono) {
                cvParamList.add(makeCvParam("MS:1001211", "parent mass type mono", psiCV));
            } else {
                cvParamList.add(makeCvParam("MS:1001212", "parent mass type average", psiCV));
            }

            cvParamList.add(makeCvParam("MS:1002494", "cross-linking search", psiCV));
            cvParamList.add(makeCvParam("MS:1003330", "noncovalently associated peptides search", psiCV));
            cvParamList.add(makeCvParam("MS:1003343", "FDR applied separately to self crosslinks and protein heteromeric crosslinks", psiCV, "true"));
            if (conf != null) {
                ArrayList<Method> fragm = conf.getFragmentMethods();
//                new String[] {"MS:1001108", "param: a ion"},
//                new String[] {"MS:1001118", "param: b ion"},
//                new String[] {"MS:1001108", "param: a ion"},
//                new String[] {"MS:1001261", "param: x ion"},
//                new String[] {"MS:1001262", "param: y ion"},
//                new String[] {"MS:1001263", "param: z ion"}
    
                for (Method m : fragm) {
                    if (m.getDeclaringClass() == AIon.class) {
                        cvParamList.add(makeCvParam("MS:1001108", "param: a ion", psiCV));
                    } else if (m.getDeclaringClass() == rappsilber.ms.sequence.ions.BIon.class) {
                        cvParamList.add(makeCvParam("MS:1001118", "param: b ion", psiCV));
                    } else if (m.getDeclaringClass() == rappsilber.ms.sequence.ions.CIon.class) {
                        cvParamList.add(makeCvParam("MS:1001119", "param: c ion", psiCV));
                    } else if (m.getDeclaringClass() == rappsilber.ms.sequence.ions.XIon.class) {
                        cvParamList.add(makeCvParam("MS:1001261", "param: x ion", psiCV));
                    } else if (m.getDeclaringClass() == rappsilber.ms.sequence.ions.YIon.class) {
                        cvParamList.add(makeCvParam("MS:1001262", "param: y ion", psiCV));
                    } else if (m.getDeclaringClass() == rappsilber.ms.sequence.ions.ZIon.class) {
                        cvParamList.add(makeCvParam("MS:1001263", "param: z ion", psiCV));
                    }
                }
            }

    //        // state that we aredoing something with cross-links
    //        cvParamList.add(makeCvParam("MS:100XXXX", "crosslinking search", psiCV));


            boolean fragmentIsMono = true;
    //        if("average".equals(inputParams.getSpectrumFragMassType())){
    //            fragmentIsMono = false;
    //            cvParamList.add(makeCvParam("MS:1001255","fragment mass type average",psiCV));
    //        }
    //        else{
            cvParamList.add(makeCvParam("MS:1001256", "fragment mass type mono", psiCV));
    //        }


            ModificationParams modParams = new ModificationParams();
            List<SearchModification> searchModList = modParams.getSearchModification();

            //residue, potential modification mass
            ArrayList<AminoModification> fixedmods = conf.getFixedModifications();
            
            for (rappsilber.ms.crosslinker.CrossLinker xl : conf.getCrossLinker()) {
                if (!(xl instanceof NonCovalentBound)){
                    List<SearchModification> searchMod = translateToSearchModification(xl, stubs, ((XiInFDR)fdr).getConfig());
                    searchModList.addAll(searchMod);
                }
            }
            convertModifications(conf.getVariableModifications(), fragmentIsMono, fdr, searchModList, conf, false);
            convertModifications(conf.getFixedModifications(), fragmentIsMono, fdr, searchModList, conf, true);
            
            // known modifications are converted as variable - even so they are not really searched that way
            convertModifications(conf.getKnownModifications(), fragmentIsMono, fdr, searchModList, conf, false);


            /*
             <ModificationParams>convertModifications
             <SearchModification fixedMod="false">
             <ModParam massDelta="15.994919" residues="M">
             <cvParam accession="UNIMOD:35" name="Oxidation" cvRef="UNIMOD" />
             </ModParam>
             </SearchModification>
             */

            /*
             <Enzymes independent="0">
             <Enzyme id="ENZ_1" CTermGain="OH" NTermGain="H" missedCleavages="1" semiSpecific="0">
             <EnzymeName>
             <cvParam accession="MS:1001251" name="Trypsin" cvRef="PSI-MS" />
             </EnzymeName>
             </Enzyme>
             </Enzymes>
             */
            //Only add this group if there are any modifications in the list:
            if (searchModList.size() > 0) {
                siProt.setModificationParams(modParams);
            }
            Enzymes enzymes = siProt.getEnzymes();

            if (enzymes == null) {
                enzymes = new Enzymes();
                siProt.setEnzymes(enzymes);
            }

            enzymes.setIndependent(false);

            List<Enzyme> enzymeList = enzymes.getEnzyme();

            for (Enzyme enzyme : getEnzyme(conf.getDigestion_method().Name(), conf.getMaxMissCleavages())) {
                enzymeList.add(enzyme);
            }

            Tolerance fragTol = new Tolerance();
            Tolerance parTol = new Tolerance();

            boolean isDaltons;

            if (conf.getFragmentTolerance().isRelative()) {
                isDaltons = false;
            } else {
                isDaltons = true;
                unimodMassError = conf.getFragmentTolerance().getValue();        //Dynamically set Unimod lookup mass error
            }

            List<CvParam> fragCvList = fragTol.getCvParam();
            CvParam fragCvPlus = getCvParamWithMassUnits(isDaltons);
            CvParam fragCvMinus = getCvParamWithMassUnits(isDaltons);


            /*
             <FragmentTolerance>
             <cvParam accession="MS:1001412" name="search tolerance plus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
             <cvParam accession="MS:1001413" name="search tolerance minus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
             </FragmentTolerance>
             <ParentTolerance>
             <cvParam accession="MS:1001412" name="search tolerance plus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
             <cvParam accession="MS:1001413" name="search tolerance minus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
             </ParentTolerance>
             */

            fragCvPlus.setAccession("MS:1001412");
            fragCvPlus.setName("search tolerance plus value");
            fragCvMinus.setAccession("MS:1001413");
            fragCvMinus.setName("search tolerance minus value");
            fragCvPlus.setValue("" + conf.getFragmentTolerance().getValue());
            fragCvMinus.setValue("" + conf.getFragmentTolerance().getValue());
            fragCvList.add(fragCvPlus);
            fragCvList.add(fragCvMinus);

            List<CvParam> parCvList = parTol.getCvParam();

            if (conf.getPrecousorTolerance().isRelative()) {
                isDaltons = false;
            } else {
                isDaltons = true;
            }

            CvParam parCvPlus = getCvParamWithMassUnits(isDaltons);
            CvParam parCvMinus = getCvParamWithMassUnits(isDaltons);

            parCvPlus.setAccession("MS:1001412");
            parCvPlus.setName("search tolerance plus value");
            parCvMinus.setAccession("MS:1001413");
            parCvMinus.setName("search tolerance minus value");
            parCvPlus.setValue("" + conf.getPrecousorTolerance().getValue());
            parCvMinus.setValue("" + conf.getPrecousorTolerance().getValue());
            parCvList.add(parCvPlus);
            parCvList.add(parCvMinus);

            siProt.setFragmentTolerance(fragTol);
            siProt.setParentTolerance(parTol);

            ParamList sip_paramList = siProt.getThreshold();
            if (sip_paramList == null) {
                sip_paramList = new ParamList();
                siProt.setThreshold(sip_paramList);
            }
            cvParamList = sip_paramList.getCvParam();
            cvParamList.add(makeCvParam("MS:1001494", "no threshold", psiCV));
            //<cvParam accession="MS:1001494" name="no threshold" cvRef="PSI-MS" />
            sipList.add(siProt);
            siProtocol.put(sID.toString(), siProt);
            if (sID.toString().trim().length()==0)
                siProtocol.put(null, siProt);
        }
        return true;
    }

    protected void convertModifications(Collection<AminoModification> modList, boolean fragmentIsMono, OfflineFDR fdr1, List<SearchModification> searchModList, RunConfig conf, boolean fixed) {
        HashMap<Double, ArrayList<AminoModification>> sameMassModification = new HashMap<>();
        HashMap<Double, ArrayList<AminoModification>> sameMassModificationProtNterm = new HashMap<>();
        HashMap<Double, ArrayList<AminoModification>> sameMassModificationProtCterm = new HashMap<>();
        HashMap<Double, ArrayList<AminoModification>> sameMassModificationPepNterm = new HashMap<>();
        HashMap<Double, ArrayList<AminoModification>> sameMassModificationPepCterm = new HashMap<>();
        
        for (AminoModification am : modList) {
            AminoAcid  base = am.BaseAminoAcid;
            Double weightDiff = am.weightDiff;
            while (base instanceof AminoModification) {
                weightDiff += ((AminoModification)base).weightDiff;
                base= ((AminoModification)base).BaseAminoAcid;
            }
            if (base != AminoAcid.X) {
                Double mass = massNormalise(weightDiff);
                ArrayList<AminoModification> mods = null;
                if (am.prot_position == AminoModification.POSITIONAL_NTERMINAL) {
                    mods = getMassModiciations(sameMassModificationProtNterm, mass);
                } else if (am.prot_position == AminoModification.POSITIONAL_CTERMINAL) {
                    mods = getMassModiciations(sameMassModificationProtCterm, mass);
                } else if (am.pep_position == AminoModification.POSITIONAL_NTERMINAL) {
                    mods = getMassModiciations(sameMassModificationPepNterm, mass);
                } else if (am.pep_position == AminoModification.POSITIONAL_CTERMINAL) {
                    mods = getMassModiciations(sameMassModificationPepCterm, mass);
                } else {
                    mods = getMassModiciations(sameMassModification, mass);
                }
                mods.add(am);
            }
        }
        
        convertTermModification(sameMassModification, fragmentIsMono, fdr1, searchModList,
                false, false, false, false);

        convertTermModification(sameMassModificationPepNterm, fragmentIsMono, fdr1, searchModList,
                true, false, false, false);

        convertTermModification(sameMassModificationPepCterm, fragmentIsMono, fdr1, searchModList,
                false, true, false, false);
        
        convertTermModification(sameMassModificationProtNterm, fragmentIsMono, fdr1, searchModList,
                false, false, true, false);

        convertTermModification(sameMassModificationProtCterm, fragmentIsMono, fdr1, searchModList,
                false, false, false, true);
        
    }

    protected void convertTermModification(
            HashMap<Double, ArrayList<AminoModification>> sameMassModification, 
            boolean fragmentIsMono, 
            OfflineFDR fdr, 
            List<SearchModification> searchModList,
            boolean pepNTerm, boolean pepCTerm,
            boolean protNTerm, boolean protCTerm) {
        for (Double mass : sameMassModification.keySet()) {
            ArrayList<AminoModification> mods = sameMassModification.get(mass);
            ArrayList<String> cResidues = new ArrayList<String>(mods.size());
            AminoModification am = null;
            AminoAcid  base = null;
            for (int i = 0; i < mods.size(); i++) {
                am = mods.get(i);
                base = am.BaseAminoAcid;
                while (base instanceof AminoModification) {
                    base= ((AminoModification)base).BaseAminoAcid;
                }                
                String residue = base.SequenceID;
                if (residue.contains("X")) {
                    cResidues.clear();
                    cResidues.add(".");
                    break;
                }
                cResidues.add(residue);
            }
            convertTermModificationAddSearcModList(cResidues, am, base, fragmentIsMono, fdr, pepNTerm, pepCTerm, protNTerm, protCTerm, searchModList);
        }
    }

    protected void convertTermModificationAddSearcModList(ArrayList<String> cResidues, AminoModification am, AminoAcid base, boolean fragmentIsMono, OfflineFDR fdr1, boolean pepNTerm, boolean pepCTerm, boolean protNTerm, boolean protCTerm, List<SearchModification> searchModList) {
        if (cResidues.size() > 0) {
            if (cResidues.size()>20) {
                cResidues.clear();
                cResidues.add(".");
            }
            SearchModification searchMod = translateToSearchModification(am, base, cResidues, fragmentIsMono, false, ((XiInFDR) fdr1).getConfig(), pepNTerm, pepCTerm, protNTerm, protCTerm);
            searchModList.add(searchMod);
        }
    }

    protected ArrayList<AminoModification> getMassModiciations(HashMap<Double, ArrayList<AminoModification>> sameMassModification, Double mass) {
        ArrayList<AminoModification> mods = sameMassModification.get(mass);
        if (mods == null) {
            mods = new ArrayList<AminoModification>();
            sameMassModification.put(mass, mods);
        }
        return mods;
    }

    protected static double massNormalise(double mass) {
        return Math.round(mass * 100000)/100000.0;
    }

    /**
     * Translated the reported modification mass to the SearchModification
     * object, complete with unimod codes, title, etc
     *
     * @param reportedMod
     * @param fragmentIsMono
     * @param isFixedMod
     * @return
     */
    private SearchModification translateToSearchModification(
            AminoModification reportedMod,
            AminoAcid base,
            Collection<String> modResidues, 
            boolean fragmentIsMono, 
            boolean isFixedMod, 
            rappsilber.config.RunConfig conf, 
            boolean pepNTerm, boolean pepCTerm,
            boolean protNTerm, boolean protCTerm) {
        ArrayList<String> residues = new ArrayList<String>();
        double monoMass = reportedMod.mass - base.mass;

        residues.addAll(modResidues);
        uk.ac.liv.unimod.ModT unimod = unimodDoc.getModByMass(monoMass, unimodMassError, fragmentIsMono, residues);

        SearchModification searchMod = new SearchModification();
        searchMod.setFixedMod(isFixedMod);

        if (unimod == null && base != AminoAcid.X) {
            String modname = reportedMod.SequenceID.replaceAll("[A-Z]", "");
            boolean[] term =new boolean[]{false};
            XLModEntry e = xlmod.guessModificationCached(new XLModQuery(monoMass, modname, new String[]{base.SequenceID}, term, term, term, term, false));
            
            if (e==null) {
                // could be a cross-linker specific modification
                for (rappsilber.ms.crosslinker.CrossLinker cl :conf.getCrossLinker()) {
                    e = xlmod.guessModificationCached(new XLModQuery(monoMass, cl.getName(), new String[]{base.SequenceID}, term, term, term, term, false));                
                }
            }
            // nothing found wih name so try without name
            if (e==null) {
                
                e = xlmod.guessModificationCached(new XLModQuery(monoMass, new String[]{base.SequenceID}, term, term, term, term, false));
            }
            if  (e == null) {
                searchMod.getCvParam().add(makeCvParam("MS:1001460", "unknown modification", psiCV));
            } else {
                searchMod.getCvParam().add(makeCvParam(e.getId(), e.getName(), xlmodCV));
            }
        } else if (base != AminoAcid.X) {

            searchMod.getCvParam().add(makeCvParam("UNIMOD:" + unimod.getRecordId(), unimod.getTitle(), unimodCV));
        }
        searchMod.setMassDelta(new Float(monoMass));
        
        if (pepNTerm) {
            searchMod.getCvParam().add(makeCvParam("MS:1001189", "modification specificity N-term", psiCV));
        } else if (pepCTerm) {
            searchMod.getCvParam().add(makeCvParam("MS:1001190", "modification specificity C-term", psiCV));
        } else if (protCTerm) {
            searchMod.getCvParam().add(makeCvParam("MS:1002058", "modification specificity protein C-term", psiCV));
        } else if (protNTerm) {
            searchMod.getCvParam().add(makeCvParam("MS:1002057", "modification specificity protein N-term", psiCV));
        }

        searchMod.getResidues().addAll(residues);

        return searchMod;
    }

    int countCrossLinker = 0;
    /**
     * Translated the reported modification mass to the SearchModification
     * object, complete with unimod codes, title, etc
     *
     * @param reportedMod
     * @param fragmentIsMono
     * @param isFixedMod
     * @return
     */
    private List<SearchModification> translateToSearchModification(rappsilber.ms.crosslinker.CrossLinker crosslinker, List<CvParam> stubs, rappsilber.config.RunConfig conf) {
        ArrayList<SearchModification> ret = new ArrayList<SearchModification>();
        //String[] temp = reportedMod.split("@");
        double monoMass = massNormalise(crosslinker.getCrossLinkedMass());


        SearchModification searchMod = new SearchModification();
        searchMod.setFixedMod(false);
        searchMod.setMassDelta((float) monoMass);
        searchMod.getCvParam().add(makeCvParam(crosslinkedDonorModAcc, crosslinkedDonorModName, psiCV, ""+countCrossLinker));

        SearchModification searchModAcceptor = new SearchModification();
        searchModAcceptor.setFixedMod(true);
        searchModAcceptor.setMassDelta(0f);
        searchModAcceptor.getCvParam().add(makeCvParam(getCrosslinkedAcceptorModAcc(), crosslinkedAcceptorModName, psiCV, ""+countCrossLinker));
        
        boolean[] term =new boolean[]{false,false};

        XLModEntry e = null;
        CvParam modParam = null;
        if (crosslinker instanceof rappsilber.ms.crosslinker.AminoAcidRestrictedCrossLinker) {
            rappsilber.ms.crosslinker.AminoAcidRestrictedCrossLinker xl = (rappsilber.ms.crosslinker.AminoAcidRestrictedCrossLinker) crosslinker;

            String[] specificity = new String[]{"*","*"};
            if (xl.getAASpecificity(0).size() >0) {
                specificity[0] = xl.getAASpecificity(0).iterator().next().SequenceID;
            }
            if (xl.getAASpecificity(1).size() >0) {
                specificity[1] = xl.getAASpecificity(1).iterator().next().SequenceID;
            }
            e = xlmod.guessModificationCached(new XLModQuery(monoMass, crosslinker.getName(), specificity, term, term, term, term, true));

            // nothing found wih name so try without name
            if (e==null) {
                e = xlmod.guessModificationCached(new XLModQuery(monoMass, specificity, term, term, term, term, false));
            }
        }
        if  (e == null) {
            modParam = makeCvParam("MS:1001460", "unknown modification", psiCV);
        } else {
            modParam = makeCvParam(e.getId(), e.getName(), xlmodCV);
        }
        searchMod.getCvParam().add(modParam);
        searchMod.getCvParam().addAll(stubs);
        searchModAcceptor.getCvParam().addAll(stubs);
        
        if (crosslinker instanceof rappsilber.ms.crosslinker.AminoAcidRestrictedCrossLinker) {
            HashSet<String> modsRes = new HashSet<>();
            boolean all = false;
            if (((rappsilber.ms.crosslinker.AminoAcidRestrictedCrossLinker)crosslinker).getAASpecificity(0).isEmpty()) {
                searchMod.getResidues().add(".");
                searchModAcceptor.getResidues().add(".");
                all = true;
            }else {
                for (AminoAcid aa: ((rappsilber.ms.crosslinker.AminoAcidRestrictedCrossLinker)crosslinker).getAASpecificity(0))
                    if (!aa.SequenceID.contains("X")) {
                        modsRes.add(aa.SequenceID.substring(0,1));
                    } else {
                        modsRes.clear();
                        modsRes.add(".");
                        break;
                    }
            }
            
            if ((!all) && ((rappsilber.ms.crosslinker.AminoAcidRestrictedCrossLinker)crosslinker).getAASpecificity(1).isEmpty()) {
                if (!all) {
                    searchMod.getResidues().add(".");
                    searchModAcceptor.getResidues().add(".");
                    all = true;
                }
            }else if (!all){
                for (AminoAcid aa: ((rappsilber.ms.crosslinker.AminoAcidRestrictedCrossLinker)crosslinker).getAASpecificity(1))
                    if (!aa.SequenceID.contains("X")) {
                        modsRes.add(aa.SequenceID.substring(0,1));
                    } else {
                        modsRes.clear();
                        modsRes.add(".");
                        break;
                    }
            }
            
            HashSet<String> existing = new HashSet<>(searchMod.getResidues());
            for (String aa : modsRes) {
                if (!existing.contains(aa)) {
                    searchMod.getResidues().add(aa);
                    searchModAcceptor.getResidues().add(aa);
                }
            }
        } else {
            allAminoAcidsToModificationResidues(conf, searchMod);
            allAminoAcidsToModificationResidues(conf, searchModAcceptor);
        }
        
        ret.add(searchMod);
        ret.add(searchModAcceptor);
        if (crosslinker.linksCTerminal(0)) {
            SearchModification searchModCTerm = new SearchModification();
            searchModCTerm.setFixedMod(true);
            searchModCTerm.setMassDelta((float) monoMass);
            searchModCTerm.getCvParam().add(modParam);
            searchModCTerm.getResidues().add(".");
            SpecificityRules sr = new SpecificityRules();
            sr.getCvParam().add(makeCvParam("MS:1002058", "modification specificity protein C-term", psiCV));
            searchModCTerm.getSpecificityRules().add(sr);
            searchModCTerm.getCvParam().add(makeCvParam(crosslinkedDonorModAcc, crosslinkedDonorModName, psiCV, ""+countCrossLinker));
            searchModCTerm.getCvParam().addAll(stubs);
            ret.add(searchModCTerm);
            SearchModification searchModCTermAcceptor = new SearchModification();
            searchModCTermAcceptor.setFixedMod(true);
            searchModCTermAcceptor.setMassDelta((float) 0);
            searchModCTermAcceptor.getCvParam().add(modParam);
            searchModCTermAcceptor.getResidues().add(".");
            searchModCTermAcceptor.getSpecificityRules().add(sr);
            searchModCTermAcceptor.getCvParam().add(makeCvParam(getCrosslinkedAcceptorModAcc(), crosslinkedAcceptorModName, psiCV, ""+countCrossLinker));
            searchModCTermAcceptor.getCvParam().addAll(stubs);
            ret.add(searchModCTermAcceptor);
        }

        if (crosslinker.linksNTerminal(0)) {
            SearchModification searchModNTerm = new SearchModification();
            searchModNTerm.setFixedMod(true);
            searchModNTerm.setMassDelta((float) monoMass);
            searchModNTerm.getCvParam().add(modParam);
            searchModNTerm.getResidues().add(".");
            SpecificityRules sr = new SpecificityRules();
            sr.getCvParam().add(makeCvParam("MS:1002057", "modification specificity protein N-term", psiCV));
            searchModNTerm.getSpecificityRules().add(sr);
            searchModNTerm.getCvParam().add(makeCvParam(crosslinkedDonorModAcc, crosslinkedDonorModName, psiCV, ""+countCrossLinker));
            searchModNTerm.getCvParam().addAll(stubs);
            ret.add(searchModNTerm);
            SearchModification searchModNTermAcceptor = new SearchModification();
            searchModNTermAcceptor.setFixedMod(true);
            searchModNTermAcceptor.setMassDelta((float) 0);
            searchModNTermAcceptor.getCvParam().add(modParam);
            searchModNTermAcceptor.getResidues().add(".");
            searchModNTermAcceptor.getSpecificityRules().add(sr);
            searchModNTermAcceptor.getCvParam().add(makeCvParam(getCrosslinkedAcceptorModAcc(), crosslinkedAcceptorModName, psiCV, ""+countCrossLinker));
            searchModNTermAcceptor.getCvParam().addAll(stubs);
            ret.add(searchModNTermAcceptor);
        }
        
        countCrossLinker ++;
        return ret;
    }

    protected void allAminoAcidsToModificationResidues(RunConfig conf, SearchModification searchMod) {
        HashSet<String> aas = new HashSet<>();
        HashSet<String> existing = new HashSet<>(searchMod.getResidues());
        for (AminoAcid aa : conf.getAllAminoAcids()) {
            String res = aa.SequenceID.substring(0,1);
            if (!(existing.contains(res) || res.contains("X")))
                aas.add(res);
        }
        for (String aa : aas) {
            searchMod.getResidues().add(aa);
        }
    }

    
    
    HashMap<String, ArrayList<SearchDatabase>> searchDatabases = new HashMap<>();
    HashMap<String, SearchDatabase> searchDatabasesByID = new HashMap<>();
    boolean showedMultipleSearchDBSWarning = false;
    public ArrayList<SearchDatabase> getDatabases(String searchID) {
        ArrayList<SearchDatabase> dbs = searchDatabases.get(searchID);
        if (dbs == null) {
            ArrayList<String> dbName = new ArrayList<String>();
            ArrayList<String> dbId = new ArrayList<String>();
            int c = ((XiInFDR)fdr).getFastas(searchID, dbId, dbName);
            
            dbs = new ArrayList<SearchDatabase>(c);
            searchDatabases.put(searchID, dbs);
            for (int i = 0; i< c; i++) {
                SearchDatabase searchDB = searchDatabasesByID.get(dbId.get(i));
                
                if (searchDB == null) {
               
                    searchDB = new SearchDatabase();
                    searchDB.setId("SDB_"+searchID+"_"+dbId.get(i));
                    //searchDB.setNumDatabaseSequences(numProts);
                    //<cvParam accession="MS:1001401" name="xtandem xml file" cvRef="PSI-MS"/>

                    UserParam param = new UserParam();
                    param.setName(dbName.get(i));
                    Param tempParam = new Param();
                    tempParam.setParam(param);
                    searchDB.setDatabaseName(tempParam);

                    //searchDB.setDatabaseName(param);
                    searchDB.setLocation(dbName.get(i));


                    FileFormat ff = new FileFormat();
                    ff.setCvParam(makeCvParam(this.databaseFileFormatID, this.databaseFileFormatName, psiCV));
                    searchDB.setFileFormat(ff);
                    inputs.getSearchDatabase().add(searchDB);
                }
                dbs.add(searchDB);
            }
            if (c>0 && showedMultipleSearchDBSWarning) {
                showedMultipleSearchDBSWarning = true;
                JOptionPane.showMessageDialog(null, "At least one search was done with multiple fasta files!\nCurrrently all found proteins will be assigned to one of these (which is likely incorrect).", "More then one fasta file", JOptionPane.WARNING_MESSAGE);
            }
        }
        return dbs;
    }
    
    
    
    public void writeMzidFile(String outputfile) {
            
        try {
            FileOutputStream fwriter = new FileOutputStream(outputfile);
            OutputStreamWriter writer = new OutputStreamWriter(fwriter, "UTF-8");
            MzIdentMLMarshaller m = new MzIdentMLMarshaller(MzIdentMLVersion.Version_1_3);

            MzIdentML mzid = new MzIdentML();
            mzid.setId(UUID.randomUUID().toString());
            
            mzid.setVersion("1.3.0");
            mzid.setCvList(cvList);
            mzid.setAnalysisSoftwareList(analysisSoftwareList);
            mzid.setProvider(provider);
            mzid.setAuditCollection(auditCollection);
            mzid.setSequenceCollection(sequenceCollection);
            mzid.setAnalysisCollection(analysisCollection);
            mzid.setAnalysisProtocolCollection(analysisProtocolCollection);
            

            //m.marshal(analysisSampleCollection, writer);     //TODO - complete this part
            //writer.write("\n");
            DataCollection dc = new DataCollection();
            dc.setInputs(inputs);
            AnalysisData ad = new AnalysisData();
            dc.setAnalysisData(ad);
            mzid.setDataCollection(dc);

            /*writer.write(m.createDataCollectionStartTag() + "\n");
            m.marshal(inputs, writer);
            writer.write("\n");


            //Inputs inputs = unmarshaller.unmarshal(MzIdentMLElement.Inputs.getXpath());
            //m.marshal(inputs, writer);
            //writer.write("\n");

            writer.write(m.createAnalysisDataStartTag() + "\n");
            */


            for (HashMap<String,SpectrumIdentificationList> run2sil : siList.values()) {
                for (SpectrumIdentificationList sil : run2sil.values()) {
                    ad.getSpectrumIdentificationList().add(sil);
                }
            }


            
            proteinDetectionList.setId("PDL_1");
            proteinDetectionList.getCvParam().add(makeCvParam("MS:1002404", "count of identified proteins",psiCV,""+ proteinDetectionList.getProteinAmbiguityGroup().size()));
            ad.setProteinDetectionList(proteinDetectionList);

//            writer.write(m.createProteinDetectionListClosingTag() + "\n");

            // creat an ID by milliseconds
            String id = "xiFDR" + LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyyMMddHHmmssSSS"));
            // make it "unique" by adding 4 digits of randomness
            id += "_" + Math.round(Math.random()*10000);

            m.marshallRoot(mzid, writer, "UTF-8", id, MzIdentMLVersion.Version_1_3);

            writer.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    /**
     * Helper method to setup a CvParam with CVRef, with either Daltons or ppm
     * as units
     *
     */
    public CvParam getCvParamWithMassUnits(boolean isDaltonUnit) {
        CvParam cvParam = new CvParam();

        //<cvParam accession="MS:1001413" name="search tolerance minus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
        cvParam.setCv(psiCV);
        cvParam.setUnitCv(unitCV);

        if (isDaltonUnit) {
            cvParam.setUnitAccession("UO:0000221");
            cvParam.setUnitName("dalton");
        } else {
            cvParam.setUnitAccession("UO:0000169");
            cvParam.setUnitName("parts per million");
        }
        return cvParam;
    }



    /**
     * Used to identify members of cross-linked PSMs
     * @return the crosslinkedSIIAcc
     */
    public static String getCrosslinkedSIIAcc() {
        return crosslinkedSIIAcc;
    }

    /**
     * Used to identify members of cross-linked PSMs
     * @param crosslinkedSIIAcc the crosslinkedSIIAcc to set
     */
    public static void setCrosslinkedSIIAcc(String Accession) {
        crosslinkedSIIAcc = Accession;
    }

    /**
     * Used to identify modifications, that span several peptides
     * @return the crosslinkedDonorModAcc
     */
    public static String getCrosslinkedModAcc() {
        return crosslinkedDonorModAcc;
    }

    /**
     * Used to identify modifications, that span several peptides
     * @param crosslinkedDonorModAcc the crosslinkedDonorModAcc to set
     */
    public static void setCrosslinkedModAcc(String crosslinkedModAcc) {
        crosslinkedDonorModAcc = crosslinkedModAcc;
    }
    /**
     * @return the crosslinkedReceptorModAcc
     */
    public static String getCrosslinkedAcceptorModAcc() {
        return crosslinkedAcceptorModAcc;
    }

    /**
     * @param aCrosslinkedReceptorModAcc the crosslinkedReceptorModAcc to set
     */
    public static void setCrosslinkedAcceptorModAcc(String aCrosslinkedReceptorModAcc) {
        crosslinkedAcceptorModAcc = aCrosslinkedReceptorModAcc;
    }

    /**
     * @return the crosslinkedReceptorModName
     */
    public static String getCrosslinkedAcceptorModName() {
        return crosslinkedAcceptorModName;
    }

    /**
     * @param aCrosslinkedReceptorModName the crosslinkedReceptorModName to set
     */
    public static void setCrosslinkedReceptorModName(String aCrosslinkedReceptorModName) {
        crosslinkedAcceptorModName = aCrosslinkedReceptorModName;
    }

    /**
     * @return the ownerFirst
     */
    public String getOwnerFirst() {
        return owner.first;
    }

    /**
     * @param ownerFirst the ownerFirst to set
     */
    public void setOwnerFirst(String ownerFirst) {
        this.owner.first = ownerFirst;
    }

    /**
     * @return the ownerLast
     */
    public String getOwnerLast() {
        return owner.last;
    }

    /**
     * @param ownerLast the ownerLast to set
     */
    public void setOwnerLast(String ownerLast) {
        this.owner.last = ownerLast;
    }

    /**
     * @return the ownerOrg
     */
    public String getOwnerOrg() {
        return owner.org;
    }

    /**
     * @param ownerOrg the ownerOrg to set
     */
    public void setOwnerOrg(String ownerOrg) {
        this.owner.org = ownerOrg;
    }

    /**
     * @return the ownerAddress
     */
    public String getOwnerAddress() {
        return owner.address;
    }

    /**
     * @param ownerAddress the ownerAddress to set
     */
    public void setOwnerAddress(String ownerAddress) {
        this.owner.address = ownerAddress;
    }

    /**
     * @return the ownerEmail
     */
    public String getOwnerEmail() {
        return owner.email;
    }

    /**
     * @param ownerEmail the ownerEmail to set
     */
    public void setOwnerEmail(String ownerEmail) {
        this.owner.email = ownerEmail;
    }

    /**
     * Used to identify modifications, that span several peptides.
     * This is denotes the donor
     * @return the crosslinkedDonorModAcc
     */
    public static String getCrosslinkedDonorModAcc() {
        return crosslinkedDonorModAcc;
    }

    /**
     * Used to identify modifications, that span several peptides.
     * This is denotes the donor
     * @param aCrosslinkedDonorModAcc the crosslinkedDonorModAcc to set
     */
    public static void setCrosslinkedDonorModAcc(String aCrosslinkedDonorModAcc) {
        crosslinkedDonorModAcc = aCrosslinkedDonorModAcc;
    }

    /**
     * @return the crosslinkedDonorModName
     */
    public static String getCrosslinkedDonorModName() {
        return crosslinkedDonorModName;
    }

    /**
     * @param aCrosslinkedDonorModName the crosslinkedDonorModName to set
     */
    public static void setCrosslinkedDonorModName(String aCrosslinkedDonorModName) {
        crosslinkedDonorModName = aCrosslinkedDonorModName;
    }

    /**
     * @return the forceExtension
     */
    public String getForceExtension() {
        return forceExtension;
    }

    /**
     * @param forceExtension the forceExtension to set
     */
    public void setForceExtension(String forceExtension) {
        this.forceExtension = forceExtension;
    }

    
    /**
     * @return the forceExtension
     */
    public String getMZMLTemplate() {
        if (mzMLscanIDTranslation != null) {
            return mzMLscanIDTranslation.getTemplate();
        }
        return null;
    }

    /**
     * @param forceExtension the forceExtension to set
     */
    public void setMZMLTemplate(String template) {
        mzMLscanIDTranslation = new MzMLScanTranslation(template);
    }
    
}
