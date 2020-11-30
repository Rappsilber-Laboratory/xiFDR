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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.lang.reflect.ParameterizedType;
import java.lang.reflect.Type;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map;
import java.util.RandomAccess;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.rappsilber.data.csv.CSVRandomAccess;
import org.rappsilber.fdr.entities.DirectionalPeptidePair;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.fdr.entities.Peptide;
import org.rappsilber.fdr.entities.PeptidePair;
import org.rappsilber.fdr.entities.Protein;
import org.rappsilber.fdr.entities.ProteinGroupDirectionalLink;
import org.rappsilber.fdr.entities.ProteinGroupDirectionalPair;
import org.rappsilber.fdr.entities.ProteinGroupLink;
import org.rappsilber.fdr.entities.ProteinGroupPair;
import org.rappsilber.fdr.entities.Site;
import org.rappsilber.fdr.entities.ProteinGroup;
import org.rappsilber.fdr.filter.PSMFilter;
import org.rappsilber.fdr.filter.SingleSubScoreFilter;
import org.rappsilber.fdr.result.FDRResult;
import org.rappsilber.fdr.result.FDRResultLevel;
import org.rappsilber.fdr.result.SubGroupFdrInfo;
import org.rappsilber.fdr.entities.AbstractFDRElement;
import org.rappsilber.fdr.entities.FDRSelfAdd;
import org.rappsilber.fdr.filter.DeltaScorePercentFilter;
import org.rappsilber.fdr.utils.HashedArrayList;
import org.rappsilber.fdr.utils.MaximisingStatus;
import org.rappsilber.fdr.utils.MaximizeLevelInfo;
import org.rappsilber.fdr.utils.MaximizeLevelInfoInteger;
import org.rappsilber.fdr.utils.MaximizingUpdate;
import org.rappsilber.fdr.utils.MiscUtils;
import org.rappsilber.utils.AutoIncrementValueMap;
import org.rappsilber.utils.CountOccurence;
import org.rappsilber.utils.DoubleArrayList;
import org.rappsilber.utils.IntArrayList;
import org.rappsilber.utils.RArrayUtils;
import org.rappsilber.utils.NullOutputStream;
import org.rappsilber.utils.SelfAdd;
import org.rappsilber.utils.SelfAddHashSet;
import org.rappsilber.utils.UpdateableInteger;
import org.rappsilber.utils.Version;
import org.rappsilber.utils.statistic.StreamingStatsEstimator;
import rappsilber.utils.MyArrayUtils;

/**
 *
 * @author lfischer
 */
public abstract class OfflineFDR {

    OfflineFDR cofilter;
    FDRResult coResult;
    String[] cofilterArgs;
    
    /**
     * store all psm
     */
    protected SelfAddHashSet<PSM> allPSMs = new SelfAddHashSet<PSM>();

    /**
     * psms that passed some form of prefilter
     */
    protected ArrayList<PSM> prefilteredPSMs = null;

    /**
     * store all peptide pairs
     */
    SelfAddHashSet<Peptide> allPeptides = new SelfAddHashSet<Peptide>();
    /**
     * store all Proteins
     */
    SelfAddHashSet<Protein> allProteins = new SelfAddHashSet<Protein>();
    /**
     * turn peptide-sequences into unique integer ids
     */
    AutoIncrementValueMap<String> m_pepIDs = new AutoIncrementValueMap<String>();
    /**
     * turn protein-accession into unique integer ids
     */
    AutoIncrementValueMap<String> m_protIDs = new AutoIncrementValueMap<String>();

    /**
     * size of the decoy independent protein pairs in terms of psms
     */
    HashMap<Integer, UpdateableInteger> protpairToSize = new HashMap<>();
    /**
     * id of a protein pair independent of target or decoy
     */
    HashMap<String, Integer> protpairToID = new HashMap<>();

    private boolean psm_directional = false;
    private boolean peptides_directional = false;
    private boolean links_directional = false;
    private boolean ppi_directional = false;
    protected int m_maximum_summed_peplength = Integer.MAX_VALUE;
    protected FDRLevel maximizeWhat = null;

    /**
     * is a higher score better than a lower score?
     */
    protected boolean PSMScoreHighBetter = true;
    /**
     * the version of xiFDR to be reported
     */
    private static Version xiFDRVersion;
    private int minPepPerProteinGroup = 1;
    private int minPepPerProteinGroupLink = 1;
    private int minPepPerProteinGroupPair = 1;

    private double targetPepDBSize = 999999999;
    private double decoyPepDBSize = 999999999;
    private double targetProtDBSize = 999999999;
    private double decoyProtDBSize = 999999999;
    private double targetLinkDBSize = 999999999;
    private double decoyLinkDBSize = 999999999;

    private double[] psmFDRSetting;
    private double[] peptidePairFDRSetting;
    private double[] ProteinGroupFDRSetting;
    private double[] linkFDRSetting;
    private double[] ppiFDRSetting;
    private double safetyFactorSetting;
    private boolean ignoreGroupsSetting;
    private boolean csvSummaryOnly = false;
    private boolean singleSummary = false;
    private PrintWriter singleSummaryOut;
    private String csvOutDirSetting;
    private String csvOutBaseSetting;
    public static final Integer MINIMUM_POSSIBLE_DECOY = 1;
    public static final Integer MINIMUM_POSSIBLE_RESULT = 1;
    private Integer m_linearPSMCount = null;
    private Integer m_XLPSMCount = null;
    private Integer m_linearPepCount = null;
    private Integer m_XLPepCount = null;

    public int m_minPepLength = 0;

    public int commandlineFDRDigits = 2;
    private boolean uniquePSMs = true;

    /**
     * We normalize psm-scores by median and MAD - but to go around some quirks
     * of our score propagation scores then get shifted so that the lowest score
     * is around one.
     */
    private double psmNormOffset = 0;

    /**
     * indicates, whether the psms went through a score normalisation
     */
    private boolean isNormalized = false;

    Boolean isNormalizedByDecoy = null;
    /**
     * what is the combinatorial limit of ambiguity, that is accepted for a
     * link.
     * <br/>0 indicates no limit
     */
    protected int m_maximumLinkAmbiguity = 0;
    /**
     * If a ProteinGroup contains more then the given number of proteins, then
     * this group is ignored.
     * <br/>0 indicates that there is no limit.
     * <br/>0 indicates no limit
     */
    protected int m_maximumProteinAmbiguity = 0;
    /**
     * If a ProteinGroupPair spans more then the given number of proteins the
     * pair is ignored
     * <br/>0 indicates that there is no limit.
     * <br/>0 indicates no limit
     */
    protected int m_maximumProteinPairAmbiguity = 0;

    /**
     * group matches by protein pairs
     */
    private boolean groupByProteinPair = false;

    /**
     * how many decoys does a fdr group need to have to be reported as result
     */
    private Integer minTDChance = 0;

    /**
     * I filter the cross-linker names through this hashmap, ensuring I have
     * only one string instance per cross-linker. That way comparison of
     * cross-linker can be reduced to A = B instead of A.equals(B)
     */
    HashMap<String, String> foundCrossLinker = new HashMap<>();
    /**
     * I filter the run names through this hashmap, ensuring I have only one
     * string instance per run. That way comparison of runs can be reduced to A
     * = B instead of A.equals(B)
     */
    HashMap<String, String> foundRuns = new HashMap<>();
    HashMap<String, Integer> runToInt = new HashMap<>();
    ArrayList<String> runs = new ArrayList<>();
    private Locale outputlocale = Locale.getDefault();
    private NumberFormat numberFormat = NumberFormat.getNumberInstance(outputlocale);
//    private String localNumberGroupingSeperator;
//    private String localNumberDecimalSeparator;
//    private String quoteDoubles;
    private boolean stopMaximizing = false;

    private ArrayList<String> extraColumns = new ArrayList<>();
    /**
     * group between links by both proteins beeing observed with self-links
     */
    private boolean groupLinksByHasInternal = false;
    /**
     * group between PPIs by both proteins beeing observed with self-links
     */
    private boolean groupPPIByHasInternal = false;
    /**
     * group between PeptidePairs by both proteins beeing observed with
     * self-links
     */
    private boolean groupPepPairByHasInternal = false;
    /**
     * group psms-by runs
     */
    private boolean groupPSMsByRun = false;

    private FDRSettings settings;

    /**
     * @return the uniquePSMs
     */
    public boolean filterUniquePSMs() {
        return uniquePSMs;
    }

    /**
     * @param uniquePSMs the uniquePSMs to set
     */
    public void setFilterUniquePSMs(boolean uniquePSMs) {
        this.uniquePSMs = uniquePSMs;
    }

//    protected HashMap<FDRLevel,HashMap<Integer,SubGroupFdrInfo>> GroupedFDRs = new HashMap<FDRLevel, HashMap<Integer, SubGroupFdrInfo>>();
    /**
     * An enum providing constants for the FDR-Levels
     */
    public static enum FDRLevel {
        PSM("PSMs"), PEPTIDE_PAIR("Peptidepairs"), PROTEINGROUP("Proteins"), PROTEINGROUPLINK("Residue Pairs"), PROTEINGROUPPAIR("Proteinpairs");
        String m_shortname;

        private FDRLevel() {
        }

        private FDRLevel(String shortName) {
            m_shortname = shortName;
        }

        @Override
        public String toString() {
            return m_shortname == null ? super.toString() : m_shortname;
        }

    }

    public OfflineFDR() {

    }

    public OfflineFDR(int[] peptideLengthGroups) {
        PeptidePair.setLenghtGroup(peptideLengthGroups);
    }

    public void normalizePSMs() {
        int count = allPSMs.size();
        int decoy = 0;
        for (PSM p : allPSMs) {
            if (p.isDecoy()) {
                decoy++;
            }
            if (decoy > count / 5) {
                normalizePSMsByDecoy();
                return;
            }
        }
        normalizePSMsAll();
    }

    public void coNormalizePSMs(OfflineFDR newData) {
        if (!this.isNormalized()) {
            this.normalizePSMsToFDR();
        }

        newData.normalizePSMsToFDR();

    }

    public String getDistributions(String what) {
        StreamingStatsEstimator sse = new StreamingStatsEstimator(0.001);
        what = what.toLowerCase();
        if (what.contentEquals("decoy")) {
            for (PSM p : allPSMs) {
                if (p.isDecoy()) {
                    sse.addValue(p.getScore());
                }
            }
        } else if (what.contentEquals("tt")) {
            for (PSM p : allPSMs) {
                if (!p.isDecoy()) {
                    sse.addValue(p.getScore());
                }
            }
        } else if (what.contentEquals("td")) {
            for (PSM p : allPSMs) {
                if (p.isTD()) {
                    sse.addValue(p.getScore());
                }
            }

        } else if (what.contentEquals("dd")) {
            for (PSM p : allPSMs) {
                if (p.isDD()) {
                    sse.addValue(p.getScore());
                }
            }
        } else {
            for (PSM p : allPSMs) {
                sse.addValue(p.getScore());
            }
        }
        return sse.dumpCSV();
    }

    public void normalizePSMsToFDR() {
        // get all psms
        ArrayList<PSM> scorePSMs = new ArrayList<>(allPSMs);
        scorePSMs.sort(new Comparator<PSM>() {
            @Override
            public int compare(PSM o1, PSM o2) {
                return Double.compare(o1.getScore(), o2.getScore());
            }
        });
        int size = scorePSMs.size();

        HashSet<String> keys = new HashSet<>(scorePSMs.size() / 2);

        boolean[] consider = new boolean[size];
        int tt = 0;
        int td = 0;
        int dd = 0;
        PSM lastTD = scorePSMs.get(size - 1);
        // do some counts and link the better PSM
        for (int p = size - 1; p >= 0; p--) {
            PSM e = scorePSMs.get(p);
            e.setLowerTD(lastTD);
            if (!keys.contains(e.getNonDirectionalUnifyingKey())) {
                keys.add(e.getNonDirectionalUnifyingKey());
                if (e.isTT()) {
                    tt++;
                } else if (e.isTD()) {
                    td++;
                    lastTD = e;
                } else {
                    dd++;
                }
                consider[p] = true;
            } else {
                consider[p] = false;
            }
        }

        double fdr = (td - dd) / (double) tt;
        lastTD = scorePSMs.get(0);
        lastTD.setFDR(fdr);
        // calculate FDR and link the higher
        for (int p = 0; p < size; p++) {
            PSM e = scorePSMs.get(p);
            e.setHigherTD(lastTD);

            // turn the FDR into a score
            //e.setScore(10*(1-e.getEstimatedFDR()));
            if (consider[p]) {
                keys.add(e.getNonDirectionalUnifyingKey());
                if (e.isTT()) {
                    tt--;
                } else if (e.isTD()) {
                    e.setFDR(fdr);
                    td--;
                    lastTD = e;
                } else {
                    dd--;
                }
            }
            double newfdr = (td - dd) / (double) tt;
            if (newfdr < fdr) {
                fdr = newfdr;
            }
        }
        for (int p = 0; p < size; p++) {
            PSM e = scorePSMs.get(p);
            // turn the FDR into a score
            e.setScore(10 * (1 - e.getEstimatedFDR()));
        }

        setOveralNormalized();
    }

    public void normalizePSMsByDecoy() {
        if (allPSMs.size() == 0) {
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Supposedly no PSMs here");
        }
        double minScore = Double.MAX_VALUE;

        StreamingStatsEstimator sse = new StreamingStatsEstimator(0.001);
        for (PSM p : allPSMs) {
            if (p.isDecoy()) {
                sse.addValue(p.getScore());
            }
            if (minScore > p.getScore()) {
                minScore = p.getScore();
            }
        }

        double mode = sse.getModeEstimation();

        // calculate MAD but as deviade from mode
        double mad = sse.getMADEstimation(mode);

        // xifdr is not particular keen on negative scores - so we define an offset
        double offset = -(minScore - mode) / mad;
        psmNormOffset = offset;
        for (PSM p : allPSMs) {
            p.setScore((p.getScore() - mode) / mad + offset);
        }

        setDecoyNormalized();
    }

    public void normalizePSMsAll() {
        if (allPSMs.size() == 0) {
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Supposedly no PSMs here");
        }

        StreamingStatsEstimator sse = new StreamingStatsEstimator(0.001);
        for (PSM p : allPSMs) {
            sse.addValue(p.getScore());
        }

        double mode = sse.getModeEstimation();

        // calculate MAD but as deviade from mode
        double mad = sse.getMADEstimation(mode);

        // xifdr is not particular keen on negative scores - so we define an offset
        double offset = -(sse.getMin() - mode) / mad;
        psmNormOffset = offset;
        for (PSM p : allPSMs) {
            p.setScore((p.getScore() - mode) / mad + offset);
        }

        setOveralNormalized();
    }

    /**
     * add a list of independently normalised psms to the current list of psms
     *
     * @param psms list of psms
     * @param offset the offset applied to this list
     */
    public void addNormalisedPsmList(SelfAddHashSet<PSM> psms, double offset) {
        double offsetdiff = psmNormOffset - offset;
        //make sure we shift the normalised median to the same place for both lists
        if (offsetdiff < 0) {
            for (PSM p : allPSMs) {
                p.setScore(p.getScore() - offsetdiff);
            }
            this.psmNormOffset -= offsetdiff;
            for (PSM p : psms) {
                p.setRun(registerRun(p.getRun()));
                if (p.getCrosslinker() != null) {
                    registerCrossLinker(p.getCrosslinker(), p);
                }
                allPSMs.add(p);
            }
        } else if (offsetdiff > 0) {
            for (PSM p : psms) {
                p.setScore(p.getScore() + offsetdiff);
                p.setRun(registerRun(p.getRun()));
                if (p.getCrosslinker() != null) {
                    registerCrossLinker(p.getCrosslinker(), p);
                }
                allPSMs.add(p);
            }
        }

//        allPSMs.addAll(psms);
    }

    /**
     * add a list of independently normalised psms to the current list of psms
     *
     * @param psms list of psms
     * @param offset the offset applied to this list
     */
    public void normaliseAndAddPsmList(OfflineFDR other) {
        coNormalizePSMs(this);
        allPSMs.addAll(other.allPSMs);
    }

    protected void levelSummary(PrintWriter summaryOut, String pepheader, FDRResultLevel level, String seperator) {

        double target_fdr = level.getTargetFDR();
        summaryOut.println("\n\"" + pepheader + "\"");
        summaryOut.print("\"Group\"");
        ArrayList<String> fdrGroups = new ArrayList<String>(level.getGroupIDs());
        java.util.Collections.sort(fdrGroups);
        for (String fg : fdrGroups) {
            summaryOut.print(seperator + "\"" + fg + "\"");
        }
        summaryOut.print("\n\"Input\"");
        for (String fg : fdrGroups) {
            summaryOut.print(seperator + (((SubGroupFdrInfo) level.getGroup(fg)).inputCount));
        }
        summaryOut.print("\n\"TT\"");
        for (String fg : fdrGroups) {
            summaryOut.print(seperator + ((SubGroupFdrInfo) level.getGroup(fg)).TT);
        }
        summaryOut.print("\n\"TD\"");
        for (String fg : fdrGroups) {
            summaryOut.print(seperator + ((SubGroupFdrInfo) level.getGroup(fg)).TD);
        }
        summaryOut.print("\n\"DD\"");
        for (String fg : fdrGroups) {
            summaryOut.print(seperator + ((SubGroupFdrInfo) level.getGroup(fg)).DD);
        }
        summaryOut.print("\n\"passing fdr (" + target_fdr + ")\"");
        for (String fg : fdrGroups) {
            summaryOut.print(seperator + ((SubGroupFdrInfo) level.getGroup(fg)).results.size());
        }
        summaryOut.print("\n\"TT\"");
        for (String fg : fdrGroups) {
            summaryOut.print(seperator + ((SubGroupFdrInfo) level.getGroup(fg)).resultTT);
        }
        summaryOut.print("\n\"TD\"");
        for (String fg : fdrGroups) {
            summaryOut.print(seperator + ((SubGroupFdrInfo) level.getGroup(fg)).resultTD);
        }
        summaryOut.print("\n\"DD\"");
        for (String fg : fdrGroups) {
            summaryOut.print(seperator + ((SubGroupFdrInfo) level.getGroup(fg)).resultDD);
        }
        summaryOut.print("\n\"last fdr > " + target_fdr + "\"");
        for (String fg : fdrGroups) {
            summaryOut.print(seperator + ((SubGroupFdrInfo) level.getGroup(fg)).firstPassingFDR);
        }
        summaryOut.print("\n\"higher fdr (> " + target_fdr + ")\"");
        for (String fg : fdrGroups) {
            summaryOut.print(seperator + ((SubGroupFdrInfo) level.getGroup(fg)).higherFDR);
        }
        summaryOut.print("\n\"lower fdr (<= " + target_fdr + ")\"");
        for (String fg : fdrGroups) {
            summaryOut.print(seperator + ((SubGroupFdrInfo) level.getGroup(fg)).lowerFDR);
        }
        summaryOut.print("\nfinal");
        for (String fg : fdrGroups) {
            summaryOut.print(seperator + ((SubGroupFdrInfo) level.getGroup(fg)).filteredResult.size());
        }
        summaryOut.print("\n");
    }

    /**
     * @return the peptides_directional
     */
    public boolean isPeptides_directional() {
        return peptides_directional;
    }

    /**
     * @param peptides_directional the peptides_directional to set
     */
    public void setPeptides_directional(boolean peptides_directional) {
        this.peptides_directional = peptides_directional;
    }

    /**
     * @return the links_directional
     */
    public boolean isLinks_directional() {
        return links_directional;
    }

    /**
     * @param links_directional the links_directional to set
     */
    public void setLinks_directional(boolean links_directional) {
        this.links_directional = links_directional;
    }

    /**
     * @return the ppi_directional
     */
    public boolean isPpi_directional() {
        return ppi_directional;
    }

    /**
     * @param ppi_directional the ppi_directional to set
     */
    public void setPpi_directional(boolean ppi_directional) {
        this.ppi_directional = ppi_directional;
    }

    /**
     * @return the psm_directional
     */
    public boolean isPsm_directional() {
        return psm_directional;
    }

    /**
     * @param psm_directional the psm_directional to set
     */
    public void setPsm_directional(boolean psm_directional) {
        this.psm_directional = psm_directional;
    }

    public static String getLongVersionString() {
        return xiFDRVersion.toLongString();
    }

    public void setLengthGroups(int[] peptideLengthGroups) {
        PeptidePair.setLenghtGroup(peptideLengthGroups);
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
     * @return
     */
    public PSM addMatch(String psmID, String pepSeq1, String pepSeq2, int peplen1, int peplen2, int site1, int site2, boolean isDecoy1, boolean isDecoy2, int charge, double score, String accession1, String description1, String accession2, String description2, int pepPosition1, int pepPosition2, double peptide1score, double peptide2score, String isSpecialCase) {
        return addMatch(psmID, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, accession1, description1, accession2, description2, pepPosition1, pepPosition2, peptide1score, peptide2score, isSpecialCase, "", "", "");
    }

    public PSM addMatch(String psmID, String pepSeq1, String pepSeq2, int peplen1, int peplen2, int site1, int site2, boolean isDecoy1, boolean isDecoy2, int charge, double score, String accession1, String description1, String accession2, String description2, int pepPosition1, int pepPosition2, double peptide1score, double peptide2score, String isSpecialCase, String crosslinker, String run, String scan) {

        long pepid1 = m_pepIDs.toIntValue(pepSeq1);
        long pepid2 = m_pepIDs.toIntValue(pepSeq2);
        long protid1 = m_protIDs.toIntValue(accession1);
        long protid2 = m_protIDs.toIntValue(accession2);

        //return addMatch(pepSeq2, pepSeq1, accession1, accession2, protid1, description2, isDecoy1, pepid1, pepPosition1, peplen1, protid2, isDecoy2, pepid2, pepPosition2, peplen2, psmID, site1, site2, charge, score, scoreRatio, isSpecialCase);
        return addMatch(psmID, pepid1, pepid2, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, protid1, accession1, description1, protid2, accession2, description2, pepPosition1, pepPosition2, "", "", peptide1score, peptide2score, isSpecialCase, crosslinker, run, scan);
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
    public PSM addMatch(String psmID, Peptide peptide1, Peptide peptide2, int peplen1, int peplen2, int site1, int site2, int charge, double score, Protein proteinId1, Protein proteinId2, int pepPosition1, int pepPosition2, double peptide1Score, double peptide2Score, String isSpecialCase, String crosslinker, String run, String Scan) {
        Peptide npepid1;
        Peptide npepid2;
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
        //double nScoreRatio = scoreRation;
        double npeptide1Score;
        double npeptide2Score;

        if (protcomp < 0 || (protcomp == 0 && pepcomp < 0) || (protcomp == 0 && pepcomp == 0 && site1 < site2)) {
            npepid1 = peptide1;
            npepid2 = peptide2;
            npeplen1 = peplen1;
            npeplen2 = peplen2;
            nsite1 = site1;
            nsite2 = site2;
            nproteinId1 = proteinId1;
            nproteinId2 = proteinId2;
            npepPosition1 = pepPosition1;
            npepPosition2 = pepPosition2;
            npeptide1Score = peptide1Score;
            npeptide2Score = peptide2Score;

        } else {
            npepid1 = peptide2;
            npepid2 = peptide1;
            npeplen1 = peplen2;
            npeplen2 = peplen1;
            nsite1 = site2;
            nsite2 = site1;
            nproteinId1 = proteinId2;
            nproteinId2 = proteinId1;
            npepPosition1 = pepPosition2;
            npepPosition2 = pepPosition1;
            npeptide1Score = peptide2Score;
            npeptide2Score = peptide1Score;
        }

        if (!PSMScoreHighBetter) {
            score = 10 - (10 * score);
        }

        PSM psm = new PSM(psmID, npepid1, npepid2, (byte) nsite1, (byte) nsite2, nproteinId1.isDecoy(), nproteinId2.isDecoy(), (byte) charge, score, npeptide1Score, npeptide2Score);
        psm.setNegativeGrouping(isSpecialCase);

        psm.setRun(registerRun(run));
        if (crosslinker == null) {
            crosslinker = "";
        }

        registerCrossLinker(crosslinker, psm);
        psm.setScan(Scan);

        PSM regpsm = getAllPSMs().register(psm);

        return regpsm;
    }

    protected void registerCrossLinker(String crosslinker, PSM psm) {
        String c = foundCrossLinker.get(crosslinker);
        if (c == null) {
            psm.setCrosslinker(crosslinker);
            foundCrossLinker.put(crosslinker, crosslinker);
        } else {
            psm.setCrosslinker(c);
        }
    }

    protected String registerRun(String run) {
        // ensure we have just a single instance of a string for each cross-linker and run
        // speeds up comparisons later
        String r = foundRuns.get(run);
        if (r == null) {
            r = run;
            foundRuns.put(run, run);
            runToInt.put(run, runs.size());
            runs.add(run);
        }
        return r;
    }

    public FDRResult calculateWriteFDR(String path, String baseName, String seperator, FDRSettings settings) throws FileNotFoundException {
        return calculateWriteFDR(path, baseName, seperator, commandlineFDRDigits, settings);
    }

    public FDRResult calculateWriteFDR(String path, String baseName, String seperator, int minDigits, FDRSettings settings) throws FileNotFoundException {
        FDRResult result = null;
        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "PATH: " + path);
        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "BaseName: " + baseName);
        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Seperator: " + seperator);

        DoubleArrayList allvalues = new DoubleArrayList();
        allvalues.addAll(getPsmFDRSetting());
        allvalues.addAll(getPeptidePairFDRSetting());
        allvalues.addAll(getProteinGroupFDRSetting());
        allvalues.addAll(getLinkFDRSetting());
        allvalues.addAll(getPpiFDRSetting());

        String format = MiscUtils.formatStringForPrettyPrintingRelatedValues(allvalues.toDoubleArray(), minDigits);

        for (double psmfdr = Math.round(getPsmFDRSetting()[0] * 1000000); psmfdr <= Math.round(getPsmFDRSetting()[1] * 1000000); psmfdr += Math.round(getPsmFDRSetting()[2] * 1000000)) {
            for (double pepfdr = Math.round(getPeptidePairFDRSetting()[0] * 1000000); pepfdr <= Math.round(getPeptidePairFDRSetting()[1] * 1000000); pepfdr += Math.round(getPeptidePairFDRSetting()[2] * 1000000)) {
                for (double pgfdr = Math.round(getProteinGroupFDRSetting()[0] * 1000000); pgfdr <= Math.round(getProteinGroupFDRSetting()[1] * 1000000); pgfdr += Math.round(getProteinGroupFDRSetting()[2] * 1000000)) {
                    for (double pglfdr = Math.round(getLinkFDRSetting()[0] * 1000000); pglfdr <= Math.round(getLinkFDRSetting()[1] * 1000000); pglfdr += Math.round(getLinkFDRSetting()[2] * 1000000)) {
                        for (double pgpfdr = Math.round(getPpiFDRSetting()[0] * 1000000); pgpfdr <= Math.round(getPpiFDRSetting()[1] * 1000000); pgpfdr += Math.round(getPpiFDRSetting()[2] * 1000000)) {

                            String fdr_basename;
                            if (getPsmFDRSetting()[0] == getPsmFDRSetting()[1]
                                    && getPeptidePairFDRSetting()[0] == getPeptidePairFDRSetting()[1]
                                    && getProteinGroupFDRSetting()[0] == getProteinGroupFDRSetting()[1]
                                    && getLinkFDRSetting()[0] == getLinkFDRSetting()[1]
                                    && getPpiFDRSetting()[0] == getPpiFDRSetting()[0]) {
                                fdr_basename = baseName;

                            } else {
                                fdr_basename = baseName + "_" + String.format(format, psmfdr / 1000000.0) + "_"
                                        + String.format(format, pepfdr / 1000000.0) + "_"
                                        + String.format(format, pgfdr / 1000000.0) + "_"
                                        + String.format(format, pglfdr / 1000000.0) + "_"
                                        + String.format(format, pgpfdr / 1000000.0) + "_"
                                        + String.format(format, getSafetyFactorSetting()) + "_"
                                        + isIgnoreGroupsSetting();
                            }

                            String testfile = path + "/" + fdr_basename + "_summary.";

                            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "next " + fdr_basename);

                            if (!(new File(testfile + "csv").exists() || new File(testfile + "txt").exists())) {
                                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Calculating:\n"
                                        + "\npsmFDR:              " + String.format(format, psmfdr / 1000000.0)
                                        + "\nPeptidePairFDR:      " + String.format(format, pepfdr / 1000000.0)
                                        + "\nProteinGroupFDR:     " + String.format(format, pgfdr / 1000000.0)
                                        + "\nProteinGroupLinkFDR: " + String.format(format, pglfdr / 1000000.0)
                                        + "\nProteinGroupPairFDR: " + String.format(format, pgpfdr / 1000000.0)
                                        + "\nReport-Factor:       " + String.format(format, getSafetyFactorSetting())
                                        + "\nIgnore Groups:       " + isIgnoreGroupsSetting());
                                FDRSettingsImpl s = new FDRSettingsImpl();
                                s.setAll(settings);
                                s.BoostingSteps = 4;
                                s.PSMFDR = psmfdr / 1000000;
                                s.PeptidePairFDR = pepfdr / 1000000;
                                s.ProteinGroupFDR = pgfdr / 1000000;
                                s.ProteinGroupLinkFDR = pglfdr / 1000000;
                                s.ProteinGroupPairFDR = pgpfdr / 1000000;
                                if (maximizeWhat == null) {
                                    result = this.calculateFDR(s, true);
                                } else {

                                    MaximisingStatus m = this.maximise(s, maximizeWhat, uniquePSMs, new MaximizingUpdate() {
                                        @Override
                                        public void setStatus(MaximisingStatus state) {

                                        }

                                        @Override
                                        public void setStatusText(String text) {
                                            System.err.println(text);
                                        }

                                        @Override
                                        public void reportError(String text, Exception ex) {
                                            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, text, ex);
                                            System.exit(-1);
                                        }
                                    });
                                    result = m.result;
                                }

                                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "PATH: " + path);
                                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "fdr_basename: " + fdr_basename);
                                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Seperator: " + seperator);

                                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Result-summary:" + summaryString(result));

                                writeFiles(path, fdr_basename, seperator, result);
                            } else {
                                Logger.getLogger(this.getClass().getName()).log(Level.INFO, fdr_basename + "skipped");

                            }

                        }

                    }

                }

            }

        }
        return result;
    }

    public void calculatePSMFDR(boolean setElementFDR, boolean ignoreGroups, FDRResult result, FDRSettings settings) {
        boolean groupByHasInternal = false;
        //boolean groupByCharge4plus = true;
        FDRResultLevel<PSM> GroupedFDRs = new FDRResultLevel<PSM>();
        GroupedFDRs.isDirectional = false;
        FDRResultLevel<PSM> GroupedFDRsS = new FDRResultLevel<PSM>();
        GroupedFDRsS.isDirectional = false;
        reset();
        result.uniquePSMs = settings.filterToUniquePSM();

        protpairToID = new HashMap<>();
        protpairToSize = new HashMap<>();
        Collection<PSM> inputPSM;
        ArrayList<PSM> allPSM = null;
        if (prefilteredPSMs == null) {
            allPSM = new ArrayList<>(getAllPSMs());
        } else {
            allPSM = new ArrayList<>(prefilteredPSMs);
        }

        if (settings.getMinPeptideFragmentsFilter() > 0) {
            PSMFilter f = new SingleSubScoreFilter("MinFragments", settings.getMinPeptideFragmentsFilter(), true);
            allPSM = f.filter(allPSM);
        }

        if (settings.getMinPeptideCoverageFilter() > 0) {
            PSMFilter f = new SingleSubScoreFilter("minPepCoverage", settings.getMinPeptideCoverageFilter(), true);
            allPSM = f.filter(allPSM);
        }

        if (settings.getMinDeltaScoreFilter() > 0) {
            PSMFilter f = new DeltaScorePercentFilter(settings.getMinDeltaScoreFilter());
            allPSM = f.filter(allPSM);
        }

        if (settings.combineScoreAndDelta()) {
            for (PSM p : allPSM) {
                p.setScore((p.getOriginalScore() + p.getDeltaScore()) / 2);
            }
        }

        if (m_maximumProteinAmbiguity > 0 && m_minPepLength == 0) {
            inputPSM = new ArrayList<PSM>(allPSM.size());
            for (PSM p : allPSM) {
                if (p.getScore() > 0 && p.getPeptide1().getProteins().size() <= m_maximumProteinAmbiguity
                        && p.getPeptide2().getProteins().size() <= m_maximumProteinAmbiguity) {
                    inputPSM.add(p);
                }
            }
        } else if (m_maximumProteinAmbiguity > 0 && m_minPepLength > 0) {
            inputPSM = new ArrayList<PSM>(allPSM.size());
            for (PSM p : allPSM) {
                Peptide pep1 = p.getPeptide1();
                Peptide pep2 = p.getPeptide2();

                if (p.getScore() > 0 && pep1.getProteins().size() <= m_maximumProteinAmbiguity
                        && pep2.getProteins().size() <= m_maximumProteinAmbiguity
                        && (pep1 == Peptide.NOPEPTIDE || pep1.length() >= m_minPepLength || (pep1.length() == 1 && pep1.getSequence().startsWith("X")))
                        && (pep2 == Peptide.NOPEPTIDE || pep2.length() >= m_minPepLength || (pep2.length() == 1 && pep2.getSequence().startsWith("X")))) {
                    inputPSM.add(p);
                }
            }
        } else if (m_minPepLength > 0) {
            inputPSM = new ArrayList<PSM>(allPSM.size());
            for (PSM p : allPSM) {
                if (p.getScore() > 0) {
                    Peptide pep1 = p.getPeptide1();
                    Peptide pep2 = p.getPeptide2();

                    if ((pep1 == Peptide.NOPEPTIDE || pep1.length() >= m_minPepLength || (pep1.length() == 1 && pep1.getSequence().startsWith("X")))
                            && (pep2 == Peptide.NOPEPTIDE || pep2.length() >= m_minPepLength || (pep2.length() == 1 && pep2.getSequence().startsWith("X")))) {
                        inputPSM.add(p);
                    }
                }
            }
        } else {
            inputPSM = new ArrayList<PSM>(allPSM.size());
            for (PSM p : allPSM) {
                if (p.getScore() > 0) {
                    inputPSM.add(p);
                }
            }
        }

//        if (filterUnique) {
//            SelfAddHashSet<PSM> uniquePSM = new SelfAddHashSet<PSM>();
//            for (PSM psm : inputPSM) 
//                uniquePSM.add(new  UniquePSM(psm));
//            inputPSM = new ArrayList<PSM>(uniquePSM);
//        }
        // filter to unique PSMs
        if (settings.filterToUniquePSM()) {
            HashMap<String, PSM> uniquePSM = filterPSMToUnique(inputPSM);
            inputPSM = new ArrayList<PSM>(uniquePSM.values());
        }

        if (settings.filterConsecutivePeptides()) {
            ArrayList<PSM> nonconsecutives = new ArrayList<>();
            for (PSM psm : inputPSM) {
                if (!psm.isConsecutive()) {
                    nonconsecutives.add(psm);
                }
            }
            inputPSM = nonconsecutives;
        }

//        for (PSM pp : inputPSM) {
//            pp.setFDRGroup(pp.getFDRGroup()+" z"+pp.getCharge());
//        }
        result.input = inputPSM;
        result.minPeptideLength = m_minPepLength;
        result.maximumProteinAmbiguity = m_maximumProteinAmbiguity;
        result.maximumLinkAmbiguity = m_maximumLinkAmbiguity;

        fdr(settings.getPSMFDR(), settings.getReportFactor(), inputPSM, GroupedFDRs, targetPepDBSize, decoyPepDBSize, 1, ignoreGroups, setElementFDR, settings.isPSMDirectional(), settings.isGroupByPSMCount(), settings.getMinTD());

        if (groupByHasInternal) {

            HashSet<ProteinGroup> internal = new HashSet<>();
            for (PSM pp : GroupedFDRs) {
                if (pp.isInternal()) {
                    internal.add(pp.getProteinGroup1());
                    internal.add(pp.getProteinGroup1().decoyComplement());
                    internal.add(pp.getProteinGroup2());
                    internal.add(pp.getProteinGroup2().decoyComplement());
                }
            }

            for (PSM pp : inputPSM) {
                if (pp.isBetween()
                        && internal.contains(pp.getProteinGroup1())
                        && internal.contains(pp.getProteinGroup2())) {
                    pp.setFDRGroup(pp.getFDRGroup() + " int_support");
                }
            }
            fdr(settings.getPSMFDR(), settings.getReportFactor(), inputPSM, GroupedFDRs, targetPepDBSize, decoyPepDBSize, 1, ignoreGroups, setElementFDR, settings.psmLocalFDR(), settings.isGroupByPSMCount(), settings.getMinTD());
            GroupedFDRs = GroupedFDRsS;
        }

        result.psmFDR = GroupedFDRs;

//        return GroupedFDRs;
    }

    protected HashMap<String, PSM> filterPSMToUnique(Collection<PSM> inputPSM) {
        HashMap<String, PSM> uniquePSM = new HashMap<String, PSM>();
        for (PSM psm : inputPSM) {
            String key = psm.getPeptide1().getSequence() + "_!xl!_"
                    + psm.getPeptide2().getSequence() + "_!xl!_"
                    + psm.getPeptideLinkSite1() + "_!xl!_"
                    + psm.getPeptideLinkSite2() + "_!xl!_"
                    + psm.getCharge();
            // do we have already something under this key?
            PSM stored = uniquePSM.get(key);
            if (stored != null) {
                // yes
                if (stored.getScore() < psm.getScore()) {
                    // new psm has higher score -> it will be the reprentative one
                    uniquePSM.put(key, psm);
                    // register the previous psm
                    psm.represents(stored);
                } else {
                    // lower score so it will be represented by the stored one
                    stored.represents(psm);
                }
            } else if (!isPsm_directional()) {
                // we need to check the inverse key as well
                key = psm.getPeptide2().getSequence() + "_!xl!_"
                        + psm.getPeptide1().getSequence() + "_!xl!_"
                        + psm.getPeptideLinkSite2() + "_!xl!_"
                        + psm.getPeptideLinkSite1() + "_!xl!_"
                        + psm.getCharge();
                stored = uniquePSM.get(key);
                if (stored != null) {
                    if (stored.getScore() < psm.getScore()) {
                        // new psm has higher score -> it will be the reprentative one
                        uniquePSM.put(key, psm);
                        // register the previous psm
                        psm.represents(stored);
                    } else {
                        // lower score so it will be represented by the stored one
                        stored.represents(psm);
                    }
                } else {
                    uniquePSM.put(key, psm);
                }
            } else {
                uniquePSM.put(key, psm);
            }
        }
        return uniquePSM;
    }

    //public void calculatePeptidePairFDR(double fdr, double safetyFactor, boolean ignoreGroups, boolean setElementFDR, FDRResult result, boolean directional) {
    public void calculatePeptidePairFDR(boolean setElementFDR, FDRResult result, FDRSettings settings, boolean ignoreGroups) {

        if (result.psmFDR != null) {
            for (PSM pp : result.psmFDR) {
                pp.setFdrPeptidePair(null);
                pp.setFdrProteinGroup(null);
            }
        }

        if (result.peptidePairFDR != null) {
            for (PeptidePair pp : result.peptidePairFDR) {
                pp.setFdrLink(null);
                pp.setFdrProteinGroup(null);
            }
            result.peptidePairFDR.clear();
        }

        if (result.proteinGroupLinkFDR != null) {
            for (ProteinGroupLink l : result.proteinGroupLinkFDR) {
                l.setFdrPPI(null);
            }
            result.proteinGroupLinkFDR.clear();
        }
        if (result.proteinGroupFDR != null) {
            result.proteinGroupFDR.clear();
        }
        if (result.proteinGroupPairFDR != null) {
            result.proteinGroupPairFDR.clear();
        }

        FDRResultLevel<PSM> psms = result.psmFDR;
        FDRResultLevel<PeptidePair> GroupedFDRs = new FDRResultLevel<PeptidePair>();
        GroupedFDRs.isDirectional = settings.isPeptidePairDirectional();

        SelfAddHashSet<PeptidePair> psmPeps = new SelfAddHashSet<PeptidePair>();
        if (!GroupedFDRs.isDirectional) {
            for (PSM psm : psms) {
                PeptidePair pp = psmPeps.register(psm.getPeptidePair());
            }
        } else {
            for (PSM psm : psms) {
                DirectionalPeptidePair dpp = new DirectionalPeptidePair(psm);
                psmPeps.register(dpp);
            }
        }
        if (groupPepPairByHasInternal) {
            HashSet<ProteinGroup> internal = new HashSet<>();
            for (PeptidePair pp : psmPeps) {
                if (pp.isInternal()) {
                    internal.add(pp.getProteinGroup1());
                    internal.add(pp.getProteinGroup1().decoyComplement());
                    internal.add(pp.getProteinGroup2());
                    internal.add(pp.getProteinGroup2().decoyComplement());
                }
            }

            for (PeptidePair pp : psmPeps) {
                if (pp.isBetween()
                        && internal.contains(pp.getProteinGroup1())
                        && internal.contains(pp.getProteinGroup2())) {
                    pp.setFDRGroup(pp.getFDRGroup() + " int_support");
                }
            }
        }

        fdr(settings.getPeptidePairFDR(), settings.getReportFactor(), psmPeps, GroupedFDRs, targetPepDBSize, decoyPepDBSize, 1, ignoreGroups, setElementFDR, settings.peppairLocalFDR(), settings.isGroupByPSMCount(), settings.getMinTD());

        result.peptidePairFDR = GroupedFDRs;
    }

    public void calculateProteinGroupFDR(boolean ignoreGroups, boolean setElementFDR, FDRSettings settings, FDRResult result) {

        FDRResultLevel<PeptidePair> peps = result.peptidePairFDR;
        FDRResultLevel<ProteinGroup> GroupedFDRs = new FDRResultLevel<ProteinGroup>();
        GroupedFDRs.isDirectional = false;

//        protFDRGroupsInput = new HashMap<Integer, Integer>();
//        nextFdrProteinGroup = new HashMap<Integer, Double>();
//        countFdrProteinGroup = new HashMap<Integer, Integer>();
        SelfAddHashSet<ProteinGroup> pepProteinGroups = new SelfAddHashSet<ProteinGroup>();
//        joinSubFDRInfos(GroupedFDRs, true);
//        SubGroupFdrInfo<PeptidePair> joined = peps.get(-1);
        CountOccurence<String> fdrgroups = new CountOccurence<String>();

        int maxAmbiguity = settings.getMaxProteinAmbiguity();
        if (maxAmbiguity == 0) {
            for (PeptidePair pp : peps) {
                Peptide p = pp.getPeptide1();
                if (p != Peptide.NOPEPTIDE) {
                    ProteinGroup pg = p.getProteinGroup();
                    pg.addPeptidePair(pp);
                    pepProteinGroups.register(pg);
                }
                p = pp.getPeptide2();

                if (p != Peptide.NOPEPTIDE) {
                    ProteinGroup pg = p.getProteinGroup();
                    pg.addPeptidePair(pp);
                    pepProteinGroups.register(pg);
                }

            }
        } else {
            for (PeptidePair pp : peps) {

                Peptide p = pp.getPeptide1();
                if (p != Peptide.NOPEPTIDE) {
                    ProteinGroup pg = p.getProteinGroup();
                    if (pg.size() <= maxAmbiguity) {
                        pg.addPeptidePair(pp);
                        pepProteinGroups.register(pg);
                    }
                }
                p = pp.getPeptide2();
                if (p != Peptide.NOPEPTIDE) {
                    ProteinGroup pg = p.getProteinGroup();
                    if (pg.size() <= maxAmbiguity) {
                        pg.addPeptidePair(pp);
                        pepProteinGroups.register(pg);
                    }
                }
            }
        }

//        double tCountMod = targetProtDBSize;
//        tCountMod = -0.5 - Math.sqrt(1 + 8 * tCountMod) / 2;
//        double dCountMod = decoyProtDBSize;
//        dCountMod = dCountMod / tCountMod;
        for (ProteinGroup pg : pepProteinGroups) {
            fdrgroups.add(pg.getFDRGroup());
        }

        if (pepProteinGroups.size() < 10 && settings.getProteinGroupFDR() < 1) {
            result.proteinGroupFDR = GroupedFDRs;
            return;
        }

        if (settings.getMinProteinPepCount() > 1) {
            filterListByPeptideSuport(pepProteinGroups, settings.getMinProteinPepCount());
        }
        //Logger.getLogger(this.getClass().getName()).log(Level.INFO, "ProteinGroup fdr " + pepProteinGroups.size() + " Groups as Input.");
        fdr(settings.getProteinGroupFDR(), settings.getReportFactor(), pepProteinGroups, GroupedFDRs, targetProtDBSize, decoyProtDBSize, settings.getMinProteinPepCount(), ignoreGroups, setElementFDR, settings.protLocalFDR(), false, settings.getMinTD());

        result.proteinGroupFDR = GroupedFDRs;
//        fdrProteinGroups = fdr(fdr, safetyFactor, pepProteinGroups, nextFdrProteinGroup, protFDRGroupsInput, countFdrProteinGroup, tCountMod, dCountMod, minPepCount, ignoreGroups, setElementFDR);

    }

    public void calculateLinkFDR(boolean ignoreGroups, boolean setElementFDR, FDRSettings settings, FDRResult result) {

        if (result.proteinGroupLinkFDR != null) {
            for (ProteinGroupLink l : result.proteinGroupLinkFDR) {
                l.setFdrPPI(null);
            }
            result.proteinGroupLinkFDR.clear();
        }
        if (result.proteinGroupPairFDR != null) {
            result.proteinGroupPairFDR.clear();
        }
//        linkFDRGroupsInput = new HashMap<Integer, Integer>();
//        nextFdrLink = new HashMap<Integer, Double>();
//        countFdrLink = new HashMap<Integer, Integer>();
        FDRResultLevel<ProteinGroupLink> GroupedFDRs = new FDRResultLevel<ProteinGroupLink>();
        GroupedFDRs.isDirectional = settings.isLinkDirectional();
        SelfAddHashSet<ProteinGroupLink> pepLinks = new SelfAddHashSet<ProteinGroupLink>();
        int maxAmbiguity = settings.getMaxLinkAmbiguity();
        if (settings.isLinkDirectional()) {
            if (maxAmbiguity == 0) {
                for (PeptidePair pp : result.peptidePairFDR) {

                    if (!(pp.isLinear() || pp.isNonCovalent()) || pp.isLoop()) {
                        ProteinGroupDirectionalLink dl = new ProteinGroupDirectionalLink(pp);
                        pp.setFdrLink(pepLinks.register(dl));
                    }

                }
            } else {

                for (PeptidePair pp : result.peptidePairFDR) {

                    if (!(pp.isLinear() || pp.isNonCovalent()) || pp.isLoop()) {
                        ProteinGroupDirectionalLink dl = new ProteinGroupDirectionalLink(pp);
                        if (dl.getAmbiguity() <= maxAmbiguity) {
                            pp.setFdrLink(pepLinks.register(dl));
                        }
                    }

                }

            }

        } else {
            if (maxAmbiguity == 0) {
                for (PeptidePair pp : result.peptidePairFDR) {

                    if (!(pp.isLinear() || pp.isNonCovalent())) {
                        pp.setFdrLink(pepLinks.register(pp.getLink()));
                    }

                }
            } else {

                for (PeptidePair pp : result.peptidePairFDR) {

                    if (!(pp.isLinear() || pp.isNonCovalent()) || pp.isLoop()) {
                        ProteinGroupLink l = pp.getLink();
                        if (l.getAmbiguity() <= maxAmbiguity) {
                            pp.setFdrLink(pepLinks.register(pp.getLink()));
                        }
                    }

                }

            }
        }

        fdr(settings.getProteinGroupLinkFDR(), settings.getReportFactor(), pepLinks, GroupedFDRs, targetLinkDBSize, decoyLinkDBSize, settings.getMinLinkPepCount(), ignoreGroups, setElementFDR, settings.linkLocalFDR(), settings.isGroupByPSMCount(), settings.getMinTD());
        if (groupLinksByHasInternal) {

            HashSet<ProteinGroup> internal = new HashSet<>();
            for (ProteinGroupLink pgl : GroupedFDRs) {
                if (pgl.isInternal()) {
                    internal.add(pgl.getProteinGroup1());
                    internal.add(pgl.getProteinGroup1().decoyComplement());
                    internal.add(pgl.getProteinGroup2());
                    internal.add(pgl.getProteinGroup2().decoyComplement());
                }
            }

            for (ProteinGroupLink pgl : pepLinks) {
                if (pgl.isBetween()
                        && internal.contains(pgl.getProteinGroup1())
                        && internal.contains(pgl.getProteinGroup2())) {
                    pgl.setFDRGroup(pgl.getFDRGroup() + " int_support");
                }
            }
            FDRResultLevel<ProteinGroupLink> GroupedFDRs2 = new FDRResultLevel<ProteinGroupLink>();
            GroupedFDRs2.isDirectional = settings.isLinkDirectional();
            fdr(settings.getProteinGroupLinkFDR(), settings.getReportFactor(), pepLinks, GroupedFDRs2, targetLinkDBSize, decoyLinkDBSize, settings.getMinLinkPepCount(), ignoreGroups, setElementFDR, settings.linkLocalFDR(), settings.isGroupByPSMCount(), settings.getMinTD());
            GroupedFDRs = GroupedFDRs2;
        }

        result.proteinGroupLinkFDR = GroupedFDRs;

//        fdrProtainGroupLinks = fdr(fdr, safetyFactor, pepLinks, nextFdrLink, linkFDRGroupsInput, countFdrLink, targetPepDBSize, decoyPepDBSize, minPepCount, ignoreGroups, setElementFDR);
    }

    protected void filterListByPeptideSuport(Collection<? extends AbstractFDRElement> list, int min) {
        ArrayList<AbstractFDRElement> remove = new ArrayList<>();
        for (AbstractFDRElement element : list) {
            if (element.getPeptidePairs().size() < min) {
                remove.add(element);
            }
        }
        list.removeAll(remove);
    }

    public void calculateProteinGroupPairFDR(boolean ignoreGroups, boolean setElementFDR, FDRSettings settings, FDRResult result) {

        int groupByMinPepSupport = 0;
        boolean directional = settings.isPPIDirectional();
        int maxAmbiguity = settings.getMaxProteinAmbiguity();
        int minPepCount = settings.getMinPPIPepCount();

        FDRResultLevel<ProteinGroupPair> GroupedFDRs = new FDRResultLevel<ProteinGroupPair>();
        GroupedFDRs.isDirectional = directional;

        SelfAddHashSet<ProteinGroupPair> linkPPIs = new SelfAddHashSet<ProteinGroupPair>();

        if (maxAmbiguity == 0) {
            if (directional) {
                for (ProteinGroupLink l : result.proteinGroupLinkFDR) {
                    ProteinGroupDirectionalPair dpp = new ProteinGroupDirectionalPair(l);
                    l.setFdrPPI(linkPPIs.register(dpp));
                }
            } else {
                for (ProteinGroupLink l : result.proteinGroupLinkFDR) {
                    l.setFdrPPI(linkPPIs.register(l.getProteinGroupPair()));
                }

            }

        } else {

            if (directional) {
                for (ProteinGroupLink l : result.proteinGroupLinkFDR) {

                    if (l.getProteins().size() - 1 <= maxAmbiguity) {
                        ProteinGroupDirectionalPair dpp = new ProteinGroupDirectionalPair(l);
                        linkPPIs.register(dpp);
                    }

                }
            } else {
                for (ProteinGroupLink l : result.proteinGroupLinkFDR) {

                    if (l.getProteins().size() - 1 <= maxAmbiguity) {
                        linkPPIs.register(l.getProteinGroupPair());
                    }

                }
            }

        }

        if (groupPPIByHasInternal) {

            HashSet<ProteinGroup> internal = new HashSet<>();
            for (ProteinGroupPair ppi : linkPPIs) {
                if (ppi.isInternal()) {
                    internal.add(ppi.getProtein1());
                    internal.add(ppi.getProtein1().decoyComplement());
                    internal.add(ppi.getProtein2());
                    internal.add(ppi.getProtein2().decoyComplement());
                }
            }
            for (ProteinGroupPair ppi : linkPPIs) {
                if (ppi.isBetween()
                        && internal.contains(ppi.getProtein1())
                        && internal.contains(ppi.getProtein2())) {
                    ppi.setFDRGroup(ppi.getFDRGroup() + " int_support");
                }
            }
        }

        fdr(settings.getProteinGroupPairFDR(), settings.getReportFactor(), linkPPIs, GroupedFDRs, targetProtDBSize, decoyProtDBSize, minPepCount, ignoreGroups, setElementFDR, settings.ppiLocalFDR(), false, settings.getMinTD());

        result.proteinGroupPairFDR = GroupedFDRs;
    }

    public void filterFDRLinksByFDRProteinGroupPairs(FDRResult result) {

        HashedArrayList<ProteinGroupLink> keep = new HashedArrayList<ProteinGroupLink>();
        for (ProteinGroupPair pp : result.proteinGroupPairFDR.filteredResults()) {
            keep.addAll(pp.getLinks());
        }
        result.proteinGroupLinkFDR.retainAll(keep);

    }

    public void filterFDRPeptidePairsByFDRProteinGroupLinks(FDRResult result) {

//        SubGroupFdrInfo<ProteinGroupLink> pgl = joinSubFDRInfos(result.proteinGroupLinkFDR, true) ;
//        SubGroupFdrInfo<PeptidePair> pps = joinSubFDRInfos(result.peptidePairFDR, true) ;
        long start = System.currentTimeMillis();
        HashedArrayList<PeptidePair> keep = new HashedArrayList<PeptidePair>();
        int count = 0;
        int total = result.proteinGroupLinkFDR.size();
        for (ProteinGroupLink l : result.proteinGroupLinkFDR.filteredResults()) {
            count++;
            if (count % 10000 == 0 && System.currentTimeMillis() - start > 5000) {
                System.err.println((count * 100f / total) + "% filterFDRPeptidePairsByFDRProteinGroupLinks");
                start = System.currentTimeMillis();
            }
            keep.addAll(l.getPeptidePairs());
        }
        for (PeptidePair pp : result.peptidePairFDR) {
            if (pp.isLinear()) {
                keep.add(pp);
            }
        }

        result.peptidePairFDR.retainAll(keep);
    }

    public void filterFDRPeptidePairsByFDRProteinGroups(FDRResult result) {
//        SubGroupFdrInfo<ProteinGroup> pg = joinSubFDRInfos(result.proteinGroupFDR, true) ;
//        SubGroupFdrInfo<PeptidePair> pps = joinSubFDRInfos(result.peptidePairFDR, true) ;

        HashedArrayList<PeptidePair> keep = new HashedArrayList<PeptidePair>();

        for (ProteinGroup pg : result.proteinGroupFDR.filteredResults()) {
            keep.addAll(pg.getPeptidePairs());
        }
//        for (PeptidePair pp : result.peptidePairFDR.filteredResults()) {
//            
//            Peptide pep1 = pp.getPeptide1();
//            ProteinGroup pg1 = result.proteinGroupFDR.filteredGet(pep1.getProteinGroup());
//            // don't throw out decoys peptides where the target protein was found
//            if (pg1==null && pep1.isDecoy()) {
//                pg1 = result.proteinGroupFDR.filteredGet(pep1.getProteinGroup().decoyComplement());
//            }
//            Peptide pep2 = pp.getPeptide1();
//            ProteinGroup pg2 = result.proteinGroupFDR.filteredGet(pep2.getProteinGroup());
//            // don't throw out decoys peptides where the target protein was found
//            if (pg2==null && pep2.isDecoy()) {
//                pg2 = result.proteinGroupFDR.filteredGet(pep2.getProteinGroup().decoyComplement());
//            }
//            
//            int fc = 0;
//            if (pg1 != null) {
//                pp.setFdrProteinGroup(pg1);
//                fc=1;
//            }
//            if (pg2 != null) {
//                pp.setFdrProteinGroup(pg2);
//                fc++;
//            }
//            if (fc == 2 || (fc == 1 && pp.isLinear()))
//                keep.add(pp);
//
//        }
        result.peptidePairFDR.retainAll(keep);
    }

    public void filterFDRProteinGroupsByFDRPeptidePairs(FDRResult result) {
//        SubGroupFdrInfo<ProteinGroup> pgs = joinSubFDRInfos(result.proteinGroupFDR, true) ;
//        SubGroupFdrInfo<PeptidePair> pps = joinSubFDRInfos(result.peptidePairFDR, true) ;
        HashedArrayList<ProteinGroup> keep = new HashedArrayList<ProteinGroup>();

        for (PeptidePair pp : result.peptidePairFDR.filteredResults()) {
            keep.add(pp.getProteinGroup1());
            keep.add(pp.getProteinGroup2());
        }

//        for (ProteinGroup pg : result.proteinGroupFDR) {
//            for (PeptidePair pp : pg.getPeptidePairs()) {
//                if (result.peptidePairFDR.filteredContains(pp)) {
//                    keep.add(pg);
//                    break;
//                }
//            }
//        }
        result.proteinGroupFDR.retainAll(keep);
    }

    public void filterFDRPSMByFDRPeptidePairs(FDRResult result) {
//        SubGroupFdrInfo<PSM> psms = joinSubFDRInfos(result.psmFDR, true) ;
//        SubGroupFdrInfo<PeptidePair> pps = joinSubFDRInfos(result.peptidePairFDR, true) ;

        HashedArrayList<PSM> keep = new HashedArrayList<PSM>();
        for (PeptidePair pp : result.peptidePairFDR.filteredResults()) {
            keep.addAll(pp.getAllPSMs());
        }

//        for (PSM psm : result.psmFDR) {
//
//            if (result.peptidePairFDR.filteredContains(psm.getPeptidePair())) {
//                keep.add(psm);
//            }
//        }
        result.psmFDR.retainAll(keep);
    }

    public FDRResult calculateFDR(FDRSettings settings, boolean setElementFDR) {
        FDRResult result = new FDRResult();
        this.settings = settings;
        boolean ignoreGroups = this.ignoreGroupsSetting;
        result.reportFactor = settings.getReportFactor();
        reset();
        if (cofilter !=null) {
            coResult = new FDRResult();
            cofilter.settings = new FDRSettingsImpl();
            cofilter.settings.setAll(settings);
            cofilter.settings.setPSMFDR(100d);
            cofilter.settings.setPeptidePairFDR(100d);
            cofilter.settings.setProteinGroupFDR(100d);
            cofilter.settings.setProteinGroupLinkFDR(100d);
            cofilter.settings.setProteinGroupPairFDR(100d);
            coResult.reportFactor = settings.getReportFactor();
            cofilter.reset();
        }

        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Input PSM :" + getAllPSMs().size() + "\n calculation psm-fdr");
        calculatePSMFDR(setElementFDR, ignoreGroups, result, settings);
        if (cofilter != null) {
            // do a 100%FDR to just apply any maybe defined filters
            cofilter.calculatePSMFDR(setElementFDR, ignoreGroups, coResult, cofilter.settings);
            
            // and filter the result by worst score
            coFilterResultLevel(result.psmFDR, coResult.psmFDR);
        }

        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "fdr PSM :" + result.psmFDR.getResultCount() + "\n calculation peptidepair-fdr");
        calculatePeptidePairFDR(setElementFDR, result, settings, ignoreGroups);
        if (cofilter != null) {
            // do a 100%FDR to just apply any maybe defined filters
            cofilter.calculatePeptidePairFDR(setElementFDR, coResult, cofilter.settings, ignoreGroups);
            
            // and filter the result by worst score
            coFilterResultLevel(result.peptidePairFDR, coResult.peptidePairFDR);
        }

        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "fdr peptide-pairs :" + result.peptidePairFDR.getResultCount() + "\n calculation protein-group-fdr");
        calculateProteinGroupFDR(ignoreGroups, setElementFDR, settings, result);
        if (cofilter != null) {
            // do a 100%FDR to just apply any maybe defined filters
            cofilter.calculateProteinGroupFDR(ignoreGroups, setElementFDR, cofilter.settings, coResult);
            
            // and filter the result by worst score
            coFilterResultLevel(result.proteinGroupFDR, coResult.proteinGroupFDR);
        }

        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "fdr protein groups :" + result.proteinGroupFDR.getResultCount() + "\n filtering peptide pairs by protein groups");
        filterFDRPeptidePairsByFDRProteinGroups(result);
        if (cofilter != null) {
            cofilter.filterFDRPeptidePairsByFDRProteinGroups(coResult);
        }

        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "fdr peptide-pairs :" + result.peptidePairFDR.getResultCount() + "\n calculation link-fdr");
        calculateLinkFDR(ignoreGroups, setElementFDR, settings, result);
        if (cofilter != null) {
            // do a 100%FDR to just apply any maybe defined filters
            cofilter.calculateLinkFDR(ignoreGroups, setElementFDR, cofilter.settings, coResult);
            
            // and filter the result by worst score
            coFilterResultLevel(result.proteinGroupLinkFDR, coResult.proteinGroupLinkFDR);
        }

        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "fdr links :" + result.proteinGroupLinkFDR.getResultCount() + "\n calculation protein-group-pair-fdr");
        calculateProteinGroupPairFDR(ignoreGroups, setElementFDR, settings, result);
        if (cofilter != null) {
            // do a 100%FDR to just apply any maybe defined filters
            cofilter.calculateProteinGroupPairFDR(ignoreGroups, setElementFDR, cofilter.settings, coResult);
            
            // and filter the result by worst score
            coFilterResultLevel(result.proteinGroupPairFDR, coResult.proteinGroupPairFDR);
        }

        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "fdr protein-group-pairs :" + result.proteinGroupPairFDR.getResultCount() + "\n filtering links by protein-group-pairs");
        filterFDRLinksByFDRProteinGroupPairs(result);
        if (cofilter != null) {
            cofilter.filterFDRLinksByFDRProteinGroupPairs(coResult);
        }

        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "fdr links :" + result.proteinGroupLinkFDR.getResultCount() + "\n filtering peptide pairs by links");
        filterFDRPeptidePairsByFDRProteinGroupLinks(result);
        if (cofilter != null) {
            cofilter.filterFDRPeptidePairsByFDRProteinGroupLinks(coResult);
        }

        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "fdr peptide-pairs :" + result.peptidePairFDR.getResultCount() + "\n filtering psm by peptide pairs");
        filterFDRPSMByFDRPeptidePairs(result);
        if (cofilter != null) {
            cofilter.filterFDRPSMByFDRPeptidePairs(coResult);
        }

        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "fdr psms :" + result.psmFDR.getResultCount() + "\n filtering ProteinGroups by peptide pairs");
        filterFDRProteinGroupsByFDRPeptidePairs(result);
        if (cofilter != null) {
            cofilter.filterFDRProteinGroupsByFDRPeptidePairs(coResult);
        }
        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "fdr protein groups :" + result.proteinGroupFDR.getResultCount());

        return result;

    }

    protected void coFilterResultLevel(FDRResultLevel resultdata, FDRResultLevel coResultdata) {
        // now filter each subgroup
        if (resultdata.getGroupIDs().containsAll(coResultdata.getGroupIDs()) &&
                coResultdata.getGroupIDs().containsAll(resultdata.getGroupIDs())) {
            for (Object gid : resultdata.getGroupIDs()) {
                SubGroupFdrInfo g = resultdata.getGroup(gid.toString());
                SubGroupFdrInfo cg = coResultdata.getGroup(gid.toString());
                double lowest_score = g.worstAcceptedScore;
                ArrayList toRemove = new ArrayList<>();
                for (Object p : cg.results) {
                    if (((FDRSelfAdd)p).getScore()<lowest_score)
                        toRemove.add(p);
                }
                cg.results.removeAll(toRemove);
            }
            
        } else {
            Type level  =((ParameterizedType)getClass().getGenericSuperclass()).getActualTypeArguments()[0];
            
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE,"Cofiltered dataset has not the same sub groups ("+level.getTypeName()+"\n ["
                    +RArrayUtils.toString(resultdata.getGroupIDs(), ",") + "] sv \n ["
                    +RArrayUtils.toString(coResultdata.getGroupIDs(), ",") + "]");
            System.exit(-1);
        }
    }

    public String summaryString(FDRResult result) {
        StringBuffer sb = new StringBuffer();
        sb.append("Input PSMs:");
        sb.append(getAllPSMs().size());
        sb.append(";  FDR PSM:");
        sb.append(result.psmFDR.getResultCount());
        sb.append(";  FDR PeptidePairs:");
        sb.append(result.peptidePairFDR.getResultCount());
        sb.append(";  FDR ProteinGroups:");
        sb.append(result.proteinGroupFDR.getResultCount());
        sb.append(";  FDR Links:");
        sb.append(result.proteinGroupLinkFDR.getResultCount());
        sb.append(";  FDR PPIs:");
        sb.append(result.proteinGroupPairFDR.getResultCount());
        return sb.toString();
    }

    public Collection<PSM> getInputPSMs() {
        return getAllPSMs();
    }

    public void writeFiles(String path, String baseName, String seperator, FDRResult result) throws FileNotFoundException {
        if (seperator.equals(",") || seperator.equals(";")) {
            writeFiles(path, baseName, ".csv", seperator, result);
        } else {
            writeFiles(path, baseName, ".tsv", seperator, result);
        }
    }

    public void writeFiles(String path, String baseName, String fileextension, String seperator, FDRResult result) throws FileNotFoundException {
        CSVRandomAccess csvFormater = new CSVRandomAccess(seperator.charAt(0), '"');
        csvFormater.setLocale(outputlocale);
        if (cofilter != null && coResult != null) {
            cofilter.writeFiles(path, "cofiltered_" + baseName, fileextension, seperator, coResult);
        }

        String extension = "_xiFDR" + OfflineFDR.xiFDRVersion + fileextension;

        CountOccurence<String> fdrPSMGroupCounts = new CountOccurence<String>();

        String outNameLinear = path + "/" + baseName + "_Linear_PSM" + extension;
        PrintWriter psmLinearOut = null;
        String outName = path + "/" + baseName + "_PSM" + extension;
        PrintWriter psmOut = null;
        if (!csvSummaryOnly) {
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Write PSM-results to " + outName);
            psmOut = new PrintWriter(outName);
            String header = csvFormater.valuesToString(getPSMHeader());
            psmOut.println(header);
            psmLinearOut = new PrintWriter(outNameLinear);
            psmLinearOut.println(header);
        } else {
            psmOut = NullOutputStream.NULLPRINTWRITER;
        }

        ArrayList<PSM> psms = new ArrayList<PSM>(result.psmFDR.getResultCount());
        for (SubGroupFdrInfo g : result.psmFDR.getGroups()) {
            psms.addAll(g.filteredResult);
        }

        java.util.Collections.sort(psms, new Comparator<PSM>() {

            public int compare(PSM o1, PSM o2) {
                return Double.compare(o2.getScore(), o1.getScore());
            }
        });
//        if (!isPSMScoreHighBetter())
//            java.util.Collections.reverse(fdrPSM);
        int psmCount = 0;
        int psmInternalTT = 0;
        int psmInternalTD = 0;
        int psmInternalDD = 0;
        int psmBetweenTT = 0;
        int psmBetweenTD = 0;
        int psmBetweenDD = 0;

        int psmLinearT = 0;
        int psmLinearD = 0;

        for (PSM pp : psms) {
            fdrPSMGroupCounts.add(pp.getFDRGroup());
            String line = csvFormater.valuesToString(getPSMOutputLine(pp));
            if (!csvSummaryOnly) {
                if (pp.isLinear()) {
                    psmLinearOut.println(line);
                } else {
                    psmOut.println(line);
                }
            }

            if (pp.getPeptide1() == Peptide.NOPEPTIDE || pp.getPeptide2() == Peptide.NOPEPTIDE) {
                if (pp.isTT()) {
                    psmLinearT++;
                } else if (pp.isTD()) {
                    psmLinearD++;
                }
            } else if (pp.isInternal()) {
                if (pp.isTT()) {
                    psmInternalTT++;
                } else if (pp.isTD()) {
                    psmInternalTD++;
                } else {
                    psmInternalDD++;
                }
            } else {
                if (pp.isTT()) {
                    psmBetweenTT++;
                } else if (pp.isTD()) {
                    psmBetweenTD++;
                } else {
                    psmBetweenDD++;
                }
            }
            psmCount++;
        }
        if (!csvSummaryOnly) {
            psmOut.flush();
            psmOut.close();
            psmLinearOut.flush();
            psmLinearOut.close();
        }

        ArrayList<PeptidePair> peps = new ArrayList<PeptidePair>(result.peptidePairFDR.getResultCount());
        for (SubGroupFdrInfo g : result.peptidePairFDR.getGroups()) {
            peps.addAll(g.filteredResult);
        }

        java.util.Collections.sort(peps, new Comparator<PeptidePair>() {

            public int compare(PeptidePair o1, PeptidePair o2) {
                return Double.compare(o2.getScore(), o1.getScore());
            }
        });

        CountOccurence<String> fdrPepPairGroupCounts = new CountOccurence<String>();

        int pepInternalTT = 0;
        int pepInternalTD = 0;
        int pepInternalDD = 0;
        int pepBetweenTT = 0;
        int pepBetweenTD = 0;
        int pepBetweenDD = 0;
        int pepCount = 0;

        int pepLinearT = 0;
        int pepLinearD = 0;

        PrintWriter pepsOut = null;
        PrintWriter pepsLinearOut = null;
        if (!csvSummaryOnly) {
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Write PeptidePairs-results to " + path + "/" + baseName + "_PSM" + extension);
            outName = path + "/" + baseName + "_PeptidePairs" + extension;
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Write linear PeptidePairs-results to " + outName);
            pepsOut = new PrintWriter(outName);
            String xlPepsHeader = csvFormater.valuesToString(getXLPepsHeader());
            pepsOut.println(xlPepsHeader);

            pepsLinearOut = new PrintWriter(path + "/" + baseName + "_Linear_Peptides" + OfflineFDR.xiFDRVersion + extension);
            String linearPepsHeader = csvFormater.valuesToString(getLinearPepsHeader());
            pepsLinearOut.println(linearPepsHeader);

        } else {
            pepsOut = NullOutputStream.NULLPRINTWRITER;
            pepsLinearOut = pepsOut;
        }

        for (PeptidePair pp : peps) {
            fdrPepPairGroupCounts.add(pp.getFDRGroup());

            if (pp.isLinear() && !pp.isLoop()) {

                if (!csvSummaryOnly) {
                    String line = csvFormater.valuesToString(getLinearPepeptideOutputLine(pp));
                    pepsLinearOut.println(line);
                }

                if (pp.isTT()) {
                    pepLinearT++;
                } else if (pp.isTD()) {
                    pepLinearD++;
                }
            } else {
                if (!csvSummaryOnly) {
                    String line = csvFormater.valuesToString(getXlPepeptideOutputLine(pp));
                    pepsOut.println(line);
//                    } catch (NullPointerException npe) {
//                        throw new Error();
//                    }
                }
                if (pp.isInternal()) {
                    if (pp.isTT()) {
                        pepInternalTT++;
                    } else if (pp.isTD()) {
                        pepInternalTD++;
                    } else {
                        pepInternalDD++;
                    }
                } else {
                    if (pp.isTT()) {
                        pepBetweenTT++;
                    } else if (pp.isTD()) {
                        pepBetweenTD++;
                    } else {
                        pepBetweenDD++;
                    }
                }
            }
            pepCount++;

        }
        if (!csvSummaryOnly) {
            pepsOut.flush();
            pepsOut.close();
            pepsLinearOut.flush();
            pepsLinearOut.close();
        }

        CountOccurence<String> fdrLinkGroupCounts = new CountOccurence<String>();
        ArrayList<ProteinGroupLink> links = new ArrayList<ProteinGroupLink>(result.proteinGroupLinkFDR.getResultCount());
        for (SubGroupFdrInfo g : result.proteinGroupLinkFDR.getGroups()) {
            links.addAll(g.filteredResult);
        }

        java.util.Collections.sort(links, new Comparator<ProteinGroupLink>() {

            public int compare(ProteinGroupLink o1, ProteinGroupLink o2) {
                return Double.compare(o2.getScore(), o1.getScore());
            }
        });

//        if (!isPSMScoreHighBetter())
//            java.util.Collections.reverse(fdrProtainGroupLinks);
        int linkInternalTT = 0;
        int linkInternalTD = 0;
        int linkInternalDD = 0;
        int linkBetweenTT = 0;
        int linkBetweenTD = 0;
        int linkBetweenDD = 0;
        int linkCount = 0;

        PrintWriter linksOut = null;
        if (!csvSummaryOnly) {
            linksOut = new PrintWriter(path + "/" + baseName + "_Links" + extension);
            String header = csvFormater.valuesToString(getLinkOutputHeader());
            linksOut.println(header);
        } else {
            linksOut = NullOutputStream.NULLPRINTWRITER;
        }

        // write out a table of all links
        for (ProteinGroupLink l : links) {
            fdrLinkGroupCounts.add(l.getFDRGroup());
            if (!csvSummaryOnly) {

                String line = csvFormater.valuesToString(getLinkOutputLine(l));

                linksOut.println(line);
            }
            if (l.isInternal) {
                if (l.isTT()) {
                    linkInternalTT++;
                } else if (l.isTD()) {
                    linkInternalTD++;
                } else {
                    linkInternalDD++;
                }
            } else {
                if (l.isTT()) {
                    linkBetweenTT++;
                } else if (l.isTD()) {
                    linkBetweenTD++;
                } else {
                    linkBetweenDD++;
                }
            }
            linkCount++;
        }
        if (!csvSummaryOnly) {
            linksOut.flush();
            linksOut.close();
        }

        CountOccurence<String> fdrPPIGroupCounts = new CountOccurence<String>();
        ArrayList<ProteinGroupPair> ppis = new ArrayList<ProteinGroupPair>(result.proteinGroupPairFDR.getResultCount());
        for (SubGroupFdrInfo g : result.proteinGroupPairFDR.getGroups()) {
            ppis.addAll(g.filteredResult);
        }

        java.util.Collections.sort(ppis, new Comparator<ProteinGroupPair>() {

            public int compare(ProteinGroupPair o1, ProteinGroupPair o2) {
                return Double.compare(o2.getScore(), o1.getScore());
            }
        });

//        if (!isPSMScoreHighBetter())
//            java.util.Collections.reverse(fdrProtainGroupPair);
        int ppiInternalTT = 0;
        int ppiInternalTD = 0;
        int ppiInternalDD = 0;
        int ppiBetweenTT = 0;
        int ppiBetweenTD = 0;
        int ppiBetweenDD = 0;
        int ppiCount = 0;

        // write out a table of all proteinpairs
        PrintWriter ppiOut = null;
        if (!csvSummaryOnly) {
            ppiOut = new PrintWriter(path + "/" + baseName + "_ppi" + extension);
            String header = csvFormater.valuesToString(getPPIOutputHeader());
            ppiOut.println(header);
        } else {
            ppiOut = NullOutputStream.NULLPRINTWRITER;
        }

        for (ProteinGroupPair pgp : ppis) {
            fdrPPIGroupCounts.add(pgp.getFDRGroup());

            if (!csvSummaryOnly) {
                String line = csvFormater.valuesToString(getPPIOutputLine(pgp));
                ppiOut.println(line);
            }

            if (pgp.isInternal()) {
                if (pgp.isTT()) {
                    ppiInternalTT++;
                } else if (pgp.isTD()) {
                    ppiInternalTD++;
                } else {
                    ppiInternalDD++;
                }
            } else {
                if (pgp.isTT()) {
                    ppiBetweenTT++;
                } else if (pgp.isTD()) {
                    ppiBetweenTD++;
                } else {
                    ppiBetweenDD++;
                }
            }
            ppiCount++;
        }
        if (!csvSummaryOnly) {
            ppiOut.flush();
            ppiOut.close();
        }

        CountOccurence<String> fdrProteinGroupCounts = new CountOccurence<String>();
        ArrayList<ProteinGroup> pgs = new ArrayList<ProteinGroup>(result.proteinGroupFDR.getResultCount());
        for (SubGroupFdrInfo g : result.proteinGroupFDR.getGroups()) {
            pgs.addAll(g.filteredResult);
        }

        java.util.Collections.sort(pgs, new Comparator<ProteinGroup>() {

            @Override
            public int compare(ProteinGroup o1, ProteinGroup o2) {
                return Double.compare(o2.getScore(), o1.getScore());
            }
        });
//        if (!isPSMScoreHighBetter())
//            java.util.Collections.reverse(fdrProteinGroups);

        int proteinGroupT = 0;
        int proteinGroupD = 0;
        int proteinCount = 0;

        // write out a table of all proteinpairs
        PrintWriter pgOut = null;
        if (!csvSummaryOnly) {
            pgOut = new PrintWriter(path + "/" + baseName + "_proteingroups" + extension);
            pgOut.println(csvFormater.valuesToString(getProteinGroupOutputHeader()));
        } else {
            pgOut = NullOutputStream.NULLPRINTWRITER;
        }

        for (ProteinGroup pg : pgs) {
            fdrProteinGroupCounts.add(pg.getFDRGroup());
            if (!csvSummaryOnly) {
                pgOut.println(csvFormater.valuesToString(getProteinGroupOutput(pg)));
            }
            if (pg.isDecoy()) {
                proteinGroupD++;
            } else {
                proteinGroupT++;
            }
            proteinCount++;

        }
        if (!csvSummaryOnly) {
            pgOut.flush();
            pgOut.close();
        }

//        if (!isPSMScoreHighBetter())
//            java.util.Collections.reverse(fdrProteinGroups);
        // write out a table of all proteinpairs
        PrintWriter summaryOut = null;
        if (singleSummary) {
            if (singleSummaryOut == null) {
                singleSummaryOut = new PrintWriter(path + "/" + baseName + "_summary" + extension);
            }

            summaryOut = singleSummaryOut;
            summaryOut.println("SummaryFile:" + baseName + "_summary" + extension);
        } else {
            summaryOut = new PrintWriter(path + "/" + baseName + "_summary" + extension);
        }

        summaryOut.println("xiFDR Version:" + seperator + OfflineFDR.xiFDRVersion);
        summaryOut.println("Source:" + seperator + csvFormater.quoteValue(getSource()));
        summaryOut.println(",\"Target FDRs:\"" + seperator + "Minimum supporting peptides" + seperator + "Directional");
        summaryOut.println("psm" + seperator + " " + result.psmFDR.getTargetFDR() + seperator + seperator + isPsm_directional());
        summaryOut.println("\"peptide pair\"" + seperator + " " + result.peptidePairFDR.getTargetFDR() + seperator + seperator + isPeptides_directional());
        summaryOut.println("\"protein group\"" + seperator + " " + result.proteinGroupFDR.getTargetFDR() + seperator + getMinPepPerProteinGroup());
        summaryOut.println("Link" + seperator + " " + result.proteinGroupLinkFDR.getTargetFDR() + seperator + getMinPepPerProteinGroupLink() + seperator + this.isLinks_directional());
        summaryOut.println("\"Protein Group Pair\"" + seperator + " " + result.proteinGroupPairFDR.getTargetFDR() + seperator + getMinPepPerProteinGroupPair() + seperator + this.isPpi_directional());
        summaryOut.println("\n\"max next level fdr factor (report-factor):\"" + seperator + result.reportFactor);
        summaryOut.println("\"minimum peptide length\"" + seperator + "" + (m_minPepLength <= 1 ? "unlimited" : m_minPepLength));
        if (result.uniquePSMs) {
            summaryOut.println("\"unique PSMs\"");
        }

        if (settings.getMinPeptideCoverageFilter() > 0) {
            summaryOut.println("\"minimum peptide coverage\"" + seperator + settings.getMinPeptideCoverageFilter());
        } else {
            summaryOut.println();
        }

        if (settings.getMinDeltaScoreFilter() > 0) {
            summaryOut.println("\"delta/score >\"" + seperator + settings.getMinDeltaScoreFilter());
        } else {
            summaryOut.println();
        }

        summaryOut.println("\n\"Accepted ambiguity:\"");
        summaryOut.println("\"Links for one peptide pair\"" + seperator + "" + (m_maximumLinkAmbiguity == 0 ? "unlimited" : m_maximumLinkAmbiguity));
        summaryOut.println("\"Protein pairs for one peptide pair\"" + seperator + "" + (m_maximumProteinAmbiguity == 0 ? "unlimited" : m_maximumProteinAmbiguity));
        summaryOut.println();

        if (ignoreGroupsSetting) {
            summaryOut.println("\"Groups Where Ignored \"");
        } else {
            summaryOut.println("\"Length-Group:\",\"" + RArrayUtils.toString(PeptidePair.getLenghtGroup(), seperator) + "\"");
        }
        summaryOut.println();

//        summaryOut.println("Input PSMs" +seperator + "fdr PSM" +seperator + "fdr peptide pairs" +seperator + "fdr links" +seperator + "fdr ppi");
        summaryOut.println("\"class\"" + seperator + "\"all\"" + seperator + "\"Internal TT\""
                + seperator + "\"Internal TD\"" + seperator + "\"Internal DD\""
                + seperator + "\"Between TT\"" + seperator + "\"Between TD\"" + seperator + "\"Between DD\""
                + seperator + "\"Linear T\"" + seperator + "\"Linear D\"");

        summaryOut.println("\"Input PSMs\"" + seperator + "" + getAllPSMs().size());

        summaryOut.println("\"fdr PSMs\"" + seperator + psmCount + seperator + psmInternalTT + seperator + psmInternalTD + seperator + psmInternalDD
                + seperator + psmBetweenTT + seperator + psmBetweenTD + seperator + psmBetweenDD
                + seperator + psmLinearT + seperator + psmLinearD);

        summaryOut.println("\"fdr Peptide Pairs\"" + seperator + pepCount + seperator + pepInternalTT + seperator + pepInternalTD + seperator + pepInternalDD
                + seperator + pepBetweenTT + seperator + pepBetweenTD + seperator + pepBetweenDD
                + seperator + pepLinearT + seperator + pepLinearD);

        summaryOut.println("\"fdr Link\"" + seperator + linkCount + seperator + linkInternalTT + seperator + linkInternalTD + seperator + linkInternalDD
                + seperator + linkBetweenTT + seperator + linkBetweenTD + seperator + linkBetweenDD);

        summaryOut.println("\"fdr Protein Group Pairs\"" + seperator + ppiCount + seperator + ppiInternalTT + seperator + ppiInternalTD + seperator + ppiInternalDD
                + seperator + ppiBetweenTT + seperator + ppiBetweenTD + seperator + ppiBetweenDD);

        summaryOut.println("\"fdr Protein Groups\"" + seperator + proteinCount + seperator + seperator + seperator
                + seperator + seperator + seperator
                + seperator + proteinGroupT + seperator + proteinGroupD);

//        summaryOut.println("\n\"linear fdr psm\"" + seperator + "" + (psmLinearT + psmLinearD) + "\n\"linear fdr peptide pairs\"" + seperator + "" + (pepLinearT + pepLinearD) + "\n\nfdr protein groups" + seperator + "" + fdrProteinGroups.size());
        String header = "Petide Spectrum Matches detailed summary";
//        HashMap<Integer,String> groups = new HashMap<Integer,String>();
//        for (String k : result.psmFDR.getGroupIDs()) {
//            groups.put(k,PSM.getFDRGroupName(k));
//        }
        FDRResultLevel level = result.psmFDR;
        levelSummary(summaryOut, header, level, seperator);

        header = "Petide Pairs detailed summary";
//        groups.clear();
//        for (Integer k : result.peptidePairFDR.getGroupIDs()) {
//            groups.put(k,PeptidePair.getFDRGroupName(k));
//        }
        level = result.peptidePairFDR;
        levelSummary(summaryOut, header, level, seperator);

        header = "Protein groups detailed summary";
//        groups.clear();
//        for (Integer k : result.proteinGroupFDR.getGroupIDs()) {
//            groups.put(k,ProteinGroup.getFDRGroupName(k));
//        }
        level = result.proteinGroupFDR;
        levelSummary(summaryOut, header, level, seperator);

        header = "Protein group links detailed summary";
        level = result.proteinGroupLinkFDR;
        levelSummary(summaryOut, header, level, seperator);

        header = "Protein group pairs detailed summary";

        level = result.proteinGroupPairFDR;
        levelSummary(summaryOut, header, level, seperator);

        summaryOut.flush();
        if (!singleSummary) {
            summaryOut.close();
        }

    }

    protected void peptidePairsToProteinGroups(ArrayList<PeptidePair> forwardPeps, SelfAddHashSet<ProteinGroup> pepProteinGroups) {
        // do we care about ambiguity
        if (m_maximumProteinAmbiguity > 0) {
            // yes we do
            // PeptidePairs to Protein groups
            peploop:
            for (PeptidePair pp : forwardPeps) {
                ArrayList<ProteinGroup> groups = new ArrayList<ProteinGroup>(2);
                for (Peptide p : pp.getPeptides()) {
                    ProteinGroup pg = p.getProteinGroup();

                    if (pg.proteinCount() > m_maximumProteinAmbiguity) {
                        continue peploop;
                    }

                    pg.addPeptidePair(pp);
                    groups.add(pg);
                }

                for (ProteinGroup pg : groups) {
                    pepProteinGroups.register(pg);
                }

            }

        } else {
            // no we report everything

            for (PeptidePair pp : forwardPeps) {
                for (Peptide p : pp.getPeptides()) {
                    ProteinGroup pg = pepProteinGroups.register(p.getProteinGroup());
                    pg.addPeptidePair(pp);
                }
            }

        }
    }

//    protected ArrayList<PeptidePair> proteinGroupFDR(double ProteinGroupFDR, SelfAddHashSet<ProteinGroup> pepProteinGroups, double safetyFactor, HashSet<ProteinGroup> fdrProteinGroupHS, ArrayList<PeptidePair> forwardPeps, int minPepCount, boolean ignoreGroups) {
//        // filter out the peptidepairs based on proteins
////        if (ProteinGroupFDR < 1) {
//
//        // we are talking about a purly linear view of things
//        // but I want to reuse the same function for fdr 
//        // therefore I need to change the database size, in a way, that 
//        // the fdr calculation is normalized for linear databses
//        double tCountMod = targetProtDBSize;
//        tCountMod = -0.5 - Math.sqrt(1 + 8 * tCountMod) / 2;
//        double dCountMod = decoyProtDBSize;
//        dCountMod = dCountMod / tCountMod;
//
//
//        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "ProteinGroup fdr " + pepProteinGroups.size() + " Groups as Input.");
//
//        fdrProteinGroups = fdr(ProteinGroupFDR, safetyFactor, pepProteinGroups, nextFdrProteinGroup, protFDRGroupsInput, countFdrProteinGroup, tCountMod, dCountMod, minPepCount, ignoreGroups, true);
//
//        fdrProteinGroupHS.addAll(fdrProteinGroups);
//
//        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Filter peptidepairs by  " + fdrProteinGroupHS.size() + " proteingroups.");
//        ArrayList<PeptidePair> pgtopep = new ArrayList<PeptidePair>();
//        //filter back to peptidepairs
//        if (ProteinGroupFDR < 1) {
//            for (PeptidePair pp : forwardPeps) {
//                boolean found = true;
//                for (Peptide p : pp.getPeptides()) {
//                    if (!fdrProteinGroupHS.contains(p.getProteinGroup())) {
//                        found = false;
//                        break;
//                    }
//                }
//                if (found) {
//                    pgtopep.add(pp);
//                }
//            }
//
//            forwardPeps = pgtopep;
//        }
//        Logger.getLogger(this.getClass().getName()).log(Level.INFO, forwardPeps.size() + " peptidepairs made the ProteinGroup fdr cutoff.");
////        } else {
////            fdrProteinGroups= new ArrayList<ProteinGroup>(pepProteinGroups);
////        }
//        return forwardPeps;
//    }
//    protected void PeptidePairsToGroupLinks(ArrayList<PeptidePair> forwardPeps, SelfAddHashSet<ProteinGroupLink> pepProteinGroupLinks, ArrayList<PeptidePair> linearPepPairs) {
//        // link level fdr
//        // turn peptide pairs to links
//        // peptide pairs to links
//        if (m_maximumLinkAmbiguity > 0) {
//            for (PeptidePair pp : forwardPeps) {
//                if (!pp.isLinear()) {
//                    ProteinGroupLink pgl = pp.getLink();
//                    if (pgl.getAmbiguity() <= m_maximumLinkAmbiguity) {
//                        pepProteinGroupLinks.register(pgl);
//                    }
//                } else {
//                    linearPepPairs.add(pp);
//                }
//            }
//        } else {
//            for (PeptidePair pp : forwardPeps) {
//                if (!pp.isLinear()) {
//                    pepProteinGroupLinks.register(pp.getLink());
//                } else {
//                    linearPepPairs.add(pp);
//                }
//            }
//        }
//    }
    /**
     * @return the m_maximumLinkAmbiguity
     */
    public int getMaximumLinkAmbiguity() {
        return m_maximumLinkAmbiguity;
    }

    /**
     * @param m_maximumLinkAmbiguity the m_maximumLinkAmbiguity to set
     */
    public void setMaximumLinkAmbiguity(int m_maximumLinkAmbiguity) {
        this.m_maximumLinkAmbiguity = m_maximumLinkAmbiguity;
    }

    /**
     *
     */
    public int getMinimumPeptideLength() {
        return m_minPepLength;
    }

    /**
     * @param m_maximumLinkAmbiguity the m_maximumLinkAmbiguity to set
     */
    public void setMinimumPeptideLength(int minimumPeptideLength) {
        this.m_minPepLength = minimumPeptideLength;
    }

    /**
     * @return the m_maximumProteinAmbiguity
     */
    public int getMaximumProteinAmbiguity() {
        return m_maximumProteinAmbiguity;
    }

    public int getMaximumProteinPairAmbiguity() {
        return m_maximumProteinPairAmbiguity;
    }

    /**
     * @param m_maximumProteinAmbiguity the m_maximumProteinAmbiguity to set
     */
    public void setMaximumProteinAmbiguity(int m_maximumProteinAmbiguity) {
        this.m_maximumProteinAmbiguity = m_maximumProteinAmbiguity;
    }

//    private <T extends FDRSelfAdd<T>> HashedArrayList<T> fdr(double fdr, double safetyfactor, Collection<T> c, HashMap<Integer, Double> nextFDR, HashMap<Integer, Integer> inputCounts, HashMap<Integer, Integer> resultCounts, double tCount, double dCount, int minPepCount, boolean ignoreGroups,boolean setElementFDR) {
//        HashMap<Integer, ArrayList<T>> groupedList = new HashMap<Integer, ArrayList<T>>(4);
//        HashMap<Integer, UpdateableInteger> gTT = new HashMap<Integer, UpdateableInteger>(8);
//        HashMap<Integer, UpdateableInteger> gTD = new HashMap<Integer, UpdateableInteger>(8);
//        HashMap<Integer, UpdateableInteger> gDD = new HashMap<Integer, UpdateableInteger>(8);
//
//        HashedArrayList<T> ret = new HashedArrayList<T>(c.size());
//
//        if (c.isEmpty()) {
//            return ret;
//        }
//
//        if (fdr == 1) {
//            fdr = 1000;
//            safetyfactor = 1000;
//        }
//
//
//        // split the data up into fdr-groups
//        for (T e : c) {
//            if (e.getPeptidePairCount() >= minPepCount) {
//                Integer fdrgroup = e.getFDRGroup();
//                if (ignoreGroups) {
//                    fdrgroup = -1;
//                } else {
//                    fdrgroup = e.getFDRGroup();
//                }
//                UpdateableInteger cTT;
//                UpdateableInteger cTD;
//                UpdateableInteger cDD;
//                ArrayList<T> gl = groupedList.get(fdrgroup);
//                if (gl == null) {
//                    gl = new ArrayList<T>();
//                    groupedList.put(fdrgroup, gl);
//                    if (e.isTT()) {
//                        gTT.put(fdrgroup, new UpdateableInteger(1));
//                        gTD.put(fdrgroup, new UpdateableInteger(0));
//                        gDD.put(fdrgroup, new UpdateableInteger(0));
//                    } else if (e.isTD()) {
//                        gTT.put(fdrgroup, new UpdateableInteger(0));
//                        gTD.put(fdrgroup, new UpdateableInteger(1));
//                        gDD.put(fdrgroup, new UpdateableInteger(0));
//                    } else {
//                        gTT.put(fdrgroup, new UpdateableInteger(0));
//                        gTD.put(fdrgroup, new UpdateableInteger(0));
//                        gDD.put(fdrgroup, new UpdateableInteger(1));
//                    }
//                } else {
//                    if (e.isTT()) {
//                        gTT.get(fdrgroup).value++;
//                    } else if (e.isTD()) {
//                        gTD.get(fdrgroup).value++;
//                    } else {
//                        gDD.get(fdrgroup).value++;
//                    }
//                }
//                gl.add(e);
//            }
//        }
//
//        for (Integer fdrgroup : groupedList.keySet()) {
//            int TT = gTT.get(fdrgroup).value;
//            int TD = gTD.get(fdrgroup).value;
//            int DD = gDD.get(fdrgroup).value;
//            ArrayList<T> group = groupedList.get(fdrgroup);
//            inputCounts.put(fdrgroup, group.size());
//            ArrayList<T> groupResult = new ArrayList<T>();
//            double prevFDR = subFDR(TD, DD, TT, group, fdr, safetyfactor, groupResult, tCount, dCount,setElementFDR);
//            nextFDR.put(fdrgroup, prevFDR);
//            ret.addAll(groupResult);
//            resultCounts.put(fdrgroup, groupResult.size());
//
//        }
//
//
//        if (ret.size() == 0 && c.size() > 100) {
//            // we didn't get any results through. try a non grouped
//            int TT = 0;
//            int TD = 0;
//            int DD = 0;
//            for (Integer fdrgroup : groupedList.keySet()) {
//                TT += gTT.get(fdrgroup).value;
//                TD += gTD.get(fdrgroup).value;
//                DD += gDD.get(fdrgroup).value;
//            }
//            ArrayList<T> all = new ArrayList<T>(c);
//            inputCounts.put(-1, c.size());
//            nextFDR.put(-1, subFDR(TD, DD, TT, all, fdr, safetyfactor, ret, tCount, dCount,setElementFDR));
//            resultCounts.put(-1, ret.size());
//        }
//
//
//        return ret;
//    }
//    private <T> SubGroupFdrInfo<T> joinSubFDRInfos(HashMap<Integer, SubGroupFdrInfo<T>> groupInfos, boolean store) {
//        SubGroupFdrInfo<T> joined = groupInfos.get(-1);
//        if (joined == null) {
//            joined = new SubGroupFdrInfo();
//            joined.results = new HashedArrayList<T>();
//        
//            for (SubGroupFdrInfo sgi : groupInfos.values()) {
//                joined.DD+=sgi.DD;
//                joined.TD+=sgi.TD;
//                joined.TT+=sgi.TT;
//                joined.higherFDR+=sgi.higherFDR;
//                joined.lowerFDR+=sgi.lowerFDR;
//                joined.resultCount+=sgi.resultCount;
//                joined.inputCount+=sgi.inputCount;
//                joined.saftyfactor+=sgi.saftyfactor;
//                joined.targteFDR+=sgi.targteFDR;
//                joined.results.addAll(sgi.results);
//                joined.filteredResult.addAll(sgi.filteredResult);
//            }
//            int groupCount = groupInfos.size();
//            joined.higherFDR/=groupCount;
//            joined.lowerFDR/=groupCount;
//            joined.saftyfactor/=groupCount;
//            joined.targteFDR/=groupCount;
//            if (store)
//                groupInfos.put(-1, joined);
//        }
//        
//        return joined;
//        
//    }
    private <T extends AbstractFDRElement<T>> void defineConnectedness(Collection<T> list) {
        double maxSupport = 0;
        SelfAddHashSet<Site> supports;
        supports = new SelfAddHashSet<Site>();

        // count how often each site was found
        for (T e : list) {
            Site s1 = supports.register(e.getLinkSite1());
            if (s1.getConnectedness() > maxSupport) {
                maxSupport = s1.getConnectedness();
            }
            Site s2 = e.getLinkSite2();
            if (s2 != null) {
                s2 = supports.register(e.getLinkSite2());
                if (s2.getConnectedness() > maxSupport) {
                    maxSupport = s2.getConnectedness();
                }
            }
        }

        // how connected are the second sites of a link
        double maxLinkedSupport = 0;
        SelfAddHashSet<Site> linkedSupports = new SelfAddHashSet<Site>();
        linkedSupports.selfAdd = false;
        for (T e : list) {
            Site s1 = supports.get(e.getLinkSite1());
            Site s2 = supports.get(e.getLinkSite2());
            Site ls1 = e.getLinkSite1();
            ls1.setConnectedness(0);
            ls1 = linkedSupports.register(ls1);
            if (s2 != null) {
                Site ls2 = e.getLinkSite2();
                ls2.setConnectedness(0);
                ls2 = linkedSupports.register(ls2);
                double lsc1 = ls1.getConnectedness() + s2.getConnectedness();
                if (maxLinkedSupport < lsc1) {
                    maxLinkedSupport = lsc1;
                }
                double lsc2 = ls2.getConnectedness() + s1.getConnectedness();
                if (maxLinkedSupport < lsc2) {
                    maxLinkedSupport = lsc2;
                }
                ls1.setConnectedness(lsc1);
                ls1.setConnectedness(lsc2);
            }
        }

        // now we know the conected ness of each link site and how connected the linked link-sites are
        // so we turn that into a metric for each link
        for (T e : list) {
            Site sup1 = supports.get(e.getLinkSite1());
            Site s1 = sup1;
            if (sup1 == null) {
                sup1 = s1;
            }
            double support = 0;
            support = s1.getConnectedness();
            double linkedSupport = linkedSupports.get(s1).getConnectedness();
            Site s2 = e.getLinkSite1();
            if (s2 != null) {
                s2 = supports.get(s2);
                support += s2.getConnectedness();
                linkedSupport += linkedSupports.get(s2).getConnectedness();
            }
            e.setLinkedSupport(support / (maxSupport) + linkedSupport / (2 * maxLinkedSupport));
        }
    }

    /**
     * test whether a result for a subgroup should be considered valid
     *
     * @param <T>
     * @param info
     * @param result
     * @return
     */
    public <T extends AbstractFDRElement<T>> boolean checkValid(SubGroupFdrInfo<T> info, ArrayList<T> result, double factor, int minTDCount) {
        // make sure we have enough targets that we could theoretically thsi number of TD 
        if (info.resultTT * info.targteFDR < (double) minTDCount) {
            return false;
        }

        if ((info.targteFDR < 1 && info.resultTT < info.resultDD) || info.resultTD < info.resultDD) {
            return false;
        }

        if ((info.resultDD + 0.00001) / (info.resultTD + 0.0001) > 1 - factor) {
            return false;
        }

        if (info.resultTT * info.targteFDR < factor * 10) {
            return false;
        }

        return true;

    }

//    private <T extends FDRSelfAdd<T>> HashedArrayList<T> fdr(double fdr, double safetyfactor, Collection<T> c, HashMap<Integer, SubGroupFdrInfo<T>> groupInfo, double tCount, double dCount, int minPepCount, boolean ignoreGroups,boolean setElementFDR) {
    /**
     * Splits up data into groups and forward each group to
     * {@link #subFDR(java.util.ArrayList, java.util.ArrayList, boolean, rappsilber.fdr.result.SubGroupFdrInfo) subFDR}.
     * results will then be collected
     *
     * @param <T> What type of Information is the FDR to be applied
     * @param fdr the target FDR that is supossed to be returned
     * @param safetyfactor don't report anything if the next higher calculatable
     * @param fdrInput the input data that should be filtered by FDR
     * @param groupInfo will hold all the information for the currently
     * calculated FDRs
     * @param tCount number of entries in the target database. This is used to
     * normalise the FDR calculation for imbalenced databases
     * @param dCount number of entries in the decoy database. This is used to
     * normalise the FDR calculation for imbalenced databases
     * @param minPepCount if an entry has less then these number of peptides
     * supporting it - it will not be considered for the FDR-calculation
     * @param ignoreGroups should we ignore groups cmpletely and just calculate
     * a joined FDR?
     * @param setElementFDR should each element be flaged up with the FDR that
     * it group has at the given score? - is used to speed up the maximation (by
     * not setting these)
     */
    private <T extends AbstractFDRElement<T>> void fdr(double fdr, double safetyfactor, Collection<T> fdrInput, FDRResultLevel<T> groupInfo, double tCount, double dCount, int minPepCount, boolean ignoreGroups, boolean setElementFDR, Boolean localFDR, boolean groupByProteinPair, int minTDChance) {
        HashMap<String, ArrayList<T>> groupedList = new HashMap<String, ArrayList<T>>(4);
        HashMap<String, UpdateableInteger> gTT = new HashMap<String, UpdateableInteger>(8);
        HashMap<String, UpdateableInteger> gTD = new HashMap<String, UpdateableInteger>(8);
        HashMap<String, UpdateableInteger> gDD = new HashMap<String, UpdateableInteger>(8);

        int resultcount = 0;
//        HashedArrayList<T> ret = new HashedArrayList<T>(c.size());
        T firstinput = null;
        boolean countmatches = false;
        if (!fdrInput.isEmpty()) {
            firstinput = fdrInput.iterator().next();
            if (firstinput instanceof PSM) {
                countmatches = true;
            }
        }

        if (fdrInput.isEmpty()) {
            return;
        }

        if (fdr == 1) {
            fdr = 1000;
            safetyfactor = 1000;
        }
        //HashMap<String,Integer> protpairToID = new HashMap<>();
        int maxID = 2;

        // split the data up into fdr-groups
        for (T e : fdrInput) {
            if (e.getPeptidePairCount() >= minPepCount) {
                String fdrgroup = null;
                if (ignoreGroups) { // don't do any grouping
                    fdrgroup = "ALL";
                } else {
                    if (groupByProteinPair) {
                        // do group by protein pairs to begin with
                        ProteinGroup pg1 = e.getProteinGroup1();
                        ProteinGroup pg2 = e.getProteinGroup2();
                        // make sure we always get the same key
                        String pgs1 = pg1.accessionsNoDecoy();
                        String pgs2 = pg2.accessionsNoDecoy();

                        String k1 = pgs1 + "_xl_" + pgs2;
                        String k2 = pgs2 + "_xl_" + pgs1;
                        Integer id = protpairToID.get(k1);
                        if (id == null) // not found with key 1 -> try reversed key
                        {
                            id = protpairToID.get(k2);
                        }
                        // new protein pair
                        if (id == null) {
                            // assign id
                            id = maxID++;
//                            if (pgs1.contentEquals(pgs2)) {
//                                id=-id;
//                            }
                            protpairToID.put(k1, id);
                            protpairToID.put(k2, id);
                            protpairToSize.put(id, new UpdateableInteger(1));
                        } else if (countmatches) {
                            protpairToSize.get(id).value++;
                        }
                        // set the id as fdr group
                        fdrgroup = "" + id;
                        // and write it to the 
                        e.setFDRGroup("" + id);
                    } else {
                        fdrgroup = e.getFDRGroup();
                    }
                }
                // get the right list of matches
                ArrayList<T> gl = groupedList.get(fdrgroup);
                if (gl == null) {
                    // does not exist yet - so make a new one
                    gl = new ArrayList<T>();
                    groupedList.put(fdrgroup, gl);
                    // start counting the types
                    if (e.isTT()) {
                        gTT.put(fdrgroup, new UpdateableInteger(1));
                        gTD.put(fdrgroup, new UpdateableInteger(0));
                        gDD.put(fdrgroup, new UpdateableInteger(0));
                    } else if (e.isTD()) {
                        gTT.put(fdrgroup, new UpdateableInteger(0));
                        gTD.put(fdrgroup, new UpdateableInteger(1));
                        gDD.put(fdrgroup, new UpdateableInteger(0));
                    } else {
                        gTT.put(fdrgroup, new UpdateableInteger(0));
                        gTD.put(fdrgroup, new UpdateableInteger(0));
                        gDD.put(fdrgroup, new UpdateableInteger(1));
                    }
                } else {
                    if (e.isTT()) {
                        gTT.get(fdrgroup).value++;
                    } else if (e.isTD()) {
                        gTD.get(fdrgroup).value++;
                    } else {
                        gDD.get(fdrgroup).value++;
                    }
                }
                gl.add(e);
            }
        }

        // join groups with similare number of PSMs
        if (groupByProteinPair) {
            UpdateableInteger max = RArrayUtils.max(protpairToSize.values());
            HashMap<String, ArrayList<T>> groupedListNew = new HashMap<String, ArrayList<T>>(4);
            HashMap<String, UpdateableInteger> gTTNew = new HashMap<String, UpdateableInteger>(8);
            HashMap<String, UpdateableInteger> gTDNew = new HashMap<String, UpdateableInteger>(8);
            HashMap<String, UpdateableInteger> gDDNew = new HashMap<String, UpdateableInteger>(8);
            for (Map.Entry<String, ArrayList<T>> glOldE : groupedList.entrySet()) {
                ArrayList<T> glOld = glOldE.getValue();
                String oldgroup = glOldE.getKey();
                UpdateableInteger protpairsize = protpairToSize.get(oldgroup);
                String sizekey = "" + (glOld.size() + 100);
                if (protpairsize != null) {
                    sizekey = "" + protpairsize.value;
                }
                int oldTT = gTT.get(oldgroup).value;
                int oldTD = gTD.get(oldgroup).value;
                int oldDD = gDD.get(oldgroup).value;
                ArrayList<T> glNew = groupedListNew.get(sizekey);
                if (glNew == null) {
                    groupedListNew.put("" + sizekey, glOld);

                    gTTNew.put(sizekey, new UpdateableInteger(oldTT));
                    gTDNew.put(sizekey, new UpdateableInteger(oldTD));
                    gDDNew.put(sizekey, new UpdateableInteger(oldDD));
                } else {
                    glNew.addAll(glOld);
                    gTTNew.get(sizekey).value += gTT.get(oldgroup).value;
                    gTDNew.get(sizekey).value += gTD.get(oldgroup).value;
                    gDDNew.get(sizekey).value += gDD.get(oldgroup).value;
                }
            }
            groupedList = groupedListNew;
            gTT = gTTNew;
            gTD = gTDNew;
            gDD = gDDNew;
        }

        // collect all subgroups that have to small number of results
        ArrayList<String> toosmall = new ArrayList<>();
        ArrayList<String> toosmallWithin = new ArrayList<>();
        int toosmallsum = 0;
        int toosmallsumWithin = 0;
        HashSet<String> collectedToSmallNames = new HashSet<String>();
        HashSet<String> collectedToSmallWithinNames = new HashSet<String>();
        SubGroupFdrInfo<T> collectedToSmall = new SubGroupFdrInfo<T>();
        SubGroupFdrInfo<T> collectedToSmallWithin = new SubGroupFdrInfo<T>();
        ArrayList<T> tooSmallGroup = new ArrayList<>();
        ArrayList<T> tooSmallGroupWithin = new ArrayList<>();

        for (String fdrgroup : groupedList.keySet()) {
            SubGroupFdrInfo<T> info = new SubGroupFdrInfo<T>();
            info.TT = gTT.get(fdrgroup).value;
            info.TD = gTD.get(fdrgroup).value;
            info.DD = gDD.get(fdrgroup).value;
            ArrayList<T> group = groupedList.get(fdrgroup);

            info.DCount = dCount;
            info.TCount = tCount;

            ArrayList<T> groupResult = new ArrayList<T>();
            info.inputCount = group.size();
            info.targteFDR = fdr;
            info.saftyfactor = safetyfactor;
            info.fdrGroup = fdrgroup;
            if (localFDR == null || localFDR) {
                subFDRLocal(group, groupResult, info);
            }
            subFDR(group, groupResult, setElementFDR, info);
            if (localFDR != null && localFDR) {
                ArrayList<T> filtered = new ArrayList<>();
                for (T e : groupResult) {
                    if (e.getPEP() <= fdr) {
                        filtered.add(e);
                    } else {
                        if (e.isTT()) {
                            info.resultTT--;
                        } else if (e.isTD()) {
                            info.resultTD--;
                        } else {
                            info.resultDD--;
                        }
                    }
                }
                groupResult = filtered;
            }

//            if (info.resultTT*fdr<getMinDecoys())  {
            if (!checkValid(info, groupResult, 0, minTDChance)) {
                Logger.getLogger(this.getClass().getName()).log(Level.FINE, "Discarded group {0}->{1}({2})", new Object[]{group.get(0).getClass(), fdrgroup, group.get(0).getFDRGroup()});
                //between
                if (fdrgroup.toLowerCase().contains("between")) {
                    toosmallsum += groupResult.size();
                    toosmall.add(fdrgroup);
                    collectedToSmallNames.add(fdrgroup);
                    collectedToSmall.TT += info.TT;
                    collectedToSmall.TD += info.TD;
                    collectedToSmall.DD += info.DD;
                    collectedToSmall.DCount += info.DCount;
                    collectedToSmall.TCount += info.TCount;
                    collectedToSmall.inputCount += info.inputCount;
                    tooSmallGroup.addAll(group);
                    groupResult.clear();
                } else {
                    // within group
                    collectedToSmallWithinNames.add(fdrgroup);
                    toosmallsumWithin += groupResult.size();
                    toosmallWithin.add(fdrgroup);
                    collectedToSmallWithin.TT += info.TT;
                    collectedToSmallWithin.TD += info.TD;
                    collectedToSmallWithin.DD += info.DD;
                    collectedToSmallWithin.DCount += info.DCount;
                    collectedToSmallWithin.TCount += info.TCount;
                    collectedToSmallWithin.inputCount += info.inputCount;
                    tooSmallGroupWithin.addAll(group);
                    groupResult.clear();

                }
                info.filteredResult = new HashedArrayList<>();
                info.didntPassCheck = true;
            } //else {
            groupInfo.addGroup(fdrgroup, info);
            //                    group, fdr, safetyfactor, groupResult, tCount, dCount,setElementFDR);
            //            nextFDR.put(fdrgroup, prevFDR);
            //            ret.addAll(groupResult);
            resultcount += groupResult.size();
            groupInfo.setLinear(groupInfo.getLinear() + info.linear);
            groupInfo.setWithin(groupInfo.getWithin() + info.within);
            groupInfo.setBetween(groupInfo.getBetween() + info.between);
            //            resultCounts.put(fdrgroup, groupResult.size());
            //}

        }

        // did we collect some subgroups with to small numbers?
        if (toosmallsum + toosmallsumWithin > getMinDecoys() / fdr) {
            Logger.getLogger(this.getClass().getName()).log(Level.FINE, "Join up discarded groups to try and get some more results");
            if (tooSmallGroup.size() > 0) {
                ArrayList<T> groupResult = new ArrayList<T>();
                collectedToSmall.targteFDR = fdr;
                collectedToSmall.saftyfactor = safetyfactor;
                collectedToSmall.fdrGroup = "CollectedSmallResults [" + RArrayUtils.toString(collectedToSmallNames, ",") + "]";
                subFDR(tooSmallGroup, groupResult, setElementFDR, collectedToSmall);
                if (localFDR == null || localFDR) {
                    subFDRLocal(tooSmallGroup, groupResult, collectedToSmall);
                    if (localFDR != null) {
                        ArrayList<T> filtered = new ArrayList<>();
                        for (T e : groupResult) {
                            if (e.getPEP() <= fdr) {
                                filtered.add(e);
                            } else {
                                if (e.isTT()) {
                                    collectedToSmall.resultTT--;
                                } else if (e.isTD()) {
                                    collectedToSmall.resultTD--;
                                } else {
                                    collectedToSmall.resultDD--;
                                }
                            }
                        }
                        groupResult = filtered;
                    }
                }
                if (checkValid(collectedToSmall, groupResult, 0, minTDChance)) {
                    groupInfo.addGroup(collectedToSmall.fdrGroup, collectedToSmall);
                    //                    group, fdr, safetyfactor, groupResult, tCount, dCount,setElementFDR);
                    //            nextFDR.put(fdrgroup, prevFDR);
                    //            ret.addAll(groupResult);
                    resultcount += groupResult.size();
                    groupInfo.setLinear(groupInfo.getLinear() + collectedToSmall.linear);
                    groupInfo.setWithin(groupInfo.getWithin() + collectedToSmall.within);
                    groupInfo.setBetween(groupInfo.getBetween() + collectedToSmall.between);
                }
            }

            if (tooSmallGroupWithin.size() > 0) {
                ArrayList<T> groupResultwithin = new ArrayList<T>();

                collectedToSmallWithin.targteFDR = fdr;
                collectedToSmallWithin.saftyfactor = safetyfactor;
                collectedToSmallWithin.fdrGroup = "CollectedSmallResultsWithin  [" + RArrayUtils.toString(collectedToSmallWithinNames, ",") + "]";
                subFDR(tooSmallGroupWithin, groupResultwithin, setElementFDR, collectedToSmallWithin);
                if (localFDR == null || localFDR) {
                    subFDRLocal(tooSmallGroupWithin, groupResultwithin, collectedToSmallWithin);
                    if (localFDR != null) {
                        ArrayList<T> filtered = new ArrayList<>();
                        for (T e : groupResultwithin) {
                            if (e.getPEP() <= fdr) {
                                filtered.add(e);
                            } else {
                                if (e.isTT()) {
                                    collectedToSmallWithin.resultTT--;
                                } else if (e.isTD()) {
                                    collectedToSmallWithin.resultTD--;
                                } else {
                                    collectedToSmallWithin.resultDD--;
                                }
                            }
                        }
                        groupResultwithin = filtered;
                    }
                }
                if (checkValid(collectedToSmallWithin, groupResultwithin, 0, minTDChance)) {
                    groupInfo.addGroup(collectedToSmallWithin.fdrGroup, collectedToSmallWithin);
                    //                    group, fdr, safetyfactor, groupResult, tCount, dCount,setElementFDR);
                    //            nextFDR.put(fdrgroup, prevFDR);
                    //            ret.addAll(groupResult);
                    resultcount += groupResultwithin.size();
                    groupInfo.setLinear(groupInfo.getLinear() + collectedToSmallWithin.linear);
                    groupInfo.setWithin(groupInfo.getWithin() + collectedToSmallWithin.within);
                    groupInfo.setBetween(groupInfo.getBetween() + collectedToSmallWithin.between);
                }
            }
        }

        if (resultcount == 0 && fdrInput.size() > 100) {
            // we didn't get any results through. try a non grouped
            SubGroupFdrInfo info = new SubGroupFdrInfo();
            info.fdrGroup = "NoSubResults";
            groupInfo.addGroup(info.fdrGroup, info);
            for (String fdrgroup : groupedList.keySet()) {
                info.TT += gTT.get(fdrgroup).value;
                info.TD += gTD.get(fdrgroup).value;
                info.DD += gDD.get(fdrgroup).value;
            }

            info.inputCount = fdrInput.size();

            ArrayList<T> allResult = new ArrayList<T>();

            ArrayList<T> all = new ArrayList<T>(fdrInput);
            subFDR(all, allResult, setElementFDR, info);
            if (localFDR == null || localFDR) {
                subFDRLocal(all, allResult, info);
                if (localFDR != null) {
                    ArrayList<T> filtered = new ArrayList<>();
                    for (T e : allResult) {
                        if (e.getPEP() <= fdr) {
                            filtered.add(e);
                        } else {

                            if (e.isTT()) {
                                info.resultTT--;
                            } else if (e.isTD()) {
                                info.resultTD--;
                            } else {
                                info.resultDD--;
                            }
                        }
                    }
                    allResult = filtered;
                }
            }
            if (allResult.size() * fdr < minTDChance) {
                allResult.clear();
            } else {
                groupInfo.setLinear(groupInfo.getLinear() + info.linear);
                groupInfo.setWithin(groupInfo.getWithin() + info.within);
                groupInfo.setBetween(groupInfo.getBetween() + info.between);
            }

        }

//        return ret;
    }

//    /**
//     * takes a sub-set of entries and does the actual FDR-estimation/cutoff.
//     * @param <T>       What type of Information is the FDR to be applied
//     * @param TD
//     * @param DD
//     * @param TT
//     * @param group     All entri
//     * @param fdr       the target FDR that is supossed to be returned
//     * @param safetyfactor  don't report anything if the next higher calculatable
//     * @param results
//     * @param TCount
//     * @param DCount
//     * @param setElementFDR
//     * @return
//     */
//    protected <T extends FDRSelfAdd<T>> double subFDR(int TD, int DD, int TT, ArrayList<T> group, double fdr, double saftyfactor, ArrayList<T> results, double TCount, double DCount,boolean setElementFDR) {
//        
//        ArrayList<T> ret = new ArrayList<T>();
//        double TTfactor = 1;
//        double DDfactor = (TCount*TCount+TCount)/(DCount*DCount+DCount);
//        double TDfactor = TCount / DCount;
//        double normTT = TT *TTfactor;
//        double normTD = TD * TDfactor;
//        double normDD = DD * DDfactor;
//
//        // total fdr rate
//        double prevFDR = 0;
//        if (normTD!=0)
//            prevFDR = ((normTD+1) + normDD*(1-2*normTD/((normTD)+Math.sqrt(normTD))))/normTT;
//        
//        if (fdr >= 1) {
//            fdr = Double.POSITIVE_INFINITY;
//        }
//
//        // check whether the fdr to reach the fdr in this dataset we would need at least the minimum number of decoys
//        // if not we assume the fdr-calculation does not make sense
//        if ((TT + TD + DD) * fdr < MINIMUM_POSSIBLE_DECOY && fdr < 1) {
//            return prevFDR;
//        }
//
//        // all the elements within this group
//        int groupSize = group.size();
//
//        // sort them by score
//        Collections.sort(group);
//
////        if (!isPSMScoreHighBetter())
////            Collections.reverse(group);
//
//        int fdrSwitch = groupSize - 1;
//        double highFDR = prevFDR;
//        // now we can just go through and find the cut-off
//        for (int i = groupSize - 1; i >= 0; i--) {
//            //double efdr = (TD / TDfactor - DD / DDfactor) * TTfactor / (double) TT;
//            normTT = TT *TTfactor;
//            normTD = TD * TDfactor;
//            normDD = DD * DDfactor;
//            
//            double efdr = 0;
//            
//            if (normTD>0)
//                efdr=(normTD + normDD*(1-2*normTD/((normTD)+Math.sqrt(normTD))))/normTT;
//            
//            if (efdr <= fdr && efdr<prevFDR &&  prevFDR / fdr < saftyfactor) {
//                
//                // as a safty meassure only accept things, if we don't have a to big a jump in fdr
////                if (prevFDR / fdr < saftyfactor) {
//                    double lastFDR = prevFDR;
//                    double setFDR = efdr;
//                    if (setElementFDR) {
//                        for (; i >= 0; i--) {
//                            T e = group.get(i);
//                            e.setFDR(setFDR);
//                            ret.add(e);
//                            if (e.isTT()) {
//                                TT--;
//                            } else if (e.isTD()) {
//                                TD--;
//                            } else if (e.isDD()) {
//                                DD--;
//                            }
//
//                            normTT = TT * TTfactor;
//                            normTD = TD * TDfactor;
//                            normDD = DD * DDfactor;
//            
//                            
//                            if (TD > 0) {
//                                setFDR = Math.min(setFDR, (normTD + normDD*(1-2*normTD/((normTD)+Math.sqrt(normTD))))/normTT);
//                            } else
//                                setFDR = 0;
//                        }
//                    } else {
//                         for (; i >= 0; i--) {
//                            T e = group.get(i);
//                            ret.add(e);                             
//                         }
//                    }
////                }
//                if (fdr >=1 || ret.size() >=MINIMUM_POSSIBLE_RESULT)
//                    results.addAll(ret);
//                return prevFDR;
//            }
//            prevFDR = Math.min(efdr, prevFDR);
//
//            T e = group.get(i);
//            if (e.isTT()) {
//                TT--;
//            } else if (e.isTD()) {
//                TD--;
//            } else if (e.isDD()) {
//                DD--;
//            } else {
//                Logger l = Logger.getLogger(FDRCalculation.class.getName());
//                l.log(Level.SEVERE, "Something is wrong here!", new Exception(""));
//                System.exit(-1);
//            }
//
//        }
//        return prevFDR;
//    }
    /**
     * Takes a sub-set of entries and does the actual FDR-estimation/cutoff.
     * <p>
     * The assumption is that we have a symmetric cross-linker. <br/>Therefore:
     * <p>
     * {@code FP(TT)=TD+DD*(1-TDdb/DDdb)}</p>
     * </p><p>
     * If TCount and DCount are of different size then the false positive
     * estimation gets scaled accordingly.</p
     *
     * @param <T> What type of Information is the FDR to be applied
     * @param TD total number of TD matches
     * @param DD total number of DD matches
     * @param TT total number of TT matches
     * @param group All entries for which the target FDR should be calculated
     * @param fdr the target FDR that is supposed to be returned
     * @param safetyfactor don't report anything if the next higher calculable
     * exceeds the target-FDR by the given factor
     * @param results the passing results will be added to this ArrayList
     * @param TCount Size of the target database
     * @param DCount size of the decoy database
     * @param isSymmetric do we calculate an FDR for a symmetric or an
     * asymmetric experiment
     * @param setElementFDR
     * @return
     */
//    protected <T extends FDRSelfAdd<T>> double subFDR(int TD, int DD, int TT, ArrayList<T> group, double fdr, double safetyfactor, ArrayList<T> results, double TCount, double DCount,boolean setElementFDR, boolean isSymmetric) {
//        
//        ArrayList<T> ret = new ArrayList<T>();
//        
//        double TTSize;
//        double TDSize;
//        double DDSize;
//        double k;
//        
//        if (isSymmetric) {
//            TDSize = TCount*DCount;
//            DDSize = (DCount*DCount+DCount)/2;
//            k=TDSize/DDSize;
//            TTSize = (TCount*TCount+TCount)/2;
//        } else {
//            TDSize = 2*TCount*DCount;
//            DDSize = DCount*DCount;
//            TTSize = TCount*TCount;
//            k=TDSize/DDSize;
//        }
//        double TTfactor = 1;
//        double DDfactor = TTSize/DDSize;
//        double TDfactor = k*TTSize/TDSize;
//
//        double normTT = TT * TTfactor;
//        double normTD = TD * TDfactor;
//        double normDD = DD * DDfactor;
//        
//
//        // total fdr rate
//        double prevFDR = 0;
//        if (normTD!=0)
//            prevFDR = (normTD + normDD*(1-k))/normTT;
//        
//        if (fdr >= 1) {
//            fdr = Double.POSITIVE_INFINITY;
//        }
//
//        // check whether the fdr to reach the fdr in this dataset we would need at least the minimum number of decoys
//        // if not we assume the fdr-calculation does not make sense
//        if ((TT + TD + DD) * fdr < MINIMUM_POSSIBLE_DECOY && fdr < 1) {
//            return prevFDR;
//        }
//
//        // all the elements within this group
//        int groupSize = group.size();
//
//        // sort them by score
//        Collections.sort(group);
//
////        if (!isPSMScoreHighBetter())
////            Collections.reverse(group);
//
//        int fdrSwitch = groupSize - 1;
//        double highFDR = prevFDR;
//        // now we can just go through and find the cut-off
//        for (int i = groupSize - 1; i >= 0; i--) {
//            //double efdr = (TD / TDfactor - DD / DDfactor) * TTfactor / (double) TT;
//            normTT = TT *TTfactor;
//            normTD = TD * TDfactor;
//            normDD = DD * DDfactor;
//            
//            double efdr = 0;
//            
//            if (normTD>0)
//                efdr=(normTD + normDD*(1-k))/normTT;
//            
//            if (efdr <= fdr && efdr<prevFDR &&  prevFDR / fdr < safetyfactor) {
//                
//                // as a safty meassure only accept things, if we don't have a to big a jump in fdr
////                if (prevFDR / fdr < safetyfactor) {
//                    double lastFDR = prevFDR;
//                    double setFDR = efdr;
//                    if (setElementFDR) {
//                        for (; i >= 0; i--) {
//                            T e = group.get(i);
//                            e.setFDR(setFDR);
//                            ret.add(e);
//                            if (e.isTT()) {
//                                TT--;
//                            } else if (e.isTD()) {
//                                TD--;
//                            } else if (e.isDD()) {
//                                DD--;
//                            }
//
//                            normTT = TT * TTfactor;
//                            normTD = TD * TDfactor;
//                            normDD = DD * DDfactor;
//            
//                            
//                            if (TD > 0) {
//                                setFDR = Math.min(setFDR, (normTD + normDD*(1-k))/normTT);
//                            } else
//                                setFDR = 0;
//                        }
//                    } else {
//                         for (; i >= 0; i--) {
//                            T e = group.get(i);
//                            ret.add(e);                             
//                         }
//                    }
////                }
//                if (fdr >=1 || ret.size() >=MINIMUM_POSSIBLE_RESULT)
//                    results.addAll(ret);
//                return prevFDR;
//            }
//            prevFDR = Math.min(efdr, prevFDR);
//
//            T e = group.get(i);
//            if (e.isTT()) {
//                TT--;
//            } else if (e.isTD()) {
//                TD--;
//            } else if (e.isDD()) {
//                DD--;
//            } else {
//                Logger l = Logger.getLogger(FDRCalculation.class.getName());
//                l.log(Level.SEVERE, "Something is wrong here!", new Exception(""));
//                System.exit(-1);
//            }
//
//        }
//        return prevFDR;
//    }
    /**
     * Takes a sub-set of entries and does the actual FDR-estimation/cutoff.
     * <p>
     * The assumption is that we have a symmetric cross-linker. <br/>Therefore:
     * <p>
     * {@code FP(TT)=TD+DD*(1-TDdb/DDdb)}</p>
     * </p><p>
     * If TCount and DCount are of different size then the false positive
     * estimation gets scaled accordingly.</p
     *
     * @param <T> What type of Information is the FDR to be applied
     * @param group All entries for which the target FDR should be calculated
     * @param fdr the target FDR that is supposed to be returned
     * @param safetyfactor don't report anything if the next higher calculable
     * exceeds the target-FDR by the given factor
     * @param results the passing results will be added to this ArrayList
     * @param TCount Size of the target database
     * @param DCount size of the decoy database
     * @param info all informations needed for this sub group
     * @param isSymmetric do we calculate an FDR for a symmetric or an
     * asymmetric experiment
     * @param setElementFDR
     * @return
     */
    protected <T extends AbstractFDRElement<T>> double subFDR(ArrayList<T> group, ArrayList<T> results, boolean setElementFDR, SubGroupFdrInfo info) {

        HashedArrayList<T> ret = new HashedArrayList<T>();
        info.results = ret;
        info.filteredResult = ret;

        double TT = info.TT;
        int TD = info.TD;
        int DD = info.DD;

        double fdr = info.targteFDR;

        if (fdr >= 1) {
            fdr = Double.POSITIVE_INFINITY;
        }

        // all the elements within this group
        int groupSize = group.size();

        // sort them by score
        Collections.sort(group, new Comparator<T>() {

            public int compare(T o1, T o2) {
                return Double.compare(o2.getScore() * o2.getLinkedSupport(), o1.getScore() * o1.getLinkedSupport());
            }
        });

        // total fdr rate
        double prevFDR = (TD - DD) / TT;
        int prevTDIndex = groupSize - 1;

//        if (!isPSMScoreHighBetter())
//            Collections.reverse(group);
        int fdrSwitch = groupSize - 1;
        double highFDR = prevFDR;
        // now we can just go through and find the cut-off
        for (int i = groupSize - 1; i >= 0; i--) {

            double efdr = (TD - DD) / TT;
            double efdr_n = ((TD - 1) - DD) / TT;
            double efdr_p = ((TD + 1) - DD) / TT;
            if (efdr_n < 0 && DD == 0) {
                efdr_n = 0;
            }

            T e = group.get(i);
            double score = e.getScore();

            // we steped below the target fdr (efdr<= fdr) 
            if (efdr <= fdr) {
                //did we just pas an fdr-step
                if ((efdr < prevFDR || prevFDR == 0.0 || prevFDR == Double.POSITIVE_INFINITY) || fdr == Double.POSITIVE_INFINITY) {
                    // does it fullfill the saftyfactor
                    if ((efdr_p / fdr < info.saftyfactor || info.saftyfactor > 1000 || info.saftyfactor == 0) || fdr == Double.POSITIVE_INFINITY) {
                        //if ((efdr_p / fdr < info.saftyfactor && efdr_n >= 0) || fdr == Double.POSITIVE_INFINITY)  {
                        if (info.firstPassingFDR == 0) {
                            info.firstPassingFDR = prevFDR;
                        }

                        info.higherFDR = prevFDR;
                        info.worstAcceptedScore
                                = info.lowerFDR = efdr;
                        info.resultCount = i;
                        int lastFDRIndex = i;

                        double lastFDR = prevFDR;
                        int lastTDIndex = prevTDIndex;
                        T lastTDElement = group.get(lastTDIndex);
                        double setFDR = efdr;
                        info.resultTT = (int) TT;
                        info.resultTD = TD;
                        info.resultDD = DD;
                        if (setElementFDR) {
                            for (; i >= 0; i--) {

                                e = group.get(i);
                                e.setFDR(setFDR);
                                e.setHigherFDR(lastFDR);
                                e.setHigherTD(lastTDElement);
                                ret.add(e);

                                if (e.isTT()) {
                                    TT--;
                                    if (e.isLinear()) {
                                        info.linear++;
                                    }

                                    if (e.isInternal()) {
                                        info.within++;
                                    } else if (e.isBetween()) {
                                        info.between++;
                                    }

                                } else if (e.isTD()) {
                                    TD--;
                                    info.resultTD++;
                                    for (int l = i + 1; l <= lastTDIndex; l++) {
                                        group.get(l).setLowerTD(e);
                                    }
                                    lastTDIndex = i;
                                    lastTDElement = e;
                                } else if (e.isDD()) {
                                    DD--;
                                    info.resultDD++;
                                }

                                double currfdr = 0;
                                if (TD > 0) {
                                    currfdr = (TD - DD) / TT;
                                }

                                if (currfdr < setFDR) {
                                    lastFDR = setFDR;
                                    // we reached a new lower fdr
                                    setFDR = currfdr;

                                    // set the lower fdr values for the previous data
                                    for (int li = lastFDRIndex; li <= i; li++) {
                                        group.get(li).setLowerFDR(currfdr);
                                    }
                                    lastFDRIndex = i - 1;
                                }

                            }
                            for (int li = lastFDRIndex; li <= i && li >= 0; li++) {
                                group.get(li).setLowerFDR(0);
                            }

                            T best = group.get(0);
                            for (int l = 0; l <= lastTDIndex; l++) {
                                group.get(l).setLowerTD(best);
                            }

                        } else {
                            for (; i >= 0; i--) {
                                e = group.get(i);
                                if (e.isTT()) {
                                    if (e.isLinear()) {
                                        info.linear++;
                                    }

                                    if (e.isInternal()) {
                                        info.within++;
                                    } else if (e.isBetween()) {
                                        info.between++;
                                    }
                                }
                                ret.add(e);
                            }
                        }
                        //                }

                        //  info.results = ret;
                        //  info.filteredResult = ret;
                        //                if (fdr >=1 || ret.size() >=MINIMUM_POSSIBLE_RESULT)
                        results.addAll(ret);
                        info.results = ret;
                        info.filteredResult = ret;
                        info.worstAcceptedScore = score;

                        return prevFDR;
                    } else {
                        info.firstPassingFDR = prevFDR;
                    }
                }
            }

            if (efdr < prevFDR) {
                prevFDR = efdr;
            }

            if (e.isTT()) {
                TT--;
            } else if (e.isTD()) {
                TD--;
                prevTDIndex = i;
            } else if (e.isDD()) {
                DD--;
            } else {
                Logger l = Logger.getLogger(this.getClass().getName());
                l.log(Level.SEVERE, "Something is wrong here!", new Exception(""));
                System.exit(-1);
            }

        }
        return prevFDR;
    }

    /**
     * Takes a sub-set of entries and does the actual FDR-estimation/cutoff.
     * <p>
     * The assumption is that we have a symmetric cross-linker. <br/>Therefore:
     * <p>
     * {@code FP(TT)=TD+DD*(1-TDdb/DDdb)}</p>
     * </p><p>
     * If TCount and DCount are of different size then the false positive
     * estimation gets scaled accordingly.</p
     *
     * @param <T> What type of Information is the FDR to be applied
     * @param group All entries for which the target FDR should be calculated
     * @param fdr the target FDR that is supposed to be returned
     * @param safetyfactor don't report anything if the next higher calculable
     * exceeds the target-FDR by the given factor
     * @param results the passing results will be added to this ArrayList
     * @param TCount Size of the target database
     * @param DCount size of the decoy database
     * @param info all informations needed for this sub group
     * @param isSymmetric do we calculate an FDR for a symmetric or an
     * asymmetric experiment
     * @param setElementFDR
     * @return
     */
    protected <T extends AbstractFDRElement<T>> void subFDRLocal(ArrayList<T> group, ArrayList<T> results, SubGroupFdrInfo info) {
        HashedArrayList<T> ret = new HashedArrayList<T>();
        info.results = ret;
        info.filteredResult = ret;

        int TT = info.TT;
        int TD = info.TD;
        int DD = info.DD;

        double TTSize;
        double TDSize;
        double DDSize;
        double k;

        // all the elements within this group
        int groupSize = group.size();

        // sort them by score
        Collections.sort(group, new Comparator<T>() {

            public int compare(T o1, T o2) {
                return Double.compare(o1.getScore(), o2.getScore());
            }
        });

        int prevTDIndex = groupSize - 1;

        int minTD = 2;
        double minWindow = Math.abs(group.get(0).getScore() - group.get(prevTDIndex).getScore()) / 10.0;

        int lastMin = 0;
        T lastMinE = group.get(0);
        double lastMinScore = lastMinE.getScore();
        int lastMax = 0;
        T lastMaxE = lastMinE;
        double lastMaxScore = lastMinScore;

        // count within the window
        int wTT = lastMinE.isTT() ? 1 : 0;
        int wTD = lastMinE.isTD() ? 1 : 0;
        int wDD = lastMinE.isDD() ? 1 : 0;
        for (T e : group) {

            double centerscore = e.getScore();
            double minscore = centerscore - minWindow;
            double maxscore = centerscore + minWindow;

            // find the highest to be included eleemnt
            while (lastMax < groupSize - 1 && maxscore > group.get(lastMax + 1).getScore()) {
                lastMax++;
                lastMaxE = group.get(lastMax);
                lastMaxScore = lastMaxE.getScore();
                if (lastMaxE.isTT()) {
                    wTT++;
                } else if (lastMaxE.isTD()) {
                    wTD++;
                } else {
                    wDD++;
                }
            }

            // shift to the lowest to be included eleemnt
            while (minscore > lastMinScore
                    && wTD >= minTD
                    && // expand the window until we have more then the required TD
                    wDD * 1.1 < wTD
                    && // also we want some more TD the DD
                    wDD * 1.1 < wTT) {  // and more TT then DD
                lastMin++;
                if (lastMinE.isTT()) {
                    wTT--;
                } else if (lastMinE.isTD()) {
                    wTD--;
                } else {
                    wDD--;
                }
                lastMinE = group.get(lastMin);
                lastMinScore = lastMinE.getScore();
            }
            // do we need to readjust the highest score?
            if (minscore > lastMinScore) {
                double targetMaxScore = lastMinScore - minscore + maxscore;
                while (lastMax < groupSize - 1) {
                    T nextMaxE = group.get(lastMax + 1);
                    double nextScore = nextMaxE.getScore();
                    if (nextScore < targetMaxScore) {
                        lastMax++;
                        lastMaxE = nextMaxE;
                        lastMaxScore = nextScore;
                        if (lastMaxE.isTT()) {
                            wTT++;
                        } else if (lastMaxE.isTD()) {
                            wTD++;
                        } else {
                            wDD++;
                        }
                    } else {
                        break;
                    }
                }
            }
            e.setPEP((wTD - wDD) / (double) wTT);

        }
    }

    /**
     * is a higher score better than a lower score?
     *
     * @return the PSMScoreHighBetter
     */
    public boolean isPSMScoreHighBetter() {
        return PSMScoreHighBetter;
    }

    /**
     * is a higher score better than a lower score?
     *
     * @param PSMScoreHighBetter the PSMScoreHighBetter to set
     */
    public void setPSMScoreHighBetter(boolean PSMScoreHighBetter) {
        this.PSMScoreHighBetter = PSMScoreHighBetter;
    }

    /**
     * @return the targetPepDBSize
     */
    public double getTargetDBSize() {
        return targetPepDBSize;
    }

    /**
     * @param targetPepDBSize the targetPepDBSize to set
     */
    public void setTargetDBSize(double targetDBSize) {
        this.targetPepDBSize = targetDBSize;
    }

    /**
     * @return the decoyPepDBSize
     */
    public double getDecoyDBSize() {
        return decoyPepDBSize;
    }

    /**
     * @param decoyPepDBSize the decoyPepDBSize to set
     */
    public void setDecoyDBSize(double decoyDBSize) {
        this.decoyPepDBSize = decoyDBSize;
    }

    /**
     * @return the targetProtDBSize
     */
    public double getTargetProtDBSize() {
        return targetProtDBSize;
    }

    /**
     * @param targetProtDBSize the targetProtDBSize to set
     */
    public void setTargetProtDBSize(double targetProtDBSize) {
        this.targetProtDBSize = targetProtDBSize;
    }

    /**
     * @return the decoyProtDBSize
     */
    public double getDecoyProtDBSize() {
        return decoyProtDBSize;
    }

    /**
     * @param targetProtDBSize the targetProtDBSize to set
     */
    public void setTargetLinkDBSize(double targetLinkDBSize) {
        this.targetLinkDBSize = targetLinkDBSize;
    }

    /**
     * @param targetProtDBSize the targetProtDBSize to set
     */
    public void setDecoyLinkDBSize(double decoyLinkDBSize) {
        this.decoyLinkDBSize = decoyLinkDBSize;
    }

    /**
     * @return the decoyProtDBSize
     */
    public double getDecoyLinkDBSize() {
        return decoyLinkDBSize;
    }

    /**
     * @return the decoyProtDBSize
     */
    public double getTargetLinkDBSize() {
        return targetLinkDBSize;
    }

    /**
     * @param decoyProtDBSize the decoyProtDBSize to set
     */
    public void setDecoyProtDBSize(double decoyProtDBSize) {
        this.decoyProtDBSize = decoyProtDBSize;
    }

    /**
     * @return the minPepPerProteinGroup
     */
    public int getMinPepPerProteinGroup() {
        return minPepPerProteinGroup;
    }

    /**
     * @param minPepPerProteinGroup the minPepPerProteinGroup to set
     */
    public void setMinPepPerProteinGroup(int minPepPerProteinGroup) {
        this.minPepPerProteinGroup = minPepPerProteinGroup;
    }

    /**
     * @return the minPepPerProteinGroupLink
     */
    public int getMinPepPerProteinGroupLink() {
        return minPepPerProteinGroupLink;
    }

    /**
     * @param minPepPerProteinGroupLink the minPepPerProteinGroupLink to set
     */
    public void setMinPepPerProteinGroupLink(int minPepPerProteinGroupLink) {
        this.minPepPerProteinGroupLink = minPepPerProteinGroupLink;
    }

    /**
     * @return the minPepPerProteinGroupPair
     */
    public int getMinPepPerProteinGroupPair() {
        return minPepPerProteinGroupPair;
    }

    /**
     * @param minPepPerProteinGroupPair the minPepPerProteinGroupPair to set
     */
    public void setMinPepPerProteinGroupPair(int minPepPerProteinGroupPair) {
        this.minPepPerProteinGroupPair = minPepPerProteinGroupPair;
    }

    public String execClass() {
        return this.getClass().getName();
    }

    public void printUsage() {
        System.out.println("java " + execClass() + " " + argList());
        System.out.println(argDescription());
    }

    public String argList() {
        return "--lenghtgroups=A,B,C "
                + "--psmfdr=X "
                + "--pepfdr=X "
                + "--proteinfdr=X "
                + "--linkfdr=X "
                + "--ppifdr=X "
                + "--reportfactor=X"
                + "--maxProteinAmbiguity=X "
                + "--maxLinkAmbiguity=X "
                + "--minPeptidesPerLink=X "
                + "--minPeptidesPerProtein=X "
                + "--minPeptidesPerPPI=X "
                + "--minPeptideLength=X "
                + "--ignoregroups "
                + "--csvOutDir=X "
                + "--csvBaseName=X "
                + "--csvSummaryOnly "
                + "--singleSummary "
                + "--uniquePSMs= "
                + "--outputlocale= "
                + "--boost= ";

    }

    public String argDescription() {
        return "--lengthgroups=A,B,C     how to group peptides by length\n"
                + "--psmfdr=X               the psm-fdr\n"
                + "                         can either be a single value or \n"
                + "                         a range (min,max,stepsize)\n"
                + "--pepfdr=X               the peptide pair fdr\n"
                + "                         can either be a single value or \n"
                + "                         a range (min,max,stepsize)\n"
                + "--proteinfdr=X           the protein pair fdr\n"
                + "                         can either be a single value or \n"
                + "                         a range (min,max,stepsize)\n"
                + "--linkfdr=X              residue pair fdr\n"
                + "                         can either be a single value or \n"
                + "                         a range (min,max,stepsize)\n"
                + "--ppifdr=X               protein pair fdr\n"
                + "                         can either be a single value or \n"
                + "                         a range (min,max,stepsize)\n"
                + "--reportfactor=X         is ignored\n"
                + "--maxProteinAmbiguity=X  any peptide that can be in more \n"
                + "                         then this number of proteins will\n"
                + "                         be ignored\n"
                + "--maxLinkAmbiguity=X     links that have more then this \n"
                + "                         number of possible residue\n"
                + "                         combinationswill be ignored\n"
                + "--minPeptidesPerLink=X   only links that have at least this \n"
                + "                         number of unique peptide pairs \n"
                + "                         supporting them will be considered\n"
                + "--minPeptidesPerProtein=X Only proteins that have at least\n"
                + "                         this number of unique peptide pairs\n"
                + "                         supporting them will be considered\n"
                + "--minPeptidesPerPPI=X    only protein pairs that have at \n"
                + "                         least this number of unique peptide\n"
                + "                         pairs supporting them will be \n"
                + "                         considered\n"
                + "--minPeptideLength=X     only accept psms where both peptides\n"
                + "                         have at least this many residues\n"
                + "--ignoregroups           don't do any grouping during FDR\n"
                + "                         calculation\n"
                + "--csvOutDir=X            where to write the output files\n"
                + "--csvBaseName=X          each file will be prepended with \n"
                + "                         this name"
                + "--csvSummaryOnly         don;t write the actuall results but\n"
                + "                         only the summary\n"
                + "--singleSummary          if fdrs where given in ranges all\n"
                + "                         summary files will be written into a\n"
                + "                         sinlge file\n"
                + "--uniquePSMs=X           filtr PSMs to unique PSMs \n"
                + "                         options are true,false,1,0\n"
                + "--outputlocale=          numbers in csv-files are writen \n"
                + "                         according to this locale\n"
                + "--boost=(pep|link|prot)  boost results on the given level\n";

    }

    public boolean setOutputLocale(String locale) {
        locale = locale.toLowerCase();
        Locale sl = null;
        for (Locale l : Locale.getAvailableLocales()) {
            if (l.toString().toLowerCase().contentEquals(locale)) {
                setOutputLocale(l);
                return true;
            }
            if (l.getDisplayName().toLowerCase().contentEquals(locale)) {
                sl = l;
            }
            if (l.getCountry().toLowerCase().contentEquals(locale)) {
                sl = l;
            }
            if (l.getDisplayScript().toLowerCase().contentEquals(locale)) {
                sl = l;
            }
            if (l.getDisplayLanguage().toLowerCase().contentEquals(locale)) {
                sl = l;
            }
        }
        if (sl == null) {
            return false;
        }
        setOutputLocale(sl);
        return true;
    }

    public void setOutputLocale(Locale locale) {
        this.outputlocale = locale;
        this.numberFormat = NumberFormat.getInstance(locale);
        numberFormat = NumberFormat.getNumberInstance(locale);
        DecimalFormat fformat = (DecimalFormat) numberFormat;
        fformat.setGroupingUsed(false);
//        DecimalFormatSymbols symbols=fformat.getDecimalFormatSymbols();
//        fformat.setMaximumFractionDigits(6);
//        localNumberDecimalSeparator= ""+symbols.getDecimalSeparator();
    }

    protected String d2s(double d) {
        return numberFormat.format(d);
    }

    protected String i2s(int i) {
        return numberFormat.format(i);
    }

    /**
     * @return the psmFDRSetting
     */
    public double[] getPsmFDRSetting() {
        return psmFDRSetting;
    }

    /**
     * @param psmFDRSetting the psmFDRSetting to set
     */
    public void setPsmFDRSetting(double from, double to, double step) {
        this.psmFDRSetting = new double[]{from, to, step};
    }

    /**
     * @return the peptidePairFDRSetting
     */
    public double[] getPeptidePairFDRSetting() {
        return peptidePairFDRSetting;
    }

    /**
     * @param peptidePairFDRSetting the peptidePairFDRSetting to set
     */
    public void setPeptidePairFDRSetting(double from, double to, double step) {
        this.peptidePairFDRSetting = new double[]{from, to, step};
    }

    /**
     * @return the ProteinGroupFDRSetting
     */
    public double[] getProteinGroupFDRSetting() {
        return ProteinGroupFDRSetting;
    }

    /**
     * @param ProteinGroupFDRSetting the ProteinGroupFDRSetting to set
     */
    public void setProteinGroupFDRSetting(double from, double to, double step) {
        this.ProteinGroupFDRSetting = new double[]{from, to, step};
    }

    /**
     * @return the linkFDRSetting
     */
    public double[] getLinkFDRSetting() {
        return linkFDRSetting;
    }

    /**
     * @param linkFDRSetting the linkFDRSetting to set
     */
    public void setLinkFDRSetting(double from, double to, double step) {
        this.linkFDRSetting = new double[]{from, to, step};;
    }

    /**
     * @return the ppiFDRSetting
     */
    public double[] getPpiFDRSetting() {
        return ppiFDRSetting;
    }

    /**
     * @param ppiFDRSetting the ppiFDRSetting to set
     */
    public void setPpiFDRSetting(double from, double to, double step) {
        this.ppiFDRSetting = new double[]{from, to, step};
    }

    /**
     * @return the safetyFactorSetting
     */
    public double getSafetyFactorSetting() {
        return safetyFactorSetting;
    }

    /**
     * @param safetyFactorSetting the safetyFactorSetting to set
     */
    public void setSafetyFactorSetting(double safetyFactorSetting) {
        this.safetyFactorSetting = safetyFactorSetting;
    }

    /**
     * @return the ignoreGroupsSetting
     */
    public boolean isIgnoreGroupsSetting() {
        return ignoreGroupsSetting;
    }

    /**
     * @param ignoreGroupsSetting the ignoreGroupsSetting to set
     */
    public void setIgnoreGroupsSetting(boolean ignoreGroupsSetting) {
        this.ignoreGroupsSetting = ignoreGroupsSetting;
    }

    /**
     * @param summaryOnly when writing out put files only write the summary file
     */
    public void setCSVSummaryOnly(boolean summaryOnly) {
        this.csvSummaryOnly = summaryOnly;
    }

    /**
     * @param summaryOnly when writing out put files only write the summary file
     */
    public boolean getCSVSummaryOnly() {
        return this.csvSummaryOnly;
    }

    /**
     * @param summaryOnly when writing out put files only write the summary file
     */
    public void setCSVSingleSummary(boolean singleSummary) {
        this.singleSummary = singleSummary;
    }

    /**
     * @return the csvOutDirSetting
     */
    public String getCsvOutDirSetting() {
        return csvOutDirSetting;
    }

    /**
     * @param csvOutDirSetting the csvOutDirSetting to set
     */
    public void setCsvOutDirSetting(String csvOutDirSetting) {
        this.csvOutDirSetting = csvOutDirSetting;
    }

    /**
     * @return the csvOutBaseSetting
     */
    public String getCsvOutBaseSetting() {
        return csvOutBaseSetting;
    }

    /**
     * @param csvOutBaseSetting the csvOutBaseSetting to set
     */
    public void setCsvOutBaseSetting(String csvOutBaseSetting) {
        this.csvOutBaseSetting = csvOutBaseSetting;
    }

    public String[] parseArgs(String[] argv, FDRSettings settings) {
        ArrayList<String> unknown = new ArrayList<String>();
        int[] lengthgroups = new int[]{4};
        double[] psmFDR = new double[]{1, 1, 1};
        double[] pepFDR = new double[]{1, 1, 1};
        double[] protFDR = new double[]{1, 1, 1};
        double[] linkFDR = new double[]{1, 1, 1};
        double[] ppiFDR = new double[]{1, 1, 1};
        int maxLinkAmbiguity = 0;
        int maxProteinGroupAmbiguity = 0;
        int minPepPerLink = 1;
        int minPepPerProtein = 1;
        int minPepPerPPI = 1;
        int minPepLength = 6;
        boolean ignoreGroups = false;
        boolean csvsummaryonly = false;
        boolean csvsinglesummary = false;
        String csvdir = null;
        String csvBase = null;
        int fdrDigits = commandlineFDRDigits;

        for (String arg : argv) {

            if (arg.startsWith("--lenghtgroups=") || arg.startsWith("--lengthgroups=")) {

                String[] slen = arg.substring(arg.indexOf("=") + 1).trim().split(",");
                lengthgroups = new int[slen.length];

                for (int i = 0; i < slen.length; i++) {
                    lengthgroups[i] = Integer.parseInt(slen[i]);
                }

            } else if (arg.toLowerCase().startsWith("--boost=")) {
                String what = arg.substring("--boost=".length()).trim().toLowerCase();
                if (what.endsWith("s")) {
                    what = what.substring(0, what.length() - 1);
                }
                if (what.contentEquals("pep")
                        || what.contentEquals("peptidepair")
                        || what.contentEquals("peppair")) {
                    this.maximizeWhat = FDRLevel.PEPTIDE_PAIR;
                } else if (what.contentEquals("link")
                        || what.contentEquals("residuepair")) {
                    this.maximizeWhat = FDRLevel.PROTEINGROUPLINK;
                } else if (what.contentEquals("prot")
                        || what.contentEquals("protein")
                        || what.contentEquals("proteinpair")
                        || what.contentEquals("ppi")) {
                    this.maximizeWhat = FDRLevel.PROTEINGROUPPAIR;
                }

            } else if (arg.startsWith("--reportfactor=")) {

            } else if (arg.startsWith("--psmfdr=")) {
                double from, to, step;
                String spsm = arg.substring(arg.indexOf("=") + 1).trim();
                String[] spsm_seq = spsm.split(",");
                if (spsm_seq.length == 1) {
                    from = to = Double.parseDouble(spsm) / 100;
                    step = 1;
                    int thisDigits = spsm.trim().replace("^[^\\.]*", "").trim().length();
                    if (thisDigits + 2 > fdrDigits) {
                        fdrDigits = thisDigits + 2;
                    }
                } else {
                    from = Double.parseDouble(spsm_seq[0]) / 100;
                    to = Double.parseDouble(spsm_seq[1]) / 100;
                    step = Double.parseDouble(spsm_seq[2]) / 100;

                    for (int i = 0; i < 3; i++) {
                        int thisDigits = spsm_seq[i].trim().replace("^[^\\.]*", "").trim().length();
                        if (thisDigits + 2 > fdrDigits) {
                            fdrDigits = thisDigits + 2;
                        }
                    }

                }
                psmFDR = new double[]{from, to, step};

            } else if (arg.startsWith("--pepfdr=")) {

                double from, to, step;
                String spsm = arg.substring(arg.indexOf("=") + 1).trim();
                String[] spsm_seq = spsm.split(",");
                if (spsm_seq.length == 1) {
                    from = to = Double.parseDouble(spsm) / 100;
                    step = 1;
                    int thisDigits = spsm.trim().replace("^[^\\.]*", "").trim().length();
                    if (thisDigits + 2 > fdrDigits) {
                        fdrDigits = thisDigits + 2;
                    }
                } else {
                    from = Double.parseDouble(spsm_seq[0]) / 100;
                    to = Double.parseDouble(spsm_seq[1]) / 100;
                    step = Double.parseDouble(spsm_seq[2]) / 100;
                    for (int i = 0; i < 3; i++) {
                        int thisDigits = spsm_seq[i].trim().replace("^[^\\.]*", "").trim().length();
                        if (thisDigits + 2 > fdrDigits) {
                            fdrDigits = thisDigits + 2;
                        }
                    }
                }
                pepFDR = new double[]{from, to, step};

            } else if (arg.startsWith("--proteinfdr=")) {

                double from, to, step;
                String spsm = arg.substring(arg.indexOf("=") + 1).trim();
                String[] spsm_seq = spsm.split(",");
                if (spsm_seq.length == 1) {
                    from = to = Double.parseDouble(spsm) / 100;
                    step = 1;
                    int thisDigits = spsm.trim().replace("^[^\\.]*", "").trim().length();
                    if (thisDigits + 2 > fdrDigits) {
                        fdrDigits = thisDigits + 2;
                    }
                } else {
                    from = Double.parseDouble(spsm_seq[0]) / 100;
                    to = Double.parseDouble(spsm_seq[1]) / 100;
                    step = Double.parseDouble(spsm_seq[2]) / 100;
                    for (int i = 0; i < 3; i++) {
                        int thisDigits = spsm_seq[i].trim().replace("^[^\\.]*", "").trim().length();
                        if (thisDigits + 2 > fdrDigits) {
                            fdrDigits = thisDigits + 2;
                        }
                    }
                }
                protFDR = new double[]{from, to, step};

            } else if (arg.startsWith("--linkfdr=")) {

                double from, to, step;
                String spsm = arg.substring(arg.indexOf("=") + 1).trim();
                String[] spsm_seq = spsm.split(",");
                if (spsm_seq.length == 1) {
                    from = to = Double.parseDouble(spsm) / 100;
                    step = 1;
                    int thisDigits = spsm.trim().replace("^[^\\.]*", "").trim().length();
                    if (thisDigits + 2 > fdrDigits) {
                        fdrDigits = thisDigits + 2;
                    }
                } else {
                    from = Double.parseDouble(spsm_seq[0]) / 100;
                    to = Double.parseDouble(spsm_seq[1]) / 100;
                    step = Double.parseDouble(spsm_seq[2]) / 100;
                    for (int i = 0; i < 3; i++) {
                        int thisDigits = spsm_seq[i].trim().replace("^[^\\.]*", "").trim().length();
                        if (thisDigits + 2 > fdrDigits) {
                            fdrDigits = thisDigits + 2;
                        }
                    }
                }
                linkFDR = new double[]{from, to, step};

            } else if (arg.startsWith("--ppifdr=")) {

                double from, to, step;
                String spsm = arg.substring(arg.indexOf("=") + 1).trim();
                String[] spsm_seq = spsm.split(",");
                if (spsm_seq.length == 1) {
                    from = to = Double.parseDouble(spsm) / 100;
                    step = 1;
                    int thisDigits = spsm.trim().replace("^[^\\.]*", "").trim().length();
                    if (thisDigits + 2 > fdrDigits) {
                        fdrDigits = thisDigits + 2;
                    }
                } else {
                    from = Double.parseDouble(spsm_seq[0]) / 100;
                    to = Double.parseDouble(spsm_seq[1]) / 100;
                    step = Double.parseDouble(spsm_seq[2]) / 100;
                    for (int i = 0; i < 3; i++) {
                        int thisDigits = spsm_seq[i].trim().replace("^[^\\.]*", "").trim().length();
                        if (thisDigits + 2 > fdrDigits) {
                            fdrDigits = thisDigits + 2;
                        }
                    }
                }
                ppiFDR = new double[]{from, to, step};

            } else if (arg.startsWith("--maxProteinAmbiguity=")) {

                String spsm = arg.substring(arg.indexOf("=") + 1).trim();
                maxProteinGroupAmbiguity = Integer.parseInt(spsm);

            } else if (arg.startsWith("--maxLinkAmbiguity=")) {

                String spsm = arg.substring(arg.indexOf("=") + 1).trim();
                maxLinkAmbiguity = Integer.parseInt(spsm);

            } else if (arg.startsWith("--minPeptidesPerLink=")) {

                String spsm = arg.substring(arg.indexOf("=") + 1).trim();
                minPepPerLink = Integer.parseInt(spsm);

            } else if (arg.startsWith("--minPeptidesPerProtein=")) {

                String spsm = arg.substring(arg.indexOf("=") + 1).trim();
                minPepPerProtein = Integer.parseInt(spsm);

            } else if (arg.startsWith("--minPeptidesPerPPI=")) {

                String spsm = arg.substring(arg.indexOf("=") + 1).trim();
                minPepPerPPI = Integer.parseInt(spsm);

            } else if (arg.startsWith("--minPeptideLength=")) {

                String spsm = arg.substring(arg.indexOf("=") + 1).trim();
                minPepLength = Integer.parseInt(spsm);

            } else if (arg.equals("--ignoregroups")) {

                ignoreGroups = true;
            } else if (arg.startsWith("--csvOutDir=")) {

                csvdir = arg.substring(arg.indexOf("=") + 1);
                if (csvBase == null) {
                    csvBase = "FDR";
                }

            } else if (arg.startsWith("--csvBaseName=")) {

                csvBase = arg.substring(arg.indexOf("=") + 1);
                if (csvdir == null) {
                    csvdir = ".";
                }

            } else if (arg.equals("--csvSummaryOnly")) {

                csvsummaryonly = true;

            } else if (arg.equals("--singleSummary")) {

                csvsummaryonly = true;
                csvsinglesummary = true;

            } else if (arg.startsWith("--outputlocale=")) {
                String l = arg.substring(arg.indexOf("=") + 1).trim();
                setOutputLocale(l);

            } else if (arg.startsWith("--uniquePSMs=")) {
                String bool = arg.substring(arg.indexOf("=") + 1).trim();
                boolean filter = bool.matches("(?i)^(T|1(\\.0*)?|-1(\\.0*)?|TRUE|Y|YES|\\+)$");
                settings.setFilterToUniquePSM(filter);
                setFilterUniquePSMs(filter);
            } else {
                unknown.add(arg);
            }

        }
        this.setLengthGroups(lengthgroups);

        commandlineFDRDigits = fdrDigits;
        setPsmFDRSetting(psmFDR[0], psmFDR[1], psmFDR[2]);
        settings.setPSMFDR(psmFDR[0]);
        setPeptidePairFDRSetting(pepFDR[0], pepFDR[1], pepFDR[2]);
        settings.setPeptidePairFDR(pepFDR[0]);
        setProteinGroupFDRSetting(protFDR[0], protFDR[1], protFDR[2]);
        settings.setProteinGroupFDR(protFDR[0]);
        setLinkFDRSetting(linkFDR[0], linkFDR[1], linkFDR[2]);
        settings.setProteinGroupLinkFDR(linkFDR[0]);
        setPpiFDRSetting(ppiFDR[0], ppiFDR[1], ppiFDR[2]);
        settings.setProteinGroupPairFDR(ppiFDR[0]);
        //setSafetyFactorSetting(reportfactor);
        //settings.setReportFactor(reportfactor);
        setIgnoreGroupsSetting(ignoreGroups);
        setCSVSummaryOnly(csvsummaryonly);
        setCSVSingleSummary(csvsinglesummary);
        setCsvOutBaseSetting(csvBase);
        setCsvOutDirSetting(csvdir);

        setMaximumLinkAmbiguity(maxLinkAmbiguity);
        settings.setMaxLinkAmbiguity(maxLinkAmbiguity);
        setMaximumProteinAmbiguity(maxProteinGroupAmbiguity);
        settings.setMaxProteinAmbiguity(maxProteinGroupAmbiguity);
        setMinPepPerProteinGroup(minPepPerProtein);
        settings.setMinProteinPepCount(minPepPerProtein);
        setMinPepPerProteinGroupLink(minPepPerLink);
        settings.setMinLinkPepCount(minPepPerLink);
        setMinPepPerProteinGroupPair(minPepPerPPI);
        settings.setMinPPIPepCount(minPepPerPPI);
        setMinimumPeptideLength(minPepLength);
        settings.setMinPeptideLength(minPepLength);

        String[] ret = new String[unknown.size()];
        ret = unknown.toArray(ret);
        return ret;
    }

    protected void reset() {
        // reset the counts
        PeptidePair.PEPTIDEPAIRCOUNT = 0;
        ProteinGroup.PROTEINGROUPCOUNT = 0;
        ProteinGroupPair.PROTEINGROUPPAIRCOUNT = 0;
        ProteinGroupLink.LINKCOUNT = 0;

        m_linearPepCount = null;
        m_XLPepCount = null;
        m_linearPSMCount = null;
        m_XLPSMCount = null;

        for (PSM psm : getAllPSMs()) {
            psm.setFDRGroup();
            psm.getAdditionalFDRGroups().remove(psm.getRun());
            psm.getAdditionalFDRGroups().remove("" + psm.getScan());
            psm.setFDR(Double.MAX_VALUE);
            psm.setFdrPeptidePair(null);
            psm.reset();
        }

        for (Protein p : allProteins) {
            p.resetFDR();
        }

        for (Peptide p : allPeptides) {
            p.resetFDR();
        }

        if (groupPSMsByRun) {
            for (PSM psm : getAllPSMs()) {
                psm.getAdditionalFDRGroups().add(psm.getRun());
            }
        } else {
            for (PSM psm : getAllPSMs()) {
                psm.getAdditionalFDRGroups().remove(psm.getRun());
            }
        }
    }

    public abstract String getSource();

    protected String getPSMOutputLine(PSM psm, String seperator) {
        return RArrayUtils.toString(getPSMOutputLine(psm), seperator);
    }

    protected ArrayList<String> getPSMOutputLine(PSM pp) {
        PeptidePair pep = pp.getFdrPeptidePair();
        ProteinGroupLink l = pep != null ? pep.getFdrLink() : null;
        ProteinGroupPair ppi = l != null ? l.getFdrPPI() : null;

        Peptide pep1 = pp.getPeptide1();
        Peptide pep2 = pp.getPeptide2();

        ProteinGroup pg1 = pep1.getProteinGroup();
        ProteinGroup pg2 = pep2.getProteinGroup();

        int pepLink1 = pp.getPeptideLinkSite1();
        int pepLink2 = pp.getPeptideLinkSite2();

        int pepLength1 = pp.getPeptideLength1();
        int pepLength2 = pp.getPeptideLength2();

        String pepSeq1 = getPeptideSequence(pep1);
        String pepSeq2 = getPeptideSequence(pep2);

        StringBuilder sbaccessions = new StringBuilder();
        StringBuilder sbdescriptions = new StringBuilder();
        StringBuilder sbPositions = new StringBuilder();
        StringBuilder sbProtLink = new StringBuilder();
        peptidePositionsToPSMOutString(pep1.getPositions(), sbaccessions, sbdescriptions, sbPositions, sbProtLink, pepLink1);
        String accessions1 = sbaccessions.toString();
        String descriptions1 = sbdescriptions.toString();
        String positons1 = sbPositions.toString();
        String proteinLinkPositons1 = pepLink1 > 0 ? sbProtLink.toString() : "";

        sbaccessions.setLength(0);
        sbdescriptions.setLength(0);
        sbPositions.setLength(0);
        sbProtLink.setLength(0);
        peptidePositionsToPSMOutString(pep2.getPositions(), sbaccessions, sbdescriptions, sbPositions, sbProtLink, pepLink2);
        String accessions2 = sbaccessions.toString();
        String descriptions2 = sbdescriptions.toString();
        String positons2 = sbPositions.toString();
        String proteinLinkPositons2 = pepLink2 > 0 ? sbProtLink.toString() : "";

        String run = pp.getRun();
        String scan = pp.getScan();
        if (run == null) {
            run = "";
        }
        if (scan == null) {
            scan = "";
        }
        ArrayList<String> ret = new ArrayList<>(37);
        ret.add(pp.getPsmID());
        ret.add(run);
        ret.add(scan);
        ret.add(pp.getPeakListName() == null ? "" : pp.getPeakListName());
        ret.add(pp.getFileScanIndex() == null ? "" : i2s(pp.getFileScanIndex()));
        ret.add(accessions1);
        ret.add(descriptions1);
        ret.add(Boolean.toString(pep1.isDecoy()));
        ret.add(accessions2);
        ret.add(descriptions2);
        ret.add(Boolean.toString(pep2.isDecoy()));
        ret.add(pepSeq1);
        ret.add(pepSeq2);
        ret.add(positons1);
        ret.add(positons2);
        ret.add((pepLength1 == 0 ? "" : i2s(pepLength1)));
        ret.add((pepLength2 == 0 ? "" : i2s(pepLength2)));
        ret.add(i2s(pepLink1));
        ret.add(i2s(pepLink2));
        ret.add(proteinLinkPositons1);
        ret.add(proteinLinkPositons2);
        ret.add(i2s(pp.getCharge()));
        ret.add(pp.getCrosslinker());
        ret.add(Double.isNaN(pp.getCrosslinkerModMass()) ? "" : d2s(pp.getCrosslinkerModMass()));
        for (String name : extraColumns) {
            Object v = pp.getOtherInfo(name);
            if (v == null) {
                ret.add("");
            } else if (v instanceof Number) {
                ret.add(d2s(((Number) v).doubleValue()));
            } else {
                ret.add(v.toString());
            }
        }
        ret.add(d2s(pp.getScore()));
        ret.add(Boolean.toString(pp.isDecoy()));
        ret.add(Boolean.toString(pp.isTT()));
        ret.add(Boolean.toString(pp.isTD()));
        ret.add(Boolean.toString(pp.isDD()));
        ret.add(pp.getFDRGroup());
        ret.add(d2s(pp.getFDR()));
        ret.add(d2s(pp.getEstimatedFDR()));
        ret.add(d2s(pp.getPEP()));
        ret.add("");
        ret.add((pep == null ? "" : d2s(pep.getFDR())));
        ret.add((pp.getFdrProteinGroup1() == null ? "" : d2s(pp.getFdrProteinGroup1().getFDR())));
        ret.add((pp.getFdrProteinGroup2() == null ? "" : d2s(pp.getFdrProteinGroup2().getFDR())));
        ret.add((l == null ? "" : d2s(l.getFDR())));
        ret.add((ppi == null ? "" : d2s(ppi.getFDR())));
        ret.add((pep == null ? "" : i2s(pep.getPeptidePairID())));
        ret.add((l == null ? "" : i2s(l.getLinkID())));
        ret.add((ppi == null ? "" : i2s(ppi.getProteinGroupPairID())));
        ret.add(pp.getInfo());

        return ret;
    }

    protected void peptidePositionsToPSMOutString(HashMap<Protein, IntArrayList> positions, StringBuilder sbaccessions, StringBuilder sbdescriptions, StringBuilder sbPositions, StringBuilder sbProtLink, int peplink) {
        for (Protein p : positions.keySet()) {
            for (int i : positions.get(p)) {
                sbaccessions.append(";").append(p.getAccession());
                sbdescriptions.append(";").append(p.getDescription());
                sbPositions.append(";").append(i2s(i));
                sbProtLink.append(";").append(i2s(i + peplink - 1));
            }
        }
        sbaccessions.deleteCharAt(0);
        sbdescriptions.deleteCharAt(0);
        sbPositions.deleteCharAt(0);
        sbProtLink.deleteCharAt(0);
    }

    protected <T extends RandomAccess & Collection> String getPSMHeader(String seperator) {
        return RArrayUtils.toString(getPSMHeader(), seperator);
    }

    protected <T extends RandomAccess & Collection> ArrayList<String> getPSMHeader() {
        ArrayList<String> ret = new ArrayList<String>(RArrayUtils.toCollection(new String[]{"PSMID", "run", "scan", "PeakListFileName", "ScanId", "Protein1", "Description1",
            "Decoy1", "Protein2", "Description2", "Decoy2",
            "PepSeq1", "PepSeq2", "PepPos1", "PepPos2",
            "PeptideLength1", "PeptideLength2", "LinkPos1", "LinkPos2",
            "ProteinLinkPos1", "ProteinLinkPos2", "Charge", "Crosslinker", "CrosslinkerModMass"}));

        if (allPSMs.size() > 0 && PSM.getOtherInfoNames().length > 0) {
            extraColumns = new ArrayList<>(MyArrayUtils.toCollection(PSM.getOtherInfoNames()));
            extraColumns.sort(new Comparator<String>() {
                @Override
                public int compare(String arg0, String arg1) {
                    return arg0.compareTo(arg1);
                }
            });
            ret.addAll(extraColumns);
        } else {
            extraColumns = new ArrayList<>(0);
        }

        ret.addAll(RArrayUtils.toCollection(new String[]{"Score",
            "isDecoy", "isTT", "isTD", "isDD", "fdrGroup", "fdr", "ifdr", "PEP",
            "", "PeptidePairFDR", "Protein1FDR", "Protein2FDR",
            "LinkFDR", "PPIFDR",
            "peptide pair id", "link id", "ppi id", "info"}));
        return ret;
    }

    protected ArrayList<String> getXLPepsHeader() {
        ArrayList<String> ret = new ArrayList<String>(RArrayUtils.toCollection(new String[]{
            "PeptidePairID",
            "PSMIDs",
            "Protein1",
            "Description1",
            "Decoy1",
            "Protein2",
            "Description2",
            "Decoy2",
            "Peptide1",
            "Peptide2",
            "Start1",
            "Start2",
            "FromSite",
            "ToSite",
            "FromProteinSite",
            "ToProteinSite",
            "psmID",
            "Crosslinker",
            "Score",
            "isDecoy",
            "isTT",
            "isTD",
            "isDD",
            "fdrGroup",
            "fdr",
            "ifdr",
            "PEP",
            "",
            "Protein1FDR",
            "Protein2FDR",
            "LinkFDR",
            "PPIFDR",
            "",
            "link id",
            "ppi id"}));
        for (String r : runs) {
            ret.add(r);
        }
        return ret;
    }

    protected String getXLPepsHeader(String seperator) {
        return RArrayUtils.toString(getXLPepsHeader(), seperator);
    }

    protected String getLinearPepsHeader(String seperator) {
        return RArrayUtils.toString(getLinearPepsHeader(), seperator);
    }

    protected ArrayList<String> getLinearPepsHeader() {
        ArrayList<String> ret = new ArrayList<>();
        ret.add("PeptidePairID");
        ret.add("PSMIDs");
        ret.add("Protein");
        ret.add("Description");
        ret.add("Decoy");
        ret.add("Peptide");
        ret.add("psmID");
        ret.add("Score");
        ret.add("isDecoy");
        ret.add("fdrGroup");
        ret.add("fdr");
        ret.add("ifdr");
        ret.add("PEP");
        ret.add("");
        ret.add("ProteinFDR");
        for (String r : runs) {
            ret.add(r);
        }

        return ret;
    }

    protected String getLinearPepeptideOutputLine(PeptidePair pp, String seperator) {
        return RArrayUtils.toString(getLinearPepeptideOutputLine(pp), seperator);
    }

    protected ArrayList<String> getLinearPepeptideOutputLine(PeptidePair pp) {
        String[] psmids = pp.getPSMids();
        ProteinGroup pg1 = pp.getPeptide1().getProteinGroup();
        ArrayList<String> ret = new ArrayList<String>();

        //           "ID" + "PSMIDs" + "Protein" + "Decoy" +  "Peptide" +  "psmID" + "Score" + "isDecoy" + "fdrGroup" + "fdr" +  + "ProteinFDR");
        ret.add(i2s(pp.getPeptidePairID()));
        ret.add(RArrayUtils.toString(psmids, ";"));
        ret.add(pg1.accessions());
        ret.add(pg1.descriptions());
        ret.add(Boolean.toString(pp.getPeptide1().isDecoy()));
        ret.add(getPeptideSequence(pp.getPeptide1()));
        ret.add(pp.getTopPSMIDs());
        ret.add(d2s(pp.getScore()));
        ret.add(Boolean.toString(pp.isDecoy()));
        ret.add(pp.getFDRGroup());
        ret.add(d2s(pp.getFDR()));
        ret.add(d2s(pp.getEstimatedFDR()));
        ret.add(d2s(pp.getPEP()));
        ret.add("");
        ret.add((pp.getFdrProteinGroup1() == null ? "" : d2s(pp.getFdrProteinGroup1().getFDR())));

        HashSet<String> linkRuns = new HashSet<>();
        HashMap<String, Double> runScore = new HashMap<>();

        for (PSM psm : pp.getAllPSMs()) {
            for (PSM upsm : psm.getRepresented()) {
                psmToRun(upsm, linkRuns, runScore);
            }
        }

        for (String r : runs) {
            Double d = runScore.get(r);
            if (d == null) {
                ret.add("");
            } else {
                ret.add(d2s(d));
            }
        }
        return ret;
    }

    protected String getXlPepeptideOutputLine(PeptidePair pp, String seperator) {
        return RArrayUtils.toString(getXlPepeptideOutputLine(pp), seperator);
    }

    protected ArrayList<String> getXlPepeptideOutputLine(PeptidePair pp) {
        ProteinGroupLink l = pp.getFdrLink();
        ProteinGroupPair ppi = l.getFdrPPI();
        String[] psmids = pp.getPSMids();
        ProteinGroup pg1 = pp.getPeptide1().getProteinGroup();
        ProteinGroup pg2 = pp.getPeptide2().getProteinGroup();
        ArrayList<String> ret = new ArrayList<String>();
        StringBuilder sbaccessions = new StringBuilder();
        StringBuilder sbdescriptions = new StringBuilder();
        StringBuilder sbPositions = new StringBuilder();
        StringBuilder sbProtLink = new StringBuilder();
        peptidePositionsToPSMOutString(pp.getPeptide1().getPositions(), sbaccessions, sbdescriptions, sbPositions, sbProtLink, pp.getPeptideLinkSite1());
        String accessions1 = sbaccessions.toString();
        String descriptions1 = sbdescriptions.toString();
        String positons1 = sbPositions.toString();
        String proteinLinkPositons1 = pp.getPeptideLinkSite1() > 0 ? sbProtLink.toString() : "";

        sbaccessions.setLength(0);
        sbdescriptions.setLength(0);
        sbPositions.setLength(0);
        sbProtLink.setLength(0);
        peptidePositionsToPSMOutString(pp.getPeptide2().getPositions(), sbaccessions, sbdescriptions, sbPositions, sbProtLink, pp.getPeptideLinkSite2());
        String accessions2 = sbaccessions.toString();
        String descriptions2 = sbdescriptions.toString();
        String positons2 = sbPositions.toString();
        String proteinLinkPositons2 = pp.getPeptideLinkSite2() > 0 ? sbProtLink.toString() : "";
        //                    try {
        ret.add(i2s(pp.getPeptidePairID()));
        ret.add(RArrayUtils.toString(psmids, ";"));
        ret.add(accessions1);
        ret.add(descriptions1);
        ret.add("" + pp.getPeptide1().isDecoy());
        ret.add(accessions2);
        ret.add(descriptions2);
        ret.add("" + pp.getPeptide2().isDecoy());
        ret.add(getPeptideSequence(pp.getPeptide1()));
        ret.add(getPeptideSequence(pp.getPeptide2()));
        ret.add(positons1);
        ret.add(pp.getPeptide2() == Peptide.NOPEPTIDE || pp.getPeptide1() == Peptide.NOPEPTIDE ? "" : positons2);
        ret.add(pp.getPeptide2() == Peptide.NOPEPTIDE || pp.getPeptide1() == Peptide.NOPEPTIDE ? "" : i2s(pp.getPeptideLinkSite1()));
        ret.add(pp.getPeptide2() == Peptide.NOPEPTIDE || pp.getPeptide1() == Peptide.NOPEPTIDE ? "" : i2s(pp.getPeptideLinkSite2()));
        ret.add(pp.getPeptide2() == Peptide.NOPEPTIDE || pp.getPeptide1() == Peptide.NOPEPTIDE ? "" : proteinLinkPositons1);
        ret.add(pp.getPeptide2() == Peptide.NOPEPTIDE || pp.getPeptide1() == Peptide.NOPEPTIDE ? "" : proteinLinkPositons2);
        ret.add(pp.getTopPSMIDs());
        ret.add(pp.getCrosslinker());
        ret.add(d2s(pp.getScore()));
        ret.add("" + pp.isDecoy());
        ret.add("" + pp.isTT());
        ret.add("" + pp.isTD());
        ret.add("" + pp.isDD());
        ret.add(pp.getFDRGroup());
        ret.add(d2s(pp.getFDR()));
        ret.add(d2s(pp.getEstimatedFDR()));
        ret.add(d2s(pp.getPEP()));
        ret.add("");
        ret.add((pp.getFdrProteinGroup1() == null ? "" : d2s(pp.getFdrProteinGroup1().getFDR())));
        ret.add((pp.getFdrProteinGroup2() == null ? "" : d2s(pp.getFdrProteinGroup2().getFDR())));
        ret.add(d2s(l.getFDR()));
        ret.add(d2s(ppi.getFDR()));
        ret.add("");
        ret.add(i2s(l.getLinkID()));
        ret.add(i2s(ppi.getProteinGroupPairID()));

        HashSet<String> linkRuns = new HashSet<>();
        HashMap<String, Double> runScore = new HashMap<>();

        for (PSM psm : pp.getAllPSMs()) {
            for (PSM upsm : psm.getRepresented()) {
                psmToRun(upsm, linkRuns, runScore);
            }
        }

        for (String r : runs) {
            Double d = runScore.get(r);
            if (d == null) {
                ret.add("");
            } else {
                ret.add(d2s(d));
            }
        }

        return ret;
    }

    protected String getLinkOutputLine(ProteinGroupLink l, String seperator) {
        return RArrayUtils.toString(getLinkOutputLine(l), seperator);
    }

    protected ArrayList<String> getLinkOutputLine(ProteinGroupLink l) {
        int[] pepids = l.getPeptidePairIDs();
        String[] psmids = l.getPSMIDs();
        HashSet<String> linkRuns = new HashSet<>();
        HashMap<String, Double> runScore = new HashMap<>();

        ProteinGroupPair ppi = l.getFdrPPI();
        double top_pepfdr = Double.MAX_VALUE;
        double top_psmfdr = Double.MAX_VALUE;
        HashSet<String> xl = new HashSet<String>();
        for (PeptidePair pp : l.getPeptidePairs()) {
            if (pp.getFDR() < top_pepfdr) {
                top_pepfdr = pp.getFDR();
            }
            xl.add(pp.getCrosslinker());
            for (PSM psm : pp.getTopPSMs()) {
                if (psm.getFDR() < top_psmfdr) {
                    top_psmfdr = psm.getFDR();
                }
            }
            for (PSM psm : pp.getAllPSMs()) {
                for (PSM upsm : psm.getRepresented()) {
                    psmToRun(upsm, linkRuns, runScore);
                }
            }
        }
        ProteinGroup pg1 = l.getProteinGroup1();
        ProteinGroup pg2 = l.getProteinGroup2();
        ArrayList<String> ret = new ArrayList<String>();

        ret.add("" + i2s(l.getLinkID()));
        ret.add(RArrayUtils.toString(pepids, ";", numberFormat));
        ret.add(RArrayUtils.toString(psmids, ";"));
        ret.add(l.site1Accessions());
        ret.add(l.site1Descriptions());
        ret.add("" + pg1.isDecoy());
        ret.add(l.site2Accessions());
        ret.add(l.site2Descriptions());
        ret.add("" + pg2.isDecoy());
        ret.add(RArrayUtils.toString(l.site1Sites(), ";", numberFormat));
        ret.add(RArrayUtils.toString(l.site2Sites(), ";", numberFormat));
        ret.add(RArrayUtils.toString(xl, ";"));
        ret.add(d2s(l.getScore()));
        ret.add("" + l.isDecoy());
        ret.add("" + l.isTT());
        ret.add("" + l.isTD());
        ret.add("" + l.isDD());
        ret.add(i2s(psmids.length));
        ret.add(i2s(pepids.length));
        ret.add(l.getFDRGroup());
        ret.add(d2s(l.getFDR()));
        ret.add(d2s(l.getEstimatedFDR()));
        ret.add(d2s(l.getPEP()));
        for (String r : runs) {
            Double d = runScore.get(r);
            if (d == null) {
                ret.add("");
            } else {
                ret.add(d2s(d));
            }
        }
        ret.add("");
        ret.add(d2s(l.getProteinGroup1().getFDR()));
        ret.add(d2s(l.getProteinGroup2().getFDR()));
        ret.add(d2s(ppi.getFDR()));
        ret.add("");
        ret.add(d2s(top_pepfdr));
        ret.add(d2s(top_psmfdr));
        ret.add("");
        ret.add(d2s(ppi.getProteinGroupPairID()));
        return ret;
    }

    private void psmToRun(PSM psm, HashSet<String> linkRuns, HashMap<String, Double> runScore) {
        String r = psm.getRun();
        linkRuns.add(r);
        Double d = runScore.get(r);
        if (d == null || d < psm.getScore()) {
            runScore.put(r, psm.getScore());
        }
    }

    protected String getLinkOutputHeader(String seperator) {
        return RArrayUtils.toString(getLinkOutputHeader(), seperator);
    }

    protected ArrayList<String> getLinkOutputHeader() {
        ArrayList<String> ret = new ArrayList<String>();
        ret.add("LinkID");
        ret.add("PeptidePairIDs");
        ret.add("PSMIDs");
        ret.add("Protein1");
        ret.add("Description1");
        ret.add("Decoy1");
        ret.add("Protein2");
        ret.add("Description2");
        ret.add("Decoy2");
        ret.add("fromSite");
        ret.add("ToSite");
        ret.add("Croslinkers");
        ret.add("Score");
        ret.add("isDecoy");
        ret.add("isTT");
        ret.add("isTD");
        ret.add("isDD");
        ret.add("count PSMs");
        ret.add("count peptide pairs");
        ret.add("fdrGroup");
        ret.add("fdr");
        ret.add("ifdr");
        ret.add("PEP");

        for (String r : runs) {
            ret.add(r);
        }

        ret.add("");
        ret.add("Protein1FDR");
        ret.add("Protein2FDR");
        ret.add("PPIFDR");
        ret.add("");
        ret.add("top pep fdr");
        ret.add("top psm fdr");
        ret.add("");
        ret.add("PPI id");
        return ret;
    }

    protected String getPPIOutputHeader(String seperator) {
        return RArrayUtils.toString(getPPIOutputHeader(), seperator);
    }

    protected ArrayList<String> getPPIOutputHeader() {
        ArrayList<String> ret = new ArrayList<String>();
        ret.add("ProteinGroupPairID");
        ret.add("LinkIDs");
        ret.add("PeptidePairIDs");
        ret.add("PSMIDs");
        ret.add("Protein1");
        ret.add("Descriptions1");
        ret.add("isDecoy1");
        ret.add("Protein2");
        ret.add("Description2");
        ret.add("isDecoy2");
        ret.add("Crosslinker");
        ret.add("Score");
        ret.add("isDecoy");
        ret.add("isTT");
        ret.add("isTD");
        ret.add("isDD");
        ret.add("count PSMs");
        ret.add("count peptide pairs");
        ret.add("count links");
        ret.add("fdrGroup");
        ret.add("fdr");
        ret.add("ifdr");
        ret.add("PEP");
        ret.add("");
        ret.add("top link fdr");
        ret.add("top peptide fdr");
        ret.add("top psm_fdr");
        return ret;
    }

    protected String getPPIOutputLine(ProteinGroupPair pgp, String seperator) {
        return RArrayUtils.toString(getPPIOutputLine(pgp), seperator);
    }

    protected ArrayList<String> getPPIOutputLine(ProteinGroupPair pgp) {
        int[] pepids = pgp.getPeptidePairIDs();
        String[] psmids = pgp.getPSMIDs();
        int[] linkids = pgp.getLinkIDs();
        double top_linkfdr = Double.MAX_VALUE;
        double top_pepfdr = Double.MAX_VALUE;
        double top_psmfdr = Double.MAX_VALUE;
        HashSet<String> xl = new HashSet<>();
        for (ProteinGroupLink l : pgp.getLinks()) {
            if (l.getFDR() < top_linkfdr) {
                top_linkfdr = l.getFDR();
            }
            for (PeptidePair pp : l.getPeptidePairs()) {
                if (pp.getFDR() < top_pepfdr) {
                    top_pepfdr = pp.getFDR();
                }
                for (PSM psm : pp.getTopPSMs()) {
                    if (psm.getFDR() < top_psmfdr) {
                        top_psmfdr = psm.getFDR();
                    }
                }
                xl.add(pp.getCrosslinker());
            }
        }
        ArrayList<String> ret = new ArrayList<String>();

        ret.add("" + pgp.getProteinGroupPairID());
        ret.add(RArrayUtils.toString(linkids, ";", numberFormat));
        ret.add(RArrayUtils.toString(pepids, ";", numberFormat));
        ret.add(RArrayUtils.toString(psmids, ";"));
        ret.add(pgp.getProtein1().accessions());
        ret.add(pgp.getProtein1().descriptions());
        ret.add("" + pgp.getProtein1().isDecoy());
        ret.add(pgp.getProtein2().accessions());
        ret.add(pgp.getProtein2().descriptions());
        ret.add("" + pgp.getProtein2().isDecoy());
        ret.add(RArrayUtils.toString(xl, ";"));
        ret.add(d2s(pgp.getScore()));
        ret.add("" + pgp.isDecoy());
        ret.add("" + pgp.isTT());
        ret.add("" + pgp.isTD());
        ret.add("" + pgp.isDD());
        ret.add(i2s(psmids.length));
        ret.add(i2s(pepids.length));
        ret.add(i2s(linkids.length));
        ret.add(pgp.getFDRGroup());
        ret.add(d2s(pgp.getFDR()));
        ret.add(d2s(pgp.getEstimatedFDR()));
        ret.add(d2s(pgp.getPEP()));
        ret.add("");
        ret.add(d2s(top_linkfdr));
        ret.add(d2s(top_pepfdr));
        ret.add(d2s(top_psmfdr));
        return ret;
    }

    protected ArrayList<String> getProteinGroupOutputHeader() {
        ArrayList<String> ret = new ArrayList<String>();
        ret.add("ProteinGroupID");
        ret.add("ProteinGroup");
        ret.add("Descriptions");
        ret.add("Sequence");
        ret.add("Crosslinker");
        ret.add("Score");
        ret.add("isDecoy");
        ret.add("isTT");
        ret.add("isTD");
        ret.add("isDD");
        ret.add("PSM IDs");
        ret.add("fdrGroup");
        ret.add("fdr");
        ret.add("ifdr");
        ret.add("PEP");
        for (String r : runs) {
            ret.add(r);
        }
        return ret;
    }

    protected ArrayList<String> getProteinGroupOutput(ProteinGroup pg) {
        ArrayList<Integer> pepids = new ArrayList<>();
        ArrayList<String> psmids = new ArrayList<>();
        HashSet<String> linkRuns = new HashSet<>();
        HashMap<String, Double> runScore = new HashMap<>();

        double top_pepfdr = Double.MAX_VALUE;
        double top_psmfdr = Double.MAX_VALUE;

        HashSet<String> crosslinker = new HashSet<>();
        for (PeptidePair pp : pg.getPeptidePairs()) {
            pepids.add(pp.getPeptidePairID());

            crosslinker.add(pp.getCrosslinker());
            if (pp.getFDR() < top_pepfdr) {
                top_pepfdr = pp.getFDR();
            }
            for (PSM psm : pp.getTopPSMs()) {
                if (psm.getFDR() < top_psmfdr) {
                    top_psmfdr = psm.getFDR();
                }
            }
            for (PSM psm : pp.getAllPSMs()) {
                for (PSM upsm : psm.getRepresented()) {
                    psmToRun(upsm, linkRuns, runScore);
                    psmids.add(upsm.getPsmID());
                }
            }
        }
        ArrayList<String> sequences = new ArrayList<>();
        for (Protein p : pg.getProteins()) {
            sequences.add(p.getSequence());
        }
        ArrayList<String> ret = new ArrayList<String>();
        ret.add(pg.ids());
        ret.add(pg.accessions());
        ret.add(pg.descriptions());
        ret.add(RArrayUtils.toString(sequences, ";"));
        ret.add(RArrayUtils.toString(crosslinker, ";"));
        ret.add(d2s(pg.getScore()));
        ret.add(pg.isDecoy() + "");
        ret.add(pg.isTT() + "");
        ret.add(pg.isTD() + "");
        ret.add(pg.isDD() + "");
        ret.add(RArrayUtils.toString(psmids, ";"));
        ret.add(pg.getFDRGroup());
        ret.add(pg.getFDR() + "");
        ret.add(pg.getEstimatedFDR() + "");
        ret.add(pg.getPEP() + "");
        for (String r : runs) {
            Double d = runScore.get(r);
            if (d == null) {
                ret.add("");
            } else {
                ret.add(d.toString());
            }
        }
        return ret;
    }

    protected String getPeptideSequence(Peptide p) {
        return p.getSequence();
    }

    public PSM addMatch(String psmID, String run, String scan, Long pepid1, Long pepid2, String pepSeq1, String pepSeq2, int peplen1, int peplen2, int site1, int site2, boolean isDecoy1, boolean isDecoy2, int charge, double score, Long protid1, String accession1, String description1, Long protid2, String accession2, String description2, int pepPosition1, int pepPosition2, double peptide1score, double peptide2score, String isSpecialCase, String crosslinker) {
        return addMatch(psmID, run, scan, pepid1, pepid2, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, protid1, accession1, description1, protid2, accession2, description2, pepPosition1, pepPosition2, "", "", peptide1score, peptide2score, isSpecialCase, crosslinker);
    }

    public PSM addMatch(String psmID, String run, String scan, Long pepid1, Long pepid2, String pepSeq1, String pepSeq2, int peplen1, int peplen2, int site1, int site2, boolean isDecoy1, boolean isDecoy2, int charge, double score, Long protid1, String accession1, String description1, Long protid2, String accession2, String description2, int pepPosition1, int pepPosition2, String Protein1Sequence, String Protein2Sequence, double peptide1score, double peptide2score, String isSpecialCase, String crosslinker) {
        PSM ret = addMatch(psmID, pepid1, pepid2, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, protid1, accession1, description1, protid2, accession2, description2, pepPosition1, pepPosition2, Protein1Sequence, Protein2Sequence, peptide1score, peptide2score, isSpecialCase, crosslinker, run, scan);
        return ret;

    }

    public PSM addMatch(String psmID, Long pepid1, Long pepid2, String pepSeq1, String pepSeq2, int peplen1, int peplen2, int site1, int site2, boolean isDecoy1, boolean isDecoy2, int charge, double score, Long protid1, String accession1, String description1, Long protid2, String accession2, String description2, int pepPosition1, int pepPosition2, String Protein1Sequence, String Protein2Sequence, double peptide1score, double peptide2score, String isSpecialCase, String crosslinker, String run, String Scan) {
//    public PSM addMatch(String pepSeq2, String pepSeq1, String accession1, String accession2, int protid1, String description2, boolean isDecoy1, int pepid1, int pepPosition1, int peplen1, int protid2, boolean isDecoy2, int pepid2, int pepPosition2, int peplen2, String psmID, int site1, int site2, int charge, double score, double scoreRatio, boolean isSpecialCase) {
        boolean linear = pepSeq2 == null || pepSeq2.isEmpty() || pepSeq1 == null || pepSeq1.isEmpty();
        boolean internal = (!linear) && (accession1.contentEquals(accession2) || ("REV_" + accession1).contentEquals(accession2) || accession1.contentEquals("REV_" + accession2));
        boolean between = !(linear || internal);
        Protein p1 = new Protein(protid1, accession1, description1, isDecoy1, linear, internal, between);

        if (Protein1Sequence != null && !Protein1Sequence.isEmpty()) {
            p1.setSequence(Protein1Sequence);
        }

        p1 = allProteins.register(p1);

        Peptide pep1 = allPeptides.register(new Peptide(pepid1, pepSeq1, isDecoy1, p1, pepPosition1, peplen1));
        Protein p2;
        Peptide pep2;

        if (linear) {
            p2 = Protein.NOPROTEIN;
            pep2 = Peptide.NOPEPTIDE;
        } else {
            p2 = new Protein(protid2, accession2, description2, isDecoy2, false, internal, between);
            if (Protein2Sequence != null && !Protein2Sequence.isEmpty()) {
                p2.setSequence(Protein2Sequence);
            }
            p2 = allProteins.register(p2);
            pep2 = allPeptides.register(new Peptide(pepid2, pepSeq2, isDecoy2, p2, pepPosition2, peplen2));
        }
        PSM psm = addMatch(psmID, pep1, pep2, peplen1, peplen2, site1, site2, charge, score, p1, p2, pepPosition1, pepPosition2, peptide1score, peptide2score, isSpecialCase, crosslinker, run, Scan);

        return psm;
    }

    /**
     * @return the allPSMs
     */
    public SelfAddHashSet<PSM> getAllPSMs() {
        return allPSMs;
    }

    /**
     * @param allPSMs the allPSMs to set
     */
    public void setAllPSMs(SelfAddHashSet<PSM> allPSMs) {
        this.allPSMs = allPSMs;
    }

    /**
     * indicates, whether the psms went through a score normalisation
     *
     * @return the isNormalized
     */
    public boolean isNormalized() {
        return isNormalized;
    }

    /**
     * indicates, whether the psms went through a score normalisation
     *
     * @param isNormalized the isNormalized to set
     */
    public void setNormalized(boolean isNormalized) {
        this.isNormalized = isNormalized;
    }

    /**
     * indicates, whether the psms went through a score normalisation
     *
     * @param isNormalized the isNormalized to set
     */
    public void setDecoyNormalized() {
        setNormalized(true);
        this.isNormalizedByDecoy = true;
        this.isNormalized = isNormalized;
    }

    /**
     * indicates, whether the psms went through a score normalisation
     *
     * @param isNormalized the isNormalized to set
     */
    public void setOveralNormalized() {
        setNormalized(true);
        this.isNormalizedByDecoy = false;
        this.isNormalized = isNormalized;
    }

    /**
     * We normalize psm-scores by median and MAD - but to go around some quirks
     * of our score propagation scores then get shifted so that the lowest score
     * is around one.
     *
     * @return the psmNormOffset
     */
    public double getPsmNormalizationOffset() {
        return psmNormOffset;
    }

    /**
     * We normalize psm-scores by median and MAD - but to go around some quirks
     * of our score propagation scores then get shifted so that the lowest score
     * is around one.
     *
     * @param psmNormOffset the psmNormOffset to set
     */
    public void setPsmNormalizationOffset(double psmNormOffset) {
        this.psmNormOffset = psmNormOffset;
    }

    /**
     * the version of xiFDR to be reported
     *
     * @return the xiFDRVersion
     */
    public static Version getXiFDRVersion() {
        if (xiFDRVersion == null) {
            xiFDRVersion = Version.parseEmbededVersion("xifdrproject.properties", "xifdr.version");
        }
        return xiFDRVersion;
    }

    /**
     * the version of xiFDR to be reported
     *
     * @param aXiFDRVersion the xiFDRVersion to set
     */
    public static void setXiFDRVersion(Version aXiFDRVersion) {
        xiFDRVersion = aXiFDRVersion;
    }

    /**
     * group matches by protein pairs
     *
     * @return the groupByProteinPair
     */
    public boolean isGroupByProteinPair() {
        return groupByProteinPair;
    }

    /**
     * group matches by protein pairs
     *
     * @param groupByProteinPair the groupByProteinPair to set
     */
    public void setGroupByProteinPair(boolean groupByProteinPair) {
        this.groupByProteinPair = groupByProteinPair;
    }

    /**
     * how many decoys does a fdr group need to have to be reported as result
     *
     * @return the minDecoys
     */
    public Integer getMinDecoys() {
        return minTDChance;
    }

    /**
     * how many decoys does a fdr group need to have to be reported as result
     *
     * @param minDecoys the minDecoys to set
     */
    public void setMinDecoys(Integer minDecoys) {
        this.minTDChance = minDecoys;
    }

    public MaximisingStatus maximise(FDRSettings fdrSettings, OfflineFDR.FDRLevel level, final boolean between, final MaximizingUpdate stateUpdate) {

        final FDRSettingsImpl settings = new FDRSettingsImpl();
        settings.setAll(fdrSettings);

        try {
            int steps = settings.getBoostingSteps();
            StringBuffer sb = new StringBuffer();
            double maxDelta = 1;
            double maxPeptideCoverage = 1;

            // get some settings, that are constant for all calculations
            boolean ignoreGroups = this.ignoreGroupsSetting;

            // if minpepcoverage is 0 and we boost on it we can probably savly start ad 0.1
            if (settings.getMinPeptideCoverageFilter() == 0 && settings.boostPepCoverage()) {
//                settings.setMinPeptideCoverageFilter(0.01);
                maxPeptideCoverage = 0.8;
            }

            // also for delta score we can probably savly start as 0.05
            if (settings.boostDeltaScore()) {
//                if (settings.getMinDeltaScoreFilter() == 0) {
//                    settings.setMinDeltaScoreFilter(0.01);
//                }
                // there is probably no point in being more restrictive then 0.7 (even that is overkill)
                if (settings.getMinDeltaScoreFilter() < 0.6) {
                    maxDelta = 0.7;
                }
            }

            final MaximizeLevelInfo absPepCover = new MaximizeLevelInfoInteger(10, true, steps);
            final MaximizeLevelInfo deltaScore = new MaximizeLevelInfo(1 - settings.getMinDeltaScoreFilter(), 1 - maxDelta, settings.boostDeltaScore(), steps);
            final MaximizeLevelInfo pepCoverage = new MaximizeLevelInfo(1 - settings.getMinPeptideCoverageFilter(), 1 - maxPeptideCoverage, settings.boostPepCoverage(), steps);
            final MaximizeLevelInfo psmFDRInfo = new MaximizeLevelInfo(settings.getPSMFDR(), settings.boostPSMs(), steps);
            final MaximizeLevelInfo pepFDRInfo = new MaximizeLevelInfo(settings.getPeptidePairFDR(), settings.boostPeptidePairs() && level.compareTo(OfflineFDR.FDRLevel.PEPTIDE_PAIR) > 0, steps);
            final MaximizeLevelInfo protFDRInfo = new MaximizeLevelInfo(settings.getProteinGroupFDR(), settings.boostProteins() && level.compareTo(OfflineFDR.FDRLevel.PROTEINGROUP) > 0, steps);
            final MaximizeLevelInfo linkFDRInfo = new MaximizeLevelInfo(settings.getProteinGroupLinkFDR(), settings.boostLinks() && level.compareTo(OfflineFDR.FDRLevel.PROTEINGROUPLINK) > 0, steps);
            final MaximizeLevelInfo ppiFDRInfo = new MaximizeLevelInfo(settings.getProteinGroupPairFDR(), false, steps);
            final ArrayList<MaximizeLevelInfo> subScoreFilter = new ArrayList<>();
            MaximizeLevelInfo targetInfoFirst;

            switch (level) {
                case PEPTIDE_PAIR:
                    targetInfoFirst = pepFDRInfo;
                    break;
                case PROTEINGROUP:
                    targetInfoFirst = protFDRInfo;
                    break;
                case PROTEINGROUPLINK:
                    targetInfoFirst = linkFDRInfo;
                    break;
                case PROTEINGROUPPAIR:
                    targetInfoFirst = ppiFDRInfo;
                    break;
                default:
                    targetInfoFirst = null;
                    break;
            }
            final MaximizeLevelInfo targetInfo = targetInfoFirst;

            boolean filterToUniquePSM = settings.filterToUniquePSM();

            double linkfdr;

            boolean optimizing = true;

            int maxCount = 0;
            int maxCountBetween = 0;

            int countDown = 5;

            int optimizingRound = 1;

            while (optimizing) {
                int lastMaxCount = maxCount;
                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Round " + optimizingRound++);

                Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                        (deltaScore.boost ? "deltaFilter from  :        " + (1 - deltaScore.fromFDR) + " to " + (1 - deltaScore.toFDR) + "\n" : "")
                        + (pepCoverage.boost ? "PepCoverage from  :        " + (1 - pepCoverage.fromFDR) + " to " + (1 - pepCoverage.toFDR) + "\n" : "")
                        + (absPepCover.boost ? "Min Pep Frags from  :        " + (10 - absPepCover.fromFDR) + " to " + (10 - absPepCover.toFDR) + "\n" : "")
                        + (psmFDRInfo.boost ? "PSM fdr from  :        " + psmFDRInfo.fromFDR + " to " + psmFDRInfo.toFDR + "\n" : "")
                        + (pepFDRInfo.boost ? "Peptide pair fdr from  " + pepFDRInfo.fromFDR + " to " + pepFDRInfo.toFDR + "\n" : "")
                        + (protFDRInfo.boost ? "Protein-groupfdr from  " + protFDRInfo.fromFDR + " to " + protFDRInfo.toFDR + "\n" : "")
                        + (linkFDRInfo.boost ? "linkfdr from  " + linkFDRInfo.fromFDR + " to " + linkFDRInfo.toFDR + "\n" : "")
                        + "Steps : " + psmFDRInfo.steps);
                FDRResult result = new FDRResult();

                // initialise subscore filter
                if (subScoreFilter.size() > 0) {
                    for (int c = 0; c < subScoreFilter.size(); c++) {
                        subScoreFilter.get(c).firstStep();
                    }
                }
                boolean optimzeSubScores = true;

                // find the combinations with the maximum number of ppis
                deltaScoreloop:
                for (deltaScore.firstStep(); deltaScore.doThisStep(); deltaScore.nextStep()) {
                    settings.setMinDeltaScoreFilter(1 - deltaScore.getCurrentFDR());
                    pepFragsloop:
                    for (absPepCover.firstStep(); absPepCover.doThisStep(); absPepCover.nextStep()) {
                        settings.setMinPeptideFragmentsFilter(10 - (int) Math.round(absPepCover.getCurrentFDR()));
                        pepCoverageloop:
                        for (pepCoverage.firstStep(); pepCoverage.doThisStep(); pepCoverage.nextStep()) {
                            settings.setMinPeptideCoverageFilter(1 - pepCoverage.getCurrentFDR());

//                        double highestPSMFDR = Double.POSITIVE_INFINITY;
//                        boolean firstPSMFDR = true;
                            psmloop:
                            for (psmFDRInfo.firstStep(); psmFDRInfo.doThisStep(); psmFDRInfo.nextStep()) {
//                            if (highestPSMFDR < psmFDRInfo.currentFDR)
//                                continue;

                                settings.setPSMFDR(psmFDRInfo.getCurrentFDR());
                                this.calculatePSMFDR(true, ignoreGroups, result, settings);
                                psmFDRInfo.setCountsPrefilter(result.psmFDR);
//                            if ( firstPSMFDR) {
//                                firstPSMFDR  =false;
//                                double maxInFDR = 0;
//                                highestPSMFDR = getHighestSubGroupInputFDR(result.psmFDR.getGroups());
//                            }

                                if (stopMaximizing) {
                                    break deltaScoreloop;
                                }
                                // if we don't get PSM - stop looking at later stages
                                if (result.psmFDR.getResultCount() == 0) {
                                    break psmloop;
                                }

//                            double highestPepFDR = Double.POSITIVE_INFINITY;
//                            boolean firstPepPairFDR = true;
                                peploop:
                                for (pepFDRInfo.firstStep(); pepFDRInfo.doThisStep(); pepFDRInfo.nextStep()) {
                                    // we can skip any steps that would not reduce the peptide FDR below what is already in the input
                                    // calculate peptide level fdr
//                                if (highestPepFDR < pepFDRInfo.currentFDR)
//                                    continue;

                                    //double highestProtFDR = Double.POSITIVE_INFINITY;
                                    //boolean firstProtFDR = true;
                                    protloop:
                                    for (protFDRInfo.firstStep(); protFDRInfo.doThisStep(); protFDRInfo.nextStep()) {
                                        // we can skip any steps that would not reduce the peptide FDR below what is already in the input
                                        // calculate peptide level fdr
//                                    if (highestProtFDR < protFDRInfo.currentFDR)
//                                        continue;

                                        settings.setPeptidePairFDR(pepFDRInfo.getCurrentFDR());
                                        this.calculatePeptidePairFDR(true, result, settings, ignoreGroups);

//                                    if ( firstPepPairFDR) {
//                                        firstPepPairFDR  =false;
//                                        double maxInFDR = 0;
//                                        highestPepFDR = getHighestSubGroupInputFDR(result.peptidePairFDR.getGroups());
//                                    }
//                                    
                                        // if we don't get peptide pairs - stop looking at later stages
                                        if (result.peptidePairFDR.getResultCount() == 0) {
                                            break peploop;
                                        }
                                        pepFDRInfo.setCountsPrefilter(result.peptidePairFDR);

                                        if (stopMaximizing) {
                                            break deltaScoreloop;
                                        }

                                        // calculate protein level fdr
                                        settings.setProteinGroupFDR(protFDRInfo.getCurrentFDR());
                                        this.calculateProteinGroupFDR(ignoreGroups, true, settings, result);

//                                    if ( firstProtFDR) {
//                                        firstProtFDR  =false;
//                                        double maxInFDR = 0;
//                                        highestProtFDR = getHighestSubGroupInputFDR(result.proteinGroupFDR.getGroups());
//                                    }
                                        if (result.proteinGroupFDR.getResultCount() == 0) {
                                            break protloop;
                                        }
                                        protFDRInfo.setCountsPrefilter(result.proteinGroupFDR);

                                        // cut down the peptides by proteins                           
                                        this.filterFDRPeptidePairsByFDRProteinGroups(result);

                                        if (stopMaximizing) {
                                            break deltaScoreloop;
                                        }

//                                    double highestLinkFDR = Double.POSITIVE_INFINITY;
//                                    boolean firstLinkFDR = true;
                                        linkloop:
                                        for (linkFDRInfo.firstStep(); linkFDRInfo.doThisStep(); linkFDRInfo.nextStep()) {
                                            // we can skip any steps that would not reduce the peptide FDR below what is already in the input
                                            // calculate peptide level fdr
//                                        if (highestLinkFDR < linkFDRInfo.currentFDR)
//                                            continue;

                                            // calculate links
                                            settings.setProteinGroupLinkFDR(linkFDRInfo.getCurrentFDR());
                                            this.calculateLinkFDR(ignoreGroups, true, settings, result);
                                            linkFDRInfo.setCountsPrefilter(result.proteinGroupLinkFDR);

//                                        if ( firstLinkFDR && linkFDRInfo.boost) {
//                                            Collection<SubGroupFdrInfo<ProteinGroupLink>> c = result.proteinGroupLinkFDR.getGroups();
//                                            firstLinkFDR  =false;
//                                            // detect the highest fdr that we can see - filtering for anything higher is pointless
//                                            highestLinkFDR = getHighestSubGroupInputFDR(c);
//                                        }
                                            if (result.proteinGroupLinkFDR.getResultCount() == 0) {
                                                break linkloop;
                                            }

                                            if (stopMaximizing) {
                                                break deltaScoreloop;
                                            }

                                            settings.setProteinGroupPairFDR(ppiFDRInfo.getCurrentFDR());
                                            this.calculateProteinGroupPairFDR(ignoreGroups, true, settings, result);

                                            if (result.proteinGroupPairFDR.getResultCount() == 0) {
                                                break linkloop;
                                            }
                                            // now we need to filter down to the required level
                                            //                                if (level.compareTo(level.PROTEINGROUPPAIR)!=0) {
                                            this.filterFDRLinksByFDRProteinGroupPairs(result);
                                            //                                }
                                            //                                if (level.compareTo(level.PROTEINGROUPLINK)!=0){
                                            this.filterFDRPeptidePairsByFDRProteinGroupLinks(result);
                                            //                                }

                                            //                                if (level.compareTo(level.PROTEINGROUP)==0){
                                            this.filterFDRProteinGroupsByFDRPeptidePairs(result);
                                            //                                }

                                            // how many links do we now have?
                                            pepCoverage.setCounts(result.psmFDR);
                                            absPepCover.setCounts(result.psmFDR);
                                            deltaScore.setCounts(result.psmFDR);
                                            psmFDRInfo.setCounts(result.psmFDR);
                                            pepFDRInfo.setCounts(result.peptidePairFDR);
                                            protFDRInfo.setCounts(result.proteinGroupFDR);
                                            linkFDRInfo.setCounts(result.proteinGroupLinkFDR);
                                            ppiFDRInfo.setCounts(result.proteinGroupPairFDR);

                                            int count = targetInfo.count;
                                            int countBetween = targetInfo.countBetween;
                                            if (count == 0 && maxCount <= targetInfo.countPreFilter) {
                                                count = targetInfo.countPreFilter;
                                            }
                                            if (countBetween == 0 && maxCountBetween <= targetInfo.countBetweenPreFilter) {
                                                countBetween = targetInfo.countBetweenPreFilter;
                                            }

                                            // is it a new best?
                                            if (((between && (countBetween > maxCountBetween || (countBetween == maxCountBetween && count > maxCount)))
                                                    || (!between && count > maxCount))) {
                                                maxCount = count;
                                                maxCountBetween = countBetween;

                                                deltaScore.setNewMaxFDR();
                                                absPepCover.setNewMaxFDR();
                                                pepCoverage.setNewMaxFDR();
                                                psmFDRInfo.setNewMaxFDR();
                                                pepFDRInfo.setNewMaxFDR();
                                                protFDRInfo.setNewMaxFDR();
                                                linkFDRInfo.setNewMaxFDR();
                                                protFDRInfo.setNewMaxFDR();

                                                // record that we found a new top
                                                String message = "link count, " + linkFDRInfo.count + "(" + linkFDRInfo.countBetween + " between), Protein Pairs, " + ppiFDRInfo.count + "(" + ppiFDRInfo.countBetween + " between)";
                                                if (linkFDRInfo.boost) {
                                                    message = "link fdr, " + (linkFDRInfo.getCurrentFDR()) + ", " + message;
                                                }
                                                if (protFDRInfo.boost) {
                                                    message = "prot fdr, " + (protFDRInfo.getCurrentFDR()) + ", " + message;
                                                }
                                                if (pepFDRInfo.boost) {
                                                    message = "pep fdr, " + (pepFDRInfo.getCurrentFDR()) + ", " + message;
                                                }
                                                if (psmFDRInfo.boost) {
                                                    message = "psm fdr, " + (psmFDRInfo.getCurrentFDR()) + ", " + message;
                                                }
                                                if (pepCoverage.boost) {
                                                    message = "PepCoverage, " + (1 - pepCoverage.getCurrentFDR()) + ", " + message;
                                                }
                                                if (absPepCover.boost) {
                                                    message = "min Fragmenst, " + (10 - absPepCover.getCurrentFDR()) + ", " + message;
                                                }
                                                if (deltaScore.boost) {
                                                    message = "deltascore, " + (1 - deltaScore.getCurrentFDR()) + ", " + message;
                                                }

                                                // String message = "psmfdr, " + psmFDRInfo.currentFDR + " , pepfdr, " + pepFDRInfo.currentFDR + " ,protfdr, " + protFDRInfo.currentFDR + ", link count, " + linkFDRInfo.count + "(" + linkFDRInfo.countBetween + " between), Protein Pairs, " + ppiFDRInfo.count + "(" + ppiFDRInfo.countBetween + " between)";
                                                sb.append(message + "\n");
                                                Logger.getLogger(this.getClass().getName()).log(Level.INFO, sb.toString());

                                                // forward the values to the gui
                                                final double showDelta = (1 - deltaScore.getCurrentFDR());
                                                final double showPepCoverage = (1 - pepCoverage.getCurrentFDR());
                                                final int showMinFrags = (int) Math.round((10 - absPepCover.getCurrentFDR()));
                                                final double showPSMFDR = psmFDRInfo.getCurrentFDR();
                                                final double showPepFDR = pepFDRInfo.getCurrentFDR();
                                                final double showProtFDR = protFDRInfo.getCurrentFDR();
                                                final double showLinkFDR = linkFDRInfo.getCurrentFDR();
                                                final double showPSMCount = psmFDRInfo.count;
                                                final double showPepCount = pepFDRInfo.count;
                                                final double showProtCount = protFDRInfo.countLinear;
                                                final double showLinkCount = linkFDRInfo.count;
                                                final double showPPICount = ppiFDRInfo.count;
                                                final double showLinkCountBetween = linkFDRInfo.countBetween;
                                                final double showPPICountBetween = ppiFDRInfo.countBetween;
                                                final double showLinkCountBetweenPreFilter = linkFDRInfo.countBetweenPreFilter;
                                                final double showLinkCountPreFilter = linkFDRInfo.countPreFilter;
                                                final MaximisingStatus state = new MaximisingStatus();
                                                state.showDelta = showDelta;
                                                state.showMinFrags = showMinFrags;
                                                state.showPepCoverage = showPepCoverage;
                                                state.showPSMFDR = showPSMFDR;
                                                state.showPepFDR = showPepFDR;
                                                state.showProtFDR = showProtFDR;
                                                state.showLinkFDR = showLinkFDR;

                                                state.showPSMCount = showPSMCount + "";
                                                state.showPepCount = showPepCount + "";
                                                state.showProtCount = showProtCount + "";
                                                state.showLinkCount = showLinkCount + (showLinkCount == 0 ? "(before PPI:" + showLinkCountPreFilter + ")" : "");
                                                state.showLinkCountBetween = showLinkCountBetween + (showLinkCountBetween == 0 ? "(before PPI:" + showLinkCountBetweenPreFilter + ")" : "");

                                                state.showPPICount = showPPICount + "";
                                                state.showPPICountBetween = showPPICountBetween + "";
                                                state.resultCount = maxCount;
                                                state.resultCountBetween = maxCountBetween;
                                                stateUpdate.setStatus(state);
                                            } else if (count == maxCount || (between && countBetween == maxCountBetween)) {
                                                deltaScore.setEqualFDR();
                                                absPepCover.setEqualFDR();
                                                pepCoverage.setEqualFDR();
                                                psmFDRInfo.setEqualFDR();
                                                pepFDRInfo.setEqualFDR();
                                                protFDRInfo.setEqualFDR();
                                                linkFDRInfo.setEqualFDR();
                                                ppiFDRInfo.setEqualFDR();
                                            }

                                            if (stopMaximizing) {
                                                break deltaScoreloop;
                                            }
                                        }

                                    }

                                }

                            }
                        }
                    }
                }

                stateUpdate.setStatusText("Max Round: " + optimizingRound + " - " + maxCount + " matches");

                // no improvement for the last few rounds?
                if ((maxCount == lastMaxCount && --countDown == 0) || stopMaximizing) {
                    optimizing = false;
                    Logger.getLogger(this.getClass().getName()).log(Level.INFO, sb.toString());
//                    FDRResult ret = this.calculateFDR(psmFDRInfo.maximumFDR, pepFDRInfo.maximumFDR, protFDRInfo.maximumFDR, linkFDRInfo.maximumFDR, ppiFDRInfo.maximumFDR, 10000, ignoreGroups, true, filterToUniquePSM);
                    settings.setMinDeltaScoreFilter(1 - deltaScore.maximumFDR);
                    settings.setMinPeptideCoverageFilter(1 - pepCoverage.maximumFDR);
                    settings.setMinPeptideFragmentsFilter(10 - (int) Math.round(absPepCover.maximumFDR));
                    settings.setPSMFDR(psmFDRInfo.maximumFDR);
                    settings.setPeptidePairFDR(pepFDRInfo.maximumFDR);
                    settings.setProteinGroupFDR(protFDRInfo.maximumFDR);
                    settings.setProteinGroupLinkFDR(linkFDRInfo.maximumFDR);
                    settings.setProteinGroupPairFDR(ppiFDRInfo.maximumFDR);
                    FDRResult ret = this.calculateFDR(settings, true);

                    final int foundCount = maxCount;
                    MaximisingStatus res = new MaximisingStatus();
                    res.showDelta = (1 - deltaScore.maximumFDR);
                    res.showMinFrags = (10 - (int) Math.round(absPepCover.maximumFDR));
                    res.showPepCoverage = (1 - pepCoverage.maximumFDR);
                    res.showPSMFDR = psmFDRInfo.maximumFDR;
                    res.showPepFDR = pepFDRInfo.maximumFDR;
                    res.showProtFDR = protFDRInfo.maximumFDR;
                    res.showLinkFDR = linkFDRInfo.maximumFDR;
                    res.result = ret;
                    res.resultCount = maxCount;
                    res.resultCountBetween = maxCountBetween;

                    if (pepCoverage.maximumFDR < 1) {
                        PSMFilter f = new SingleSubScoreFilter("minPepCoverage", 1 - pepCoverage.maximumFDR, true);
                        prefilteredPSMs = f.filter(getAllPSMs());
                        if (prefilteredPSMs.size() < 50) {
                            break;
                        }
                    } else {
                        prefilteredPSMs = null;
                    }

                    stopMaximizing = false;

                    stateUpdate.setStatus(res);

                    return res;

                } else {
                    if (maxCount > lastMaxCount) {
                        // yes we improved
                        countDown = 5;
                    } else {
                        stateUpdate.setStatusText("Max Link Round: " + optimizingRound + " - " + maxCount + " links - count down: " + countDown);

                    }
                    // so see if we make the resoltuion finer
                    // can we get a better result?
                    boolean startLow = false;
                    int stepChange = 0;
                    if (countDown == 2) {
                        startLow = true;
                        stepChange = 1;
                    }
                    if (countDown == 1) {
                        stepChange = 1;
                    }
                    deltaScore.calcNextFDRRange(startLow, stepChange);
                    absPepCover.calcNextFDRRange(startLow, stepChange);
                    pepCoverage.calcNextFDRRange(startLow, stepChange);
                    psmFDRInfo.calcNextFDRRange(startLow, stepChange);
                    pepFDRInfo.calcNextFDRRange(startLow, stepChange);
                    linkFDRInfo.calcNextFDRRange(startLow, stepChange);
                    protFDRInfo.calcNextFDRRange(startLow, stepChange);
                    ppiFDRInfo.calcNextFDRRange(startLow, stepChange);

                }

            }
        } catch (Exception ex) {
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Error maximizing links", ex);
            stateUpdate.reportError("Error maximizing", ex);
        }
        return null;

    }

    public MaximisingStatus maximiseSequential(FDRSettings fdrSettings, OfflineFDR.FDRLevel level, final boolean between, final MaximizingUpdate stateUpdate) {

        final FDRSettingsImpl settings = new FDRSettingsImpl();
        settings.setAll(fdrSettings);

        try {
            int steps = settings.getBoostingSteps();
            StringBuffer sb = new StringBuffer();
            double maxDelta = 1;
            double maxPeptideCoverage = 1;

            // get some settings, that are constant for all calculations
            boolean ignoreGroups = this.ignoreGroupsSetting;

            // if minpepcoverage is 0 and we boost on it we can probably savly start ad 0.1
            if (settings.getMinPeptideCoverageFilter() == 0 && settings.boostPepCoverage()) {
//                settings.setMinPeptideCoverageFilter(0.01);
                maxPeptideCoverage = 0.8;
            }

            // also for delta score we can probably savly start as 0.05
            if (settings.boostDeltaScore()) {
//                if (settings.getMinDeltaScoreFilter() == 0) {
//                    settings.setMinDeltaScoreFilter(0.01);
//                }
                // there is probably no point in being more restrictive then 0.7 (even that is overkill)
                if (settings.getMinDeltaScoreFilter() < 0.6) {
                    maxDelta = 0.7;
                }
            }

            final MaximizeLevelInfo absPepCover = new MaximizeLevelInfoInteger(10, true, steps);
            final MaximizeLevelInfo deltaScore = new MaximizeLevelInfo(1 - settings.getMinDeltaScoreFilter(), 1 - maxDelta, settings.boostDeltaScore(), steps);
            final MaximizeLevelInfo pepCoverage = new MaximizeLevelInfo(1 - settings.getMinPeptideCoverageFilter(), 1 - maxPeptideCoverage, settings.boostPepCoverage(), steps);
            final MaximizeLevelInfo psmFDRInfo = new MaximizeLevelInfo(settings.getPSMFDR(), settings.boostPSMs(), steps);
            final MaximizeLevelInfo pepFDRInfo = new MaximizeLevelInfo(settings.getPeptidePairFDR(), settings.boostPeptidePairs() && level.compareTo(OfflineFDR.FDRLevel.PEPTIDE_PAIR) > 0, steps);
            final MaximizeLevelInfo protFDRInfo = new MaximizeLevelInfo(settings.getProteinGroupFDR(), settings.boostProteins() && level.compareTo(OfflineFDR.FDRLevel.PROTEINGROUP) > 0, steps);
            final MaximizeLevelInfo linkFDRInfo = new MaximizeLevelInfo(settings.getProteinGroupLinkFDR(), settings.boostLinks() && level.compareTo(OfflineFDR.FDRLevel.PROTEINGROUPLINK) > 0, steps);
            final MaximizeLevelInfo ppiFDRInfo = new MaximizeLevelInfo(settings.getProteinGroupPairFDR(), false, steps);
            absPepCover.name = "min frags";
            deltaScore.name = "delta";
            pepCoverage.name = "pep coverage";
            psmFDRInfo.name = "psm";
            pepFDRInfo.name = "pep";
            protFDRInfo.name = "prot";
            linkFDRInfo.name = "link";
            ppiFDRInfo.name = "ppi";

            final ArrayList<MaximizeLevelInfo> subScoreFilter = new ArrayList<>();
            MaximizeLevelInfo targetInfoFirst;

            ArrayList<MaximizeLevelInfo> boostLevels = new ArrayList<>();
            if (settings.boostMinFragments()) {
                boostLevels.add(absPepCover);
            }
            if (settings.boostPepCoverage()) {
                boostLevels.add(pepCoverage);
            }
            if (settings.boostDeltaScore()) {
                boostLevels.add(deltaScore);
            }
            switch (level) {
                case PROTEINGROUPPAIR:
                    if (settings.boostLinks()) {
                        boostLevels.add(linkFDRInfo);
                    }
                case PROTEINGROUPLINK:
                    if (settings.boostProteins()) {
                        boostLevels.add(protFDRInfo);
                    }
                    if (settings.boostPeptidePairs()) {
                        boostLevels.add(pepFDRInfo);
                    }
                case PEPTIDE_PAIR:
                    if (settings.boostPSMs()) {
                        boostLevels.add(psmFDRInfo);
                    }
            }

            switch (level) {
                case PSM:
                    targetInfoFirst = psmFDRInfo;
                    break;
                case PEPTIDE_PAIR:
                    targetInfoFirst = pepFDRInfo;
                    break;
                case PROTEINGROUP:
                    targetInfoFirst = protFDRInfo;
                    break;
                case PROTEINGROUPLINK:
                    targetInfoFirst = linkFDRInfo;
                    break;
                case PROTEINGROUPPAIR:
                    targetInfoFirst = ppiFDRInfo;
                    break;
                default:
                    targetInfoFirst = null;
                    break;
            }
            final MaximizeLevelInfo targetInfo = targetInfoFirst;

            boolean filterToUniquePSM = settings.filterToUniquePSM();

            MaximisingStatus res = new MaximisingStatus();

            int maxCount = 0;
            int maxCountBetween = 0;


            for (int i = 2; i>0;i--) {
                for (MaximizeLevelInfo ml : boostLevels) {
                    ml.fromFDR = ml.minimumFDR;
                    boolean optimizing = true;

                    int countDown = 5;

                    int optimizingRound = 0;
                    if (ml == absPepCover) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                "Min Pep Frags from  :        " + (10 - ml.fromFDR) + " to " + (10 - ml.toFDR));
                    } else if (ml == pepCoverage) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                "Min Pep Coverage from  :        " + (1 - ml.fromFDR) + " to " + (1 - ml.toFDR));
                    } else if (ml == deltaScore) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                "delta score from  :        " + (1 - ml.fromFDR) + "*score to " + (1 - ml.toFDR) + "*score");
                    } else if (ml == psmFDRInfo) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                "psmFDR from  :        " + ml.fromFDR + " to " + ml.toFDR);
                    } else if (ml == pepFDRInfo) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                "pep FDR from  :        " + ml.fromFDR + " to " + ml.toFDR);
                    } else if (ml == protFDRInfo) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                "prot FDR from  :        " + ml.fromFDR + " to " + ml.toFDR);
                    } else if (ml == linkFDRInfo) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                "link FDR from  :        " + ml.fromFDR + " to " + ml.toFDR);
                    }
                    while (optimizing && !stopMaximizing) {
                        optimizingRound++;
                        int lastMaxCount = maxCount;
                        for (ml.firstStep(); ml.doThisStep(); ml.nextStep()) {
                            FDRResult result = new FDRResult();
                            settings.setMinDeltaScoreFilter(1 - deltaScore.getCurrentFDR());
                            settings.setMinPeptideFragmentsFilter(10 - (int) Math.round(absPepCover.getCurrentFDR()));
                            settings.setMinPeptideCoverageFilter(1 - pepCoverage.getCurrentFDR());
                            settings.setPeptidePairFDR(pepFDRInfo.getCurrentFDR());
                            settings.setProteinGroupFDR(protFDRInfo.getCurrentFDR());
                            settings.setProteinGroupLinkFDR(linkFDRInfo.getCurrentFDR());
                            this.calculatePSMFDR(true, ignoreGroups, result, settings);
                            if (result.psmFDR.getResultCount() == 0 || stopMaximizing) {
                                break;
                            }
                            psmFDRInfo.setCountsPrefilter(result.psmFDR);

                            this.calculatePeptidePairFDR(true, result, settings, ignoreGroups);
                            if (result.peptidePairFDR.getResultCount() == 0 || stopMaximizing) {
                                break;
                            }
                            pepFDRInfo.setCountsPrefilter(result.peptidePairFDR);
                            this.calculateProteinGroupFDR(ignoreGroups, true, settings, result);
                            if (result.proteinGroupFDR.getResultCount() == 0 || stopMaximizing) {
                                break;
                            }
                            protFDRInfo.setCountsPrefilter(result.proteinGroupFDR);
                            // cut down the peptides by proteins                           
                            this.filterFDRPeptidePairsByFDRProteinGroups(result);
                            this.calculateLinkFDR(ignoreGroups, true, settings, result);
                            linkFDRInfo.setCountsPrefilter(result.proteinGroupLinkFDR);
                            if (result.proteinGroupLinkFDR.getResultCount() == 0 || stopMaximizing) {
                                break;
                            }
                            settings.setProteinGroupPairFDR(ppiFDRInfo.getCurrentFDR());
                            this.calculateProteinGroupPairFDR(ignoreGroups, true, settings, result);

                            if (result.proteinGroupPairFDR.getResultCount() == 0 || stopMaximizing) {
                                break;
                            }
                            // now we need to filter down 
                            this.filterFDRLinksByFDRProteinGroupPairs(result);
                            if (stopMaximizing) {
                                break;
                            }
                            this.filterFDRPeptidePairsByFDRProteinGroupLinks(result);
                            if (stopMaximizing) {
                                break;
                            }
                            this.filterFDRProteinGroupsByFDRPeptidePairs(result);
                            if (stopMaximizing) {
                                break;
                            }
                            this.filterFDRPSMByFDRPeptidePairs(result);

                            // how many links do we now have?
                            pepCoverage.setCounts(result.psmFDR);
                            absPepCover.setCounts(result.psmFDR);
                            deltaScore.setCounts(result.psmFDR);
                            psmFDRInfo.setCounts(result.psmFDR);
                            pepFDRInfo.setCounts(result.peptidePairFDR);
                            protFDRInfo.setCounts(result.proteinGroupFDR);
                            linkFDRInfo.setCounts(result.proteinGroupLinkFDR);
                            ppiFDRInfo.setCounts(result.proteinGroupPairFDR);

                            int count = targetInfo.count;
                            int countBetween = targetInfo.countBetween;
                            if (count == 0 && maxCount <= targetInfo.countPreFilter) {
                                count = targetInfo.countPreFilter;
                            }
                            if (countBetween == 0 && maxCountBetween <= targetInfo.countBetweenPreFilter) {
                                countBetween = targetInfo.countBetweenPreFilter;
                            }

                            // is it a new best?
                            if (((between && (countBetween > maxCountBetween || (countBetween == maxCountBetween && count > maxCount)))
                                    || (!between && count > maxCount))) {
                                maxCount = count;
                                maxCountBetween = countBetween;
                                ml.setNewMaxFDR();
                                String message = " , counts: " + count + ", of wich between : " + countBetween;
                                if (ml == absPepCover) {
                                    message = "Min Pep Frags : " + (10 - ml.getCurrentFDR()) + " , counts: " + count + ", of wich between : " + countBetween;
                                } else if (ml == pepCoverage) {
                                    message = "Min Pep Coverage : " + (1 - ml.getCurrentFDR()) + " , counts: " + count + ", of wich between : " + countBetween;
                                } else if (ml == deltaScore) {
                                    message = "min delta score : " + (1 - ml.getCurrentFDR()) + "*score , counts: " + count + ", of wich between : " + countBetween;
                                } else if (ml == psmFDRInfo) {
                                    message = "psmFDR : " + ml.getCurrentFDR() + "*score , counts: " + count + ", of wich between : " + countBetween;
                                } else if (ml == pepFDRInfo) {
                                    message = "peptide pair FDR : " + ml.getCurrentFDR() + "*score , counts: " + count + ", of wich between : " + countBetween;
                                    Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                            "pep FDR from  :        " + ml.fromFDR + " to " + ml.toFDR);
                                } else if (ml == protFDRInfo) {
                                    message = "protein FDR : " + ml.getCurrentFDR() + "*score , counts: " + count + ", of wich between : " + countBetween;
                                } else if (ml == linkFDRInfo) {
                                    message = "link FDR : " + ml.getCurrentFDR() + "*score , counts: " + count + ", of wich between : " + countBetween;
                                }
                                Logger.getLogger(this.getClass().getName()).log(Level.INFO, sb.toString());

                                // forward the values to the gui
                                final double showDelta = (1 - deltaScore.getCurrentFDR());
                                final double showPepCoverage = (1 - pepCoverage.getCurrentFDR());
                                final int showMinFrags = (int) Math.round((10 - absPepCover.getCurrentFDR()));
                                final double showPSMFDR = psmFDRInfo.getCurrentFDR();
                                final double showPepFDR = pepFDRInfo.getCurrentFDR();
                                final double showProtFDR = protFDRInfo.getCurrentFDR();
                                final double showLinkFDR = linkFDRInfo.getCurrentFDR();
                                final double showPSMCount = psmFDRInfo.count;
                                final double showPepCount = pepFDRInfo.count;
                                final double showProtCount = protFDRInfo.countLinear;
                                final double showLinkCount = linkFDRInfo.count;
                                final double showPPICount = ppiFDRInfo.count;
                                final double showLinkCountBetween = linkFDRInfo.countBetween;
                                final double showPPICountBetween = ppiFDRInfo.countBetween;
                                final double showLinkCountBetweenPreFilter = linkFDRInfo.countBetweenPreFilter;
                                final double showLinkCountPreFilter = linkFDRInfo.countPreFilter;
                                final MaximisingStatus state = new MaximisingStatus();
                                state.showDelta = showDelta;
                                state.showMinFrags = showMinFrags;
                                state.showPepCoverage = showPepCoverage;
                                state.showPSMFDR = showPSMFDR;
                                state.showPepFDR = showPepFDR;
                                state.showProtFDR = showProtFDR;
                                state.showLinkFDR = showLinkFDR;

                                state.showPSMCount = showPSMCount + "";
                                state.showPepCount = showPepCount + "";
                                state.showProtCount = showProtCount + "";
                                state.showLinkCount = showLinkCount + (showLinkCount == 0 ? "(before PPI:" + showLinkCountPreFilter + ")" : "");
                                state.showLinkCountBetween = showLinkCountBetween + (showLinkCountBetween == 0 ? "(before PPI:" + showLinkCountBetweenPreFilter + ")" : "");

                                state.showPPICount = showPPICount + "";
                                state.showPPICountBetween = showPPICountBetween + "";
                                state.resultCount = maxCount;
                                state.resultCountBetween = maxCountBetween;
                                stateUpdate.setStatus(state);
                            } else if (count == maxCount || (between && countBetween == maxCountBetween)) {
                                ml.setEqualFDR();
                            }

                        }

                        stateUpdate.setStatusText("Max Round: " + optimizingRound + " ("+ ml.name+") - " + maxCount + " matches");

                        // no improvement for the last few rounds?
                        if ((maxCount == lastMaxCount && --countDown == 0) || stopMaximizing) {
                            optimizing = false;
                            Logger.getLogger(this.getClass().getName()).log(Level.INFO, sb.toString());
                            //                    FDRResult ret = this.calculateFDR(psmFDRInfo.maximumFDR, pepFDRInfo.maximumFDR, protFDRInfo.maximumFDR, linkFDRInfo.maximumFDR, ppiFDRInfo.maximumFDR, 10000, ignoreGroups, true, filterToUniquePSM);
                            settings.setMinDeltaScoreFilter(1 - deltaScore.maximumFDR);
                            settings.setMinPeptideCoverageFilter(1 - pepCoverage.maximumFDR);
                            settings.setMinPeptideFragmentsFilter(10 - (int) Math.round(absPepCover.maximumFDR));
                            settings.setPSMFDR(psmFDRInfo.maximumFDR);
                            settings.setPeptidePairFDR(pepFDRInfo.maximumFDR);
                            settings.setProteinGroupFDR(protFDRInfo.maximumFDR);
                            settings.setProteinGroupLinkFDR(linkFDRInfo.maximumFDR);
                            //settings.setProteinGroupPairFDR(settings.getProteinGroupPairFDR());
                            res.showDelta = (1 - deltaScore.maximumFDR);
                            res.showMinFrags = (10 - (int) Math.round(absPepCover.maximumFDR));
                            res.showPepCoverage = (1 - pepCoverage.maximumFDR);
                            res.showPSMFDR = psmFDRInfo.maximumFDR;
                            res.showPepFDR = pepFDRInfo.maximumFDR;
                            res.showProtFDR = protFDRInfo.maximumFDR;
                            res.showLinkFDR = linkFDRInfo.maximumFDR;
                            res.resultCount = maxCount;
                            res.resultCountBetween = maxCountBetween;

                            if (pepCoverage.maximumFDR < 1) {
                                PSMFilter f = new SingleSubScoreFilter("minPepCoverage", 1 - pepCoverage.maximumFDR, true);
                                prefilteredPSMs = f.filter(getAllPSMs());
                                if (prefilteredPSMs.size() < 50) {
                                    break;
                                }
                            } else {
                                prefilteredPSMs = null;
                            }

                            stateUpdate.setStatus(res);
                            optimizing = false;
                        } else {
                            if (maxCount > lastMaxCount) {
                                // yes we improved
                                countDown = 5;
                            } else {
                                stateUpdate.setStatusText("Max Link Round: " + optimizingRound + " - " + maxCount + " links - count down: " + countDown);

                            }
                            // so see if we make the resoltuion finer
                            // can we get a better result?
                            boolean startLow = false;
                            int stepChange = 0;
                            if (countDown == 2) {
                                startLow = true;
                                stepChange = 1;
                            }
                            if (countDown == 1) {
                                stepChange = 1;
                            }
                            ml.calcNextFDRRange(startLow, stepChange);

                        }

                    }
                    ml.setCurrentFDR(ml.maximumFDR);

                }
            }
            FDRResult ret = this.calculateFDR(settings, true);
            res.result = ret;

            final int foundCount = maxCount;
            
            stopMaximizing = false;
            return res;
        } catch (Exception ex) {
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Error maximizing links", ex);
            stateUpdate.reportError("Error maximizing", ex);
        }
        return null;

    }

    protected <T extends FDRSelfAdd> double getHighestSubGroupInputFDR(Collection<SubGroupFdrInfo<T>> c) {
        double maxInFDR = 0;
        for (SubGroupFdrInfo g : c) {
            double ginfdr = (g.TD - g.DD) / (double) g.TT;
            if (maxInFDR < ginfdr) {
                maxInFDR = ginfdr;
            }
        }
        return maxInFDR;
    }

    /**
     * @return the numberFormat
     */
    public NumberFormat getNumberFormat() {
        return numberFormat;
    }

    /**
     * @param numberFormat the numberFormat to set
     */
    public void setNumberFormat(NumberFormat numberFormat) {
        this.numberFormat = numberFormat;
    }

    public void stopMaximizing() {
        stopMaximizing = true;
    }

    public void cleanup() {

    }

    /**
     * group between links by both proteins beeing observed with self-links
     *
     * @return the setGroupLinksByHasInternal
     */
    public boolean groupLinksByHasInternal() {
        return groupLinksByHasInternal;
    }

    /**
     * group between links by both proteins beeing observed with self-links
     *
     * @param groupLinksByHasInternal the setGroupLinksByHasInternal to set
     */
    public void setGroupLinksByHasInternal(boolean groupLinksByHasInternal) {
        this.groupLinksByHasInternal = groupLinksByHasInternal;
    }

    /**
     * group between PPIs by both proteins beeing observed with self-links
     *
     * @return the groupPPIByHasInternal
     */
    public boolean groupPPIByHasInternal() {
        return groupPPIByHasInternal;
    }

    /**
     * group between PPIs by both proteins beeing observed with self-links
     *
     * @param groupPPIByHasInternal the groupPPIByHasInternal to set
     */
    public void setGroupPPIByHasInternal(boolean groupPPIByHasInternal) {
        this.groupPPIByHasInternal = groupPPIByHasInternal;
    }

    /**
     * group between PeptidePairs by both proteins beeing observed with
     * self-links
     *
     * @return the groupPepPairByHasInternal
     */
    public boolean groupPepPairByHasInternal() {
        return groupPepPairByHasInternal;
    }

    /**
     * group between PeptidePairs by both proteins beeing observed with
     * self-links
     *
     * @param groupPepPairByHasInternal the groupPepPairByHasInternal to set
     */
    public void setGroupPepPairByHasInternal(boolean groupPepPairByHasInternal) {
        this.groupPepPairByHasInternal = groupPepPairByHasInternal;
    }

    /**
     * group between PeptidePairs by both proteins beeing observed with
     * self-links
     *
     * @param groupByrun
     */
    public void setGroupPSMsByRun(boolean groupByrun) {
        this.groupPSMsByRun = groupByrun;
    }

    /**
     * group between PeptidePairs by both proteins beeing observed with
     * self-links
     *
     * @param groupByrun
     */
    public boolean groupPSMsByRun() {
        return this.groupPSMsByRun;
    }

}
