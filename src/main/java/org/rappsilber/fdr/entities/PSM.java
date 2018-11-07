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
package org.rappsilber.fdr.entities;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import org.rappsilber.fdr.groups.ProteinGroup;
import org.rappsilber.fdr.utils.AbstractFDRElement;
import org.rappsilber.fdr.utils.FDRGroupNames;
import org.rappsilber.utils.SelfAddHashSet;


/**
 * represents a single PSM
 */
public class PSM extends AbstractFDRElement<PSM> { 

    /**
     * unique id for the PSM.
     */
    private String psmID;
    /**
     * Run name
     */
    private String run;
    /**
     * scan-number
     */
    private String scan;
    /**
     * the peptide pair represented by this PSM
     */
    private PeptidePair m_peppair;
    

    /** link-site in peptide1 */
    private byte pepsite1;
    /** link-site in peptide2 */
    private byte pepsite2;
    /** positions of peptide1 */
    Peptide peptide1;
    /** positions of peptide2 */
    Peptide peptide2;
    /** overall score of this pair */
    private double score;
    /** top-score per charge state */
    private byte charge = 0;
    /**
     * M/Z value of the spectrum
     */
    private double experimentalMZ;
    /**
     * Mass value of the matched peptide(pair)
     */
    private double calcMass;
    /**
     * Charge value of the matched peptide(pair)
     */
    private byte exp_charge;
    
    double peptide1Score;
    double peptide2Score;
    /** hash code for fast access in hashmaps */
    int hashcode;
    /** first peptide comes from a decoy sequence */
    private boolean isDecoy1;
    /** second peptide comes from a decoy sequence */
    private boolean isDecoy2;
    /** short cut for testing of tt (0) td(1) or dd(2) matches */
    private byte isTT_TD_DD;
    /** it could be a protein internal link */
    private boolean isInternal = false;
    /** is this a linear "peptide pair" - meaning only one peptide */
    private boolean isLinear = false;
    /** is this a loop-linked peptide */
    private boolean isLoop = false;
//    /** are all peptides from the target database */
//    private boolean isTT = false;
//    /** is one peptide from the decoy database */
//    private boolean isTD = false;
//    /** are all peptides from the decoy database */
//    private boolean isDD = false;
    /** an id for the fdr-group of this PSM */
    private String fdrGroup;
    /** what is the calculated FDR for the score of this PSM */
    public double m_fdr = -1;
    /** protein group for the first peptide*/
    public ProteinGroup fdrProteinGroup1 = null;
    /** protein group for the second peptide */
    public ProteinGroup fdrProteinGroup2 = null;
    /**
     * arbitrary additional informations - will be written out to csv-files unquoted
     */
    private String   info=null;
    /** 
     * Is this a special case. E.g. PSMs that where matched with unknown 
     * charge-state have a inherently higher chance to be wrong. Therefore it
     * makes sense to calculate the FDR for these separately
     */
    //private boolean specialcase = false;

//    private HashSet<String> negativeGrouping = null;
    
    private PSM partOfUniquePSM;
    
    private final static String NOLINKER = "";
    
    /**
     * The cross-linker for this match
     */
    public String crosslinker = NOLINKER;

    /**
     * the mass of the cross-linker
     */
    private double crosslinkerModMass = Double.NaN;
    
    /**
     * I filter all cross-linker names through this HashMap.
     * This way all PSMs that have the same cross-linker will actually refer to 
     * the same string for the cross-linker and I can compare them via 
     * a.crosslinker == b.crosslinker instead of 
     * a.crosslinker.equals(b.crosslinker).
     * Should be a wee bit faster.
     */
    private static HashMap<String,String> allLinker = new HashMap<String, String>();

    /**
     * what was the score prior normalisation
     */
    private double preNormalisedScore;
    
    /**
     * If we filter to unique PSMs then the top-scoring will represent a set of 
     * psms. As we don't want to lose the info which ones are represented we 
     * list them here.
     */
    private ArrayList<PSM> represents;
    /**
     * This indicates that the PSM is not of a cross-linked peptide pair two
     * peptide that are not cross-linked but only stayed together do to 
     * non-covalent interactions 
     */
    private boolean isNonCovalent = false;
    
    /**
     * the rank of the psm for a the referenced spectrum
     */
    private int rank = 1;
    
    /**
     * for the mzIdentML export we need the index in the referenced files.
     */
    private Integer fileScanIndex;
    /**
     * for the mzIdentML export we need the actual name of the peaklist that was 
     * the spectra originated from.
     */
    private String peakListName;
    
    /**
     * will be set to true if a psm gets assigned a peakListName
     */
    private static boolean peakListNameFound= false;
    

    /**
     * creates a new instance of a PSM.
     * @param psmID a unique ID for the PSM (e.g. used when a CSV-file contains 
     * several rows for a single PSM that represent each possible origin of that 
     * PSM).
     * @param run run-name
     * @param scan scan-id
     * @param peptide1 first peptide
     * @param peptide2 second peptide
     * @param site1 link-site in the first peptide
     * @param site2 link site in the second peptide
     * @param isDecoy1 is the first peptide a decoy?
     * @param isDecoy2 is the second peptide a decoy?
     * @param charge precursor charge state
     * @param score score for this PSM
     * @param scoreRatio how to split the score between the peptides - this is 
     * only used, if a protein FDR (each single protein) is to be calculated.
     */
    public PSM(String psmID, String run, String scan, Peptide peptide1, Peptide peptide2, byte site1, byte site2, boolean isDecoy1, boolean isDecoy2, byte charge, double score, double peptide1Score, double peptide2Score) {
        this(psmID, peptide1, peptide2, (byte)site1, (byte)site2, isDecoy1, isDecoy2, (byte)charge, score, peptide1Score, peptide2Score);
        this.run = run;
        this.scan = scan;
    }
    
    
    /**
     * creates a new instance of a PSM.
     * @param psmID a unique ID for the PSM (e.g. used when a CSV-file contains 
     * several rows for a single PSM that represent each possible origin of that 
     * PSM).
     * @param peptide1 first peptide
     * @param peptide2 second peptide
     * @param site1 link-site in the first peptide
     * @param site2 link site in the second peptide
     * @param isDecoy1 is the first peptide a decoy?
     * @param isDecoy2 is the second peptide a decoy?
     * @param charge precursor charge state
     * @param score score for this PSM
     * @param scoreRatio how to split the score between the peptides - this is 
     * only used, if a protein FDR (each single protein) is to be calculated.
     */
    public PSM(String psmID, Peptide peptide1, Peptide peptide2, byte site1, byte site2, boolean isDecoy1, boolean isDecoy2, byte charge, double score, double peptide1Score, double peptide2Score) {
        this.psmID = psmID;
//        this.id1 = id1;
//        this.id2 = id2;
        
        this.pepsite1 = site1;
        this.pepsite2 = site2;
        this.isDecoy1 = isDecoy1;
        this.isDecoy2 = isDecoy2;
        if (isDecoy1) {
            isTT_TD_DD++;
        } 
        if (isDecoy2) {
            isTT_TD_DD++;
        }
        this.charge = charge;
        this.score = score;
        this.preNormalisedScore = score;
        this.peptide1Score = peptide1Score;
        this.peptide2Score = peptide2Score;
        this.peptide1 = peptide1;
        this.peptide2 = peptide2;
        
        this.hashcode = (peptide1.hashCode() + peptide2.hashCode()) % 100000 * 31 + (pepsite1 + pepsite2) % 10;
        
        // first distinction is whether it could be a protein internal link
        this.isLinear = peptide1 == Peptide.NOPEPTIDE || peptide2 == Peptide.NOPEPTIDE;
        this.isLoop = this.isLinear && site1 >= 0 && site2 >= 0;
        this.isInternal = peptide1.sameProtein(peptide2);
        represents = new ArrayList<PSM>();
        represents.add(this);
        isNonCovalent = pepsite1<=0 && pepsite2<=0 && !isLinear;
        
    }

    @Override
    public Site getLinkSite1() {
        if (peptide1 == Peptide.NOPEPTIDE)
            return null;
        return new PeptideSite(peptide1,pepsite1);
    }

    @Override
    public Site getLinkSite2() {
        if (peptide2 == Peptide.NOPEPTIDE)
            return null;
        return new PeptideSite(peptide2,pepsite2);
    }

    /**
     * hash code of this PSM
     * @return 
     */
    @Override
    public int hashCode() {
        return hashcode;
    }

    /**
     * tests, whether to PSMs are the same.
     * Used to join up PSM entries that where generated for different source of 
     * the peptides (ambiguous peptides/non prototypic peptides)
     * @param l
     * @return 
     */
    @Override
    public boolean equals(Object l) {
        PSM c = (PSM) l;
        if (isNonCovalent != c.isNonCovalent)
            return false;
//        return this.score == c.score && this.charge == c.charge &&  this.psmID.contentEquals(c.psmID);
        return this.score == c.score && this.crosslinker == c.crosslinker && this.charge == c.charge &&  this.psmID.contentEquals(c.psmID) &&
                ((((c.peptide1Score == this.peptide1Score &&c.peptide2Score == this.peptide2Score)) && c.peptide1.equals(this.peptide1) && c.peptide2.equals(this.peptide2) && c.pepsite1 == pepsite1 && c.pepsite2 == pepsite2) 
                || /* to be safe from binary inaccuracy we make an integer-comparison*/
                ((c.peptide1Score == this.peptide2Score &&c.peptide2Score == this.peptide1Score) && c.peptide2.equals(this.peptide1) && c.peptide1.equals(this.peptide2) && c.pepsite2 == pepsite1 && c.pepsite1 == pepsite2));
    }

    /**
     * get a peptide pair represented by this PSM
     * @return 
     */
    public PeptidePair getPeptidePair() {
                             //psmID, id1, id2, site1,    site2,    isDecoy1, isDecoy2, peptide1Positions, peptide2Positions, charge, double score, double scoreRatio
        return new PeptidePair(this);
    }
    
    /**
     * get all peptidepairs represented by this PSM - it's exactly one :)
     * @return 
     */    
    public ArrayList<PeptidePair> getPeptidePairs() {
        ArrayList<PeptidePair> ret =  new ArrayList<PeptidePair>(1);
        ret.add(getPeptidePair());
        return ret;
    }


    /**
     * get all proteins involved with this PSM 
     * @return 
     */
    public ArrayList<Protein> getProteins() {
        HashSet<Protein> prots = new HashSet<Protein>();
        prots.addAll(peptide1.getProteins());
        prots.addAll(peptide2.getProteins());
        ArrayList<Protein> ret = new ArrayList<Protein>(prots);
        return ret;
    }
    
    /**
     * get the proteingroups for both peptides
     * @return 
     */
    public ArrayList<ProteinGroup> getProteinGroupss() {
        HashSet<ProteinGroup> prots = new HashSet<ProteinGroup>();
        for (Peptide p : getPeptides()) {
            if (p != Peptide.NOPEPTIDE) {
                prots.add(p.getProteinGroup());
            }
        }
        ArrayList<ProteinGroup> ret = new ArrayList<ProteinGroup>(prots);
        return ret;
    }
    
    /**
     * Return the supporting proteins, but filter them first through the 
     * given {@link SelfAddHashSet} to make sure that all PSMs return the same
     * instance of the proteins
     * @param allProteins
     * @return 
     */
    public ArrayList<Protein> getProteins(SelfAddHashSet<Protein> allProteins) {
        return allProteins.registerAll(getProteins());
    }


//    public int compareTo(PSM o) {
//        return Double.compare(o.getScore(), this.getScore());
//    }

    /**
     * adds the information of another PSM to this one.
     * @param p 
     */
    @Override
    public void add(PSM p) {
        boolean invertScoreRatio = false;
        if (p == this)
            return;
        if (!equals(p))
            throw new UnsupportedOperationException("Can only add the same PSM together (e.g. for assigning other protein positions)");
        
        
        if (p.isInternal && !isInternal) {
            isInternal = true;
          //  setFDRGroup();
        }
    }

    /**
     * @return the psmID
     */
    public String getPsmID() {
        return psmID;
    }

    /**
     * @return the id of the first peptide
     */
    public Peptide getPeptide1() {
        return peptide1;
    }

    /**
     * @return the id of the second peptide
     */
    public Peptide getPeptide2() {
        return peptide2;
    }
    
    
    

    /**
     * @return the link site of the first peptide
     */
    public byte getPeptideLinkSite1() {
        return pepsite1;
    }

    /**
     * @return the link site of the second peptide
     */
    public byte getPeptideLinkSite2() {
        return pepsite2;
    }

    /**
     * @param peptide
     * @return the link site of the given peptide
     */
    public byte getPeptideLinkSite(int peptide) {
        return peptide == 0 ? pepsite1 : pepsite2;
    }
    

    /**
     * @return the score
     */
    @Override
    public double getScore() {
        return score;
    }

    /**
     * @return the score
     */
    public void setScore(double score) {
        this.score = score;
    }
    
    /**
     * @return the precursor charge state of the PSM
     */
    public byte getCharge() {
        return charge;
    }

//    /**
//     * Score ratio is used to split the score for the as support for the matched 
//     * protein groups
//     * @return the scoreRatio
//     */
//    public double getScoreRatio() {
//        return scoreRatio;
//    }

    public double getPeptide1Score() {
        return peptide1Score;
    }

    public double getPeptide2Score() {
        return peptide1Score;
    }
    
    /**
     * @return is the first peptide a decoy
     */
    public boolean isDecoy1() {
        return isDecoy1;
    }

    /**
     * @return is the second peptide a decoy
     */
    public boolean isDecoy2() {
        return isDecoy2;
    }

    /**
     * This returns whether the PSM could be a protein-internal match. Meaning 
     * is there an overlap between the Protein-groups for each peptide.
     * @return whether the PSM could be a protein-internal match.
     * 
     */
    @Override
    public boolean isInternal() {
        return isInternal;
    }

    /**
     * Are all peptides target peptides
     * @return the isTT
     */
    @Override
    public boolean isTT() {
        return isTT_TD_DD==0;
    }

    /**
     * is one peptide from the decoy database?
     * @return TD
     */
    @Override
    public boolean isTD() {
        return isTT_TD_DD==1;
    }

    /**
     * are all peptides decoy peptides
     * @return the isDD
     */
    @Override
    public boolean isDD() {
        return isTT_TD_DD==2;
    }

    /**
     * numeric id for the FDR-group of this PSM
     * @return 
     */
    @Override
    public String getFDRGroup() {
        return fdrGroup;
    }

//    /**
//     * a name for the FDR-group
//     * @return 
//     */
//    @Override
//    public String getFDRGroupName() {
//        return PeptidePair.fdrGroupNames.get(fdrGroup);
//    }
    
//    /**
//     * Returns the name for the given FDR-group-id
//     * @param fdrgroup
//     * @return 
//     */
//    public static String getFDRGroupName(int fdrgroup) {
//        return PeptidePair.getFDRGroupName(fdrgroup);
//    }    
//    
    /**
     * Returns a list of links supported by this PSM.
     * As the ambiguity is handled in the Protein groups this will report 
     * exactly one link.
     * @return 
     */
    public Collection<ProteinGroupLink> getLinks() {
        return this.getPeptidePair().getLinks();
    }

    
    /**
     * return all peptides involved in this PSM
     * @return 
     */
    public Collection<Peptide> getPeptides() {
        ArrayList<Peptide> ret = new ArrayList<Peptide>(2);
        if (getPeptide2() != Peptide.NOPEPTIDE) {
            ret.add(getPeptide2());
        }
        ret.add(getPeptide1());
        return ret;
    }
    
    /**
     * A textual representation of this PSM
     * @return 
     */
    @Override
    public String toString() {
        
        if (isTT())
            return "TT - " + peptide1.toString() + " xl " + peptide2.toString();
        if (isTD())
            return "TD - " + peptide1.toString() + " xl " + peptide2.toString();
        
        return "DD - " + peptide1.toString() + " xl " + peptide2.toString();
        
    }
    
    /**
     * is at least one of the peptides a decoy
     * @return 
     */
    @Override
    public boolean isDecoy() {
        return isDecoy1 || isDecoy2;
    }
    
    
    /**
     * set the FDR group according to the information on this PSM
     */
    public void setFDRGroup() {
        
        fdrGroup = PeptidePair.getFDRGroup(peptide1, peptide2, isLinear(), isInternal, this.getNegativeGrouping(), getPositiveGrouping(),isNonCovalent ? "NonCovalent":"");
       
    }     

    /**
     * set the FDR group according to the information on this PSM
     */
    public void setFDRGroup(String fdrGroup) {
        
        this.fdrGroup = fdrGroup;
       
    }     

    /**
     * store the FDR assigned to the score of this PSM
     * @param fdr 
     */
    @Override
    public void setFDR(double fdr) {
        m_fdr = fdr;
    }

    
    /**
     * The FDR assigned to the score of this PSM
     * @return the score FDR
     */
    @Override
    public double getFDR() {
        return m_fdr;
    }


    /**
     * @return the length of peptide 1
     */
    public int getPeptideLength1() {
        return peptide1.length();
    }

    /**
     * @return the length of peptide 2
     */
    public int getPeptideLength2() {
        return peptide2.length();
    }

    /**
     * @param peptide
     * @return the length of the given peptide
     */
    public int getPeptideLength(int peptide) {
        return peptide ==0 ? peptide1.length() : peptide2.length() ;
    }

    

    /**
     * Set the peptide pair that was supported by this PSM and passed the 
     * peptide pair FDR
     * @param pp
     */
    public void setFdrPeptidePair(PeptidePair pp) {
        this.m_peppair = pp;
    }

    /**
     * Get the peptide pair that was supported by this PSM and passed the 
     * peptide pair FDR
     * @return 
     */
    public PeptidePair getFdrPeptidePair() {
        return this.m_peppair;
    }
    
    /**
     * reset the stored information, that might be different for a new FDR-calculation
     */
    public void reset() {
        resetFdrProteinGroup();
        setFdrPeptidePair(null);
        setLinkedSupport(1);
        setPartOfUniquePSM(null);
        represents.clear();
        represents.add(this);
    }
    
    /**
     * unset the fdrproteingroups
     */
    public void resetFdrProteinGroup() {
            this.fdrProteinGroup1 = null;
            this.fdrProteinGroup2 = null;
        
    }
    
    
    /**
     * set the ProteinGroup that passed the fdr and is supported by this PSM
     * @param pg
    */
    public void setFdrProteinGroup(ProteinGroup pg) {
        if (pg == null) {
            this.fdrProteinGroup1 = null;
            this.fdrProteinGroup2 = null;
        } else {
            if (pg.equals(peptide1.getProteinGroup()))
                this.fdrProteinGroup1 = pg;
            if (pg.equals(peptide2.getProteinGroup()))
                this.fdrProteinGroup2 = pg;
        }
    }
    
    /**
     * get the ProteinGroup that passed the fdr and is supported by the first
     * peptide of this PSM 
     * @return 
     */
    public ProteinGroup getFdrProteinGroup1() {
        return this.fdrProteinGroup1;
    }
    

    /**
     * get the ProteinGroup that passed the fdr and is supported by the second
     * peptide of this PSM 
     * @return 
     */
    public ProteinGroup getFdrProteinGroup2() {
        return this.fdrProteinGroup2;
    }

    /** 
     * number of supporting peptidepairs - is for a psm always 1.
     * @return  
     */
    @Override
    public int getPeptidePairCount() {
        return 1;
    }

//    /**
//     * Is this a special case. E.g. PSMs that where matched with unknown 
//     * charge-state have a inherently higher chance to be wrong. Therefore it
//     * makes sense to calculate the FDR for these separately
//     * @return the specialcase
//     */
//    public boolean hasNegativeGrouping() {
//        return negativeGrouping != null;
//    }
//
//    public HashSet<String> getNegativeGrouping() {
//        return negativeGrouping;
//    }
//
//    /**
//     * Is this a special case. E.g. PSMs that where matched with unknown 
//     * charge-state have a inherently higher chance to be wrong. Therefore it
//     * makes sense to calculate the FDR for these separately
//     * @param specialcase the specialcase to set
//     */
//    public void setNegativeGrouping(boolean specialcase) {
//        this.negativeGrouping = specialcase?"Special":null;
//    }

//    /**
//    /**
//     * Is this a special case. E.g. PSMs that where matched with unknown 
//     * charge-state have a inherently higher chance to be wrong. Therefore it
//     * makes sense to calculate the FDR for these separately
//     * @param specialcase the reason specialcase to set
//     */
//    public void setNegativeGrouping(String specialcase) {
//        if (specialcase == null) {
//            this.negativeGrouping = null;
//        } else {
//            this.negativeGrouping = new HashSet<String>();
//            this.negativeGrouping.add(specialcase);
//        }
//    }
    
    /**
     * Is this a non-cross-linked PSM?
     * @return the isLinear
     */
    @Override
    public boolean isLinear() {
        return isLinear;
    }

    /**
     * Is this a non-cross-linked PSM?
     * @param isLinear 
     */
    public void setLinear(boolean isLinear) {
        this.isLinear = isLinear;
    }

    /**
     * the run-name for the PSM
     * @return the run
     */
    public String getRun() {
        return run;
    }

    /**
     * the run-name for the PSM
     * @param run the run to set
     */
    public void setRun(String run) {
        this.run = run;
    }

    /**
     * the scan-id for the PSM
     * @return the scan
     */
    public String getScan() {
        return scan;
    }

    /**
     * the scan-id for the PSM
     * @param scan the scan to set
     */
    public void setScan(String scan) {
        this.scan = scan;
    }

    /**
     * Get the Peptide for the given site.
     * @param n
     * @return 
     */
    @Override
    public Object getSite(int n) {
        return n==0 ? getPeptide1() : getPeptide2();
    }

    /**
     * how many sites does this PSM have.
     * Is either 2 for cross-linked PSMs or 1 for linear PSMs
     * @return 
     */
    @Override
    public int getSites() {
        if (isLinear)
            return 1;
        else 
            return 2;
    }

    /**
     * is this a cross-linked PSM where both peptides come definitely from two 
     * different Proteins
     * @return 
     */
    @Override
    public boolean isBetween() {
        return !isLinear && !isInternal;
    }

    /**
     * is this a linear loop-linked peptide PSM
     * @return the isLoop
     */
    public boolean isLoop() {
        return isLoop;
    }

    /**
     * is this a linear loop-linked peptide PSM
     * @param isLoop the isLoop to set
     */
    public void setLoop(boolean isLoop) {
        this.isLoop = isLoop;
    }

    /**
     * M/Z value of the spectrum
     * @return the experimentalMZ
     */
    public double getExperimentalMZ() {
        return experimentalMZ;
    }

    /**
     * M/Z value of the spectrum
     * @param experimentalMZ the experimentalMZ to set
     */
    public void setExperimentalMZ(double experimentalMZ) {
        this.experimentalMZ = experimentalMZ;
    }

    /**
     * Mass value of the matched peptide(pair)
     * @return the calcMZ
     */
    public double getCalcMass() {
        return calcMass;
    }

    /**
     * Mass value of the matched peptide(pair)
     * @param calcMass
     */
    public void setCalcMass(double calcMass) {
        this.calcMass = calcMass;
    }

    /**
     * Charge value of the matched peptide(pair)
     * @return the match_charge
     */
    public int getExpCharge() {
        return exp_charge;
    }

    /**
     * Charge value of the matched peptide(pair)
     * @param match_charge the match_charge to set
     */
    public void setExpCharge(byte match_charge) {
        this.exp_charge = match_charge;
    }

    /**
     * @return the crosslinker
     */
    public String getCrosslinker() {
        return crosslinker;
    }

    /**
     * @param crosslinker the crosslinker to set
     */
    public void setCrosslinker(String crosslinker) {
        
        // doing this makes sure I can just make a comparison via = instead of String.equal
        // and I expect that the number of cros-linker in anygiven settup will be rather limited
        String prevXL = allLinker.get(crosslinker);
        if (prevXL == null) {
            prevXL = crosslinker;
            allLinker.put(crosslinker, crosslinker);
        }

        this.crosslinker = prevXL;
    }

    /**
     * @return the info
     */
    public String getInfo() {
        return info;
    }

    /**
     * @param info the info to set
     */
    public void setInfo(String info) {
        this.info = info;
    }

    /**
     * @return the normScore
     */
    public double getNonnormalizedScore() {
        return preNormalisedScore;
    }

    /**
     * @param normScore the normScore to set
     */
    public void setNonnormalizedScore(double score) {
        this.preNormalisedScore = score;
    }

    /**
     * @return the partOfUniquePSM
     */
    public PSM getPartOfUniquePSM() {
        return partOfUniquePSM;
    }

    /**
     * @param partOfUniquePSM the partOfUniquePSM to set
     */
    public void setPartOfUniquePSM(PSM partOfUniquePSM) {
        this.partOfUniquePSM = partOfUniquePSM;
    }
    
    
    public void represents(PSM psm) {
        this.represents.add(psm);
        if (psm.represents.size() >0) {
            this.represents.addAll(psm.represents);
        }
    }
    
    public ArrayList<PSM> getRepresented() {
        return represents;
    }

    @Override
    public ProteinGroup getProteinGroup1() {
        return getPeptide1().getProteinGroup();
    }

    @Override
    public ProteinGroup getProteinGroup2() {
        return getPeptide2().getProteinGroup();
    }

    /**
     * @return the isNonCovalent
     */
    public boolean isNonCovalent() {
        return isNonCovalent;
    }

    /**
     * @param isNonCovalent the isNonCovalent to set
     */
    public void setNonCovalent(boolean isNonCovalent) {
        this.isNonCovalent = isNonCovalent;
    }

//    @Override
//    public boolean hasPositiveGrouping() {
//        return this.positiveGroups != null;
//    }
//    
//    @Override
//    public void setPositiveGrouping(String av) {
//        if (av == null)
//            this.positiveGroups = null;
//        else {
//            this.positiveGroups = new HashSet<String>();
//            this.positiveGroups.add(av);
//        }
//    }
//
//    @Override
//    public HashSet<String> getPositiveGrouping() {
//        return this.positiveGroups;
//    }

    /**
     * the rank of the psm for a the referenced spectrum
     * @return the rank
     */
    public int getRank() {
        return rank;
    }

    /**
     * the rank of the psm for a the referenced spectrum
     * @param rank the rank to set
     */
    public void setRank(int rank) {
        this.rank = rank;
    }

    
    /**
     * @return the fileScanIndex
     */
    public Integer getFileScanIndex() {
        return fileScanIndex;
    }

    /**
     * @param fileScanIndex the fileScanIndex to set
     */
    public void setFileScanIndex(int fileScanIndex) {
        this.fileScanIndex = fileScanIndex;
    }    
    
    
    
    /**
     * @return the fileScanIndex
     */
    public String getPeakListName() {
        return peakListName;
    }

    /**
     * @param fileScanIndex the fileScanIndex to set
     */
    public void setPeakListName(String peakListName) {
        this.peakListName = peakListName;
        this.peakListNameFound = true;
    }    

    /**
     * @returwas a peakfilename provided
     */
    public static boolean getPeakListNameFound() {
        return peakListNameFound;
    }

    /**
     * the mass of the cross-linker
     * @return the crosslinkerModMass
     */
    public double getCrosslinkerModMass() {
        return crosslinkerModMass;
    }

    /**
     * the mass of the cross-linker
     * @param crosslinkerModMass the crosslinkerModMass to set
     */
    public void setCrosslinkerModMass(double crosslinkerModMass) {
        this.crosslinkerModMass = crosslinkerModMass;
    }
    
}
