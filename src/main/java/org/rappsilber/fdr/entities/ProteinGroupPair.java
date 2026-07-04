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
import java.util.Comparator;
import java.util.HashSet;
import java.util.TreeSet;
import org.rappsilber.fdr.utils.FDRGroupNames;
import org.rappsilber.utils.IntArrayList;
import org.rappsilber.utils.RArrayUtils;

/**
 *
 * @author lfischer
 */
public class ProteinGroupPair extends AbstractFDRElement<ProteinGroupPair> { //implements Comparable<ProteinGroupPair> {

    /**
     * counts all instances of this class.
     * used for assigning unique IDs
     */
    public static int PROTEINGROUPPAIRCOUNT = 0;
    
    
    protected int proteinGroupPairID = PROTEINGROUPPAIRCOUNT++;
    private HashSet<ProteinGroupLink> links = new HashSet<ProteinGroupLink>();
    private TreeSet<ProteinGroupLink> sortedLinks = new TreeSet<>(new Comparator<ProteinGroupLink>() {
        @Override
        public int compare(ProteinGroupLink o1, ProteinGroupLink o2) {
            return Double.compare(o2.getScore(), o1.getScore());
        }
    });
    private ProteinGroup protein1;
    private ProteinGroup protein2;
    private double score;
    private Double scoreTopN;
    private int lastTopN = 0;
    int hashcode;
    /** first peptide comes from a decoy sequence */
    private boolean isDecoy1;
    /** second peptide comes from a decoy sequence */
    private boolean isDecoy2;
    /** it could be a protein internal link */
    private boolean isInternal = false;
    /** self-link supported only by overlapping peptide pairs */
    private boolean isOverlapping = false;
    private boolean isTT = false;
    private boolean isTD = false;
    private boolean isDD = false;
    private String fdrGroup;
    public double m_fdr = -1;
//    private boolean m_specialOnly = true;
    private boolean isNonCovalent = false;
//    private HashSet<String> positiveGroups;
//    private HashSet<String> m_NegativeGrouping;

    
    protected ProteinGroupPair(ProteinGroup Prot1, ProteinGroup Prot2, double score, boolean isSpecialOnly) {
        this.protein1 = Prot1;
        this.protein2 = Prot2;
        
        hashcode = (Prot1.hashCode() + Prot2.hashCode()) % 10000;
        
        this.isDecoy1 = !Prot1.isTT();;
        this.isDecoy2 = !Prot2.isTT();;
        this.isInternal = Prot1.hasOverlap(Prot2);
        this.score = score;
        if (isDecoy1) {
            if (isDecoy2) {
                isDD = true;
            } else {
                isTD = true;
            }
        } else if (isDecoy2) {
            isTD = true;
        } else {
            isTT = true;
        }
        if (isSpecialOnly)
            setNegativeGrouping("Special");
        
    }

    public void setFDRGroup() {
        if (protein1 == protein2 || protein1.hasOverlap(protein2))
            fdrGroup = isOverlapping ? "Self Overlapping" : "Self";
        else
            fdrGroup = "Between";

        if (m_negativeGroups != null)
            fdrGroup += " [n" + RArrayUtils.toString(m_negativeGroups,", n") +"]";

        if (isNonCovalent) 
            fdrGroup += " NonCovalent";

        if (m_positiveGroups!= null)
            fdrGroup += " has [p" + RArrayUtils.toString(m_positiveGroups,", p") + "]";
        
        fdrGroup = FDRGroupNames.get(fdrGroup);
    }
    
    
    public ProteinGroupPair(PeptidePair pp) {
        this(pp.getPeptide1().getProteinGroup(), pp.getPeptide2().getProteinGroup(), pp.getScore(), pp.hasNegativeGrouping());
        isNonCovalent = pp.isNonCovalent();
        this.m_positiveGroups = pp.getPositiveGrouping();
        this.m_negativeGroups = pp.getNegativeGrouping();
        this.isOverlapping = pp.isOverlapping();
    }

    public ProteinGroupPair(ProteinGroupLink l) {
        this(l.getProteinGroup1(), l.getProteinGroup2(),  l.getScore(), l.hasNegativeGrouping());
        this.links.add(l);
        isNonCovalent = l.isNonCovalent();
        this.m_positiveGroups = l.getPositiveGrouping();
        this.m_negativeGroups = l.getNegativeGrouping();
        this.sortedLinks.add(l);
        this.isOverlapping = l.isOverlapping();
    }

    @Override
    public Site getLinkSite1() {
        return new ProteinGroupSite(protein1);
    }

    @Override
    public Site getLinkSite2() {
        return new ProteinGroupSite(protein2);
    }

    @Override
    public int hashCode() {
        return isOverlapping ? hashcode ^ 0x8000 : hashcode;
    }

    @Override
    public boolean equals(Object l) {
        if (l ==  this)
            return true;
        ProteinGroupPair c = (ProteinGroupPair) l;
        if (isNonCovalent != c.isNonCovalent)
            return false;
        if (isInternal && isOverlapping != c.isOverlapping)
            return false;
        return (c.protein1.equals(protein1) && c.protein2.equals(protein2)) || (c.protein1.equals(protein2) && c.protein2.equals(protein1));
    }

    public int compareTo(ProteinGroupPair o) {
        return Double.compare(o.getScore(), this.getScore());
    }

    public void add(ProteinGroupPair pp) {
        if (pp == this)
            return;
        this.links.addAll(pp.links);
        this.sortedLinks.addAll(pp.links);
        this.score = Math.sqrt(this.score * this.score + pp.score * pp.score);
        this.lastTopN=0;
        addFDRGroups(pp);
        
    }

    /**
     * @return the links
     */
    public HashSet<ProteinGroupLink> getLinks() {
        return links;
    }

    /**
     * @return the protein1
     */
    public ProteinGroup getProtein1() {
        return protein1;
    }

    /**
     * @return the protein2
     */
    public ProteinGroup getProtein2() {
        return protein2;
    }

    /**
     * @return the score
     */
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
     * @return the score
     */
    public double getScore(int topN) {
        if (topN == this.lastTopN)
            return this.scoreTopN;
        int i=0;
        for (ProteinGroupLink l : this.sortedLinks) {
            if (i++>topN)
                break;
            score+= l.getScore()*l.getScore();
        }
        this.lastTopN = topN;
        this.score=Math.sqrt(score);
        return this.score;
    }
    
//    /**
//     * @return the fdrGroup
//     */
//    public int getFdrGroup() {
//        return fdrGroup;
//    }


    /**
     * @return the isDecoy1
     */
    public boolean isDecoy1() {
        return isDecoy1;
    }

    /**
     * @return the isDecoy2
     */
    public boolean isDecoy2() {
        return isDecoy2;
    }

    /**
     * @return the isInternal
     */
    public boolean isInternal() {
        return isInternal;
    }

    public boolean isOverlapping() {
        return isOverlapping;
    }

    public void setOverlapping(boolean isOverlapping) {
        this.isOverlapping = isOverlapping;
        this.fdrGroup = null; // force recomputation via setFDRGroup()
    }

    /**
     * @return the isTT
     */
    public boolean isTT() {
        return isTT;
    }

    /**
     * @return the isTD
     */
    public boolean isTD() {
        return isTD;
    }

    /**
     * @return the isDD
     */
    public boolean isDD() {
        return isDD;
    }

    public String getFDRGroup() {
        if (fdrGroup == null) {
            setFDRGroup();
        }
        return fdrGroup;
    }


    public static String getFDRGroupName(int fdrGroup) {
        int g=fdrGroup;
        if (fdrGroup < -1)
            fdrGroup = (-fdrGroup)-2;
        String n;
        switch (fdrGroup) {
            case -1 : n= "all combined"; break;
            case 0 : n= "Between"; break;
            case 1 : n= "Within"; break;
            case 2 : n= "Special Between"; break;
            case 3 : n= "Special Within"; break;
            default : return "Unknown"; 
        }
        if (g<-2) {
            return "NonCovalent " + n;
        }
        return n;
        
    }
    
    
    public Collection<PeptidePair> getPeptidePairs() {
        HashSet<PeptidePair> ret = new HashSet<PeptidePair>();
        for (ProteinGroupLink l : links)
            ret.addAll(l.getPeptidePairs());
        return ret;
    }

    public int[] getPeptidePairID() {
        IntArrayList retal = new IntArrayList();
        for (ProteinGroupLink l : links)
            retal.addAll(l.getPeptidePairIDs());
        return retal.toIntArray();
    }
    
    public Collection<ProteinGroupPair> getProteinPairs() {
        ArrayList<ProteinGroupPair> ret = new ArrayList<ProteinGroupPair>();
        ret.add(this);
        return ret;
        
    }

    public Collection<Protein> getProteins() {
        ArrayList<Protein> ret = new ArrayList<Protein>(2);
        ret.addAll(getProtein1().getProteins());
        ret.addAll(getProtein2().getProteins());
        return ret;
    }

    
    
    public String toString() {
        if (isTT()) 
            return "TT - " + getProtein1().toString() + " xl " + getProtein2().toString();
        if (isTD()) 
            return "TD - " + getProtein1().toString() + " xl " + getProtein2().toString();
        return "DD - " + getProtein1().toString() + " xl " + getProtein2().toString();
    }
    
    public boolean isDecoy() {
        return isDecoy1 || isDecoy2;
    }
 
    
    public void setFDR(double fdr) {
        m_fdr = fdr;
        // flag up Links
        for (ProteinGroupLink pgl : links) {
            pgl.setFdrPPI(this);
        }
    }

    public double getFDR() {
        return m_fdr;
    }

    
    /**
     * @return the peptidePairID
     */
    public int getProteinGroupPairID() {
        return proteinGroupPairID;
    }
    
    /**
     * @return the peptidePairID
     */
    public int[] getLinkIDs() {
        int[] ret = new int[links.size()];
        int c =0;
        for (ProteinGroupLink l : links) {
            ret[c++] = l.getLinkID();
        }
        return ret;
    }

    /**
     * @return the peptidePairID
     */
    public int[] getPeptidePairIDs() {
        IntArrayList ret = new IntArrayList();
        for (ProteinGroupLink l : links) {
            ret.addAll(l.getPeptidePairIDs());
        }
        return ret.toIntArray();
    }
    
    public String[] getPSMIDs() {

        ArrayList<String> ret = new ArrayList<String>();
        for (ProteinGroupLink l : links) {
            for (String id : l.getPSMIDs())
                ret.add(id);
        }
        String[] reta = new String[ret.size()];
        return ret.toArray(reta);
        
    }

    public int getPeptidePairCount() {
        int ret =0;
        for (ProteinGroupLink l : links)
            ret+=l.getPeptidePairCount();
        return ret;
    }

    public Object getSite(int n) {
        return n==0 ? getProtein1() : getProtein2();
    }

    public int getSites() {
        return 2;
    }    

    public boolean isLinear() {
        return false;
    }

    public boolean isBetween() {
        return !isInternal;
    }

    @Override
    public ProteinGroup getProteinGroup1() {
        return protein1;
    }

    @Override
    public ProteinGroup getProteinGroup2() {
        return protein2;
    }

    /**
     * @param fdrGroup the fdrGroup to set
     */
    public void setFDRGroup(String fdrGroup) {
        this.fdrGroup = FDRGroupNames.get(fdrGroup);
    }
    
    
//    @Override
//    public boolean hasPositiveGrouping() {
//        return this.positiveGroups!=null;
//    }
//    
//    @Override
//    public void setPositiveGrouping(String av) {
//        if (av == null) {
//            this.positiveGroups = null;
//        }else {
//            this.positiveGroups = new HashSet<String>(1);
//            this.positiveGroups.add(av);
//        }
//    }
//
//    @Override
//    public HashSet<String> getPositiveGrouping() {
//        return positiveGroups;
//    }
        
//    /**
//     * are all supporting PSMs "special" cases?
//     * @return the specialcase
//     */
//    public boolean hasNegativeGrouping() {
//        return m_NegativeGrouping!=null;
//    }
//
////    /**
////     * are all supporting PSMs "special" cases?
////     * @param specialcase 
////     */
////    public void setNegativeGrouping(boolean specialcase) {
////        if (specialcase) {
////            this.m_NegativeGrouping = "Special";
////        } else {
////            this.m_NegativeGrouping = null;
////        }
////    }
//
//    /**
//     * are all supporting PSMs "special" cases?
//     * @param specialcase 
//     */
//    @Override
//    public void setNegativeGrouping(String cause) {
//        if (cause == null) {
//            this.m_NegativeGrouping = null;
//        } else {
//            this.m_NegativeGrouping = new HashSet<>(1);
//            this.m_NegativeGrouping.add(cause);
//        }
//    }
//    
//    @Override
//    public HashSet<String> getNegativeGrouping() {
//        return this.m_NegativeGrouping;
//    }
 
    public boolean isNonCovalent() {
        return isNonCovalent;
    }

    
}
