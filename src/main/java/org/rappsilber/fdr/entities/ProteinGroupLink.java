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
import java.util.Iterator;
import org.rappsilber.fdr.groups.ProteinGroup;
import org.rappsilber.fdr.utils.AbstractFDRElement;
import org.rappsilber.utils.IntArrayList;
import org.rappsilber.utils.MapUtils;

/**
 *
 * @author lfischer
 */
public class ProteinGroupLink extends AbstractFDRElement<ProteinGroupLink> { //implements Comparable<ProteinGroupLink> {

    public static int LINKCOUNT = 0;
    
    protected int linkID = LINKCOUNT++;
    
    private double score = Double.POSITIVE_INFINITY;
    public boolean isInternal = false;
    private ProteinGroup pg1;
    private HashMap<Protein, IntArrayList> position1 = new HashMap<Protein, IntArrayList>(1);
    private ProteinGroup pg2;
    private HashMap<Protein, IntArrayList> position2 = new HashMap<Protein, IntArrayList>(1);;
    int hashcode = 5;
    private HashSet<PeptidePair> support = new HashSet<PeptidePair>(1);
    private PeptidePair m_reference;
    public double m_fdr = -1;
    protected double m_PPIfdr = -1;
    protected ProteinGroupPair m_ppi;
    protected boolean m_linkssorted = false;
    private boolean m_specialOnly = true;
    private HashMap<String, HashSet<String>> runtoScan = null;  
    
    public static int MIN_DISTANCE_FOR_LONG = 0;
    


    public ProteinGroupLink(PeptidePair pp) {
        Peptide peptide1 =pp.getPeptide1();
        Peptide peptide2 = pp.getPeptide2();
        int pepSite1 = pp.getPeptideLinkSite1();
        int pepSite2 = pp.getPeptideLinkSite2();
        position1=new HashMap<Protein, IntArrayList>(peptide1.getPositions().size());
        for (Protein p : peptide1.getPositions().keySet()) {
            IntArrayList peppos = peptide1.getPositions().get(p).clone();
            peppos.addToAll(pepSite1-1);
            position1.put(p, peppos);
        }

        pg1 = pp.getFdrProteinGroup1();
        if (pg1==null)
            pg1=peptide1.getProteinGroup();
        
        
        if (pp.isLoop()) {
            position2 = position1;
            pg2 = pg1;
            isInternal = true;
        } else {
            position2=new HashMap<Protein, IntArrayList>(peptide2.getPositions().size());
            for (Protein p : peptide2.getPositions().keySet()) {
                IntArrayList peppos = peptide2.getPositions().get(p).clone();
                peppos.addToAll(pepSite2-1);
                position2.put(p, peppos);
            }
            pg2 = pp.getFdrProteinGroup2();
            if (pg2 == null)
                pg2 = peptide2.getProteinGroup();
            isInternal = pg1.hasOverlap(pg2);
        }
        
        
        //if we have run and scan we can implement some additional checks
        Iterator<PSM> psmIter = pp.getTopPSMs().iterator();
        PSM firstpsm = psmIter.next();
        if (firstpsm.getRun() != null  && firstpsm.getScan() != null && !firstpsm.getScan().isEmpty()) {
            runtoScan = new HashMap<String, HashSet<String>>();
            HashSet<String> scans = new HashSet<String>(1);
            scans.add(firstpsm.getScan());
            runtoScan.put(firstpsm.getRun(), scans);
            while (psmIter.hasNext()) {
                PSM psm = psmIter.next();
                String  run = psm.getRun();
                String scan = psm.getScan();
                scans = runtoScan.get(run);
                if (scans == null) {
                    scans = new HashSet<String>(1);
                    runtoScan.put(run, scans);
                }
                scans.add(scan);
            }
        }
        
        
        
        support.add(pp);
        m_reference = pp;
        score = pp.getScore();
        hashcode += (pg1.hashCode() + pg2.hashCode())*100;
        for (IntArrayList ia : position1.values()) {
            for (int p : ia) 
                hashcode += p;
        }
        for (IntArrayList ia : position2.values()) {
            for (int p : ia) 
                hashcode += p;
        }
        m_specialOnly = pp.isSpecialcase();
    }

    @Override
    public Site getLinkSite1() {
        return new LinkSite(this,0);
    }

    @Override
    public Site getLinkSite2() {
        return new LinkSite(this,1);
    }


    @Override
    public int hashCode() {
        return hashcode;
    }

    /**
     * a link group is considered equal, if it contains all the same links
     * @param o
     * @return
     */
    @Override
    public boolean equals(Object o) {
        // if it is actaually the same object
        if (this == o)
            return true;
        
        if (o instanceof ProteinGroupLink) {
            
            ProteinGroupLink pgl = (ProteinGroupLink) o;

            // is it a complete internal link? - meaning a link within the same protein-group (as opposed links between groups that contain a comon protein)
            if (pg1.equals(pg2) && pg1.equals(pgl.pg1) && pg1.equals(pgl.pg2))
                return ((MapUtils.sameMapCollection(position1, pgl.position1) && MapUtils.sameMapCollection(position2, pgl.position2)) ||
                        (MapUtils.sameMapCollection(position2, pgl.position1) && MapUtils.sameMapCollection(position1, pgl.position2)));
            // same groups
            else if (pg1.equals(pgl.pg1) && pg2.equals(pgl.pg2))
                return (MapUtils.sameMapCollection(position1, pgl.position1) && MapUtils.sameMapCollection(position2, pgl.position2));
            // same groups but reversed
            else if (pg2.equals(pgl.pg1) && pg1.equals(pgl.pg2))
                return (MapUtils.sameMapCollection(position2, pgl.position1) && MapUtils.sameMapCollection(position1, pgl.position2));
        }
        return false;

    }

    public double getScore() {
        if (score != Double.POSITIVE_INFINITY) {
            return score;
        }
        score = 0;
        for (PeptidePair pp : support) {
            score+=(pp.getScore() * pp.getScore());
        }
        score= Math.sqrt(score);
        return score;
    }

 
    public int compareTo(ProteinGroupLink o) {
        return Double.compare(o.getScore(), this.getScore());
    }

    public void add(ProteinGroupLink o) {
        double addScore = 0;
        if (o == this)
            return;
        if (runtoScan != null) {
            for (PeptidePair pp : o.getPeptidePairs()) {
                support.add(pp);
                for (PSM psm : pp.getTopPSMs()) {
                    String run = psm.getRun();
                    String scan = psm.getScan();
                    HashSet<String> scans = runtoScan.get(run);
                    if (scans == null) {
                        scans = new HashSet<String>(1);
                        scans.add(scan);
                        runtoScan.put(run, scans);
                        addScore += psm.getScore() * psm.getScore();
                    } else {
                        if (!scans.contains(scan)) {
                            // this is only the second match for the same spectra supporting the same peptide.
                            // this usally happens if there is a modification with two equally possible positions
                            addScore += psm.getScore() * psm.getScore();
                            scans.add(scan);
                        }
                    }

                }
            }        
            this.score = Math.sqrt(addScore + score* score);
        } else {
            support.addAll(o.getPeptidePairs());
            this.score = Math.sqrt(score*score + o.getScore() * o.getScore());
        }
        this.m_specialOnly &= o.m_specialOnly;
    }


    public int getFDRGroup() {
        int group = (isInternal?1:0) + (m_specialOnly?2:0);
        if (MIN_DISTANCE_FOR_LONG > 0 && isInternal && getProteinGroup1().size() + getProteinGroup1().size() ==2) {
            for (int f : getPosition1().values().iterator().next()) {
                for (int t : getPosition2().values().iterator().next()) {
                    if (Math.abs(f-t) > MIN_DISTANCE_FOR_LONG) {
                        group+=4;
                    }
                }
            }
        }
        return group;
    }

    public String getFDRGroupName() {
        return getFDRGroupName(getFDRGroup());
    }

    public static String getFDRGroupName(int group) {
        String name = "";
        switch(group % 4) {
            case 0 : name = "Between";
                break;
            case 1 : name =  "Within";
                break;
            case 2 : name = "Special Between";
                break;
            case 3 : name = "Special Within";
                break;
            case -1 : name = "ALL";
        }
        if (group / 4 == 1 && !name.isEmpty()) {
            name += " long";
        }
        if (name.isEmpty())
            return "Unknown";
        return name;
    }
    
//    public static String getFDRGroupName(int group) {
//        return group == 1? "Internal" : "Between";
//    }
    
    public int getAmbiguity () {
        int c1 = 0;
        for (IntArrayList ia : position1.values()) {
            c1+=ia.size();
        }
        
        int c2 = 0;
        for (IntArrayList ia : position2.values()) {
            c2+=ia.size();
        }
        
        return c1*c2;

    }
    
    public boolean isTT() {
        return m_reference.isTT();
    }
    public boolean isTD() {
        return m_reference.isTD();
    }
    public boolean isDD() {
        return m_reference.isDD();
    }

    public ProteinGroup getProteinGroup1() {
        return pg1;
    }
    public ProteinGroup getProteinGroup2() {
        return pg2;
    }

    public Collection<PeptidePair> getPeptidePairs() {
        return support;
    }


    public int[] getPeptidePairIDs() {
        int[] ret = new int[support.size()];
        int c =0;
        for (PeptidePair pp : support)
            ret[c++] = pp.getPeptidePairID();
        return ret;
    }

    public String[] getPSMIDs() {
        ArrayList<String> ret = new ArrayList<String>();
        for (PeptidePair pp : support) {
            for (String id : pp.getPSMids())
                ret.add(id);
        }
        String[] reta = new String[ret.size()];
        return ret.toArray(reta);
    }
    

    public Collection<Protein> getProteins() {
        HashSet<Protein> ret = new HashSet<Protein>();
        ret.addAll(pg1.getProteins());
        ret.addAll(pg2.getProteins());
        return ret;
    }

    public Collection<Peptide> getPeptides() {
        HashSet<Peptide> ret = new HashSet<Peptide>();
        for  (PeptidePair pp : support)
            ret.addAll(pp.getPeptides());
        return ret;
    }
    
    
    public ProteinGroupPair getProteinGroupPair() {
        return new ProteinGroupPair(this);
    }
    
    public String site1ToString() {
        StringBuffer sb = new StringBuffer(); 
        for (Protein p : position1.keySet()) {
            IntArrayList pos = position1.get(p);
            for (Integer i : pos) {
                sb.append(";");
                sb.append(p.getAccession());
                sb.append("(");
                sb.append(i);
                sb.append(")");
            }
        }
        return sb.substring(1);
    }

    public String site2ToString() {
        StringBuffer sb = new StringBuffer(); 
        for (Protein p : position2.keySet()) {
            IntArrayList pos = position2.get(p);
            for (Integer i : pos) {
                sb.append(";");
                sb.append(p.getAccession());
                sb.append("(");
                sb.append(i);
                sb.append(")");
            }
        }
        return sb.substring(1);
    }

    public String site1Accessions() {
        return getAccessions(position1);
    }

    public String site1Descriptions() {
        return getDescriptions(position1);
    }

    public String site2Accessions() {
        return getAccessions(position2);
    }    

    public String site2Descriptions() {
        return getDescriptions(position2);
    }
    
    public String site1Sites() {
        StringBuffer sb = new StringBuffer(); 
        for (Protein p : position1.keySet()) {
            IntArrayList pos = position1.get(p);
            for (Integer i : pos) {
                sb.append(";");
                sb.append(i);
            }
        }
        return sb.substring(1);
    }

    public String site2Site() {
        StringBuffer sb = new StringBuffer(); 
        for (Protein p : position2.keySet()) {
            IntArrayList pos = position2.get(p);
            for (Integer i : pos) {
                sb.append(";");
                sb.append(i);
            }
        }
        return sb.substring(1);
    }    
    

    public String toString() {
        if (isTT()) 
            return "TT - " + site1ToString() + " xl " + site2ToString();
        if (isTD()) 
            return "TD - " + site1ToString() + " xl " + site2ToString();
        return "DD - " + site1ToString() + " xl " + site2ToString();
    }
    
    public boolean isDecoy() {
        return pg1.isDecoy() || pg2.isDecoy();
    }    
    
    public void setFDR(double fdr) {
        m_fdr = fdr;
        // flag up the peptide pairs
        for (PeptidePair pp : support) {
            pp.setFdrLink(this);
        }
    }

    public double getFDR() {
        return m_fdr;
    }

//    /**
//     * @return the m_PPIfdr
//     */
//    public double getPPIfdr() {
//        return m_PPIfdr;
//    }

//    /**
//     * @param m_PPIfdr the m_PPIfdr to set
//     */
//    public void setPPIfdr(double fdr) {
//        this.m_PPIfdr = fdr;
//        // flag up the peptide pairs
//        for (PeptidePair pp : support) {
//            pp.setPPIFdr(fdr);
//        }
//    }

        /**
     * @param m_PPIfdr the m_PPIfdr to set
     */
    public void setFdrPPI(ProteinGroupPair ppi) {
        this.m_ppi = ppi;
        // flag up the peptide pairs
        for (PeptidePair pp : support) {
            pp.setFdrLink(this);
        }
    }

        /**
     * @param m_PPIfdr the m_PPIfdr to set
     */
    public ProteinGroupPair getFdrPPI() {
        return this.m_ppi;

    }    

    /**
     * @return the linkID
     */
    public int getLinkID() {
        return linkID;
    }

    public int getPeptidePairCount() {
        return support.size();
    }

    protected String getAccessions(HashMap<Protein, IntArrayList> positions) {
        StringBuffer sb = new StringBuffer(); 
        boolean isdecoy = false;
        for (Protein p : positions.keySet()) {
            isdecoy|= p.isDecoy();
            IntArrayList pos = positions.get(p);
            for (Integer i : pos) {
                sb.append(";");
                sb.append(p.getAccession());
            }
        }
        if (isdecoy)
            return "DECOY:" + sb.substring(1);
        return sb.substring(1);
    }

    protected String getDescriptions(HashMap<Protein, IntArrayList> positions) {
        StringBuffer sb = new StringBuffer(); 
        boolean isdecoy = false;
        for (Protein p : positions.keySet()) {
            isdecoy|= p.isDecoy();
            IntArrayList pos = positions.get(p);
            for (Integer i : pos) {
                sb.append(";");
                sb.append(p.getDescription());
            }
        }
        if (isdecoy)
            return "DECOY:" + sb.substring(1);
        return sb.substring(1);
    }
    
    /**
     * @return the m_specialOnly
     */
    public boolean isSpecialOnly() {
        return m_specialOnly;
    }

    /**
     * @return the position1
     */
    public HashMap<Protein, IntArrayList> getPosition1() {
        return position1;
    }

    /**
     * @return the position2
     */
    public HashMap<Protein, IntArrayList> getPosition2() {
        return position2;
    }

    public Object getSite(int n) {
        return n==0 ? getProteinGroup1() : getProteinGroup2();
    }

    public int getSites() {
        return 2;
    }    

    public boolean isLinear() {
        return false;
    }

    public boolean isInternal() {
        return isInternal;
    }

    public boolean isBetween() {
        return !isInternal;
    }
    
}
