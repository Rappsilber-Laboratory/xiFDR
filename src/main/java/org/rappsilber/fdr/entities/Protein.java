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
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.rappsilber.fdr.utils.FDRGroupNames;
import org.rappsilber.utils.DoubleArrayList;

/**
 * Represents a single protein. 
 * @author lfischer
 */
public class Protein extends AbstractFDRElement<Protein> {//implements Comparable<Protein> {
    private long id;
    private String searchId;
    private String accession;
    private String name;
    private String description;
    private String sequence="";
    //private HashSet<PeptidePair> peps = new HashSet<PeptidePair>();
    private HashSet<Peptide> peps = new HashSet<Peptide>();
    private HashSet<PeptidePair> peppairs = new HashSet<PeptidePair>();
    private double score = 0;
    private DoubleArrayList scores = new DoubleArrayList(2);
    private boolean scoresSorted = true;
    private boolean isDecoy;
    private boolean linearSupport = false;
    private boolean internalSupport = false;
    private boolean betweenSupport = false;
    private boolean crosslinkerModSupport = false;
    //private boolean isAmbigious = true;
    private double m_fdr = -1;
    private int size = -1;
    private Protein decoy_complement = null;
    
    public static String DECOY_PREFIX = null;
    public static Pattern DEFAULT_DECOY_PATTERN = Pattern.compile("^(REV|RAN|RANDOM|SHUFFLE|REVERSE)[:_]", Pattern.CASE_INSENSITIVE);
    String zero = ""+(char)0;
    private Pattern zerosplit = Pattern.compile(zero);
    private Pattern spacesplit = Pattern.compile("^\\s*(\"[^\"]*\"|'[^']*'|\\([^)]*\\)|[^\\s]+)");
    
    
    private String fdrgroup = null;
    
    /**
     * Used as the protein of origin for the NOPEPTIDE peptide
     */
    public static Protein NOPROTEIN = new Protein(Integer.MAX_VALUE, "", "", false,true,true, true);
    private String validated;

    /**
     * Constructor
     * @param id
     * @param accession
     * @param description
     * @param isDecoy
     * @param linear
     * @param internal
     * @param between 
     */
    public Protein(long id, String accession, String name, String description, boolean isDecoy, boolean linear, boolean internal, boolean between, String decoyprefix) {
        this(id, accession, description, isDecoy, linear, internal, between, decoyprefix);
        this.name = name;
    }

    /**
     * Constructor
     * @param id
     * @param accession
     * @param description
     * @param isDecoy
     * @param linear
     * @param internal
     * @param between 
     */
    public Protein(long id, String accession, String name, String description, boolean isDecoy, boolean linear, boolean internal, boolean between) {
        this(id, accession, name, description, isDecoy, linear, internal, between, DECOY_PREFIX);
    }    

    /**
     * Constructor
     * @param id
     * @param accession
     * @param description
     * @param isDecoy
     * @param linear
     * @param internal
     * @param between 
     */
    public Protein(long id, String accession, String description, boolean isDecoy, boolean linear, boolean internal, boolean between) {
        this(id, accession, description, isDecoy, linear, internal, between, DECOY_PREFIX);
    }    
    /**
     * Constructor
     * @param id
     * @param accession
     * @param description
     * @param isDecoy
     * @param linear
     * @param internal
     * @param between 
     */
    public Protein(long id, String accession, String description, boolean isDecoy, boolean linear, boolean internal, boolean between, String decoyPrefix) {
        this.id = id;
        setAccession(accession, decoyPrefix);
        
        this.description = description;
        this.isDecoy = isDecoy;
        this.linearSupport = linear;
        this.internalSupport = internal;
        this.betweenSupport = between;
        
        if (description != null && description.length() > 0) {
            if (description.indexOf(zero)>=0) {
                this.name = zerosplit.split(description)[1];
                this.description = zerosplit.split(description)[0];
            } else {
                // try parsing a name from description
                Matcher m = spacesplit.matcher(description);
                if (m.find()) {
                    this.name = m.group(1);
                } else {
                    this.name = this.description;
                }
            }
        }
        
    }
    
     
    @Override
    public int hashCode() {
        return accession.hashCode();
    }

    @Override
    public boolean equals(Object p) {
        return ((Protein) p).accession.contentEquals(accession) && 
                isDecoy == ((Protein) p).isDecoy && 
                ( 
                    (!isDecoy) ||  // unluckily we have different ways to do decoy at times - so we have to ignore the sequence here
                    sequence.contentEquals(((Protein) p).sequence) || 
                    sequence.isEmpty() || 
                    ((Protein) p).sequence.isEmpty()
                ); // && ((Protein) p).sequence.contentEquals(sequence);
    }

    /**
     * Is p the same protein as this.
     * Decoy is ignored
     * @param p
     * @return 
     */
    public boolean equalsDecoysUnaware(Protein p) {
        return (p.accession.toLowerCase().contentEquals(accession.toLowerCase()) || p.accession.toLowerCase().contentEquals("rev_" + accession.toLowerCase()) ||
            accession.toLowerCase().contentEquals("rev_" + p.accession.toLowerCase())) && (p.sequence.length() == sequence.length() || p.sequence.isEmpty() || sequence.isEmpty()); // && ((Protein) p).sequence.contentEquals(sequence);
    }

    /**
     * Add the information (supporting peptides) from the given protein instance
     * to this instance
     * @param o 
     */
    @Override
    public void add(Protein o) {
        if (o == this)
            return;
        this.score = Math.sqrt(this.score*this.score + o.score*o.score);
        scores.add(o.score);
        scoresSorted = false;
        this.scoresSorted = false;
        peps.addAll(o.peps);
        linearSupport |= o.linearSupport;
        internalSupport |= o.internalSupport;
        this.betweenSupport |= o.betweenSupport;
        if (this.description == null || this.description.length()==0)
            this.description = o.description;
        //setFDRGroup();
       //isAmbigious &= o.isAmbigious;        
        //            score+= o.score;
//        if (o.hasPositiveGrouping()) {
//            if (this.validated == null ) {
//                this.validated = o.getPositiveGrouping();
//            } else if (!this.validated.contentEquals(o.getPositiveGrouping())) {
//                this.validated += " " + o.getPositiveGrouping();
//            }
//        }
    }

    /**
     * Adds the peptide pair as support for this protein.
     * @param pp
     * @param score
     */
    public void add(PeptidePair pp, double score) {
        //            double s = pp.score/(pp.peptide1Positions.size() * pp.peptide1Positions.size());
        //            score += s*s;
        peppairs.add(pp);
        
        this.score = Math.sqrt(this.score*this.score + score*score);
        scoresSorted = false;
        this.internalSupport |= pp.isInternal();
        this.betweenSupport |= !(pp.isInternal() || pp.isLinear);
        this.linearSupport |= pp.isLinear;
        //setFDRGroup();
//        if (pp.hasPositiveGrouping()) {
//            if (this.validated == null) {
//                this.validated = pp.getPositiveGrouping();
//            } else if (!this.validated.contentEquals(pp.getPositiveGrouping())) {
//                this.validated += " " + pp.getPositiveGrouping();
//            }
//        }
    }


    @Override
    public double getScore() {
        return score;
    }

    @Override
    public double getScore(int topN) {
        if (!scoresSorted) {
            java.util.Collections.sort(scores);
            scoresSorted = true;
        }
        
        int i = 0;
        double score = 0;
        double lastScore = Double.NaN;
        for (double s : this.scores) {
            if (lastScore != s) {
                lastScore = s;
                if (i++>topN)
                    break;
            }
            score+= s;
        }
        return score;
    }    

    public int compare(Protein p) {
        if (p.equals(this)) {
            return 0;
        }
        if (this.id < p.id) {
            return -1;
        }
        return 1;
    }

    public int compareDecoyUnAware(Protein p) {
        if (p.equalsDecoysUnaware(this)) {
            return 0;
        }
        if (this.id < p.id) {
            return -1;
        }
        return 1;
    }

    /**
     * @return the id
     */
    public long getId() {
        return id;
    }

    /**
     * @return the search id
     */
    public String getSearchId() {
        return searchId;
    }
    
    /**
     * @return the accession
     */
    public String getAccession() {
        return accession;
    }

    /**
     * @return the description
     */
    public String getDescription() {
        return description;
    }

    /**
     * @return the description
     */
    public String getName() {
        return this.name;
    }
    
    Pattern sptr_regex = Pattern.compile("^(?:sp|tr)\\|([^\\|]*)\\|([^ ]*)?.*");
    
    /**
     * @return the accession
     */
    public void setAccession(String a, String decoy_prefix) {
        accession = a;
        if (decoy_prefix == null) {
            Matcher ddm = DEFAULT_DECOY_PATTERN.matcher(accession);
            if (ddm.matches()) {
                accession = accession.substring(ddm.end());
                this.isDecoy = true;
            }
        } else {
            if (accession.startsWith(decoy_prefix)) {
                this.isDecoy = true;
                accession = accession.substring(decoy_prefix.length());
            }
                
        }
        Matcher m = sptr_regex.matcher(a);
        if (m.matches()) {
            accession = m.group(1);
            //if (m.groupCount()>1 && m.group(2) != null && m.group(2).trim().length()>0)
            //    this.setName(m.group(2));
        }

        if (decoy_prefix == null) {
            Matcher ddm = DEFAULT_DECOY_PATTERN.matcher(accession);
            if (ddm.matches()) {
                accession = accession.substring(ddm.end());
                this.isDecoy = true;
            }
        } else {
            if (accession.startsWith(decoy_prefix)) {
                this.isDecoy = true;
                accession = accession.substring(decoy_prefix.length());
            }
                
        }
        this.decoy_complement = null;        
    }

    /**
     * @return the accession
     */
    public void setAccession(String a) {
        setAccession(a,null);
    }
    
    /**
     * @return the accession
     */
    public void setName(String n) {
        name = n;
        this.decoy_complement = null;

    }
    
    /**
     * @return the description
     */
    public void setDescription(String d) {
        description = d;
        if (description != null && description.length() > 0) {
            if (description.indexOf(zero)>=0) {
                this.name = zerosplit.split(description)[1];
                this.description = zerosplit.split(description)[0];
            }
        }
        this.decoy_complement = null;
    }

    /**
     * @return the peps
     */
    public HashSet<PeptidePair> getPeptidePairs() {
        return peppairs;
    }

    /**
     * @return the isDecoy
     */
    @Override
    public boolean isDecoy() {
        return isDecoy;
    }

    /**
     * @return the linearSupport
     */
    public boolean hasLinearSupport() {
        return linearSupport;
    }

    /**
     * @return the internalSupport
     */
    public boolean hasInternalSupport() {
        return internalSupport;
    }

    public int compareTo(Protein o) {
        return Double.compare(o.getScore(), this.getScore());
    }

    private void setFDRGroup() {
        fdrgroup = "";
        if (hasLinearSupport()) {
            fdrgroup = "Linear";
        } 
        if (hasInternalSupport()) {
            fdrgroup += "Internal";
        } 
        if (fdrgroup.isEmpty()){
            fdrgroup = "Between";
        }
        fdrgroup=FDRGroupNames.get(fdrgroup);
        
    }
    
    @Override
    public String getFDRGroup() {
        if (fdrgroup == null)
            setFDRGroup();
        return fdrgroup;
    }
    
    
//    @Override
//    public String getFDRGroupName() {
//        return getFDRGroupName(fdrgroup);
//    }
//    
//    public static String getFDRGroupName(int group) {
//        switch (group) {
//            case -1 : return "all combined";
//            case 0 : return "Only PPI";
//            case 1 : return "Linear"    ;
//            case 2 : return "Within";
//            case 3 : return "Linear + Within";
//            case 4 : return "between only";
//                
//            default : return "?no support?";
//        }
//    }
    
    @Override
    public boolean isTT() {
        return !isDecoy();
    }

    @Override
    public boolean isTD() {
        return isDecoy();
    }

    @Override
    public boolean isDD() {
        return false;
    }



    public Collection<Protein> getProteins() {
        ArrayList<Protein> ret = new ArrayList<Protein>();
        ret.add(this);
        return ret;
    }

    public Collection<Peptide> getPeptides() {
        HashSet<Peptide> ret = new HashSet<Peptide>();
        for (PeptidePair pp : peppairs) {
            ret.addAll(pp.getPeptides());
        }
        return ret;
    }
    
    @Override
    public String toString() {
        
        if (isDecoy)
            return "D - " + accession;
        return "T - " + accession;
        
    }    
    
    public Protein decoyComplement() {
        if (this.decoy_complement == null) {
            this.decoy_complement = new Protein(id, accession, name, description, !isDecoy, linearSupport, internalSupport, betweenSupport);
        }
        return this.decoy_complement;
    }
    
    public int support() {
        return peppairs.size();
    }
    
    @Override
    public void setFDR(double fdr) {
        m_fdr = fdr;
    }

    @Override
    public double getFDR() {
        return m_fdr;
    }

    public void resetFDR() {
        score = 0;
        peps = new HashSet<Peptide>();
        peppairs = new HashSet<PeptidePair>();
        m_fdr = Double.MAX_VALUE;
        setLinkedSupport(1);
    }
    

    
    @Override
    public int getPeptidePairCount() {
        return peppairs.size();
    }

    /**
     * @return the betweenSupport
     */
    public boolean hasBetweenSupport() {
        return betweenSupport;
    }
    
    /**
     * @param betweenSupport the betweenSupport to set
     */
    public void setBetweenSupport(boolean betweenSupport) {
        this.betweenSupport = betweenSupport;
    }
    
    @Override
    public Object getSite(int n) {
        return this;
    }

    @Override
    public int getSites() {
        return 1;
    }    

    /**
     * @return the sequence
     */
    public String getSequence() {
        return sequence;
    }

    /**
     * @param sequence the sequence to set
     */
    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    @Override
    public boolean isLinear() {
        return true;
    }

    @Override
    public boolean isInternal() {
        return false;
    }

    @Override
    public boolean isBetween() {
        return false;
    }
    
    @Override
    public Site getLinkSite1() {
        return new ProteinSite(this);
    }

    @Override
    public Site getLinkSite2() {
        return null;
    }

    @Override
    public ProteinGroup getProteinGroup1() {
        ArrayList<Protein> pg = new ArrayList<>(1);
        pg.add(this);
        return new ProteinGroup(pg, peppairs);
    }

    @Override
    public ProteinGroup getProteinGroup2() {
        return ProteinGroup.NOPROTEINGROUP;
    }

    @Override
    public void setFDRGroup(String fdrGroup) {
        this.fdrgroup = FDRGroupNames.get(fdrGroup);
    }
//    
//    @Override
//    public boolean hasPositiveGrouping() {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }
//    
//    @Override
//    public void setPositiveGrouping(String av) {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }
//
//    @Override
//    public HashSet<String> getPositiveGrouping() {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }
//
//    @Override
//    public boolean hasNegativeGrouping() {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }
//
//    @Override
//    public void setNegativeGrouping(String v) {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }
//
//    @Override
//    public HashSet<String> getNegativeGrouping() {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }
//        

    /**
     * @return the size
     */
    public int getSize() {
        if (size < 0) {
            this.size = this.sequence.replaceAll("[^A-Z]", "").length();
        }
        return size;
    }

    /**
     * @param size the size to set
     */
    public void setSize(int size) {
        this.size = size;
    }
    
    public boolean isNonCovalent() {
        return false;
    }
    
}
