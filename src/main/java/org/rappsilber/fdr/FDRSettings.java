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

import java.awt.event.ActionListener;

/**
 *
 * @author lfischer
 */
public interface FDRSettings {
    public static int DEFAULT_MIN_TD_COUNT=2;

    public enum BoostMode {
        NONE("No Boosting", "Don't try to boost"),
        ALL("All", "Boost all crosslinked items as one list"),
        BETWEEN("Between", "Boost to increase protein heteromeric matches"),
        SELFBETWEEN("Self+Between", "Bost first self links and then protein heteromeric matches");

        private final String name;
        private final String description;

        BoostMode(String name, String description) {
            this.name = name;
            this.description = description;
        }


        public String getName() {
            return name;
        }

        public String getDescription() {
            return description;
        }

        @Override
        public String toString() {
            return name;
        }
        
        
        
    }
    
    /**
     * should the boosting be done in two steps (lower-FDR and non-FDR filters).
     * @return 
     */
    boolean twoStepOptimization();
    void twoStepOptimization(boolean stepped);
    
    void addCalcListener(ActionListener al);

    OfflineFDR.FDRLevel doOptimize();
    void doOptimize(OfflineFDR.FDRLevel level);

    int getMaxLinkAmbiguity();


    int getMaxProteinAmbiguity();

    int getMinLinkPepCount();

    int getMinPPIPepCount();

    int getMinPeptideLength();

    int getMinProteinPepCount();

    double getPSMFDR();

    double getPeptidePairFDR();

    double getProteinGroupFDR();

    double getProteinGroupLinkFDR();

    double getProteinGroupPairFDR();

    void setPSMFDR(Double fdr);
    
    void setPeptidePairFDR(Double fdr);

    void setProteinGroupFDR(Double fdr);

    void setProteinGroupLinkFDR(Double fdr);
    
    
    int getBoostingSteps();
    public void setBoostingSteps(int steps);

    BoostMode getBoostMode();
    void setBoostMode(BoostMode mode);

    /** @deprecated Use {@link #getBoostMode()} instead. */
    default boolean getBoostBetween() { return getBoostMode() == BoostMode.BETWEEN; }
    /** @deprecated Use {@link #setBoostMode(BoostMode)} instead. */
    default void setBoostBetween(boolean between) { if (between) setBoostMode(BoostMode.BETWEEN); }
    /** @deprecated Use {@link #getBoostMode()} instead. */
    default boolean getBoostSelfAndBetween() { return getBoostMode() == BoostMode.SELFBETWEEN; }
    /** @deprecated Use {@link #setBoostMode(BoostMode)} instead. */
    default void setBoostSelfAndBetween(boolean selfAndBetween) { if (selfAndBetween) setBoostMode(BoostMode.SELFBETWEEN); }

    enum CrosslinkType { ALL, SELF, BETWEEN }

    default CrosslinkType getCrosslinkTypeFilter() { return CrosslinkType.ALL; }
    default void setCrosslinkTypeFilter(CrosslinkType t) { }

    boolean isLinkDirectional();

    boolean isPPIDirectional();

    boolean isPSMDirectional();

    boolean isPeptidePairDirectional();

    void removeCalcListener(ActionListener al);

    void setLinkDirectional(boolean directional);

    void setMaxLinkAmbiguity(Integer maxAmbiguity);


    void setMaxProteinAmbiguity(Integer maxAmbiguity);

    void setMinLinkPepCount(Integer minPep);

    void setMinPPIPepCount(Integer minPep);

    void setMinPeptideLength(Integer minLength);

    void setMinProteinPepCount(Integer minPep);

    void setPPIDirectional(boolean directional);

    void setPSMDirectional(boolean directional);

    void setPeptidePairDirectional(boolean directional);


    void setProteinGroupPairFDR(Double fdr);
    
    double getReportFactor();
    
    public void setReportFactor(double factor);       
    
    boolean filterToUniquePSM();
    
    void setFilterToUniquePSM(boolean filterToUnique);
    
    boolean filterConsecutivePeptides();
    
    void setFilterConsecutivePeptides(boolean filterConsecutive);

    boolean filterBySelfAndMono();
    
    void setfilterBySelfAndMono(boolean filter);
    
    public void setAll(FDRSettings settings);
    
    public boolean boostProteins();
    public void boostProteins(boolean boost);

    public boolean boostPepCoverage();    
    public void boostPepCoverage(boolean boost);    

    public boolean boostMinFragments();
    public void boostMinFragments(boolean boost);
    
    public boolean boostDeltaScore();
    public void boostDeltaScore(boolean boost);

    public boolean boostPeptideStubs();
    public void boostPeptideStubs(boolean boost);

    public boolean boostPeptideDoublets();
    public void boostPeptideDoublets(boolean boost);

    public boolean boostMinScore();
    public void boostMinScore(boolean boost);

    public Double minScore();
    public void minScore(Double minScore);
    
    public boolean boostPSMs();
    public void boostPSMs(boolean boost);
//    public boolean boostSubScores();
//    public void boostSubScores(boolean boost);
    public boolean boostPeptidePairs();
    public void boostPeptidePairs(boolean boost);
    public boolean boostLinks();

    public void boostLinks(boolean boost);
    public boolean isGroupByPSMCount();
    public void setGroupByPSMCount(boolean groupByPPI);
    
    public double getMinPeptideCoverageFilter();
    public void setMinPeptideCoverageFilter(double d);
    
    public int getMinPeptideFragmentsFilter();
    public void setMinPeptideFragmentsFilter(int frags);

    public double getMinDeltaScoreFilter();
    public void setMinDeltaScoreFilter(double d);
    
    public double getMinPeptideStubFilter();    
    public void setMinPeptideStubFilter(double d);
    
    public double getMinPeptideDoubletFilter();
    public void setMinPeptideDoubletFilter(double d);
    
    public void setMinTD(Integer c) ;

    public int getMinTD();

    public Boolean psmLocalFDR();
    public Boolean peppairLocalFDR();
    public Boolean protLocalFDR();
    public Boolean linkLocalFDR();
    public Boolean ppiLocalFDR();
    public void psmLocalFDR(Boolean local);
    public void peppairLocalFDR(Boolean local);
    public void protLocalFDR(Boolean local);
    public void linkLocalFDR(Boolean local);
    public void ppiLocalFDR(Boolean local);
    public boolean combineScoreAndDelta();
    public void combineScoreAndDelta(boolean c);
    public boolean ignoreValidityChecks();
    public void ignoreValidityChecks(boolean ignore);
//    public double getSubScoreCutOff();
//    public void setSubScoreCutOff(double localfdr);
    public void setGroupByCrosslinkerStubs(boolean group);
    public boolean getGroupByCrosslinkerStubs();
    
    public Integer getScoreTopNAggregate();
    public void setScoreTopNAggregate(Integer n);
    
}
