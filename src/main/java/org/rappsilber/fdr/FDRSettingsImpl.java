/*
 * Copyright 2016 lfischer.
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


public class FDRSettingsImpl implements FDRSettings {

    OfflineFDR.FDRLevel doOptimize = OfflineFDR.FDRLevel.PROTEINGROUPLINK;
    int MaxLinkAmbiguity = 0;
    int MaxProteinAmbiguity = 0;
    int MinLinkPepCount = 0;
    int MinPPIPepCount = 0;
    int MinPeptideLength = 5;
    int MinProteinPepCount = 0;
    double PSMFDR = 1;
    double PeptidePairFDR = 1;
    double ProteinGroupFDR = 1;
    double ProteinGroupLinkFDR = 0.05;
    double ProteinGroupPairFDR = 1;
    int BoostingSteps = 4;
    FDRSettings.BoostMode boostMode = FDRSettings.BoostMode.NONE;
    FDRSettings.CrosslinkType crosslinkTypeFilter = FDRSettings.CrosslinkType.ALL;
    boolean LinkDirectional;
    boolean PPIDirectional;
    boolean PSMDirectional;
    boolean PeptidePairDirectional;
    double ReportFactor = 1000000;
    boolean filterToUniquePSM = true;
    protected boolean boostPSM = true;
    protected boolean boostPeptidePairs = true;
    protected boolean boostProteins = true;
    protected boolean boostLinks = true;
//    private boolean boostSubScores = false;
    protected boolean groupByPPI = false;
    int minTD = DEFAULT_MIN_TD_COUNT;
    private boolean filterConsectutivePeptides = false;
    private boolean filterBySelfAndMono = false;
    private double subScoreCutOff = 1;
    private boolean boostPepCoverage = true;
    private Boolean psmLocalFDR = null;
    private Boolean peppairLocalFDR = null;
    private Boolean protLocalFDR = null;
    private Boolean linkLocalFDR = null;
    private Boolean ppiLocalFDR = null;
    private double minPeptideCoverageFilter = 0;
    private double minDeltaScoreFilter = 0;
    private double minPeptideStubs = 0;
    private double minDoublets = 0;
    private boolean boostPepDeltaScore =true;
    private boolean combineScoreAndDelta = false;
    private int minFragments =  0;
    private boolean boostMinFragments = true;
    private boolean ignoreValidityChecks = true;
    private boolean groubByCrosslinkerStubs = false;
    private boolean twoStepBoost = true;
    private boolean boostPeptideStubs;
    private boolean boostPeptideDoublets;
    private boolean boostMinScore;
    private double minMinPeptideDoubpletFilter;
    private double minPeptideStubFilter;
    private double minScore = 0d;
    private Integer scoreTopNAggregate = null;

    public FDRSettingsImpl() {
    }
    
    
    public FDRSettingsImpl(FDRSettings settings) {
        transferSettings(settings, this);
    }
    
    
    @Override
    public void addCalcListener(ActionListener al) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public OfflineFDR.FDRLevel doOptimize() {
        return doOptimize;
    }

    @Override
    public void doOptimize(OfflineFDR.FDRLevel level) {
        doOptimize = level;
    }

    @Override
    public int getMaxLinkAmbiguity() {
        return MaxLinkAmbiguity;
    }

    @Override
    public int getMaxProteinAmbiguity() {
        return MaxProteinAmbiguity;
    }

    @Override
    public int getMinLinkPepCount() {
        return MinLinkPepCount;
    }

    @Override
    public int getMinPPIPepCount() {
        return MinPPIPepCount;
    }

    @Override
    public int getMinPeptideLength() {
        return MinPeptideLength;
    }

    @Override
    public int getMinProteinPepCount() {
        return MinProteinPepCount;
    }

    @Override
    public double getPSMFDR() {
        return PSMFDR; 
    }

    @Override
    public double getPeptidePairFDR() {
        return PeptidePairFDR;
    }

    @Override
    public double getProteinGroupFDR() {
        return ProteinGroupFDR;
    }

    @Override
    public double getProteinGroupLinkFDR() {
        return ProteinGroupLinkFDR;
    }

    @Override
    public double getProteinGroupPairFDR() {
        return ProteinGroupPairFDR;
    }

    @Override
    public void setPSMFDR(Double fdr) {
        PSMFDR = fdr;
    }

    @Override
    public void setPeptidePairFDR(Double fdr) {
        PeptidePairFDR=fdr;
    }

    @Override
    public void setProteinGroupFDR(Double fdr) {
        ProteinGroupFDR=fdr;
    }

    @Override
    public void setProteinGroupLinkFDR(Double fdr) {
        ProteinGroupLinkFDR=fdr;
    }

    @Override
    public int getBoostingSteps() {
        return BoostingSteps;
    }

    @Override
    public void setBoostingSteps(int steps) {
        BoostingSteps=steps;
    }

    @Override
    public FDRSettings.BoostMode getBoostMode() {
        return boostMode;
    }

    @Override
    public void setBoostMode(FDRSettings.BoostMode mode) {
        boostMode = mode;
    }

    @Override
    public FDRSettings.CrosslinkType getCrosslinkTypeFilter() {
        return crosslinkTypeFilter;
    }

    @Override
    public void setCrosslinkTypeFilter(FDRSettings.CrosslinkType t) {
        this.crosslinkTypeFilter = t;
    }

    @Override
    public boolean isLinkDirectional() {
        return LinkDirectional;
    }

    @Override
    public boolean isPPIDirectional() {
        return PPIDirectional;
    }

    @Override
    public boolean isPSMDirectional() {
        return PSMDirectional;
    }

    @Override
    public boolean isPeptidePairDirectional() {
        return PeptidePairDirectional;
    }

    @Override
    public void removeCalcListener(ActionListener al) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setLinkDirectional(boolean directional) {
        LinkDirectional = directional;
    }

    @Override
    public void setMaxLinkAmbiguity(Integer maxAmbiguity) {
        MaxLinkAmbiguity = maxAmbiguity;
    }

    @Override
    public void setMaxProteinAmbiguity(Integer maxAmbiguity) {
        MaxProteinAmbiguity = maxAmbiguity;
    }

    @Override
    public void setMinLinkPepCount(Integer minPep) {
        MinLinkPepCount = minPep;
    }

    @Override
    public void setMinPPIPepCount(Integer minPep) {
        MinPPIPepCount = minPep;
    }

    @Override
    public void setMinPeptideLength(Integer minLength) {
        MinPeptideLength = minLength;
    }

    @Override
    public void setMinProteinPepCount(Integer minPep) {
        MinProteinPepCount = minPep;
    }

    @Override
    public void setPPIDirectional(boolean directional) {
        PPIDirectional=directional;
    }

    @Override
    public void setPSMDirectional(boolean directional) {
        PSMDirectional = directional;
    }

    @Override
    public void setPeptidePairDirectional(boolean directional) {
        PeptidePairDirectional = directional;
    }

    @Override
    public void setProteinGroupPairFDR(Double fdr) {
        ProteinGroupPairFDR = fdr;
    }

    @Override
    public double getReportFactor() {
        return ReportFactor;
    }

    @Override
    public void setReportFactor(double factor) {
        ReportFactor =factor;
    }

    @Override
    public boolean filterToUniquePSM() {
        return filterToUniquePSM;
    }

    @Override
    public void setFilterToUniquePSM(boolean filterToUnique) {
        filterToUniquePSM=filterToUnique;
    }

    @Override
    public void setAll(FDRSettings settings) {
        transferSettings(settings, this);
    }
    
    public static void transferSettings(FDRSettings from, FDRSettings to) {
        to.setBoostingSteps(from.getBoostingSteps());
        to.setMaxLinkAmbiguity(from.getMaxLinkAmbiguity());
        to.setMaxProteinAmbiguity(from.getMaxProteinAmbiguity());
        to.setMinLinkPepCount(from.getMinLinkPepCount());
        to.setMinPPIPepCount(from.getMinPPIPepCount());
        to.setMinPeptideLength(from.getMinPeptideLength());
        to.setMinProteinPepCount(from.getMinProteinPepCount());
        to.setMinTD(from.getMinTD());
        to.setPSMFDR(from.getPSMFDR());
        to.setPSMDirectional(from.isPSMDirectional());
        to.setPPIDirectional(from.isPPIDirectional());
        to.setLinkDirectional(from.isLinkDirectional());
        to.setPeptidePairDirectional(from.isPeptidePairDirectional());
        to.setPeptidePairFDR(from.getPeptidePairFDR());
        to.setProteinGroupFDR(from.getProteinGroupFDR());
        to.setProteinGroupLinkFDR(from.getProteinGroupLinkFDR());
        to.setProteinGroupPairFDR(from.getProteinGroupPairFDR());
        to.setReportFactor(from.getReportFactor());
        to.doOptimize(from.doOptimize());
        to.setFilterToUniquePSM(from.filterToUniquePSM());
        to.setBoostMode(from.getBoostMode());
        to.setCrosslinkTypeFilter(from.getCrosslinkTypeFilter());
        to.boostLinks(from.boostLinks());
        to.boostPSMs(from.boostPSMs());
        to.boostPeptidePairs(from.boostPeptidePairs());
        to.boostProteins(from.boostProteins());
        to.boostMinFragments(from.boostMinFragments());
        to.setFilterConsecutivePeptides(from.filterConsecutivePeptides());
        to.psmLocalFDR(from.psmLocalFDR());    
        to.peppairLocalFDR(from.peppairLocalFDR());    
        to.protLocalFDR(from.protLocalFDR());    
        to.linkLocalFDR(from.linkLocalFDR());    
        to.ppiLocalFDR(from.ppiLocalFDR());    
        to.setMinDeltaScoreFilter(from.getMinDeltaScoreFilter());    
        to.setMinPeptideFragmentsFilter(from.getMinPeptideFragmentsFilter()); 
        to.setMinPeptideCoverageFilter(from.getMinPeptideCoverageFilter()); 
        to.setMinPeptideStubFilter(from.getMinPeptideStubFilter());    
        to.setMinPeptideDoubletFilter(from.getMinPeptideDoubletFilter());    
        to.boostPepCoverage(from.boostPepCoverage());
        to.boostDeltaScore(from.boostDeltaScore());
        to.boostPeptideStubs(from.boostPeptideStubs());
        to.boostPeptideDoublets(from.boostPeptideDoublets());
        to.combineScoreAndDelta(from.combineScoreAndDelta());
        to.ignoreValidityChecks(from.ignoreValidityChecks());
        to.setGroupByCrosslinkerStubs(from.getGroupByCrosslinkerStubs());
        to.twoStepOptimization(from.twoStepOptimization());
        to.setfilterBySelfAndMono(from.filterBySelfAndMono());
        to.boostMinScore(from.boostMinScore());
        to.minScore(from.minScore());
        to.setScoreTopNAggregate(from.getScoreTopNAggregate());
    }

    public boolean combineScoreAndDelta() {
        return this.combineScoreAndDelta;
    }
    public void combineScoreAndDelta(boolean c) {
        this.combineScoreAndDelta = c;
    }


    public boolean boostProteins() {
        return boostProteins;
    }

    public void boostProteins(boolean boost) {
        boostProteins = boost;

    }

    
//    public boolean boostSubScores() {
//        return boostSubScores;
//    }
//
//    public void boostSubScores(boolean boost) {
//        boostSubScores = boost;
//    }
    
    public boolean boostPSMs() {
        return boostPSM;
    }

    public void boostPSMs(boolean boost) {
        boostPSM = boost;

    }

    public boolean boostPeptidePairs() {
        return boostPeptidePairs;
    }

    public void boostPeptidePairs(boolean boost) {
        boostPeptidePairs = boost;
    }

    public boolean boostLinks() {
        return boostLinks;
    }

    public void boostLinks(boolean boost) {
        boostLinks = boost;
    }    

    @Override
    public boolean isGroupByPSMCount() {
        return groupByPPI;
    }

    @Override
    public void setGroupByPSMCount(boolean groupByPPI) {
        this.groupByPPI=groupByPPI;
    }

    @Override
    public void setMinTD(Integer c) {
        minTD =c;
    }

    @Override
    public int getMinTD() {
        return minTD;
    }

    @Override
    public boolean filterConsecutivePeptides() {
        return filterConsectutivePeptides;
    }

    @Override
    public void setFilterConsecutivePeptides(boolean filterConsecutive) {
        this.filterConsectutivePeptides=filterConsecutive;
    }
    
    public double getSubScoreCutOff() {
        return this.subScoreCutOff;
    }
    public void setSubScoreCutOff(double localfdr) {
        this.subScoreCutOff = localfdr;
    }
    
    @Override
    public boolean boostDeltaScore(){
        return this.boostPepDeltaScore;
    }

    @Override
    public void boostDeltaScore(boolean boost){
        this.boostPepDeltaScore = boost;
    }
    
    @Override
    public boolean boostPeptideStubs(){
        return this.boostPeptideStubs;
    }

    @Override
    public void boostPeptideStubs(boolean boost){
        this.boostPeptideStubs = boost;
    }
    
    @Override
    public boolean boostMinScore(){
        return this.boostMinScore;
    }

    @Override
    public void boostMinScore(boolean boost){
        this.boostMinScore = boost;
    }
    
    @Override
    public Double minScore(){
        return this.minScore;
    }

    @Override
    public void minScore(Double minScore){
        this.minScore = minScore;
    }

    @Override
    public boolean boostPeptideDoublets(){
        return this.boostPeptideDoublets;
    }

    @Override
    public void boostPeptideDoublets(boolean boost){
        this.boostPeptideDoublets = boost;
    }

    @Override
    public boolean boostPepCoverage(){
        return this.boostPepCoverage;
    }

    @Override
    public void boostPepCoverage(boolean boost){
        this.boostPepCoverage = boost;
    }

    @Override
    public boolean boostMinFragments(){
        return this.boostMinFragments;
    }

    @Override
    public void boostMinFragments(boolean boost){
        this.boostMinFragments = boost;
    }

    @Override
    public Boolean psmLocalFDR() {
        return this.psmLocalFDR;
    }

    @Override
    public Boolean peppairLocalFDR() {
        return this.peppairLocalFDR;
    }

    @Override
    public Boolean protLocalFDR() {
        return this.protLocalFDR;
    }

    @Override
    public Boolean linkLocalFDR() {
        return this.linkLocalFDR;
    }

    @Override
    public Boolean ppiLocalFDR() {
        return this.ppiLocalFDR;
    }

    @Override
    public void psmLocalFDR(Boolean local) {
       this.psmLocalFDR = local;
    }

    @Override
    public void peppairLocalFDR(Boolean local) {
        this.peppairLocalFDR = local;
    }

    @Override
    public void protLocalFDR(Boolean local) {
        this.protLocalFDR = local;
    }

    @Override
    public void linkLocalFDR(Boolean local) {
        this.linkLocalFDR = local;
    }

    @Override
    public void ppiLocalFDR(Boolean local) {
        this.ppiLocalFDR = local;
    }
    
    public double getMinPeptideCoverageFilter() {
        return this.minPeptideCoverageFilter;
    }
    public void setMinPeptideCoverageFilter(double d) {
        this.minPeptideCoverageFilter = d;
    }

    public double getMinDeltaScoreFilter() {
        return this.minDeltaScoreFilter;
    }
    
    public void setMinDeltaScoreFilter(double d) {
        this.minDeltaScoreFilter = d;
    }

    public double getMinPeptideStubFilter() {
        return this.minPeptideStubFilter;
    }
    
    public void setMinPeptideStubFilter(double d) {
        this.minPeptideStubFilter = d;
    }
    
    public double getMinPeptideDoubletFilter() {
        return this.minMinPeptideDoubpletFilter;
    }
    
    public void setMinPeptideDoubletFilter(double d) {
        this.minMinPeptideDoubpletFilter = d;
    }
    
    
    void setMinAbsolutePepCoverFilter(double d) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int getMinPeptideFragmentsFilter() {
        return this.minFragments;
    }

    @Override
    public void setMinPeptideFragmentsFilter(int frags) {
        this.minFragments = frags;
    }

    @Override
    public boolean ignoreValidityChecks() {
        return this.ignoreValidityChecks;
    }

    @Override
    public void ignoreValidityChecks(boolean ignore) {
        this.ignoreValidityChecks = ignore;
    }

    @Override
    public void setGroupByCrosslinkerStubs(boolean group) {
        this.groubByCrosslinkerStubs = group;
    }

    @Override
    public boolean getGroupByCrosslinkerStubs() {
        return this.groubByCrosslinkerStubs;
    }

    @Override
    public boolean twoStepOptimization() {
        return this.twoStepBoost ;
    }

    @Override
    public void twoStepOptimization(boolean stepped) {
        this.twoStepBoost = stepped;
    }

    @Override
    public boolean filterBySelfAndMono() {
        return filterBySelfAndMono;
    }

    @Override
    public void setfilterBySelfAndMono(boolean filter) {
        filterBySelfAndMono = filter;
        boostMinScore(filter);
    }

    @Override
    public Integer getScoreTopNAggregate() {
        return this.scoreTopNAggregate;
    }

    @Override
    public void setScoreTopNAggregate(Integer n) {
        this.scoreTopNAggregate = n;
    }
    
    
}
