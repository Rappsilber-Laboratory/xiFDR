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
package org.rappsilber.fdr.result;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.fdr.entities.PeptidePair;
import org.rappsilber.fdr.entities.ProteinGroupLink;
import org.rappsilber.fdr.entities.ProteinGroupPair;
import org.rappsilber.fdr.entities.ProteinGroup;
import org.rappsilber.fdr.utils.MaximisingStatus;

/**
 *
 * @author lfischer
 */
public  class FDRResult {
    private Collection<PSM> input;
    private FDRResultLevel<PSM> psmFDR;
    private FDRResultLevel<PeptidePair> peptidePairFDR;
    private FDRResultLevel<ProteinGroup> proteinGroupFDR;
    private FDRResultLevel<ProteinGroupLink> proteinGroupLinkFDR;
    private FDRResultLevel<ProteinGroupPair> proteinGroupPairFDR;

    private int minPeptideLength = 0;
    private int maximumProteinAmbiguity = 0;
    private int maximumLinkAmbiguity = 0;
    private double reportFactor = 0;
    private boolean uniquePSMs = false;
    private ArrayList<String> excludedGroups = new ArrayList();

    /** Non-null when this result came from a sequential self/between boost: converged settings for round 1 (Self). */
    private MaximisingStatus selfBoostStatus = null;
    /** Non-null when this result came from a sequential self/between boost: converged settings for round 2 (Between). */
    private MaximisingStatus betweenBoostStatus = null;

    

    /**
     * @return the input
     */
    public Collection<PSM> getInput() {
        return input;
    }

    /**
     * @param input the input to set
     */
    public void setInput(Collection<PSM> input) {
        this.input = input;
    }

    /**
     * @return the psmFDR
     */
    public FDRResultLevel<PSM> getPsmFDR() {
        return psmFDR;
    }

    /**
     * @param psmFDR the psmFDR to set
     */
    public void setPsmFDR(FDRResultLevel<PSM> psmFDR) {
        this.psmFDR = psmFDR;
    }

    /**
     * @return the peptidePairFDR
     */
    public FDRResultLevel<PeptidePair> getPeptidePairFDR() {
        return peptidePairFDR;
    }

    /**
     * @param peptidePairFDR the peptidePairFDR to set
     */
    public void setPeptidePairFDR(FDRResultLevel<PeptidePair> peptidePairFDR) {
        this.peptidePairFDR = peptidePairFDR;
    }

    /**
     * @return the proteinGroupFDR
     */
    public FDRResultLevel<ProteinGroup> getProteinGroupFDR() {
        return proteinGroupFDR;
    }

    /**
     * @param proteinGroupFDR the proteinGroupFDR to set
     */
    public void setProteinGroupFDR(FDRResultLevel<ProteinGroup> proteinGroupFDR) {
        this.proteinGroupFDR = proteinGroupFDR;
    }

    /**
     * @return the proteinGroupLinkFDR
     */
    public FDRResultLevel<ProteinGroupLink> getProteinGroupLinkFDR() {
        return proteinGroupLinkFDR;
    }

    /**
     * @param proteinGroupLinkFDR the proteinGroupLinkFDR to set
     */
    public void setProteinGroupLinkFDR(FDRResultLevel<ProteinGroupLink> proteinGroupLinkFDR) {
        this.proteinGroupLinkFDR = proteinGroupLinkFDR;
    }

    /**
     * @return the proteinGroupPairFDR
     */
    public FDRResultLevel<ProteinGroupPair> getProteinGroupPairFDR() {
        return proteinGroupPairFDR;
    }

    /**
     * @param proteinGroupPairFDR the proteinGroupPairFDR to set
     */
    public void setProteinGroupPairFDR(FDRResultLevel<ProteinGroupPair> proteinGroupPairFDR) {
        this.proteinGroupPairFDR = proteinGroupPairFDR;
    }

    /**
     * @return the minPeptideLength
     */
    public int getMinPeptideLength() {
        return minPeptideLength;
    }

    /**
     * @param minPeptideLength the minPeptideLength to set
     */
    public void setMinPeptideLength(int minPeptideLength) {
        this.minPeptideLength = minPeptideLength;
    }

    /**
     * @return the maximumProteinAmbiguity
     */
    public int getMaximumProteinAmbiguity() {
        return maximumProteinAmbiguity;
    }

    /**
     * @param maximumProteinAmbiguity the maximumProteinAmbiguity to set
     */
    public void setMaximumProteinAmbiguity(int maximumProteinAmbiguity) {
        this.maximumProteinAmbiguity = maximumProteinAmbiguity;
    }

    /**
     * @return the maximumLinkAmbiguity
     */
    public int getMaximumLinkAmbiguity() {
        return maximumLinkAmbiguity;
    }

    /**
     * @param maximumLinkAmbiguity the maximumLinkAmbiguity to set
     */
    public void setMaximumLinkAmbiguity(int maximumLinkAmbiguity) {
        this.maximumLinkAmbiguity = maximumLinkAmbiguity;
    }

    /**
     * @return the reportFactor
     */
    public double getReportFactor() {
        return reportFactor;
    }

    /**
     * @param reportFactor the reportFactor to set
     */
    public void setReportFactor(double reportFactor) {
        this.reportFactor = reportFactor;
    }

    /**
     * @return the uniquePSMs
     */
    public boolean isUniquePSMs() {
        return uniquePSMs;
    }

    /**
     * @param uniquePSMs the uniquePSMs to set
     */
    public void setUniquePSMs(boolean uniquePSMs) {
        this.uniquePSMs = uniquePSMs;
    }

    /**
     * @return the excludedGroups
     */
    public ArrayList<String> getExcludedGroups() {
        return excludedGroups;
    }

    /**
     * @param excludedGroups the excludedGroups to set
     */
    public void setExcludedGroups(ArrayList<String> excludedGroups) {
        this.excludedGroups = excludedGroups;
    }

    /**
     * @return the selfBoostStatus
     */
    public MaximisingStatus getSelfBoostStatus() {
        return selfBoostStatus;
    }

    /**
     * @param selfBoostStatus the selfBoostStatus to set
     */
    public void setSelfBoostStatus(MaximisingStatus selfBoostStatus) {
        this.selfBoostStatus = selfBoostStatus;
    }

    /**
     * @return the betweenBoostStatus
     */
    public MaximisingStatus getBetweenBoostStatus() {
        return betweenBoostStatus;
    }

    /**
     * @param betweenBoostStatus the betweenBoostStatus to set
     */
    public void setBetweenBoostStatus(MaximisingStatus betweenBoostStatus) {
        this.betweenBoostStatus = betweenBoostStatus;
    }
    
    public void clear() {
        if (psmFDR != null)
            psmFDR.clear();
        if (peptidePairFDR != null)
            peptidePairFDR.clear();
        if (proteinGroupFDR != null)
            proteinGroupFDR.clear();
        if (proteinGroupLinkFDR != null)
            proteinGroupLinkFDR.clear();
        if (proteinGroupPairFDR != null)
            proteinGroupPairFDR.clear();
    }
    

    
}
