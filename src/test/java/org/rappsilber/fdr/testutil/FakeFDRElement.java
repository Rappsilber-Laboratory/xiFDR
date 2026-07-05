/*
 * Copyright 2026 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
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
package org.rappsilber.fdr.testutil;

import java.util.Collection;
import java.util.Collections;
import org.rappsilber.fdr.entities.AbstractFDRElement;
import org.rappsilber.fdr.entities.PeptidePair;
import org.rappsilber.fdr.entities.ProteinGroup;
import org.rappsilber.fdr.entities.Site;

/**
 * Minimal {@link AbstractFDRElement} test double for unit-testing
 * {@link org.rappsilber.fdr.calculation.FDRImplement} and
 * {@link org.rappsilber.fdr.utils.MaximizeLevelInfo} without needing real
 * PSM/PeptidePair/ProteinGroup objects.
 */
public class FakeFDRElement extends AbstractFDRElement<FakeFDRElement> {

    public enum MatchType {TT, TD, DD}

    private double score;
    private MatchType type = MatchType.TT;
    private boolean linear;
    private boolean internal;
    private boolean between;
    private int peptidePairCount = 1;
    private String fdrGroup;
    private double fdr;

    public FakeFDRElement score(double score) {
        this.score = score;
        return this;
    }

    public FakeFDRElement type(MatchType type) {
        this.type = type;
        return this;
    }

    public FakeFDRElement linear(boolean linear) {
        this.linear = linear;
        return this;
    }

    public FakeFDRElement internal(boolean internal) {
        this.internal = internal;
        return this;
    }

    public FakeFDRElement between(boolean between) {
        this.between = between;
        return this;
    }

    public FakeFDRElement peptidePairCount(int peptidePairCount) {
        this.peptidePairCount = peptidePairCount;
        return this;
    }

    public FakeFDRElement fdrGroup(String fdrGroup) {
        this.fdrGroup = fdrGroup;
        return this;
    }

    @Override
    public double getScore(int topN) {
        return score;
    }

    @Override
    public double getScore() {
        return score;
    }

    @Override
    public String getFDRGroup() {
        return fdrGroup;
    }

    @Override
    public void setFDRGroup(String fdrGroup) {
        this.fdrGroup = fdrGroup;
    }

    @Override
    public boolean isTT() {
        return type == MatchType.TT;
    }

    @Override
    public boolean isTD() {
        return type == MatchType.TD;
    }

    @Override
    public boolean isDD() {
        return type == MatchType.DD;
    }

    @Override
    public boolean isDecoy() {
        return type != MatchType.TT;
    }

    @Override
    public boolean isLinear() {
        return linear;
    }

    @Override
    public boolean isInternal() {
        return internal;
    }

    @Override
    public boolean isBetween() {
        return between;
    }

    @Override
    public void setFDR(double fdr) {
        this.fdr = fdr;
    }

    @Override
    public double getFDR() {
        return fdr;
    }

    @Override
    public int getPeptidePairCount() {
        return peptidePairCount;
    }

    @Override
    public Object getSite(int n) {
        return null;
    }

    @Override
    public int getSites() {
        return 0;
    }

    @Override
    public boolean isNonCovalent() {
        return false;
    }

    @Override
    public void add(FakeFDRElement o) {
    }

    @Override
    public Site getLinkSite1() {
        return null;
    }

    @Override
    public Site getLinkSite2() {
        return null;
    }

    @Override
    public ProteinGroup getProteinGroup1() {
        return null;
    }

    @Override
    public ProteinGroup getProteinGroup2() {
        return null;
    }

    @Override
    public Collection<PeptidePair> getPeptidePairs() {
        return Collections.emptyList();
    }
}
