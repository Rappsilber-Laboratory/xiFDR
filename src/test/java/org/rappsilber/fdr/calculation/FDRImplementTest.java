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
package org.rappsilber.fdr.calculation;

import java.util.ArrayList;
import org.junit.jupiter.api.Test;
import org.rappsilber.fdr.FDRSettingsImpl;
import org.rappsilber.fdr.result.FDRResultLevel;
import org.rappsilber.fdr.result.SubGroupFdrInfo;
import org.rappsilber.fdr.testutil.FakeFDRElement;
import org.rappsilber.fdr.testutil.FakeFDRElement.MatchType;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class FDRImplementTest {

    private static final double DELTA = 1e-12;

    private FDRImplement newFdrImplement() {
        return new FDRImplement(new ValidityCheckImplement(0, 2));
    }

    /**
     * 10 elements, scores 10 down to 1 (descending), ranks 1-8 = TT, ranks 9 and 10
     * have the given types.
     */
    private ArrayList<FakeFDRElement> tenElementGroup(MatchType rank9, MatchType rank10) {
        ArrayList<FakeFDRElement> group = new ArrayList<>();
        for (int score = 10; score >= 3; score--) {
            group.add(new FakeFDRElement().score(score).type(MatchType.TT).internal(true));
        }
        group.add(new FakeFDRElement().score(2).type(rank9).internal(true));
        group.add(new FakeFDRElement().score(1).type(rank10).internal(true));
        return group;
    }

    @Test
    void subFDR_basicCutoff_landsBetweenRank9AndRank10() {
        ArrayList<FakeFDRElement> group = tenElementGroup(MatchType.TD, MatchType.TD);
        ArrayList<FakeFDRElement> results = new ArrayList<>();
        SubGroupFdrInfo<FakeFDRElement> info = new SubGroupFdrInfo<>();
        info.TT = 8;
        info.TD = 2;
        info.DD = 0;
        info.targetFDR = 0.15;
        info.saftyfactor = 1_000_000;

        double returned = newFdrImplement().subFDR(group, results, true, info);

        // top-9 (excluding rank 10): TT=8,TD=1,DD=0 -> FDR=1/8=0.125 (passes)
        // top-10 (all):               TT=8,TD=2,DD=0 -> FDR=2/8=0.25 (fails)
        assertEquals(0.25, returned, DELTA);
        assertEquals(8, info.resultCount);
        assertEquals(9, results.size());
        assertEquals(8, info.resultTT);
        assertEquals(1, info.resultTD);
        assertEquals(0, info.resultDD);
        assertEquals(0.125, info.lowerFDR, DELTA);
        assertEquals(0.25, info.higherFDR, DELTA);
        assertEquals(8, info.within);

        FakeFDRElement rank1 = group.get(0);
        FakeFDRElement rank9 = group.get(8);
        FakeFDRElement rank10 = group.get(9);

        assertEquals(0.125, rank9.getFDR(), DELTA);
        for (int i = 0; i <= 7; i++) {
            assertEquals(0.0, group.get(i).getFDR(), DELTA);
        }

        // rank 9 (the lowest-scored kept TD) sees rank 10 as its "higher FDR" neighbour
        assertSame(rank10, rank9.getHigherFDRTD());
        // ranks 1-8 see rank 9 as their "higher FDR" neighbour
        assertSame(rank9, rank1.getHigherFDRTD());
    }

    @Test
    void subFDR_ddInsteadOfTd_lowersEstimatedFdrAndKeepsAllElements() {
        // Same 10-element layout as above, but rank 10 is a DD instead of a second TD.
        ArrayList<FakeFDRElement> group = tenElementGroup(MatchType.TD, MatchType.DD);
        ArrayList<FakeFDRElement> results = new ArrayList<>();
        SubGroupFdrInfo<FakeFDRElement> info = new SubGroupFdrInfo<>();
        info.TT = 8;
        info.TD = 1;
        info.DD = 1;
        info.targetFDR = 0.15;
        info.saftyfactor = 1_000_000;

        double returned = newFdrImplement().subFDR(group, results, true, info);

        // top-10 (all): TT=8,TD=1,DD=1 -> FDR=(1-1)/8=0 (passes), one more element than the TD-only case
        assertEquals(0.125, returned, DELTA);
        assertEquals(9, info.resultCount);
        assertEquals(10, results.size());
        assertEquals(8, info.resultTT);
        assertEquals(1, info.resultTD);
        assertEquals(1, info.resultDD);
        assertEquals(0.0, info.lowerFDR, DELTA);
        assertEquals(0.125, info.higherFDR, DELTA);

        for (FakeFDRElement e : group) {
            assertEquals(0.0, e.getFDR(), DELTA);
        }
    }

    @Test
    void subFDR_safetyFactorSkipsFirstPassingStep() {
        // 4 TT + 1 TD. Naive cutoff at the TD (efdr=0.25<=0.3) would pass, but
        // efdr_p/fdr = 0.5/0.3 >= saftyfactor(1.0), so that step is skipped.
        ArrayList<FakeFDRElement> group = new ArrayList<>();
        for (int score = 5; score >= 2; score--) {
            group.add(new FakeFDRElement().score(score).type(MatchType.TT).internal(true));
        }
        group.add(new FakeFDRElement().score(1).type(MatchType.TD).internal(true));

        ArrayList<FakeFDRElement> results = new ArrayList<>();
        SubGroupFdrInfo<FakeFDRElement> info = new SubGroupFdrInfo<>();
        info.TT = 4;
        info.TD = 1;
        info.DD = 0;
        info.targetFDR = 0.3;
        info.saftyfactor = 1.0;

        double returned = newFdrImplement().subFDR(group, results, true, info);

        // the skipped step (at the TD) is recorded as the first passing FDR ...
        assertEquals(0.5, info.firstPassingFDR, DELTA);
        // ... but the search continues and finds the all-TT cutoff instead
        assertEquals(0.25, returned, DELTA);
        assertEquals(3, info.resultCount);
        assertEquals(4, results.size());
        assertEquals(4, info.resultTT);
        assertEquals(0, info.resultTD);
        assertEquals(0, info.resultDD);
        assertEquals(0.0, info.lowerFDR, DELTA);
        assertEquals(0.25, info.higherFDR, DELTA);
    }

    @Test
    void subFDR_noStepPassesSafetyFactor_resultsStayEmpty() {
        // all-TT group with target FDR 0: efdr==0<=fdr always, but saftyfactor=1.0
        // makes efdr_p/fdr (a division by zero -> +Infinity) fail the safety check
        // at every step, so nothing is ever accepted.
        ArrayList<FakeFDRElement> group = new ArrayList<>();
        group.add(new FakeFDRElement().score(3).type(MatchType.TT).internal(true));
        group.add(new FakeFDRElement().score(2).type(MatchType.TT).internal(true));
        group.add(new FakeFDRElement().score(1).type(MatchType.TT).internal(true));
        group.add(new FakeFDRElement().score(1).type(MatchType.TD).internal(true));

        ArrayList<FakeFDRElement> results = new ArrayList<>();
        SubGroupFdrInfo<FakeFDRElement> info = new SubGroupFdrInfo<>();
        info.TT = 3;
        info.TD = 1;
        info.DD = 0;
        info.targetFDR = 0;
        info.saftyfactor = 1.0;

        double returned = newFdrImplement().subFDR(group, results, true, info);

        assertEquals(0.0, returned, DELTA);
        assertEquals(0, info.resultCount);
        assertTrue(results.isEmpty());
        assertTrue(info.filteredResult.isEmpty());
        assertEquals(0, info.firstPassingFDR, DELTA);
    }

    private ArrayList<FakeFDRElement> twentyElementInput() {
        ArrayList<FakeFDRElement> input = new ArrayList<>();
        for (int score = 20; score >= 5; score--) {
            input.add(new FakeFDRElement().score(score).type(MatchType.TT).internal(true));
        }
        input.add(new FakeFDRElement().score(4).type(MatchType.TD).internal(true));
        input.add(new FakeFDRElement().score(3).type(MatchType.TD).internal(true));
        input.add(new FakeFDRElement().score(2).type(MatchType.TD).internal(true));
        input.add(new FakeFDRElement().score(1).type(MatchType.DD).internal(true));
        return input;
    }

    @Test
    void fdr_ignoreGroups_matchesStandaloneSubFDR() {
        ArrayList<FakeFDRElement> input = twentyElementInput();

        FDRSettingsImpl settings = new FDRSettingsImpl();
        settings.ignoreValidityChecks(false);
        settings.setMinTD(2);

        FDRImplement fdrImpl = newFdrImplement();
        FDRResultLevel<FakeFDRElement> groupInfo = new FDRResultLevel<>();
        fdrImpl.fdr(0.2, settings, input, groupInfo, input.size(), input.size(), 0, true, true, false, false);

        SubGroupFdrInfo<FakeFDRElement> all = groupInfo.getGroup("ALL");
        assertNotNull(all);
        assertNull(all.didntPassCheck);

        // standalone subFDR on a fresh copy of the same data
        ArrayList<FakeFDRElement> standaloneGroup = twentyElementInput();
        SubGroupFdrInfo<FakeFDRElement> standaloneInfo = new SubGroupFdrInfo<>();
        standaloneInfo.TT = 16;
        standaloneInfo.TD = 3;
        standaloneInfo.DD = 1;
        standaloneInfo.targetFDR = 0.2;
        standaloneInfo.saftyfactor = settings.getReportFactor();
        ArrayList<FakeFDRElement> standaloneResults = new ArrayList<>();
        fdrImpl.subFDR(standaloneGroup, standaloneResults, true, standaloneInfo);

        assertEquals(standaloneInfo.resultCount, all.resultCount);
        assertEquals(standaloneInfo.resultTT, all.resultTT);
        assertEquals(standaloneInfo.resultTD, all.resultTD);
        assertEquals(standaloneInfo.resultDD, all.resultDD);

        assertEquals(16, all.within);
        assertEquals(all.within, groupInfo.getWithin());
        assertEquals(all.between, groupInfo.getBetween());
        assertEquals(all.linear, groupInfo.getLinear());
    }

    @Test
    void fdr_discardsGroupThatFailsValidityCheck() {
        // 2 TT, 5 TD, 5 DD: resultTT can be at most 2, so resultTT*targetFDR(0.05) < minTD(2) always.
        ArrayList<FakeFDRElement> input = new ArrayList<>();
        input.add(new FakeFDRElement().score(12).type(MatchType.TT).internal(true));
        input.add(new FakeFDRElement().score(11).type(MatchType.TT).internal(true));
        for (int score = 10; score >= 6; score--) {
            input.add(new FakeFDRElement().score(score).type(MatchType.TD).internal(true));
        }
        for (int score = 5; score >= 1; score--) {
            input.add(new FakeFDRElement().score(score).type(MatchType.DD).internal(true));
        }

        FDRSettingsImpl settings = new FDRSettingsImpl();
        settings.ignoreValidityChecks(false);
        settings.setMinTD(2);

        FDRImplement fdrImpl = newFdrImplement();
        FDRResultLevel<FakeFDRElement> groupInfo = new FDRResultLevel<>();
        fdrImpl.fdr(0.05, settings, input, groupInfo, input.size(), input.size(), 0, true, true, false, false);

        SubGroupFdrInfo<FakeFDRElement> all = groupInfo.getGroup("ALL");
        assertNotNull(all);
        assertEquals("not enough TT", all.didntPassCheck);
        assertTrue(all.filteredResult.isEmpty());
    }
}
