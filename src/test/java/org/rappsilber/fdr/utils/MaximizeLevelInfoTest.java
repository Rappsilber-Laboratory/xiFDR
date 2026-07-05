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
package org.rappsilber.fdr.utils;

import org.junit.jupiter.api.Test;
import org.rappsilber.fdr.result.FDRResultLevel;
import org.rappsilber.fdr.testutil.FakeFDRElement;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class MaximizeLevelInfoTest {

    private static final double DELTA = 1e-12;

    @Test
    void constructor_boost_picksMinimumFromFdrAndComputesStepWidth() {
        // minimumFDR (0) is below the 0.0005 floor, and toFDR/steps (0.005) is above it,
        // so fromFDR ends up at the 0.0005 floor.
        MaximizeLevelInfo info = new MaximizeLevelInfo(0.05, 0, true, 10);

        assertEquals(0.05, info.toFDR, DELTA);
        assertEquals(0.05, info.firstFDR, DELTA);
        assertEquals(0.05, info.maximumFDR, DELTA);
        assertEquals(0.05, info.getCurrentFDR(), DELTA);
        assertEquals(0.0005, info.fromFDR, DELTA);
        assertEquals(0.00495, info.stepWidth, DELTA);
    }

    @Test
    void constructor_noBoost_fromEqualsToAndStepWidthIs100() {
        MaximizeLevelInfo info = new MaximizeLevelInfo(0.05, 0, false, 10);

        assertEquals(0.05, info.toFDR, DELTA);
        assertEquals(0.05, info.fromFDR, DELTA);
        assertEquals(100, info.stepWidth, DELTA);
    }

    @Test
    void stepping_iteratesFromToFdrDownToFromFdrInclusive() {
        // toFDR=1.0, fromFDR=0.5, stepWidth=0.25 -> currentFDR visits 1.0, 0.75, 0.5, then 0.25 stops.
        MaximizeLevelInfo info = new MaximizeLevelInfo(1.0, 0.5, true, 2);

        info.firstStep();
        assertEquals(1.0, info.getCurrentFDR(), DELTA);

        int steps = 0;
        while (info.doThisStep()) {
            steps++;
            info.nextStep();
        }

        assertEquals(3, steps);
        assertEquals(0.25, info.getCurrentFDR(), DELTA);
        assertFalse(info.doThisStep());
    }

    @Test
    void stepping_noBoost_doThisStepTrueOnlyForFirstStep() {
        MaximizeLevelInfo info = new MaximizeLevelInfo(5.0, 0, false, 10);

        info.firstStep();
        assertTrue(info.doThisStep());

        info.nextStep();
        assertFalse(info.doThisStep());
    }

    @Test
    void setNewMaxFDR_setsMaximumAndEqualBoundsToCurrentFDR() {
        MaximizeLevelInfo info = new MaximizeLevelInfo(1.0, 0.5, true, 2);

        info.setCurrentFDR(0.75);
        info.setNewMaxFDR();

        assertEquals(0.75, info.maximumFDR, DELTA);
        assertEquals(0.75, info.smallestEqualFDR, DELTA);
        assertEquals(0.75, info.largestEqualFDR, DELTA);
    }

    @Test
    void setEqualFDR_widensSmallestAndLargestBounds() {
        MaximizeLevelInfo info = new MaximizeLevelInfo(1.0, 0.5, true, 2);
        info.smallestEqualFDR = 0.5;
        info.largestEqualFDR = 0.5;

        info.setCurrentFDR(0.75);
        info.setEqualFDR();
        assertEquals(0.5, info.smallestEqualFDR, DELTA);
        assertEquals(0.75, info.largestEqualFDR, DELTA);

        info.setCurrentFDR(0.25);
        info.setEqualFDR();
        assertEquals(0.25, info.smallestEqualFDR, DELTA);
        assertEquals(0.75, info.largestEqualFDR, DELTA);
    }

    @Test
    void calcNextFDRRange_recomputesBoundsAndStepWidth() {
        // fromFDR=0.5, toFDR=1.0, stepWidth=0.25, firstFDR=1.0, minimumFDR=0.5
        MaximizeLevelInfo info = new MaximizeLevelInfo(1.0, 0.5, true, 2);
        info.smallestEqualFDR = 0.5;
        info.largestEqualFDR = 0.75;

        info.calcNextFDRRange(false, 0);

        // toFDR = min(largestEqualFDR + stepWidth*3/4, firstFDR) = min(0.75+0.1875, 1.0) = 0.9375
        assertEquals(0.9375, info.toFDR, DELTA);
        // newFrom = max(smallestEqualFDR - stepWidth*3/4, 0) = max(0.5-0.1875, 0) = 0.3125 (!= old fromFDR 0.5)
        assertEquals(0.3125, info.fromFDR, DELTA);
        assertEquals(0, info.stepChange, DELTA);
        // stepWidth = (toFDR-fromFDR)/(steps+stepChange) = (0.9375-0.3125)/2 = 0.3125
        assertEquals(0.3125, info.stepWidth, DELTA);
    }

    @Test
    void setCountsAndSetCountsPrefilter_readFromResultLevel() {
        MaximizeLevelInfo info = new MaximizeLevelInfo(0.05, 0, false, 1);
        FDRResultLevel<FakeFDRElement> level = new FDRResultLevel<>();
        level.setWithin(5);
        level.setBetween(3);
        level.setLinear(2);

        info.setCounts(level);
        assertEquals(8, info.count);
        assertEquals(3, info.countBetween);
        assertEquals(2, info.countLinear);

        info.setCountsPrefilter(level);
        assertEquals(8, info.countPreFilter);
        assertEquals(3, info.countBetweenPreFilter);
        assertEquals(2, info.countLinearPreFilter);
    }
}
