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

import static org.junit.jupiter.api.Assertions.assertEquals;

public class MaximizeLevelInfoIntegerTest {

    private static final double DELTA = 1e-12;

    @Test
    void constructor_floorsStepWidthBelowOneToOne() {
        // raw stepWidth = (1.0-0.0005)/10 = 0.09995 -> rounds to 0 -> floored to 1
        MaximizeLevelInfoInteger info = new MaximizeLevelInfoInteger(1.0, true, 10);

        assertEquals(1, info.stepWidth, DELTA);
    }

    @Test
    void constructor_roundsStepWidthToNearestInt() {
        // raw stepWidth = (10.0-0.0005)/4 = 2.499875 -> rounds to 2
        MaximizeLevelInfoInteger info = new MaximizeLevelInfoInteger(10.0, true, 4);

        assertEquals(2, info.stepWidth, DELTA);
    }

    @Test
    void getCurrentFDR_roundsUnderlyingDouble() {
        MaximizeLevelInfoInteger info = new MaximizeLevelInfoInteger(10.0, true, 4);

        info.setCurrentFDR(7.6);
        assertEquals(8, info.getCurrentFDR(), DELTA);

        info.setCurrentFDR(7.4);
        assertEquals(7, info.getCurrentFDR(), DELTA);

        info.setCurrentFDR(7.5);
        assertEquals(8, info.getCurrentFDR(), DELTA);
    }

    @Test
    void setNewMaxFDR_usesRoundedCurrentFDR() {
        MaximizeLevelInfoInteger info = new MaximizeLevelInfoInteger(10.0, true, 4);

        info.setCurrentFDR(7.6);
        info.setNewMaxFDR();

        assertEquals(8, info.maximumFDR, DELTA);
        assertEquals(8, info.smallestEqualFDR, DELTA);
        assertEquals(8, info.largestEqualFDR, DELTA);
    }

    @Test
    void calcNextFDRRange_recomputesAndRoundsStepWidth() {
        // fromFDR=0.0005, toFDR=10.0, stepWidth=2 (after construction-time rounding), firstFDR=10.0
        MaximizeLevelInfoInteger info = new MaximizeLevelInfoInteger(10.0, true, 4);
        info.smallestEqualFDR = 8.0;
        info.largestEqualFDR = 8.0;

        info.calcNextFDRRange(false, 0);

        // toFDR = min(largestEqualFDR + stepWidth*3/4, firstFDR) = min(8 + 1.5, 10) = 9.5
        assertEquals(9.5, info.toFDR, DELTA);
        // newFrom = max(smallestEqualFDR - stepWidth*3/4, 0) = max(8 - 1.5, 0) = 6.5 (!= old fromFDR 0.0005)
        assertEquals(6.5, info.fromFDR, DELTA);
        // raw stepWidth = (9.5-6.5)/4 = 0.75 -> rounds to 1
        assertEquals(1, info.stepWidth, DELTA);
    }
}
