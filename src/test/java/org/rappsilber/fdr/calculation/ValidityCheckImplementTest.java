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

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.rappsilber.fdr.result.SubGroupFdrInfo;
import org.rappsilber.fdr.testutil.FakeFDRElement;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;

public class ValidityCheckImplementTest {

    @ParameterizedTest(name = "{index}: targetFDR={0} TT={1} TD={2} DD={3} minTD={4} factor={5} -> {6}")
    @CsvSource(value = {
        "1.0,    0,    0, 100, 2, 0,   NULL",
        "0.05,  10,    0,   0, 2, 0,   not enough TT",
        "0.5,    2,   10,   5, 0, 0,   to many DD",
        "0.5,  100,    2,   5, 2, 0,   to many DD",
        "0.5,   10,    2,   2, 2, 0.5, resolution to bad (TD vs DD)",
        "0.1,    5,   10,   1, 0, 0.1, resolution to bad (TT count)",
        "0.05, 1000,  10,   1, 2, 0,   NULL"
    }, nullValues = "NULL")
    void checkValid_evaluatesEachBranch(double targetFDR, int resultTT, int resultTD, int resultDD,
            int minTD, double factor, String expected) {
        SubGroupFdrInfo<FakeFDRElement> info = new SubGroupFdrInfo<>();
        info.targteFDR = targetFDR;
        info.resultTT = resultTT;
        info.resultTD = resultTD;
        info.resultDD = resultDD;

        ValidityCheckImplement valid = new ValidityCheckImplement(factor, minTD);

        assertEquals(expected, valid.checkValid(info, minTD, factor));
    }

    @Test
    void checkValid_oneArgDelegatesToConstructorFactorAndMinTD() {
        SubGroupFdrInfo<FakeFDRElement> info = new SubGroupFdrInfo<>();
        info.targteFDR = 0.05;
        info.resultTT = 10;
        info.resultTD = 0;
        info.resultDD = 0;

        ValidityCheckImplement valid = new ValidityCheckImplement(0, 2);

        assertEquals("not enough TT", valid.checkValid(info));

        SubGroupFdrInfo<FakeFDRElement> passing = new SubGroupFdrInfo<>();
        passing.targteFDR = 0.05;
        passing.resultTT = 1000;
        passing.resultTD = 10;
        passing.resultDD = 1;

        assertNull(valid.checkValid(passing));
    }
}
