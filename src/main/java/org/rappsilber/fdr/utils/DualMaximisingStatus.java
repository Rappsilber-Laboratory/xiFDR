/*
 * Copyright 2018 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
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

import org.rappsilber.fdr.result.FDRResult;

/**
 * Result of a dual boosting run that searches for the link/PPI FDR cutoffs
 * that maximise the Self counts and the Between counts independently, while
 * sharing all upstream (PSM/peptide-pair/protein-group FDR and prefilter)
 * settings between both searches.
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class DualMaximisingStatus {
    /** status of the search optimised for Self link/PPI counts */
    public MaximisingStatus self;
    /** status of the search optimised for Between link/PPI counts */
    public MaximisingStatus between;
    /** merged result: Self/linear groups from {@link #self}, Between groups from {@link #between} */
    public FDRResult result;
}
