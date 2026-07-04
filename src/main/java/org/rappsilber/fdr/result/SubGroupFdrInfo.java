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

import java.util.Iterator;
import org.rappsilber.fdr.entities.FDRSelfAdd;
import org.rappsilber.fdr.utils.HashedArrayList;

/**
 *
 * @author lfischer
 */
public class SubGroupFdrInfo<T extends FDRSelfAdd> implements Iterable<T> {
    public String fdrGroup;
    /**  how many TD cases are there in total */
    public int TD;
    /**  how many DD cases are there in total */
    public int DD;
    /**  how many TT cases are there in total */
    public int TT;
    /**  how many TD cases are there in total */
    public int resultTD;
    /**  how many DD cases are there in total */
    public int resultDD;
    /**  how many TT cases are there in total */
    public int resultTT;

    /**  target fdr */
    public double targteFDR;
    /**  some safety factor for reporting */
    public double saftyfactor;

    public int inputCount;

    /** how many passed the fdr */
    public int resultCount;

    public HashedArrayList<T> results;

    public HashedArrayList<T> filteredResult;

    /** what was the next higher fdr above the target fdr */
    public double firstPassingFDR;
    /** 
     * what was the next higher or equal fdr then the accepted fdr 
     * unless the saftyfactor comes into play it should be the same as firstPassingFDR
     */
    public double higherFDR;
    /** what was the next lower fdr then the target fdr */
    public double lowerFDR;

    /** normalisation factor for the target database */
    public double TCount=1;
    public double DCount = 1;

    public double worstAcceptedScore = Double.NaN;
    
    public int within;
    public int between;
    public int linear;
    
    public String didntPassCheck = null;

    /**
     * when true the FDR is calculated as (TD+DD)/TT instead of (TD-DD)/TT
     * (used for self-links with overlapping peptides)
     */
    public boolean useAdditiveFDR = false;
    
    public Iterator<T> filteredIterator() {
        return filteredResult.iterator();
    }

    public Iterator<T> iterator() {
        return results.iterator();
    }
    
    public boolean filteredContains(T e)  {
        return filteredResult.contains(e);
    }

    public T filteredGet(T e)  {
        return filteredResult.get(e);
    }
    
    /**
     * @return the within
     */
    public int getWithin() {
        return within;
    }

    /**
     * @param within the within to set
     */
    public void setWithin(int within) {
        this.within = within;
    }

    /**
     * @return the between
     */
    public int getBetween() {
        return between;
    }

    /**
     * @param between the between to set
     */
    public void setBetween(int between) {
        this.between = between;
    }

    /**
     * @return the linear
     */
    public int getLinear() {
        return linear;
    }

    /**
     * @param linear the linear to set
     */
    public void setLinear(int linear) {
        this.linear = linear;
    }

    void clear() {
        results.clear();
        filteredResult.clear();
    }


}

