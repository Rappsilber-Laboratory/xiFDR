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

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import org.rappsilber.fdr.entities.FDRSelfAdd;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.fdr.entities.PeptidePair;
import org.rappsilber.fdr.entities.ProteinGroup;
import org.rappsilber.fdr.entities.ProteinGroupLink;
import org.rappsilber.fdr.entities.ProteinGroupPair;
import org.rappsilber.fdr.utils.HashedArrayList;

/**
 *
 * @author lfischer
 */
public class FDRResultLevel<T extends FDRSelfAdd>  implements Iterable<T> {
    HashMap<String,SubGroupFdrInfo<T>> groups = new HashMap<>();
    public boolean isDirectional = false;
    private boolean localFDR = false;
    private int within;
    private int between;
    private int linear;
    
    public PSM contains_psm_id(String id) {
        T f = groups.values().iterator().next().results.iterator().next();
        if (f instanceof PSM) {
            for (SubGroupFdrInfo<T>  sg : groups.values()) {
                for (T e: sg.results) {
                    if (((PSM)e).getPsmID().contentEquals(id))
                        return (PSM)e;
                }
                
            }
            return null;
        }else if (f instanceof PeptidePair) {
            for (SubGroupFdrInfo<T>  sg : groups.values()) {
                for (T e: sg.results) {
                    for (PSM p : ((PeptidePair) e).getAllPSMs())
                        if (p.getPsmID().contentEquals(id))
                            return p;
                }
            }
            return null;
        }else if (f instanceof ProteinGroup) {
            for (SubGroupFdrInfo<T>  sg : groups.values()) {
                for (T e: sg.results) {
                    for (PeptidePair pp : ((ProteinGroup) e).getPeptidePairs())
                        for (PSM p : pp.getAllPSMs())
                            if (p.getPsmID().contentEquals(id))
                                return p;
                }
            }
            return null;
        }else if (f instanceof ProteinGroupLink) {
            for (SubGroupFdrInfo<T>  sg : groups.values()) {
                for (T e: sg.results) {
                    for (PeptidePair pp : ((ProteinGroupLink) e).getPeptidePairs())
                        for (PSM p : pp.getAllPSMs())
                            if (p.getPsmID().contentEquals(id))
                                return p;
                }
            }
            return null;
        }else if (f instanceof ProteinGroupPair) {
            for (SubGroupFdrInfo<T>  sg : groups.values()) {
                for (T e: sg.results) {
                    for (PeptidePair pp : ((ProteinGroupPair) e).getPeptidePairs())
                        for (PSM p : pp.getAllPSMs())
                            if (p.getPsmID().contentEquals(id))
                                return p;
                }
            }
            return null;
        }        
        return null;
    }

    public int getInputCount() {
        int c=0;
        for (SubGroupFdrInfo<T> g: groups.values())
            c+=g.inputCount;
        return c;
    }

    public int getResultCount() {
        int c=0;
        for (SubGroupFdrInfo<T> g: groups.values())
            c+=g.filteredResult.size();
        return c;
    }

    public double getLowerFDR() {
        int c=0;
        double fdr = 0d;
        for (SubGroupFdrInfo<T> g: groups.values()) {
            fdr+=g.results.size()*g.lowerFDR;
            c+=g.results.size();
        }
        return fdr/c;
    }

    public double getHigherFDR() {
        int c=0;
        double fdr = 0d;
        for (SubGroupFdrInfo<T> g: groups.values()) {
            fdr+=g.results.size()*g.higherFDR;
            c+=g.results.size();
        }
        return fdr/c;
    }

    public double geFirstParsingFDR() {
        int c=0;
        double fdr = 0d;
        for (SubGroupFdrInfo<T> g: groups.values()) {
            fdr+=g.results.size()*g.firstPassingFDR;
            c+=g.results.size();
        }
        return fdr/c;
    }

    public double getTargetFDR() {
        int c=0;
        double fdr = 0d;
        for (SubGroupFdrInfo<T> g: groups.values()) {
            fdr+=g.results.size()*g.targetFDR;
            c+=g.results.size();
        }
        return fdr/c;
    }

    public Iterator<T>iterator() {
        final Iterator<SubGroupFdrInfo<T>> gi = groups.values().iterator();

        return new Iterator<T>() {
            SubGroupFdrInfo<T> ng = null;
            Iterator<T> ngi = null;
            public boolean hasNext() {
                if (ngi!= null  && ngi.hasNext())
                    return true;

                while (gi.hasNext()) {
                    ng=gi.next();
                    ngi = ng.iterator();
                    if (ngi.hasNext())
                        return true;
                }
                return false;
            }

            public T next() {
                return ngi.next();
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
            }
        };

    }
    
    public Iterable<T> filteredResults() {
        return new Iterable<T>() {

            @Override
            public Iterator<T> iterator() {
                return filteredIterator();
            }
            
        };
    }

    public Iterator<T> filteredIterator() {
        final Iterator<SubGroupFdrInfo<T>> gi = groups.values().iterator();

        return new Iterator<T>() {
            SubGroupFdrInfo<T> ng = null;
            Iterator<T> ngi = null;
            public boolean hasNext() {
                if (ngi!= null  && ngi.hasNext())
                    return true;

                while (gi.hasNext()) {
                    ng=gi.next();
                    ngi = ng.filteredResult.iterator();
                    if (ngi.hasNext())
                        return true;
                }
                return false;
            }

            public T next() {
                return ngi.next();
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
            }
        };

    }
    
    public boolean filteredContains(T e)  {
        for (SubGroupFdrInfo<T> g : groups.values()) {
            if (g.filteredContains(e))
                return true;
        }
        return false;
    }

    public T filteredGet(T e)  {
        for (SubGroupFdrInfo<T> g : groups.values()) {
            T ret = g.filteredGet(e);
            if (ret!=null) {
                return ret;
            }
        }
        return null;
    }    
    
    public void retainAll(Collection<T> k) {
        between =0;
        linear =0;
        within = 0;
        for (SubGroupFdrInfo<T> g: groups.values()) {
            g.filteredResult = new HashedArrayList<T>(g.filteredResult);
            g.filteredResult.retainAll(k);
  
            for (T e : g.filteredResult) {
                if (e.isTT()) {
                    if (e.isLinear()) 
                        linear++;
                    if (e.isBetween()) 
                        between++;
                    else if (e.isInternal()) 
                        within++;
                }
            }
        }
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

    public Set<String> getGroupIDs() {
        return groups.keySet();
    }

    public Collection<SubGroupFdrInfo<T>> getGroups() {
        return groups.values();
    }

    public SubGroupFdrInfo<T> getGroup(String id) {
        return groups.get(id);
    }

    public void addGroup(String id, SubGroupFdrInfo<T> group) {
        groups.put(id, group);
    }    

    public void addGroup(SubGroupFdrInfo<T> group) {
        groups.put(group.fdrGroup, group);
    }    
    
    public int size() {
        int ret = 0;
        for (SubGroupFdrInfo<T> sg : getGroups()) {
            ret += sg.within + sg.between + sg.linear;
        }
        return ret;
    }
    
    public void clear() {
        for (SubGroupFdrInfo si : groups.values()) {
            si.clear();
        }        
    }

    /**
     * @return the localFDR
     */
    public boolean isLocalFDR() {
        return localFDR;
    }

    /**
     * @param localFDR the localFDR to set
     */
    public void setLocalFDR(boolean localFDR) {
        this.localFDR = localFDR;
    }
}
