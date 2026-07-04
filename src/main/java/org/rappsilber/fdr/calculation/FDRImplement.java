/*
 * Copyright 2020 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
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
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.rappsilber.fdr.FDRSettings;
import org.rappsilber.fdr.entities.AbstractFDRElement;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.fdr.entities.ProteinGroup;
import org.rappsilber.fdr.result.FDRResultLevel;
import org.rappsilber.fdr.result.SubGroupFdrInfo;
import org.rappsilber.fdr.utils.HashedArrayList;
import org.rappsilber.utils.RArrayUtils;
import org.rappsilber.utils.UpdateableInteger;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class FDRImplement implements FDR {
    
    CheckValid valid;
    
    /**
     * size of the decoy independent protein pairs in terms of psms
     */
    HashMap<Integer, UpdateableInteger> protpairToSize = new HashMap<>();
    /**
     * id of a protein pair independent of target or decoy
     */
    HashMap<String, Integer> protpairToID = new HashMap<>();

    public FDRImplement(CheckValid valid) {
        this.valid = valid;
    }

    
    public void setValidityCheck(CheckValid c) {
        this.valid = c;
    }
    
//    private <T extends FDRSelfAdd<T>> HashedArrayList<T> fdr(double fdr, double safetyfactor, Collection<T> c, HashMap<Integer, SubGroupFdrInfo<T>> groupInfo, double tCount, double dCount, int minPepCount, boolean ignoreGroups,boolean setElementFDR) {
    /**
     * Splits up data into groups and forward each group to
     * {@link #subFDR(java.util.ArrayList, java.util.ArrayList, boolean, rappsilber.fdr.result.SubGroupFdrInfo) subFDR}.
     * results will then be collected
     *
     * @param <T> What type of Information is the FDRImplement to be applied
     * @param fdr the target FDRImplement that is supossed to be returned
     * @param safetyfactor don't report anything if the next higher calculatable
     * @param fdrInput the input data that should be filtered by FDRImplement
     * @param groupInfo will hold all the information for the currently
     * calculated FDRs
     * @param tCount number of entries in the target database. This is used to
 normalise the FDRImplement calculation for imbalenced databases
     * @param dCount number of entries in the decoy database. This is used to
 normalise the FDRImplement calculation for imbalenced databases
     * @param minPepCount if an entry has less then these number of peptides
 supporting it - it will not be considered for the FDRImplement-calculation
     * @param ignoreGroups should we ignore groups cmpletely and just calculate
 a joined FDRImplement?
     * @param setElementFDR should each element be flaged up with the FDRImplement that
 it group has at the given score? - is used to speed up the maximation (by
 not setting these)
     */
    @Override
    public <T extends AbstractFDRElement<T>> void fdr(double fdr, FDRSettings settings, Collection<T> fdrInput, FDRResultLevel<T> groupInfo, double tCount, double dCount, int minPepCount, boolean ignoreGroups, boolean setElementFDR, Boolean localFDR, boolean groupByProteinPair) {
        HashMap<String, ArrayList<T>> groupedList = new HashMap<String, ArrayList<T>>(4);
        HashMap<String, UpdateableInteger> gTT = new HashMap<String, UpdateableInteger>(8);
        HashMap<String, UpdateableInteger> gTD = new HashMap<String, UpdateableInteger>(8);
        HashMap<String, UpdateableInteger> gDD = new HashMap<String, UpdateableInteger>(8);
        int discardedGroupsBetween = 0;
        int discardedGroupsOther = 0;
        double safetyfactor = settings.getReportFactor();
        int minTDChance = settings.getMinTD();

        int resultcount = 0;

        T firstinput = null;
        // we count the PSMs 
        boolean countmatches = false;
        // if we group by protein pairs we do that by number of PSMs
        // so we only count the size of a protein pair on level of PSMs
        if (groupByProteinPair && !fdrInput.isEmpty()) {
            firstinput = fdrInput.iterator().next();
            if (firstinput instanceof PSM) {
                countmatches = true;
            }
        }

        if (fdrInput.isEmpty()) {
            return;
        }

        if (fdr == 1) {
            fdr = 1000;
            safetyfactor = 1000;
        }

        int maxID = 2;

        // split the data up into fdr-groups
        for (T e : fdrInput) {
            if (e.getPeptidePairCount() >= minPepCount) {
                String fdrgroup = null;
                if (ignoreGroups) { // don't do any grouping
                    fdrgroup = "ALL";
                } else {
                    if (groupByProteinPair) {
                        fdrgroup = fdrGroupByProteinPair(e, maxID, countmatches);
                    } else {
                        fdrgroup = e.getFDRGroup();
                    }
                }
                addElementToGroup(groupedList, fdrgroup, e, gTT, gTD, gDD);
            }
        }

        // join groups with similare number of PSMs
        if (groupByProteinPair) {
            HashMap<String, ArrayList<T>> groupedListNew = new HashMap<String, ArrayList<T>>(4);
            HashMap<String, UpdateableInteger> gTTNew = new HashMap<String, UpdateableInteger>(8);
            HashMap<String, UpdateableInteger> gTDNew = new HashMap<String, UpdateableInteger>(8);
            HashMap<String, UpdateableInteger> gDDNew = new HashMap<String, UpdateableInteger>(8);
            groupedList = joinProteinGroupsBySize(groupedList, gTT, gTD, gDD, groupedListNew, gTTNew, gTDNew, gDDNew);
            gTT = gTTNew;
            gTD = gTDNew;
            gDD = gDDNew;
        }

        // subgroups that do not check as valid result get collected 
        SubGroupFdrInfo<T> collectedBetween = new SubGroupFdrInfo<T>();
        SubGroupFdrInfo<T> collectedOthers = new SubGroupFdrInfo<T>();
        ArrayList<T> collectedElementsBetween = new ArrayList<>();
        ArrayList<T> collectedElementsOthers = new ArrayList<>();

        for (String fdrgroup : groupedList.keySet()) {
            SubGroupFdrInfo<T> info = new SubGroupFdrInfo<T>();
            info.TT = gTT.get(fdrgroup).value;
            info.TD = gTD.get(fdrgroup).value;
            info.DD = gDD.get(fdrgroup).value;
            ArrayList<T> group = groupedList.get(fdrgroup);

            info.DCount = dCount;
            info.TCount = tCount;

            ArrayList<T> groupResult = new ArrayList<T>();
            info.inputCount = group.size();
            info.targteFDR = fdr;
            info.saftyfactor = safetyfactor;
            info.fdrGroup = fdrgroup;
            info.useAdditiveFDR = fdrgroup.contains("Overlapping");
            // do a local fdr
            if (localFDR == null || localFDR) {
                subFDRLocal(group, info);
            }
            subFDR(group, groupResult, setElementFDR, info);
            if (localFDR != null && localFDR) {
                groupResult = filterByPEP(groupResult, fdr, info);
            }

            String valid = this.valid.checkValid(info, settings.getMinTD(), 0);
            if (valid != null) {
                if ((!settings.ignoreValidityChecks())) {
                    Logger.getLogger(this.getClass().getName()).log(Level.FINE, "Discarded group {0}->{1}({2})", new Object[]{group.get(0).getClass(), fdrgroup, group.get(0).getFDRGroup()});
                    //between
                    if (fdrgroup.toLowerCase().contains("between")) {
                        discardedGroupsBetween ++;
                        addGroupToGroup(collectedElementsBetween, collectedBetween, info, group);
                        groupResult.clear();

                        
                    } else {
                        // within group
                        discardedGroupsOther ++;
                        addGroupToGroup(collectedElementsOthers, collectedOthers, info, group);
                        groupResult.clear();
                        groupResult.clear();

                    }
                    info.filteredResult = new HashedArrayList<>();
                } else {
                    Logger.getLogger(this.getClass().getName()).log(Level.FINE, "FDR-Group {0}->{1}({2}) deemed unreliable - be carefull", new Object[]{group.get(0).getClass(), fdrgroup, group.get(0).getFDRGroup()});
                }
                info.didntPassCheck = valid;
            } //else {
            groupInfo.addGroup(info);

            resultcount += groupResult.size();
            groupInfo.setLinear(groupInfo.getLinear() + info.linear);
            groupInfo.setWithin(groupInfo.getWithin() + info.within);
            groupInfo.setBetween(groupInfo.getBetween() + info.between);
            //            resultCounts.put(fdrgroup, groupResult.size());
            //}

        }

        // did we collect some subgroups with to small numbers?
        if (discardedGroupsBetween + discardedGroupsOther > 1) {
            Logger.getLogger(this.getClass().getName()).log(Level.FINE, "Join up discarded groups to try and get some more results");
            if (discardedGroupsBetween >1) {
                ArrayList<T> groupResult = new ArrayList<T>();
                collectedBetween.targteFDR = fdr;
                collectedBetween.saftyfactor = safetyfactor;
                collectedBetween.fdrGroup = "CollectedResultsBetween [" + collectedBetween.fdrGroup + "]";
                if (localFDR == null || localFDR) {
                    subFDRLocal(collectedElementsBetween, collectedBetween);
                }
                subFDR(collectedElementsBetween, groupResult, setElementFDR, collectedBetween);
                if (localFDR != null && localFDR) {
                    groupResult = filterByPEP(groupResult, fdr, collectedBetween);
                }
                if (settings.ignoreValidityChecks() || this.valid.checkValid(collectedBetween, settings.getMinTD(), 0) == null) {
                    groupInfo.addGroup(collectedBetween);
                    //                    group, fdr, safetyfactor, groupResult, tCount, dCount,setElementFDR);
                    //            nextFDR.put(fdrgroup, prevFDR);
                    //            ret.addAll(groupResult);
                    resultcount += groupResult.size();
                    groupInfo.setLinear(groupInfo.getLinear() + collectedBetween.linear);
                    groupInfo.setWithin(groupInfo.getWithin() + collectedBetween.within);
                    groupInfo.setBetween(groupInfo.getBetween() + collectedBetween.between);
                }
            } else if (collectedElementsOthers.size() > 0 && discardedGroupsBetween == 1) {
                // we have only one failed between group, but could try if we join with the other groups.
                addGroupToGroup(collectedElementsOthers, collectedOthers, collectedBetween, collectedElementsBetween);
            }

            if (collectedElementsOthers.size() > 1) {
                ArrayList<T> groupResultwithin = new ArrayList<T>();

                collectedOthers.targteFDR = fdr;
                collectedOthers.saftyfactor = safetyfactor;
                collectedOthers.fdrGroup = "CollectedSmallResultsOthers  [" + collectedOthers.fdrGroup + "]";
                if (localFDR == null || localFDR) {
                    subFDRLocal(collectedElementsOthers, collectedOthers);
                }
                subFDR(collectedElementsOthers, groupResultwithin, setElementFDR, collectedOthers);
                if (localFDR != null && localFDR) {
                    groupResultwithin = filterByPEP(groupResultwithin, fdr, collectedOthers);
                }
                if (settings.ignoreValidityChecks() || this.valid.checkValid(collectedOthers, settings.getMinTD(), 0) == null) {
                    groupInfo.addGroup(collectedOthers);
                    //                    group, fdr, safetyfactor, groupResult, tCount, dCount,setElementFDR);
                    //            nextFDR.put(fdrgroup, prevFDR);
                    //            ret.addAll(groupResult);
                    resultcount += groupResultwithin.size();
                    groupInfo.setLinear(groupInfo.getLinear() + collectedOthers.linear);
                    groupInfo.setWithin(groupInfo.getWithin() + collectedOthers.within);
                    groupInfo.setBetween(groupInfo.getBetween() + collectedOthers.between);
                }
            }
        }

        if (resultcount == 0 && fdrInput.size() > 100) {
            // we didn't get any results through. try a non grouped
            SubGroupFdrInfo info = new SubGroupFdrInfo();
            info.fdrGroup = "NoSubResults";
            ArrayList<T> nothing = new ArrayList<>(0);
            for (String fdrgroup : groupedList.keySet()) {
                SubGroupFdrInfo gi =  groupInfo.getGroup(fdrgroup);
                addGroupToGroup(nothing, info, gi, nothing);
            }

            info.inputCount = fdrInput.size();

            ArrayList<T> allResult = new ArrayList<T>();

            ArrayList<T> all = new ArrayList<T>(fdrInput);
            if (localFDR == null || localFDR) {
                subFDRLocal(all, info);
            }
            
            subFDR(all, allResult, setElementFDR, info);
            if (localFDR != null && localFDR) {
                allResult = filterByPEP(allResult, fdr, info);
            }
            if (settings.ignoreValidityChecks() || valid.checkValid(info, settings.getMinTD(), 0) == null) {
                groupInfo.addGroup(info);
                groupInfo.setLinear(groupInfo.getLinear() + info.linear);
                groupInfo.setWithin(groupInfo.getWithin() + info.within);
                groupInfo.setBetween(groupInfo.getBetween() + info.between);
            } else {
                allResult.clear();
            } 
        }

//        return ret;
    }

    protected <T extends AbstractFDRElement<T>> ArrayList<T> filterByPEP(ArrayList<T> groupResultwithin, double fdr, SubGroupFdrInfo<T> collectedToSmallWithin) {
        ArrayList<T> filtered = new ArrayList<>();
        for (T e : groupResultwithin) {
            if (e.getPEP() <= fdr) {
                filtered.add(e);
            } else {
                if (e.isTT()) {
                    collectedToSmallWithin.resultTT--;
                } else if (e.isTD()) {
                    collectedToSmallWithin.resultTD--;
                } else {
                    collectedToSmallWithin.resultDD--;
                }
            }
        }
        groupResultwithin = filtered;
        return groupResultwithin;
    }

    protected <T extends AbstractFDRElement<T>> void addGroupToGroup(ArrayList<T> collectedElements, SubGroupFdrInfo<T> collectedInfo, SubGroupFdrInfo<T> addInfo, ArrayList<T> addElements) {
        collectedInfo.TT += addInfo.TT;
        collectedInfo.TD += addInfo.TD;
        collectedInfo.DD += addInfo.DD;
        collectedInfo.DCount += addInfo.DCount;
        collectedInfo.TCount += addInfo.TCount;
        collectedInfo.inputCount += addInfo.inputCount;
        collectedInfo.useAdditiveFDR = collectedInfo.useAdditiveFDR || addInfo.useAdditiveFDR;
        if (collectedInfo.fdrGroup == null || collectedInfo.fdrGroup.isEmpty())
            collectedInfo.fdrGroup = addInfo.fdrGroup;
        else {
            collectedInfo.fdrGroup += ", " + addInfo.fdrGroup;
        }
        collectedElements.addAll(addElements);
    }

    /**
     * This method will join different groups by amount of matches in the group.
     * @param <T>
     * @param groupedList
     * @param gTT
     * @param gTD
     * @param gDD
     * @param groupedListNew
     * @param gTTNew
     * @param gTDNew
     * @param gDDNew
     * @return 
     */
    protected <T extends AbstractFDRElement<T>> HashMap<String, ArrayList<T>> joinProteinGroupsBySize(HashMap<String, ArrayList<T>> groupedList, HashMap<String, UpdateableInteger> gTT, HashMap<String, UpdateableInteger> gTD, HashMap<String, UpdateableInteger> gDD, HashMap<String, ArrayList<T>> groupedListNew, HashMap<String, UpdateableInteger> gTTNew, HashMap<String, UpdateableInteger> gTDNew, HashMap<String, UpdateableInteger> gDDNew) {
        UpdateableInteger max = RArrayUtils.max(protpairToSize.values());
        for (Map.Entry<String, ArrayList<T>> glOldE : groupedList.entrySet()) {
            ArrayList<T> glOld = glOldE.getValue();
            String oldgroup = glOldE.getKey();
            UpdateableInteger protpairsize = protpairToSize.get(oldgroup);
            String sizekey = "" + (glOld.size() + 100);
            if (protpairsize != null) {
                sizekey = "" + protpairsize.value;
            }
            int oldTT = gTT.get(oldgroup).value;
            int oldTD = gTD.get(oldgroup).value;
            int oldDD = gDD.get(oldgroup).value;
            ArrayList<T> glNew = groupedListNew.get(sizekey);
            if (glNew == null) {
                groupedListNew.put("" + sizekey, glOld);
                
                gTTNew.put(sizekey, new UpdateableInteger(oldTT));
                gTDNew.put(sizekey, new UpdateableInteger(oldTD));
                gDDNew.put(sizekey, new UpdateableInteger(oldDD));
            } else {
                glNew.addAll(glOld);
                gTTNew.get(sizekey).value += gTT.get(oldgroup).value;
                gTDNew.get(sizekey).value += gTD.get(oldgroup).value;
                gDDNew.get(sizekey).value += gDD.get(oldgroup).value;
            }
        }
        groupedList = groupedListNew;
        return groupedList;
    }

    protected <T extends AbstractFDRElement<T>> void addElementToGroup(HashMap<String, ArrayList<T>> groupedList, String fdrgroup, T e, HashMap<String, UpdateableInteger> gTT, HashMap<String, UpdateableInteger> gTD, HashMap<String, UpdateableInteger> gDD) {
        // get the right list of matches
        ArrayList<T> gl = groupedList.get(fdrgroup);
        if (gl == null) {
            // does not exist yet - so make a new one
            gl = new ArrayList<T>();
            groupedList.put(fdrgroup, gl);
            // start counting the types
            if (e.isTT()) {
                gTT.put(fdrgroup, new UpdateableInteger(1));
                gTD.put(fdrgroup, new UpdateableInteger(0));
                gDD.put(fdrgroup, new UpdateableInteger(0));
            } else if (e.isTD()) {
                gTT.put(fdrgroup, new UpdateableInteger(0));
                gTD.put(fdrgroup, new UpdateableInteger(1));
                gDD.put(fdrgroup, new UpdateableInteger(0));
            } else {
                gTT.put(fdrgroup, new UpdateableInteger(0));
                gTD.put(fdrgroup, new UpdateableInteger(0));
                gDD.put(fdrgroup, new UpdateableInteger(1));
            }
        } else {
            if (e.isTT()) {
                gTT.get(fdrgroup).value++;
            } else if (e.isTD()) {
                gTD.get(fdrgroup).value++;
            } else {
                gDD.get(fdrgroup).value++;
            }
        }
        gl.add(e);
    }

    protected <T extends AbstractFDRElement<T>> String fdrGroupByProteinPair(T e, int maxID, boolean countmatches) {
        String fdrgroup;
        // do group by protein pairs to begin with
        ProteinGroup pg1 = e.getProteinGroup1();
        ProteinGroup pg2 = e.getProteinGroup2();
        // make sure we always get the same key
        String pgs1 = pg1.accessionsNoDecoy();
        String pgs2 = pg2.accessionsNoDecoy();
        String k1 = pgs1 + "_xl_" + pgs2;
        String k2 = pgs2 + "_xl_" + pgs1;
        Integer id = protpairToID.get(k1);
        if (id == null) // not found with key 1 -> try reversed key
        {
            id = protpairToID.get(k2);
        }
        // new protein pair
        if (id == null) {
            // assign id
            id = maxID++;
            protpairToID.put(k1, id);
            protpairToID.put(k2, id);
            protpairToSize.put(id, new UpdateableInteger(1));
        } else if (countmatches) {
            protpairToSize.get(id).value++;
        }
        // set the id as fdr group
        fdrgroup = "" + id;
        // and write it to the
        e.setFDRGroup("" + id);
        return fdrgroup;
    }

    /**
     * Takes a sub-set of entries and does the actual FDR-estimation/cutoff.
     * <p>
     * The assumption is that we have a symmetric cross-linker. <br/>Therefore:
     * <p>
     * {@code FP(TT)=TD+DD*(1-TDdb/DDdb)}</p>
     * </p><p>
     * If TCount and DCount are of different size then the false positive
     * estimation gets scaled accordingly.</p
     *
     * @param <T> What type of Information is the FDR to be applied
     * @param TD total number of TD matches
     * @param DD total number of DD matches
     * @param TT total number of TT matches
     * @param group All entries for which the target FDR should be calculated
     * @param fdr the target FDR that is supposed to be returned
     * @param safetyfactor don't report anything if the next higher calculable
     * exceeds the target-FDR by the given factor
     * @param results the passing results will be added to this ArrayList
     * @param TCount Size of the target database
     * @param DCount size of the decoy database
     * @param isSymmetric do we calculate an FDR for a symmetric or an
     * asymmetric experiment
     * @param setElementFDR
     * @return
     */

    /**
     * Takes a sub-set of entries and does the actual FDR-estimation/cutoff.
     * <p>
     * The assumption is that we have a symmetric cross-linker. <br/>Therefore:
     * <p>
     * {@code FP(TT)=TD+DD*(1-TDdb/DDdb)}</p>
     * </p><p>
     * If TCount and DCount are of different size then the false positive
     * estimation gets scaled accordingly.</p
     *
     * @param <T> What type of Information is the FDRImplement to be applied
     * @param group All entries for which the target FDRImplement should be calculated
     * @param fdr the target FDRImplement that is supposed to be returned
     * @param safetyfactor don't report anything if the next higher calculable
 exceeds the target-FDRImplement by the given factor
     * @param results the passing results will be added to this ArrayList
     * @param TCount Size of the target database
     * @param DCount size of the decoy database
     * @param info all informations needed for this sub group
     * @param isSymmetric do we calculate an FDRImplement for a symmetric or an
 asymmetric experiment
     * @param setElementFDR
     * @return
     */
    protected <T extends AbstractFDRElement<T>> double subFDR(ArrayList<T> group, ArrayList<T> results, boolean setElementFDR, SubGroupFdrInfo info) {

        HashedArrayList<T> ret = new HashedArrayList<T>();
        info.results = ret;
        info.filteredResult = ret;

        double TT = info.TT;
        int TD = info.TD;
        int DD = info.DD;

        double fdr = info.targteFDR;

        if (fdr >= 1) {
            fdr = Double.POSITIVE_INFINITY;
        }

        // all the elements within this group
        int groupSize = group.size();

        // sort them by score
        Collections.sort(group, new Comparator<T>() {

            public int compare(T o1, T o2) {
                return Double.compare(o2.getScore() * o2.getLinkedSupport(), o1.getScore() * o1.getLinkedSupport());
            }
        });

        // total fdr rate
        // for self-links with overlapping peptides TD should not occur so FDR = (TD+DD)/TT
        final boolean addDD = info.useAdditiveFDR;
        double prevFDR = addDD ? ((TD+1) + DD) / TT : ((TD+1) - DD) / TT;
        int prevTDIndex = groupSize - 1;

//        if (!isPSMScoreHighBetter())
//            Collections.reverse(group);
        int fdrSwitch = groupSize - 1;
        double highFDR = prevFDR;
        // now we can just go through and find the cut-off
        for (int i = groupSize - 1; i >= 0; i--) {

            double efdr = addDD ? (TD + DD) / TT : (TD - DD) / TT;
            double efdr_n = addDD ? ((TD - 1) + DD) / TT : ((TD - 1) - DD) / TT;
            double efdr_p = addDD ? ((TD + 1) + DD) / TT : ((TD + 1) - DD) / TT;
            if (efdr_n < 0 && DD == 0) {
                efdr_n = 0;
            }

            T e = group.get(i);
            double score = e.getScore();

            // we steped below the target fdr (efdr<= fdr) 
            if (efdr <= fdr) {
                //did we just pas an fdr-step
                if ((efdr < prevFDR || prevFDR <= 0.0 || prevFDR == Double.POSITIVE_INFINITY) || fdr == Double.POSITIVE_INFINITY) {
                    // does it fullfill the saftyfactor
                    if ((efdr_p / fdr < info.saftyfactor || info.saftyfactor > 1000 || info.saftyfactor == 0) || fdr == Double.POSITIVE_INFINITY) {
                        //if ((efdr_p / fdr < info.saftyfactor && efdr_n >= 0) || fdr == Double.POSITIVE_INFINITY)  {
                        if (info.firstPassingFDR == 0) {
                            info.firstPassingFDR = prevFDR;
                        }

                        info.higherFDR = prevFDR;
                        info.lowerFDR = efdr;
                        info.resultCount = i;
                        int lastFDRIndex = i;

                        double lastFDR = prevFDR;
                        int lastTDIndex = prevTDIndex;
                        T lastTDElement = group.get(lastTDIndex);
                        double setFDR = efdr;
                        info.resultTT = (int) TT;
                        info.resultTD = TD;
                        info.resultDD = DD;
                        if (setElementFDR) {
                            for (; i >= 0; i--) {

                                e = group.get(i);
                                e.setFDR(setFDR);
                                e.setHigherFDR(lastFDR);
                                e.setHigherFDRTD(lastTDElement);
                                ret.add(e);

                                if (e.isTT()) {
                                    TT--;
                                    if (e.isLinear()) {
                                        info.linear++;
                                    }

                                    if (e.isInternal()) {
                                        info.within++;
                                    } else if (e.isBetween()) {
                                        info.between++;
                                    }

                                } else if (e.isTD()) {
                                    TD--;
                                    for (int l = i + 1; l <= lastTDIndex; l++) {
                                        group.get(l).setLowerFDRTD(e);
                                    }
                                    lastTDIndex = i;
                                    lastTDElement = e;
                                } else if (e.isDD()) {
                                    DD--;
                                }

                                double currfdr = 0;
                                if (TD > 0) {
                                    currfdr = addDD ? (TD + DD) / TT : (TD - DD) / TT;
                                }

                                if (currfdr < setFDR) {
                                    lastFDR = setFDR;
                                    // we reached a new lower fdr
                                    setFDR = currfdr;

                                    // set the lower fdr values for the previous data
                                    for (int li = lastFDRIndex; li <= i; li++) {
                                        group.get(li).setLowerFDR(currfdr);
                                    }
                                    lastFDRIndex = i - 1;
                                    }

                            }
                            for (int li = lastFDRIndex; li <= i && li >= 0; li++) {
                                group.get(li).setLowerFDR(0);
                            }

                            T best = group.get(0);
                            for (int l = 0; l <= lastTDIndex; l++) {
                                group.get(l).setLowerFDRTD(best);
                            }

                        } else {
                            for (; i >= 0; i--) {
                                e = group.get(i);
                                if (e.isTT()) {
                                    if (e.isLinear()) {
                                        info.linear++;
                                    }

                                    if (e.isInternal()) {
                                        info.within++;
                                    } else if (e.isBetween()) {
                                        info.between++;
                                    }
                                }
                                ret.add(e);
                            }
                        }
                        //                }

                        //  info.results = ret;
                        //  info.filteredResult = ret;
                        //                if (fdr >=1 || ret.size() >=MINIMUM_POSSIBLE_RESULT)
                        results.addAll(ret);
                        info.results = ret;
                        info.filteredResult = ret;
                        info.worstAcceptedScore = score;

                        return prevFDR;
                    } else {
                        info.firstPassingFDR = prevFDR;
                    }
                }
            }

            if (efdr < prevFDR) {
                prevFDR = efdr;
            }

            if (e.isTT()) {
                TT--;
            } else if (e.isTD()) {
                TD--;
                prevTDIndex = i;
            } else if (e.isDD()) {
                DD--;
            } else {
                Logger l = Logger.getLogger(this.getClass().getName());
                l.log(Level.SEVERE, "Something is wrong here!", new Exception(""));
                System.exit(-1);
            }

        }
        return prevFDR;
    }

    /**
     * Takes a sub-set of entries and does the actual FDR-estimation/cutoff.
     * <p>
     * The assumption is that we have a symmetric cross-linker. <br/>Therefore:
     * <p>
     * {@code FP(TT)=TD+DD*(1-TDdb/DDdb)}</p>
     * </p><p>
     * If TCount and DCount are of different size then the false positive
     * estimation gets scaled accordingly.</p
     *
     * @param <T> What type of Information is the FDRImplement to be applied
     * @param group All entries for which the target FDRImplement should be calculated
     * @param fdr the target FDRImplement that is supposed to be returned
     * @param safetyfactor don't report anything if the next higher calculable
 exceeds the target-FDRImplement by the given factor
     * @param results the passing results will be added to this ArrayList
     * @param TCount Size of the target database
     * @param DCount size of the decoy database
     * @param info all informations needed for this sub group
     * @param isSymmetric do we calculate an FDRImplement for a symmetric or an
 asymmetric experiment
     * @param setElementFDR
     * @return
     */
    protected <T extends AbstractFDRElement<T>> void subFDRLocal(ArrayList<T> group, SubGroupFdrInfo info) {

        final boolean addDD = info.useAdditiveFDR;
        int TT = info.TT;
        int TD = info.TD;
        int DD = info.DD;

        double TTSize;
        double TDSize;
        double DDSize;
        double k;
        //@TODO remove
        boolean printOut  =false;

        // all the elements within this group
        int groupSize = group.size();

        // sort them by score
        Collections.sort(group, new Comparator<T>() {

            public int compare(T o1, T o2) {
                return Double.compare(o1.getScore(), o2.getScore());
            }
        });

        int prevTDIndex = groupSize - 1;

        int minTD = 2;
        int minTotal = 20;
        int id5 = (int)(groupSize*0.5);
        int id90 = (int)(groupSize*0.90);
        double minWindow = Math.abs(group.get(id5).getScore() - group.get(id90).getScore()) / 10.0;

        int lastMin = 0;
        T lastMinE = group.get(0);
        double lastMinScore = lastMinE.getScore();
        int lastMax = 0;
        T lastMaxE = lastMinE;
        double lastMaxScore = lastMinScore;

        // count within the window
        int wTT = lastMinE.isTT() ? 1 : 0;
        int wTD = lastMinE.isTD() ? 1 : 0;
        int wDD = lastMinE.isDD() ? 1 : 0;
        for (T e : group) {

            double centerscore = e.getScore();
            double minscore = centerscore - minWindow;
            double maxscore = centerscore + minWindow;

            // find the highest to be included element
            while (lastMax < groupSize - 1 && maxscore > group.get(lastMax + 1).getScore()) {
                lastMax++;
                lastMaxE = group.get(lastMax);
                lastMaxScore = lastMaxE.getScore();

                if (lastMaxE.isTT()) {
                    wTT++;
                } else if (lastMaxE.isTD()) {
                    wTD++;
                } else {
                    wDD++;
                }
                
            }

            // shift to the lowest to be included eleemnt
            while (minscore > lastMinScore
                    && wTD >= minTD
                    && // expand the window until we have more then the required TD
                    wDD * 1.1 < wTD
                    && // also we want some more TD the DD
                    wDD * 1.1 < wTT &&
                    wTT+wTD+wDD >=minTotal) {  // and more TT then DD
                lastMin++;
                if (lastMinE.isTT()) {
                    wTT--;
                } else if (lastMinE.isTD()) {
                    wTD--;
                } else {
                    wDD--;
                }
                lastMinE = group.get(lastMin);
                lastMinScore = lastMinE.getScore();
            }
            
            // do we need to readjust the highest score?
            if (minscore > lastMinScore) {
                double targetMaxScore = lastMinScore - minscore + maxscore;
                while (lastMax < groupSize - 1) {
                    T nextMaxE = group.get(lastMax + 1);
                    double nextScore = nextMaxE.getScore();
                    if (nextScore < targetMaxScore) {
                        lastMax++;
                        lastMaxE = nextMaxE;
                        lastMaxScore = nextScore;
                        if (lastMaxE.isTT()) {
                            wTT++;
                        } else if (lastMaxE.isTD()) {
                            wTD++;
                        } else {
                            wDD++;
                        }
                    } else {
                        break;
                    }
                }
            }
            e.setPEP(addDD ? (wTD + wDD) / (double) wTT : (wTD - wDD) / (double) wTT);
            if (printOut) {
                System.out.println(e.getScore() + ", " + e.getPEP() + ", "  + e.isTT() + ", " + e.isTD() + ", " + e.isDD() + ", " + minscore +", " + maxscore);
            }
        }
    }


    
}
