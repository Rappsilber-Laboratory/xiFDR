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
package org.rappsilber.fdr.gui.components.settings;

import javax.swing.AbstractButton;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JCheckBox;
import javax.swing.JRadioButton;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SwingUtilities;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import org.rappsilber.fdr.FDRSettings;
import org.rappsilber.fdr.FDRSettingsImpl;
import org.rappsilber.fdr.OfflineFDR;
import org.rappsilber.fdr.OfflineFDR.FDRLevel;
import org.rappsilber.fdr.gui.components.FDRSpinnerModel;
import static org.rappsilber.fdr.gui.components.settings.FDRSettingsComplete.setSpinnerModel;
import rappsilber.ms.statistics.utils.UpdateableDouble;

/**
 *
 * @author lfischer
 */
public class FDRSettingsMedium extends FDRSettingsPanel  {

//    private ArrayList<java.awt.event.ActionListener> m_calc_listener = new ArrayList<ActionListener>();
//    
//    private OfflineFDR.FDRLevel m_optimizeWhat; 
    
    private int m_minPepLength = 0;
    private int m_boostingSteps = 4;
    private int m_minTD = DEFAULT_MIN_TD_COUNT;

    
    private boolean m_filterToUniquePSM = true;
    
    private boolean scaleByContectedness = false;
    private Boolean ppiLocalFDR;
    private Boolean linkLocalFDR;
    private Boolean protLocalFDR;
    private Boolean peppairLocalFDR;
    private Boolean psmLocalFDR;
    private double minPeptideCoverageFilter;
    private double minDeltaScoreFilter;
    private boolean combineScoreAndDelta;
    private int minFragments;
    private boolean ignoreValidityChecks = true;
    private boolean psmDirectional;
    private boolean peptidePairDirectional;
    private boolean linkDirectional;
    private boolean ppiDirectional;
    private double reportFactor;
    UpdateableDouble protFDR = new UpdateableDouble(0);

    @Override
    public FDRSettings.BoostMode getBoostMode() {
        FDRSettings.BoostMode mode = (FDRSettings.BoostMode) cmbBoostMode.getSelectedItem();
        return mode != null ? mode : FDRSettings.BoostMode.NONE;
    }

    @Override
    public void setBoostMode(final FDRSettings.BoostMode mode) {
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                cmbBoostMode.setSelectedItem(mode);
            }
        });
    }


    /**
     * Creates new form FDRSettingsComplete
     */
    public FDRSettingsMedium() {
        initComponents();
        
        ChangeListener max100Listener = new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                SpinnerModel sp = (SpinnerModel) e.getSource();
                if (((Double)sp.getValue()) >100)
                    sp.setValue(100d);
            }
        };        
        setSpinnerModel(spPsmFDR, 100);
        setSpinnerModel(spPepFDR, 100);
        setSpinnerModel(spLinkFDR, 5);
        setSpinnerModel(spPPIFdr, 100);
        cbBoostOther.setVisible(false);
        rbBoostOther.setVisible(false);
        this.setAll(new FDRSettingsImpl());
        
    }
    
    private void doCalc() {
 //       if (ckMaximize.isSelected()) {
        protFDR.value = 1;

        this.psmLocalFDR = tsLocalFDR.getSelectionState();
        this.peppairLocalFDR = tsLocalFDR.getSelectionState();
        this.protLocalFDR = tsLocalFDR.getSelectionState();
        this.linkLocalFDR = tsLocalFDR.getSelectionState();
        this.ppiLocalFDR = tsLocalFDR.getSelectionState();

        btnStopBoost.setEnabled(getBoostMode() != FDRSettings.BoostMode.NONE);
        raiseStartCalc(getBoostMode() != FDRSettings.BoostMode.NONE);

    }
    
 
    
    @Override
    public double getPSMFDR() {
        return ((FDRSpinnerModel)(spPsmFDR.getModel())).getFDR();
    }

    
    @Override
    public double getPeptidePairFDR() {
        return ((FDRSpinnerModel)(spPepFDR.getModel())).getFDR();
    }

    @Override
    public double getProteinGroupFDR() {
        return protFDR.value;
    }

    @Override
    public double getProteinGroupLinkFDR() {
        return ((FDRSpinnerModel)(spLinkFDR.getModel())).getFDR();
    }

    @Override
    public double getProteinGroupPairFDR() {
        return ((FDRSpinnerModel)(spPPIFdr.getModel())).getFDR();
    }

    @Override
    public void setPSMFDR( Double fdr) {
        setFDRLater(spPsmFDR,fdr);
    }

    
    @Override
    public void setPeptidePairFDR( Double fdr) {
        setFDRLater(spPepFDR,fdr);
    }

    @Override
    public void setProteinGroupLinkFDR( Double fdr) {
        setFDRLater(spLinkFDR,fdr);
    }

    @Override
    public void setProteinGroupPairFDR( Double fdr) {
        setFDRLater(spPPIFdr,fdr);
    }

    @Override
    public void setProteinGroupFDR( Double fdr) {
        protFDR.value = fdr;
    }

    @Override
    public int getMinProteinPepCount() {
        return 1;
    }

    @Override
    public int getMinLinkPepCount() {
        return 1;
    }

    @Override
    public int getMinPPIPepCount() {
        return 1;
    }

    @Override
    public void setMinProteinPepCount(Integer minPep) {
    }

    @Override
    public void setMinLinkPepCount(Integer minPep) {
    }

    @Override
    public void setMinPPIPepCount(Integer minPep) {
    }

    
    

    @Override
    public int getMinPeptideLength() {
        return m_minPepLength;
    }

    @Override
    public void setMinPeptideLength(Integer minLength) {
        m_minPepLength = minLength;
    }

    
    @Override
    public int getMaxLinkAmbiguity() {
        return 0;
    }

    @Override
    public void setMaxLinkAmbiguity(Integer maxAmbiguity) {
        
    }

    @Override
    public int getMaxProteinAmbiguity() {
        return 0;
    }

    @Override
    public void setMaxProteinAmbiguity(Integer maxAmbiguity) {
    }
    
    @Override
    public boolean isPSMDirectional() {
        return this.psmDirectional;
    }
    
    @Override
    public boolean isPeptidePairDirectional() {
        return this.peptidePairDirectional;
    }

    @Override
    public boolean isLinkDirectional() {
        return this.linkDirectional;
    }

    @Override
    public boolean isPPIDirectional() {
        return this.ppiDirectional;
    }

    @Override
    public void setPSMDirectional(boolean directional) {
        psmDirectional=directional;
    }
    
    @Override
    public void setPeptidePairDirectional(boolean directional) {
        peptidePairDirectional = directional;
    }

    @Override
    public void setLinkDirectional(boolean directional) {
        linkDirectional = directional;
    }

    @Override
    public void setPPIDirectional(boolean directional) {
        ppiDirectional = directional;
    }
    
    
    @Override
    public OfflineFDR.FDRLevel doOptimize() {
        if (getBoostMode() == FDRSettings.BoostMode.NONE)
            return null;
        if (rbBoostPeptidePairs.isSelected()) {
            return OfflineFDR.FDRLevel.PEPTIDE_PAIR;
        }
        if (rbBoostLinks.isSelected()) {
            return OfflineFDR.FDRLevel.PROTEINGROUPLINK;
        }
        if (rbBoostPPI.isSelected()) {
            return OfflineFDR.FDRLevel.PROTEINGROUPPAIR;
        } 
        return (OfflineFDR.FDRLevel) cbBoostOther.getSelectedItem();
    }

    @Override
    public void doOptimize(OfflineFDR.FDRLevel level) {
        if (level == null) {
            setBoostMode(FDRSettings.BoostMode.NONE);
        } else {
            switch(level) {
                case PEPTIDE_PAIR:
                    bgBoost.setSelected(rbBoostPeptidePairs.getModel(), true);
                    cbBoostOther.setVisible(false);
                    rbBoostOther.setVisible(false);
                    break;
                case PROTEINGROUPLINK:
                    bgBoost.setSelected(rbBoostLinks.getModel(), true);
                    cbBoostOther.setVisible(false);
                    rbBoostOther.setVisible(false);
                    break;
                case PROTEINGROUPPAIR:
                    bgBoost.setSelected(rbBoostPPI.getModel(), true);
                    cbBoostOther.setVisible(false);
                    rbBoostOther.setVisible(false);
                    break;
                default:
                    bgBoost.getSelection().setSelected(false);
                    cbBoostOther.setSelectedItem(level);
                    cbBoostOther.setVisible(true);
                    rbBoostOther.setVisible(true);
                    bgBoost.setSelected(rbBoostOther.getModel(), true);
            }
        }
    }
    
    public int getBoostingSteps() {
        return m_boostingSteps;
    }
    
    public void setBoostingSteps(int steps) {
        m_boostingSteps = steps;
    }
    
    public double getReportFactor() {
        return this.reportFactor;
    }    
    
    public void setReportFactor(double factor) {
        this.reportFactor = factor;
    }       

    
    public boolean filterToUniquePSM() {
        return m_filterToUniquePSM;
    }

    public void setFilterToUniquePSM(boolean filterToUniquePSM) {
        m_filterToUniquePSM = filterToUniquePSM;
    }
    
    public boolean combineScoreAndDelta() {
        return this.combineScoreAndDelta;
    }
    public void combineScoreAndDelta(boolean c) {
        this.combineScoreAndDelta = c;
    }

    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        bgBoost = new javax.swing.ButtonGroup();
        cmbBoostMode = new javax.swing.JComboBox<>(FDRSettings.BoostMode.values());
        btnCalc = new javax.swing.JButton();
        btnStopBoost = new javax.swing.JButton();
        jLabel5 = new javax.swing.JLabel();
        spPsmFDR = new javax.swing.JSpinner();
        jLabel9 = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        spPepFDR = new javax.swing.JSpinner();
        jLabel3 = new javax.swing.JLabel();
        spLinkFDR = new javax.swing.JSpinner();
        jLabel4 = new javax.swing.JLabel();
        spPPIFdr = new javax.swing.JSpinner();
        rbBoostLinks = new javax.swing.JRadioButton();
        rbBoostPPI = new javax.swing.JRadioButton();
        rbBoostPeptidePairs = new javax.swing.JRadioButton();
        cbBoostOther = new javax.swing.JComboBox<>();
        rbBoostOther = new javax.swing.JRadioButton();
        tsLocalFDR = new org.rappsilber.gui.components.TriStateCheckBox();
        jLabel6 = new javax.swing.JLabel();
        btnBoostHelp = new javax.swing.JButton("?");

        cmbBoostMode.setToolTipText("Select what the boost should optimise");
        cmbBoostMode.setSelectedItem(FDRSettings.BoostMode.NONE);

        btnBoostHelp.setMargin(new java.awt.Insets(0, 4, 0, 4));
        btnBoostHelp.setToolTipText("Show boost mode descriptions");
        btnBoostHelp.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnBoostHelpActionPerformed(evt);
            }
        });

        btnCalc.setText("Calculate");
        btnCalc.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnCalcActionPerformed(evt);
            }
        });

        btnStopBoost.setText("stop boost");
        btnStopBoost.setEnabled(false);
        btnStopBoost.setMargin(new java.awt.Insets(2, 7, 2, 7));
        btnStopBoost.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnStopBoostActionPerformed(evt);
            }
        });

        jLabel5.setText("PSM");

        spPsmFDR.setToolTipText("FDR value accepted for PSMs");

        jLabel9.setText("Max FDRs");

        jLabel2.setText("Peptide Pair");

        spPepFDR.setToolTipText("FDR value accepted for Peptide Pairs (including linksite within the peptide)");

        jLabel3.setText("Residue Pairs");

        spLinkFDR.setToolTipText("FDR value accepted for Links (Protein-group-links)");

        jLabel4.setText("Protein Pairs");

        spPPIFdr.setToolTipText("FDR value accepted for protein-pairs");

        bgBoost.add(rbBoostLinks);
        rbBoostLinks.setSelected(true);
        rbBoostLinks.setToolTipText("try to maximize the number of residue pairs");
        rbBoostLinks.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                setBoostingLevel(evt);
            }
        });

        bgBoost.add(rbBoostPPI);
        rbBoostPPI.setToolTipText("try to maximize the number of protein pairs");
        rbBoostPPI.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                setBoostingLevel(evt);
            }
        });

        bgBoost.add(rbBoostPeptidePairs);
        rbBoostPeptidePairs.setToolTipText("try to maximize the number of peptide pairs");
        rbBoostPeptidePairs.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                setBoostingLevel(evt);
            }
        });

        cmbBoostMode.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                boolean boost = getBoostMode() != FDRSettings.BoostMode.NONE;
                rbBoostLinks.setEnabled(boost);
                rbBoostPPI.setEnabled(boost);
                rbBoostPeptidePairs.setEnabled(boost);
                rbBoostOther.setEnabled(boost);
                cbBoostOther.setEnabled(boost);
            }
        });

        cbBoostOther.setModel(new DefaultComboBoxModel<OfflineFDR.FDRLevel>(new OfflineFDR.FDRLevel[]{OfflineFDR.FDRLevel.PSM, OfflineFDR.FDRLevel.PROTEINGROUP}));
        cbBoostOther.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cbBoostOtherActionPerformed(evt);
            }
        });

        bgBoost.add(rbBoostOther);
        rbBoostOther.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbBoostOthersetBoostingLevel(evt);
            }
        });

        tsLocalFDR.setToolTipText("calculate Local FDR(=PEP)  (square); calculate and filter by Local FDR (tick); or do not calculate Local FDR(empty)");
        tsLocalFDR.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                tsLocalFDRActionPerformed(evt);
            }
        });

        jLabel6.setText("Local FDR");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel5)
                    .addComponent(jLabel2)
                    .addComponent(jLabel3)
                    .addComponent(jLabel4)
                    .addComponent(jLabel6))
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(16, 16, 16)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                            .addComponent(spPPIFdr, javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spLinkFDR, javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spPepFDR, javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spPsmFDR, javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel9, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.PREFERRED_SIZE, 92, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(18, 18, 18)
                        .addComponent(tsLocalFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.CENTER)
                    .addComponent(rbBoostPeptidePairs)
                    .addComponent(rbBoostLinks)
                    .addComponent(rbBoostPPI)
                    .addComponent(rbBoostOther))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.CENTER)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(cmbBoostMode, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(btnBoostHelp))
                    .addComponent(cbBoostOther, javax.swing.GroupLayout.PREFERRED_SIZE, 95, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 86, Short.MAX_VALUE)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(btnCalc)
                    .addComponent(btnStopBoost))
                .addContainerGap())
        );

        layout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {btnCalc, btnStopBoost});

        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel9)
                    .addComponent(cmbBoostMode, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(btnBoostHelp))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(spPsmFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel5))
                        .addGap(11, 11, 11)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                .addComponent(spPepFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(jLabel2))
                            .addComponent(rbBoostPeptidePairs))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                .addComponent(spLinkFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(jLabel3))
                            .addComponent(rbBoostLinks, javax.swing.GroupLayout.Alignment.TRAILING))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                .addComponent(spPPIFdr, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(jLabel4))
                            .addComponent(rbBoostPPI, javax.swing.GroupLayout.Alignment.TRAILING))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(cbBoostOther, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                .addComponent(jLabel6)
                                .addComponent(tsLocalFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))))
                    .addComponent(rbBoostOther))
                .addGap(0, 57, Short.MAX_VALUE))
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(btnStopBoost)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(btnCalc)
                .addContainerGap())
        );
    }// </editor-fold>//GEN-END:initComponents

    private void btnCalcActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnCalcActionPerformed
        // TODO add your handling code here:
        doCalc();
        
    }//GEN-LAST:event_btnCalcActionPerformed

    private void btnStopBoostActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnStopBoostActionPerformed
        raiseStopMaximizing();
    }//GEN-LAST:event_btnStopBoostActionPerformed


    private void setBoostingLevel(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_setBoostingLevel
        if (evt.getSource() instanceof JRadioButton) {
            cbBoostOther.setVisible(false);
            rbBoostOther.setVisible(false);
        }
    }//GEN-LAST:event_setBoostingLevel

    private void rbBoostOthersetBoostingLevel(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbBoostOthersetBoostingLevel
        // TODO add your handling code here:
    }//GEN-LAST:event_rbBoostOthersetBoostingLevel

    private void cbBoostOtherActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cbBoostOtherActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_cbBoostOtherActionPerformed

    private void tsLocalFDRActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_tsLocalFDRActionPerformed


    }//GEN-LAST:event_tsLocalFDRActionPerformed

    private void btnBoostHelpActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnBoostHelpActionPerformed
        javax.swing.JPopupMenu popup = new javax.swing.JPopupMenu();
        StringBuilder sb = new StringBuilder("<html>");
        for (FDRSettings.BoostMode mode : FDRSettings.BoostMode.values()) {
            sb.append("<b>").append(mode.getName()).append("</b>: ")
              .append(mode.getDescription()).append("<br>");
        }
        sb.append("</html>");
        popup.add(new javax.swing.JLabel(sb.toString()));
        popup.show(btnBoostHelp, 0, btnBoostHelp.getHeight());
    }//GEN-LAST:event_btnBoostHelpActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.ButtonGroup bgBoost;
    private javax.swing.JButton btnBoostHelp;
    private javax.swing.JButton btnCalc;
    public javax.swing.JButton btnStopBoost;
    private javax.swing.JComboBox<FDRLevel> cbBoostOther;
    private javax.swing.JComboBox<FDRSettings.BoostMode> cmbBoostMode;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel9;
    private javax.swing.JRadioButton rbBoostLinks;
    private javax.swing.JRadioButton rbBoostOther;
    private javax.swing.JRadioButton rbBoostPPI;
    private javax.swing.JRadioButton rbBoostPeptidePairs;
    private javax.swing.JSpinner spLinkFDR;
    private javax.swing.JSpinner spPPIFdr;
    private javax.swing.JSpinner spPepFDR;
    private javax.swing.JSpinner spPsmFDR;
    private org.rappsilber.gui.components.TriStateCheckBox tsLocalFDR;
    // End of variables declaration//GEN-END:variables

    @Override
    public void setEnabled(boolean e) {
        super.setEnabled(e);
        btnCalc.setEnabled(e);
        if (e) 
            btnStopBoost.setEnabled(false);
    }

//    @Override
//    public boolean isEnabeled(boolean e) {
//    }

    @Override
    public void setMinTD(Integer c) {
        m_minTD = c;
    }

    @Override
    public int getMinTD() {
        return m_minTD;
    }

    @Override
    public Boolean psmLocalFDR() {
        return this.psmLocalFDR;
    }

    @Override
    public Boolean peppairLocalFDR() {
        return this.peppairLocalFDR;
    }

    @Override
    public Boolean protLocalFDR() {
        return this.protLocalFDR;
    }

    @Override
    public Boolean linkLocalFDR() {
        return this.linkLocalFDR;
    }

    @Override
    public Boolean ppiLocalFDR() {
        return this.ppiLocalFDR;
    }

    @Override
    public void psmLocalFDR(Boolean local) {
        this.psmLocalFDR = local;
        setLocalFDR();
    }

    protected void setLocalFDR() {
        Boolean hasTrue = false;
        Boolean hasFalse = false;
        if (this.peppairLocalFDR!= null) {
            hasTrue |= this.peppairLocalFDR;
            hasFalse &= !this.peppairLocalFDR;
        }
        if (this.psmLocalFDR!= null) {
            hasTrue |= this.psmLocalFDR;
            hasFalse &= !this.psmLocalFDR;
        }
        if (this.protLocalFDR!= null) {
            hasTrue |= this.protLocalFDR;
            hasFalse &= !this.protLocalFDR;
        }
        if (this.linkLocalFDR!= null) {
            hasTrue |= this.linkLocalFDR;
            hasFalse &= !this.linkLocalFDR;
        }
        if (this.ppiLocalFDR!= null) {
            hasTrue |= this.ppiLocalFDR;
            hasFalse &= !this.ppiLocalFDR;
        }
        if (hasTrue) {
            this.tsLocalFDR.setSelectionState(true);
        } else if (hasFalse) {
            this.tsLocalFDR.setSelectionState(false);
        } else {
            this.tsLocalFDR.setSelectionState(null);
        }
    }

    @Override
    public void peppairLocalFDR(Boolean local) {
        this.peppairLocalFDR = local;
        setLocalFDR();
    }

    @Override
    public void protLocalFDR(Boolean local) {
        this.protLocalFDR = local;
        setLocalFDR();
    }

    @Override
    public void linkLocalFDR(Boolean local) {
        this.linkLocalFDR = local;
        setLocalFDR();
    }

    @Override
    public void ppiLocalFDR(Boolean local) {
        this.ppiLocalFDR = local;
        setLocalFDR();
    }    

    
    @Override
    public double getMinPeptideCoverageFilter() {
        return this.minPeptideCoverageFilter;
    }
    @Override
    public void setMinPeptideCoverageFilter(double d) {
        this.minPeptideCoverageFilter = d;
    }

    @Override
    public double getMinDeltaScoreFilter() {
        return this.minDeltaScoreFilter;
    }
    
    @Override
    public void setMinDeltaScoreFilter(double d) {
        this.minDeltaScoreFilter = d;
    }

    @Override
    public int getMinPeptideFragmentsFilter() {
        return this.minFragments;
    }

    @Override
    public void setMinPeptideFragmentsFilter(int frags) {
        this.minFragments = frags;
    }

    @Override
    public boolean ignoreValidityChecks() {
        return this.ignoreValidityChecks;
    }

    @Override
    public void ignoreValidityChecks(boolean ignore) {
        this.ignoreValidityChecks = ignore;
    }
    
    
    
}
