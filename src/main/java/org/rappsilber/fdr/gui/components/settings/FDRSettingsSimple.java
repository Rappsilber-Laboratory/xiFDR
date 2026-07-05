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

import javax.swing.JCheckBox;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SwingUtilities;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import org.rappsilber.fdr.FDRSettings;
import org.rappsilber.fdr.FDRSettingsImpl;
import org.rappsilber.fdr.OfflineFDR;

/**
 *
 * @author lfischer
 */
public class FDRSettingsSimple extends FDRSettingsPanel  {

//    private ArrayList<java.awt.event.ActionListener> m_calc_listener = new ArrayList<ActionListener>();
//    
//    private OfflineFDR.FDRLevel m_optimizeWhat; 
    
    private int m_minPepLength = 0;
    private int m_boostingSteps = 4;
    private double m_psmfdr = 100;
    private double m_pepfdr = 100;
    private double m_protfdr = 100;
    private double m_linkfdr = 100;
    private double m_ppifdr = 100;
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
    public FDRSettingsSimple() {
        initComponents();
        
        ChangeListener max100Listener = new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                SpinnerModel sp = (SpinnerModel) e.getSource();
                if (((Double)sp.getValue()) >100)
                    sp.setValue(100d);
            }
        };        
        spFDR.getModel().addChangeListener(max100Listener);
        spReportFactor.setVisible(false);
        lblReportFactor.setVisible(false);
        this.setAll(new FDRSettingsImpl());
    }
    
    private void doCalc() {
 //       if (ckMaximize.isSelected()) {
            switch ((OfflineFDR.FDRLevel) cbFDRLevel.getSelectedItem()) {
                case PSM: 
                    m_pepfdr = m_linkfdr = m_protfdr = m_ppifdr = 1;
                    break;
                case PEPTIDE_PAIR: 
                    m_psmfdr = m_linkfdr = m_protfdr = m_ppifdr = 1;
                    break;
                case PROTEINGROUP: 
                    m_pepfdr = m_linkfdr = m_psmfdr = m_ppifdr = 1;
                    break;
                case PROTEINGROUPLINK: 
                    m_pepfdr = m_psmfdr = m_protfdr = m_ppifdr = 1;
                    break;
                case PROTEINGROUPPAIR: 
                    m_pepfdr = m_linkfdr = m_protfdr = m_psmfdr = 1;
                    break;
            }
            btnStopBoost.setEnabled(getBoostMode() != FDRSettings.BoostMode.NONE);
            raiseStartCalc(getBoostMode() != FDRSettings.BoostMode.NONE);

    }
    
 
    
    @Override
    public double getPSMFDR() {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PSM)
            return ((Double) spFDR.getValue())/100.0;
        else 
            return m_psmfdr;
    }

    
    @Override
    public double getPeptidePairFDR() {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PEPTIDE_PAIR)
            return ((Double) spFDR.getValue())/100.0;
        else 
            return m_pepfdr;
    }

    @Override
    public double getProteinGroupFDR() {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PROTEINGROUP)
            return ((Double) spFDR.getValue()) /100;
        else 
            return m_protfdr;
    }

    @Override
    public double getProteinGroupLinkFDR() { 
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PROTEINGROUPLINK)
            return ((Double) spFDR.getValue())/100.0;
        else 
            return m_linkfdr;
    }

    @Override
    public double getProteinGroupPairFDR() {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PROTEINGROUPPAIR)
            return ((Double) spFDR.getValue())/100.0;
        else 
            return m_ppifdr;
    }

    @Override
    public void setPSMFDR( Double fdr) {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PSM)
            setValueLater(spFDR, fdr*100.0);
        m_psmfdr = fdr;
    }

    
    @Override
    public void setPeptidePairFDR( Double fdr) {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PEPTIDE_PAIR)
            setValueLater(spFDR, fdr*100.0);
        m_pepfdr = fdr;
    }

    @Override
    public void setProteinGroupFDR( Double fdr) {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PROTEINGROUP)
            setValueLater(spFDR, fdr*100.0);
        m_protfdr = fdr;
    }

    @Override
    public void setProteinGroupLinkFDR( Double fdr) {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PROTEINGROUPLINK)
            setValueLater(spFDR, fdr*100.0);
        m_linkfdr = fdr;
    }

    @Override
    public void setProteinGroupPairFDR( Double fdr) {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PROTEINGROUPPAIR)
            setValueLater(spFDR, fdr*100.0);
        m_ppifdr = fdr;
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
        else return (OfflineFDR.FDRLevel) cbFDRLevel.getSelectedItem();
    }

    @Override
    public void doOptimize(OfflineFDR.FDRLevel level) {
        if (level == null)
            setBoostMode(FDRSettings.BoostMode.NONE);
        else
            cbFDRLevel.getModel().setSelectedItem(level);
    }
    
    public int getBoostingSteps() {
        return m_boostingSteps;
    }
    
    public void setBoostingSteps(int steps) {
        m_boostingSteps = steps;
    }
    
    public double getReportFactor() {
        return (Double) spReportFactor.getValue();
    }    
    
    public void setReportFactor(double factor) {
        spReportFactor.setValue(factor);
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

        lblReportFactor = new javax.swing.JLabel();
        spReportFactor = new javax.swing.JSpinner();
        jLabel1 = new javax.swing.JLabel();
        btnCalc = new javax.swing.JButton();
        lblBoost = new javax.swing.JLabel();
        spFDR = new javax.swing.JSpinner();
        cbFDRLevel = new org.rappsilber.fdr.gui.components.FDRLevelComboBox();
        btnStopBoost = new javax.swing.JButton();
        cmbBoostMode = new javax.swing.JComboBox<>(FDRSettings.BoostMode.values());
        btnBoostHelp = new javax.swing.JButton("?");

        lblReportFactor.setText("Report Factor");
        lblReportFactor.setToolTipText("maximum factor the next step in fdr is permited to exced the target fdr");

        spReportFactor.setModel(new javax.swing.SpinnerNumberModel(1.1d, 1.0d, null, 0.1d));
        spReportFactor.setToolTipText("maximum factor the next step in fdr is permited to exced the target fdr");

        jLabel1.setText("FDR");

        btnCalc.setText("Calculate");
        btnCalc.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnCalcActionPerformed(evt);
            }
        });

        lblBoost.setText("Boost result");
        lblBoost.setToolTipText("Use prefiltering on lower level to boost the results on the chossen level of information");

        spFDR.setModel(new javax.swing.SpinnerNumberModel(100.0d, 0.0d, null, 1.0d));
        spFDR.setToolTipText("What is the maximal acceptable FDR");

        cbFDRLevel.setSelectedIndex(3);
        cbFDRLevel.setToolTipText("To what to apply the FDR");

        btnStopBoost.setText("stop boost");
        btnStopBoost.setEnabled(false);
        btnStopBoost.setMargin(new java.awt.Insets(2, 7, 2, 7));
        btnStopBoost.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnStopBoostActionPerformed(evt);
            }
        });

        cmbBoostMode.setToolTipText("Select the boost mode");
        cmbBoostMode.setSelectedItem(FDRSettings.BoostMode.NONE);

        btnBoostHelp.setMargin(new java.awt.Insets(0, 4, 0, 4));
        btnBoostHelp.setToolTipText("Show boost mode descriptions");
        btnBoostHelp.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnBoostHelpActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(btnStopBoost, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jLabel1)
                    .addComponent(lblBoost)
                    .addComponent(lblReportFactor, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(spReportFactor, javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spFDR)
                            .addGroup(javax.swing.GroupLayout.Alignment.LEADING, layout.createSequentialGroup()
                                .addComponent(cmbBoostMode, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(btnBoostHelp)))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(cbFDRLevel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addGap(0, 0, Short.MAX_VALUE)
                        .addComponent(btnCalc)))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(spFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(cbFDRLevel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(lblBoost)
                    .addComponent(cmbBoostMode, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(btnBoostHelp))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(lblReportFactor)
                    .addComponent(spReportFactor, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 32, Short.MAX_VALUE)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(btnCalc)
                    .addComponent(btnStopBoost))
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
    private javax.swing.JButton btnBoostHelp;
    private javax.swing.JButton btnCalc;
    public javax.swing.JButton btnStopBoost;
    private org.rappsilber.fdr.gui.components.FDRLevelComboBox cbFDRLevel;
    private javax.swing.JComboBox<FDRSettings.BoostMode> cmbBoostMode;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel lblBoost;
    public javax.swing.JLabel lblReportFactor;
    private javax.swing.JSpinner spFDR;
    private javax.swing.JSpinner spReportFactor;
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
    }

    @Override
    public void peppairLocalFDR(Boolean local) {
        this.peppairLocalFDR = local;
    }

    @Override
    public void protLocalFDR(Boolean local) {
        this.protLocalFDR = local;
    }

    @Override
    public void linkLocalFDR(Boolean local) {
        this.linkLocalFDR = local;
    }

    @Override
    public void ppiLocalFDR(Boolean local) {
        this.ppiLocalFDR = local;
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
