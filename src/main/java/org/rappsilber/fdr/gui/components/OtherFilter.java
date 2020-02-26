/*
 * Copyright 2019 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
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
package org.rappsilber.fdr.gui.components;

import java.awt.Container;
import javax.swing.JFrame;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class OtherFilter extends javax.swing.JPanel {

    /**
     * Creates new form OtherFilter
     */
    public OtherFilter() {
        initComponents();
        ckCombineScoreAndDelta.setVisible(false);
        txtMinPepFrags.setHorizontalAlignment(SwingConstants.RIGHT);
        txtDeltaScore.setHorizontalAlignment(SwingConstants.RIGHT);
        txtMinPepCoverage.setHorizontalAlignment(SwingConstants.RIGHT);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        txtMinPepCoverage = new javax.swing.JTextField();
        jLabel1 = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        txtDeltaScore = new javax.swing.JTextField();
        ckCombineScoreAndDelta = new javax.swing.JCheckBox();
        txtMinPepFrags = new javax.swing.JTextField();
        jLabel3 = new javax.swing.JLabel();
        ckTwoStepBoost = new javax.swing.JCheckBox();

        txtMinPepCoverage.setText("0");
        txtMinPepCoverage.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                higherThen1Test(evt);
            }
        });

        jLabel1.setText("Min Peptide Coverage");

        jLabel2.setText("Minimum Delta/Score");

        txtDeltaScore.setText("0");
        txtDeltaScore.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                higherThen1Test(evt);
            }
        });

        ckCombineScoreAndDelta.setText("PSM Score=(Score+Delta)/2");

        txtMinPepFrags.setText("0");
        txtMinPepFrags.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                txtMinPepFragshigherThen1Test(evt);
            }
        });

        jLabel3.setText("Min Peptide Fragments");

        ckTwoStepBoost.setText("Boost Separately");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(ckCombineScoreAndDelta, javax.swing.GroupLayout.DEFAULT_SIZE, 309, Short.MAX_VALUE)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(javax.swing.GroupLayout.Alignment.LEADING, layout.createSequentialGroup()
                                .addComponent(jLabel3)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(txtMinPepFrags))
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(jLabel1)
                                    .addComponent(jLabel2))
                                .addGap(28, 28, 28)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(txtMinPepCoverage)
                                    .addComponent(txtDeltaScore))))
                        .addGap(28, 28, 28))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(ckTwoStepBoost)
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(txtMinPepFrags, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel3))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(txtMinPepCoverage, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel1))
                .addGap(12, 12, 12)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel2)
                    .addComponent(txtDeltaScore, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(ckCombineScoreAndDelta)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(ckTwoStepBoost)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
    }// </editor-fold>//GEN-END:initComponents

    private void higherThen1Test(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_higherThen1Test
        JTextField txt  = (JTextField) (evt.getSource());
        String sValue = txt.getText();
        sValue = sValue.replaceAll("[^0-9+\\-\\.,Ee]", "");
        if (sValue.length()< txt.getText().length())
            txt.setText(sValue);
        Double v = Double.parseDouble(sValue);
        if (v>1) {
            v=v/100.0;
            txt.setText(v.toString());
        }
    }//GEN-LAST:event_higherThen1Test

    private void txtMinPepFragshigherThen1Test(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_txtMinPepFragshigherThen1Test
        JTextField txt  = (JTextField) (evt.getSource());
        String sValue = txt.getText();
        sValue = sValue.replaceAll("[^0-9+\\-\\.Ee]", "");
        if (sValue.length()< txt.getText().length())
            txt.setText(sValue);
    }//GEN-LAST:event_txtMinPepFragshigherThen1Test

    public double getMinDeltaScore() {
        return Double.parseDouble(txtDeltaScore.getText());
    }

    public double getMinPepCoverage() {
        return Double.parseDouble(txtMinPepCoverage.getText());
    }


    public void setMinDeltaScore(Double d) {
        txtDeltaScore.setText(d.toString());
    }

    public void setMinPepCoverage(Double d) {
        txtMinPepCoverage.setText(d.toString());
    }
    
    public void combineScoreAndDelta(boolean c) {
        ckCombineScoreAndDelta.setSelected(c);
    }
    
    public boolean combineScoreAndDelta() {
        return ckCombineScoreAndDelta.isSelected();
    }

    public void twoStepBoost(boolean stepped) {
        ckTwoStepBoost.setSelected(stepped);
    }

    public boolean twoStepBoost() {
        return ckTwoStepBoost.isSelected();
    }

    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JCheckBox ckCombineScoreAndDelta;
    private javax.swing.JCheckBox ckTwoStepBoost;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    public javax.swing.JTextField txtDeltaScore;
    public javax.swing.JTextField txtMinPepCoverage;
    public javax.swing.JTextField txtMinPepFrags;
    // End of variables declaration//GEN-END:variables
}
