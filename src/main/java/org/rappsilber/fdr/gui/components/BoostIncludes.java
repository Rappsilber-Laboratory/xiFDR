/*
 * Copyright 2016 lfischer.
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

import org.rappsilber.fdr.FDRSettings;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import javax.swing.Action;
import javax.swing.event.EventListenerList;

/**
 *
 * @author lfischer
 */
public class BoostIncludes extends javax.swing.JPanel {

//    protected EventListenerList listenerList = new EventListenerList();
    private FDRSettings settings;

    /**
     * Creates new form BoostIgnores
     */
    public BoostIncludes() {
        initComponents();
    }

    public BoostIncludes(FDRSettings settings) {
        this();
        this.settings = settings;
        ckPSM.setSelected(settings.boostPSMs());
        ckPepPair.setSelected(settings.boostPeptidePairs());
        ckProt.setSelected(settings.boostProteins());
        ckLinks.setSelected(settings.boostLinks());
        ckMinFragments.setSelected(settings.boostMinFragments());
        ckPeptideCoverage.setSelected(settings.boostPepCoverage());
        ckDeltaScore.setSelected(settings.boostDeltaScore());
//        ckSubScores.setSelected(settings.boostSubScores());
    }

    /**
     * Adds an <code>ActionListener</code> to the button.
     *
     * @param l the <code>ActionListener</code> to be added
     */
    public void addActionListener(ActionListener l) {
        listenerList.add(ActionListener.class, l);
    }

    /**
     * Removes an <code>ActionListener</code> from the button.
     *
     * @param l the listener to be removed
     */
    public void removeActionListener(ActionListener l) {
        listenerList.remove(ActionListener.class, l);
    }

    /**
     * Notifies all listeners that have registered interest for notification on
     * this event type. The event instance is lazily created using the
     * <code>event</code> parameter.
     *
     * @param event the <code>ActionEvent</code> object
     * @see EventListenerList
     */
    protected void fireActionPerformed(ActionEvent event) {
        // Guaranteed to return a non-null array
        Object[] listeners = listenerList.getListenerList();
        String actionCommand = "Ignores Set";
        ActionEvent e = new ActionEvent(this,
                ActionEvent.ACTION_PERFORMED,
                actionCommand,
                event.getWhen(),
                event.getModifiers());
        for (int i = listeners.length - 2; i >= 0; i -= 2) {
            if (listeners[i] == ActionListener.class) {
                ((ActionListener) listeners[i + 1]).actionPerformed(e);
            }
        }
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        ckPSM = new javax.swing.JCheckBox();
        ckPepPair = new javax.swing.JCheckBox();
        ckProt = new javax.swing.JCheckBox();
        ckLinks = new javax.swing.JCheckBox();
        jLabel1 = new javax.swing.JLabel();
        jButton1 = new javax.swing.JButton();
        ckPeptideCoverage = new javax.swing.JCheckBox();
        ckDeltaScore = new javax.swing.JCheckBox();
        ckMinFragments = new javax.swing.JCheckBox();

        ckPSM.setText("PSM");

        ckPepPair.setText("Peptide Pairs");

        ckProt.setText("Protein Groups");

        ckLinks.setText("Links");

        jLabel1.setText("Use for boosting:");

        jButton1.setText("OK");
        jButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton1ActionPerformed(evt);
            }
        });

        ckPeptideCoverage.setText("Peptide Coverage");

        ckDeltaScore.setText("Delta Score");

        ckMinFragments.setText("Minimum Fragmnets");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addComponent(ckDeltaScore)
                        .addComponent(ckPSM)
                        .addComponent(ckPepPair)
                        .addComponent(ckProt)
                        .addComponent(ckLinks)
                        .addComponent(ckPeptideCoverage)
                        .addComponent(ckMinFragments))
                    .addComponent(jLabel1))
                .addContainerGap(76, Short.MAX_VALUE))
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(jButton1)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(ckPeptideCoverage)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(ckMinFragments)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(ckDeltaScore)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(ckPSM)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(ckPepPair)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(ckProt)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(ckLinks)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(jButton1)
                .addContainerGap())
        );
    }// </editor-fold>//GEN-END:initComponents

    private void jButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton1ActionPerformed
        settings.boostDeltaScore(ckDeltaScore.isSelected());
        settings.boostPepCoverage(ckPeptideCoverage.isSelected());
        settings.boostPSMs(ckPSM.isSelected());
        settings.boostPeptidePairs(ckPepPair.isSelected());
        settings.boostProteins(ckProt.isSelected());
        settings.boostLinks(ckLinks.isSelected());
        settings.boostMinFragments(ckMinFragments.isSelected());
        //settings.boostSubScores(ckSubScores.isSelected());
        fireActionPerformed(evt);
    }//GEN-LAST:event_jButton1ActionPerformed


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JCheckBox ckDeltaScore;
    private javax.swing.JCheckBox ckLinks;
    private javax.swing.JCheckBox ckMinFragments;
    private javax.swing.JCheckBox ckPSM;
    private javax.swing.JCheckBox ckPepPair;
    private javax.swing.JCheckBox ckPeptideCoverage;
    private javax.swing.JCheckBox ckProt;
    private javax.swing.JButton jButton1;
    private javax.swing.JLabel jLabel1;
    // End of variables declaration//GEN-END:variables
}
