/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr.gui.components;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class AdditionalFDRGroups extends javax.swing.JPanel {
    final JFrame window = new JFrame("FDR Groups");

    /**
     * Creates new form AdditionalFDRGroups
     */
    public AdditionalFDRGroups() {
        initComponents();
        JButton ok = new JButton("OK");
        window.getContentPane().setLayout(new BoxLayout(window.getContentPane(), BoxLayout.Y_AXIS));
        window.getContentPane().add(this);
        window.add(ok);
        ok.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                window.setVisible(false);
            }
        });
        window.pack();

    }
    
    public void showWindow() {
        window.setVisible(true);
    }
    
    public boolean getFlagSearchID() {
        return ckSearchID.isSelected();
    }

    public boolean getFlagRun() {
        return ckRun.isSelected();
    }
    
    public boolean getFlagModified() {
        return ckModified.isSelected();
    }

    public boolean getFlagAutoValidated() {
        return ckAutoValidated.isSelected();
    }
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        ckAutoValidated = new javax.swing.JCheckBox();
        ckModified = new javax.swing.JCheckBox();
        ckSearchID = new javax.swing.JCheckBox();
        ckRun = new javax.swing.JCheckBox();

        ckAutoValidated.setText("Flag auto-validated");

        ckModified.setText("Flag Modified peptides");

        ckSearchID.setText("Group PSM by Search");

        ckRun.setText("Group PSM by Run");
        ckRun.setEnabled(false);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(ckAutoValidated)
                    .addComponent(ckModified)
                    .addComponent(ckSearchID)
                    .addComponent(ckRun))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(ckAutoValidated)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(ckModified)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(ckSearchID)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(ckRun)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
    }// </editor-fold>//GEN-END:initComponents


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JCheckBox ckAutoValidated;
    private javax.swing.JCheckBox ckModified;
    private javax.swing.JCheckBox ckRun;
    private javax.swing.JCheckBox ckSearchID;
    // End of variables declaration//GEN-END:variables
}