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
package org.rappsilber.fdr.utils;

import java.sql.Connection;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import org.rappsilber.utils.RArrayUtils;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class ClearValidation extends javax.swing.JFrame {

    /**
     * Creates new form ClearValidation
     */
    public ClearValidation() {
        initComponents();
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        getSearch = new rappsilber.gui.components.db.GetSearch();
        btnResValidateClearAll = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        btnResValidateClearAll.setText("clear all validation");
        btnResValidateClearAll.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnResValidateClearAllActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap(277, Short.MAX_VALUE)
                .addComponent(btnResValidateClearAll)
                .addContainerGap())
            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addComponent(getSearch, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 453, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap(281, Short.MAX_VALUE)
                .addComponent(btnResValidateClearAll)
                .addContainerGap())
            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(layout.createSequentialGroup()
                    .addComponent(getSearch, javax.swing.GroupLayout.DEFAULT_SIZE, 270, Short.MAX_VALUE)
                    .addGap(48, 48, 48)))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void btnResValidateClearAllActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnResValidateClearAllActionPerformed
        try {
            int[] ids = getSearch.getSelectedSearchIds();
            if (JOptionPane.showConfirmDialog(this, "Clear validation on searches " + RArrayUtils.toString(ids, ",") + "?", "Clear Validation?", JOptionPane.OK_CANCEL_OPTION) == JOptionPane.OK_OPTION) {
                String query = "UPDATE spectrum_match set validated=null where search_id in (" + RArrayUtils.toString(ids, ",") + ") and validated is not null";
                Connection c = getSearch.getConnection();
                c.setAutoCommit(false);
                Statement s =  c.createStatement();
                int ret = s.executeUpdate(query);
                if (ret >0) {
                    if (JOptionPane.showConfirmDialog(this, "This will clear the validation state on " + ret + " matches ?", "Confirm?", JOptionPane.OK_CANCEL_OPTION) == JOptionPane.OK_OPTION) {
                        c.commit();
                    } else {
                        c.rollback();
                    }
                }  else {
                    JOptionPane.showMessageDialog(this, "Search had nothing validated", "Nothing changed?", JOptionPane.INFORMATION_MESSAGE);
                }

                s.close();
                c.close();
            }
        } catch (SQLException ex) {
            Logger.getLogger(ClearValidation.class.getName()).log(Level.SEVERE, null, ex);
        }
    }//GEN-LAST:event_btnResValidateClearAllActionPerformed

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(ClearValidation.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(ClearValidation.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(ClearValidation.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(ClearValidation.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        shoWDialog(JFrame.EXIT_ON_CLOSE);
    }

    public static  ClearValidation shoWDialog() {
        return shoWDialog(JFrame.DISPOSE_ON_CLOSE);
    }
    
    public static  ClearValidation shoWDialog(int closeOperation) {
        ClearValidation ret = new ClearValidation();
        ret.setDefaultCloseOperation(closeOperation);
        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new ClearValidation().setVisible(true);
            }
        });
        return ret;
    }
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton btnResValidateClearAll;
    public rappsilber.gui.components.db.GetSearch getSearch;
    // End of variables declaration//GEN-END:variables
}
