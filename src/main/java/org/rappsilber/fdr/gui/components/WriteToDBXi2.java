/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr.gui.components;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.UUID;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JOptionPane;
import org.rappsilber.config.LocalProperties;
import org.rappsilber.fdr.DB2inFDR;
import org.rappsilber.fdr.DBinFDR;
import org.rappsilber.fdr.gui.FDRGUI;
import rappsilber.gui.components.db.DatabaseProvider;

/**
 *
 * @author lfischer
 */
public class WriteToDBXi2 extends javax.swing.JPanel {

    FDRGUI m_gui;
    
    class DBUser {
        long id;
        String username;

        public DBUser(long id, String username) {
            this.id = id;
            this.username = username;
        }

        @Override
        public String toString() {
            return username; //To change body of generated methods, choose Tools | Templates.
        }
        
    }
    
    /**
     * Creates new form WriteToDB
     */
    public WriteToDBXi2() {
        initComponents();
        
        this.txtNotes.setText("Add your notes here \n" + DB2inFDR.summary_marker_long);
        ArrayList<String> scores = new ArrayList<>();
        for (DB2inFDR.FDRScores s  : DB2inFDR.FDRScores.values()) {
            scores.add(s.getName());
        }
        cmbMainScore.setModel(new javax.swing.DefaultComboBoxModel<>(scores.toArray(new String[scores.size()])));
        cmbMainScore.setSelectedIndex(2);
        
    }

    public void setDBProvider(DatabaseProvider provider) throws SQLException {
        ArrayList<DBUser> users = new ArrayList<>();
        String user_querry = "select id, username from core_user where is_active OR username = 'xiFDR' order by username <> 'xiFDR', username";
        Statement st  = (((GetSearch) provider).getConnectionAdmin()).createStatement();
        String default_user = LocalProperties.getProperty("xi2db_default_user", "xiFDR");
        if (cbUser.getSelectedItem()!=null)  {
            default_user = ((DBUser)cbUser.getSelectedItem()).username;
        }
        ResultSet rs = st.executeQuery(user_querry);
        int default_user_id = -1;
        while (rs.next()) {
            DBUser u = new DBUser(rs.getLong(1), rs.getString(2));
            if (u.username.contentEquals(default_user))
                default_user_id = users.size();
            users.add(u);
        }
        this.cbUser.setModel(new DefaultComboBoxModel(users.toArray(new DBUser[users.size()])));
        if (default_user_id>=0) {
            this.cbUser.setSelectedIndex(default_user_id);
        }
        rs.close();
        st.close();
        
    }
    
    private void writeDB() {
        final DB2inFDR ofdr = (DB2inFDR) m_gui.getFdr();
        
        
        m_gui.setEnableCalc(false);
        this.setEnableWrite(false);
        m_gui.setEnableRead(false);
        final boolean within = rbWithin.isSelected() || rbAll.isSelected();
        final boolean between = rbBetween.isSelected()  || rbAll.isSelected();
        final long userid = ((DBUser)this.cbUser.getSelectedItem()).id;
        

        Runnable runnable = new Runnable() {
            public void run() {
                try {
                    m_gui.setStatus("Writing to db: " + ofdr.summaryString(m_gui.getResult()));
                    ofdr.setMainScore(cmbMainScore.getSelectedItem().toString());
                    ofdr.writeResult(txtName.getText(), txtNotes.getText(), (UUID)null, m_gui.getResult(),within, between, userid);
                    m_gui.setStatus("finished writing: " + ofdr.summaryString(m_gui.getResult()));

                } catch (SQLException sex) {
                    Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, sex);
                    m_gui.setStatus(sex.toString());
                    SQLException ex = sex;
                    int chain = 1;
                    while ((ex = ex.getNextException()) != null && chain++ <10) {
                        Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "next exception", ex);
                    }

                } catch (Exception e) {
                    Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, "Error writing Results", e);
                    m_gui.setStatus(e.toString());
                }
                m_gui.setEnableCalc(true);
                WriteToDBXi2.this.setEnableWrite(true);
                m_gui.setEnableRead(true);
            }
        };

        new Thread(runnable).start();


    }
 
    
    public void setEnableWrite(final boolean enable) {
        Runnable setModel = new Runnable() {
            public void run() {
                btnWriteDB.setEnabled(enable);
            }
        };
        javax.swing.SwingUtilities.invokeLater(setModel);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        buttonGroup1 = new javax.swing.ButtonGroup();
        btnWriteDB = new javax.swing.JButton();
        rbAll = new javax.swing.JRadioButton();
        rbWithin = new javax.swing.JRadioButton();
        rbBetween = new javax.swing.JRadioButton();
        jLabel1 = new javax.swing.JLabel();
        txtName = new javax.swing.JTextField();
        jLabel2 = new javax.swing.JLabel();
        jScrollPane1 = new javax.swing.JScrollPane();
        txtNotes = new javax.swing.JTextArea();
        cbUser = new javax.swing.JComboBox<>();
        jLabel3 = new javax.swing.JLabel();
        cmbMainScore = new javax.swing.JComboBox<>();
        jLabel4 = new javax.swing.JLabel();

        btnWriteDB.setText("Write");
        btnWriteDB.setEnabled(false);
        btnWriteDB.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnWriteDBActionPerformed(evt);
            }
        });

        buttonGroup1.add(rbAll);
        rbAll.setSelected(true);
        rbAll.setText("All");

        buttonGroup1.add(rbWithin);
        rbWithin.setText("Within");

        buttonGroup1.add(rbBetween);
        rbBetween.setText("Between");
        rbBetween.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbBetweenActionPerformed(evt);
            }
        });

        jLabel1.setText("Name");

        jLabel2.setText("Notes");

        txtNotes.setColumns(20);
        txtNotes.setRows(5);
        jScrollPane1.setViewportView(txtNotes);

        cbUser.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cbUserActionPerformed(evt);
            }
        });

        jLabel3.setText("As User");

        cmbMainScore.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "CSM FDR", "CSM PEP", "CSM Score", "PeptidePair FDR", "PeptidePair PEP", "PeptdePair Score", "ResiduePair FDR", "ResiduePair PEP", "ResiduePair Score", "ProteinPair FDR", "ProteinPair PEP", "ProteinPair Score", "Protein1 FDR", "Protein1 PEP", "Protein1 Score", "Protein2 FDR", "Protein2 PEP", "Protein2 Score" }));
        cmbMainScore.setSelectedIndex(2);

        jLabel4.setText("Main Score");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addComponent(jLabel1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(txtName))
            .addComponent(jScrollPane1)
            .addGroup(layout.createSequentialGroup()
                .addComponent(jLabel2)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(jLabel4)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(cmbMainScore, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addGap(0, 0, Short.MAX_VALUE)
                .addComponent(btnWriteDB))
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addComponent(jLabel3)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(cbUser, javax.swing.GroupLayout.PREFERRED_SIZE, 126, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 19, Short.MAX_VALUE)
                .addComponent(rbAll)
                .addGap(18, 18, 18)
                .addComponent(rbWithin)
                .addGap(18, 18, 18)
                .addComponent(rbBetween))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(txtName, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(jLabel2)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(cmbMainScore, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jLabel4)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 130, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(rbAll)
                    .addComponent(rbWithin)
                    .addComponent(rbBetween)
                    .addComponent(cbUser, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel3))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(btnWriteDB)
                .addContainerGap())
        );
    }// </editor-fold>//GEN-END:initComponents

    private void btnWriteDBActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnWriteDBActionPerformed

        writeDB();

    }//GEN-LAST:event_btnWriteDBActionPerformed

    private void rbBetweenActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbBetweenActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_rbBetweenActionPerformed

    private void cbUserActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cbUserActionPerformed
        LocalProperties.setProperty("xi2db_default_user", ((DBUser)cbUser.getSelectedItem()).username);
    }//GEN-LAST:event_cbUserActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton btnWriteDB;
    private javax.swing.ButtonGroup buttonGroup1;
    private javax.swing.JComboBox<String> cbUser;
    private javax.swing.JComboBox<String> cmbMainScore;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JRadioButton rbAll;
    private javax.swing.JRadioButton rbBetween;
    private javax.swing.JRadioButton rbWithin;
    private javax.swing.JTextField txtName;
    private javax.swing.JTextArea txtNotes;
    // End of variables declaration//GEN-END:variables

    public void setFdrgui(FDRGUI gui) {
        this.m_gui= gui;
    }
}
