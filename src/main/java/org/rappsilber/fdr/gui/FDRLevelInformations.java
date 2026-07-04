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
package org.rappsilber.fdr.gui;

import java.util.Set;
import org.rappsilber.data.csv.CSVRandomAccess;
import org.rappsilber.fdr.result.FDRResultLevel;
import org.rappsilber.fdr.result.SubGroupFdrInfo;
import org.rappsilber.fdr.utils.MiscUtils;

/**
 *
 * @author lfischer
 */
public class FDRLevelInformations extends javax.swing.JFrame {

    /**
     * Creates new form FDRLevelInformations
     */
    public FDRLevelInformations(FDRResultLevel level, String title) {
        initComponents();
        this.setTitle(title);

        CSVRandomAccess csv = new CSVRandomAccess(',', '\"');
        csv.setHeader(new String[]{"Group","Input","TargetFDR","Next FDR","Accepted FDR", "Lower FDR","Passed","Filtered Result","Worst Score", "passed check"});
        csvLevelInfo.setCSV(csv);
        Set<String> ids = level.getGroupIDs();
        for (String fg: ids) {
            String[] rowString = new String[10];
            SubGroupFdrInfo sg = level.getGroup(fg);
            
            
            String formatString = MiscUtils.formatStringForPrettyPrintingRelatedValues(
                    new double[] {
                        sg.higherFDR*100, 
                        sg.lowerFDR*100, 
                        sg.targteFDR*100, 
                        sg.firstPassingFDR*100
                    },1);
            

            rowString[0] = sg.fdrGroup;
            rowString[1] = ""+sg.inputCount;
            
            rowString[2] = sg.targteFDR >=+ 1 ? "unrestricted" : ""+(sg.targteFDR * 100);
            rowString[3] = String.format(formatString + "%% FP",sg.firstPassingFDR*100);
            rowString[4] = String.format("<"+ formatString + "%% H", sg.higherFDR*100);
            rowString[5] = String.format(">"+ formatString + "%% L", sg.lowerFDR*100);
            
            rowString[6] = "" + sg.results.size()  + "(" + sg.resultTT + " TT)";
            rowString[7] = "" + sg.filteredResult.size() ;
            rowString[8] = "" + sg.worstAcceptedScore;

            rowString[9] = (sg.didntPassCheck == null ? "true" :  sg.didntPassCheck);
            
            csv.insertLine(csv.getRowCount(), rowString);
        }
    }

    /**
     * Creates new form FDRLevelInformations
     */
    public FDRLevelInformations(FDRResultLevel[] levels, String[] levelNames, String title) {
        initComponents();
        this.setTitle(title);

        CSVRandomAccess csv = new CSVRandomAccess(',', '\"');
        csv.setHeader(new String[]{"","Group","Input","TargetFDR","Next FDR","Accepted FDR", "Lower FDR","Passed","Filtered Result","Worst Score", "passed check"});
        csvLevelInfo.setCSV(csv);
        for (int l =0; l<levels.length; l++) {
            FDRResultLevel level = levels[l];
            String ln = levelNames[l];
            Set<String> ids = level.getGroupIDs();
            for (String fg: ids) {
                String[] rowString = new String[11];
                SubGroupFdrInfo sg = level.getGroup(fg);


                String formatString = MiscUtils.formatStringForPrettyPrintingRelatedValues(
                        new double[] {
                            sg.higherFDR*100, 
                            sg.lowerFDR*100, 
                            sg.targteFDR*100, 
                            sg.firstPassingFDR*100
                        },1);


                rowString[0] = ln;
                rowString[1] = sg.fdrGroup;
                rowString[2] = ""+sg.inputCount;

                rowString[3] = sg.targteFDR >=+ 1 ? "unrestricted" : ""+(sg.targteFDR * 100);
                rowString[4] = String.format(formatString + "%% FP",sg.firstPassingFDR*100);
                rowString[5] = String.format("<"+ formatString + "%% H", sg.higherFDR*100);
                rowString[6] = String.format(">"+ formatString + "%% L", sg.lowerFDR*100);

                rowString[7] = "" + sg.results.size()  + "(" + sg.resultTT + " TT)";
                rowString[8] = "" + sg.filteredResult.size() ;
                rowString[9] = "" + sg.worstAcceptedScore;

                rowString[10] = (sg.didntPassCheck == null ? "true" :  sg.didntPassCheck);

                csv.insertLine(csv.getRowCount(), rowString);
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

        btnOK = new javax.swing.JButton();
        jLabel1 = new javax.swing.JLabel();
        jTextField1 = new javax.swing.JTextField();
        jLabel2 = new javax.swing.JLabel();
        jTextField2 = new javax.swing.JTextField();
        jLabel3 = new javax.swing.JLabel();
        jTextField3 = new javax.swing.JTextField();
        jTextField4 = new javax.swing.JTextField();
        jLabel4 = new javax.swing.JLabel();
        jLabel5 = new javax.swing.JLabel();
        jTextField5 = new javax.swing.JTextField();
        jLabel6 = new javax.swing.JLabel();
        jTextField6 = new javax.swing.JTextField();
        csvLevelInfo = new org.rappsilber.data.csv.gui.CSVPanel();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);

        btnOK.setText("ok");
        btnOK.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnOKActionPerformed(evt);
            }
        });

        jLabel1.setText("Group:");

        jTextField1.setEditable(false);
        jTextField1.setText("subgroup, that was used to calculate the fdr");
        jTextField1.setBorder(javax.swing.BorderFactory.createEtchedBorder());

        jLabel2.setText("Target FDR:");

        jTextField2.setEditable(false);
        jTextField2.setText("What was the requested FDR");
        jTextField2.setBorder(javax.swing.BorderFactory.createEtchedBorder());

        jLabel3.setText("Next FDR:");

        jTextField3.setEditable(false);
        jTextField3.setText("The next fdr-level, that was found to higher the the targed fdr");
        jTextField3.setBorder(javax.swing.BorderFactory.createEtchedBorder());

        jTextField4.setEditable(false);
        jTextField4.setText("Taking the report factor into account, what was the actually accepted FDR");
        jTextField4.setBorder(javax.swing.BorderFactory.createEtchedBorder());

        jLabel4.setText("Accepted FDR:");

        jLabel5.setText("Passed:");

        jTextField5.setEditable(false);
        jTextField5.setText("How many data passed this fdr");
        jTextField5.setBorder(javax.swing.BorderFactory.createEtchedBorder());

        jLabel6.setText("Filtered :");

        jTextField6.setEditable(false);
        jTextField6.setText("how many of the results passed also the higher level fdr");
        jTextField6.setBorder(javax.swing.BorderFactory.createEtchedBorder());

        csvLevelInfo.setCSV(null);
        csvLevelInfo.setShowLoadPanel(false);
        csvLevelInfo.setShowSavePanel(false);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(csvLevelInfo, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(0, 661, Short.MAX_VALUE)
                        .addComponent(btnOK))
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel1)
                            .addComponent(jLabel2)
                            .addComponent(jLabel3)
                            .addComponent(jLabel4)
                            .addComponent(jLabel5)
                            .addComponent(jLabel6))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jTextField1, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE)
                            .addComponent(jTextField2, javax.swing.GroupLayout.PREFERRED_SIZE, 1, Short.MAX_VALUE)
                            .addComponent(jTextField3, javax.swing.GroupLayout.PREFERRED_SIZE, 1, Short.MAX_VALUE)
                            .addComponent(jTextField4, javax.swing.GroupLayout.PREFERRED_SIZE, 1, Short.MAX_VALUE)
                            .addComponent(jTextField5, javax.swing.GroupLayout.PREFERRED_SIZE, 1, Short.MAX_VALUE)
                            .addComponent(jTextField6, javax.swing.GroupLayout.PREFERRED_SIZE, 1, Short.MAX_VALUE))))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(csvLevelInfo, javax.swing.GroupLayout.DEFAULT_SIZE, 138, Short.MAX_VALUE)
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jTextField1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(jTextField2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel2))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(jTextField3, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel3))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(jTextField4, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel4))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(jTextField5, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel5))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(jTextField6, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel6)))
                    .addComponent(jLabel1))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(btnOK)
                .addContainerGap())
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void btnOKActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnOKActionPerformed
        this.dispose();
    }//GEN-LAST:event_btnOKActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton btnOK;
    private org.rappsilber.data.csv.gui.CSVPanel csvLevelInfo;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JTextField jTextField1;
    private javax.swing.JTextField jTextField2;
    private javax.swing.JTextField jTextField3;
    private javax.swing.JTextField jTextField4;
    private javax.swing.JTextField jTextField5;
    private javax.swing.JTextField jTextField6;
    // End of variables declaration//GEN-END:variables
}
