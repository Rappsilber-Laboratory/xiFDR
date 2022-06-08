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

import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.EventObject;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Properties;
import java.util.Set;
import java.util.logging.Filter;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;
import javax.swing.DefaultCellEditor;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.table.TableModel;
import org.rappsilber.config.LocalProperties;
import org.rappsilber.data.csv.CsvParser;
import org.rappsilber.data.csv.condition.CsvCondition;
import org.rappsilber.fdr.CSVinFDR;
import org.rappsilber.fdr.DBinFDR;
import org.rappsilber.fdr.FDRSettings;
import org.rappsilber.fdr.FDRSettingsImpl;
import org.rappsilber.fdr.result.FDRResult;
import org.rappsilber.fdr.result.FDRResultLevel;
import org.rappsilber.fdr.MZIdentXLFDR;
import org.rappsilber.fdr.OfflineFDR;
import org.rappsilber.fdr.XiCSVinFDR;
import org.rappsilber.fdr.XiInFDR;
import org.rappsilber.fdr.result.SubGroupFdrInfo;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.fdr.entities.PeptidePair;
import org.rappsilber.fdr.entities.ProteinGroupLink;
import org.rappsilber.fdr.entities.ProteinGroupPair;
import org.rappsilber.fdr.gui.components.FDRSettingsPanel;
import org.rappsilber.fdr.utils.MZIdentMLExport;
import org.rappsilber.fdr.utils.MaximisingStatus;
import org.rappsilber.fdr.utils.MaximizingUpdate;
import org.rappsilber.fdr.utils.MiscUtils;
import org.rappsilber.fdr.utils.MZIdentMLOwner;
import org.rappsilber.gui.components.AutoAddTableModelListener;
import org.rappsilber.gui.GenericTextPopUpMenu;
import org.rappsilber.gui.components.JoinedThreadedTextOuput;
import org.rappsilber.utils.RArrayUtils;
import org.rappsilber.gui.logging.JTextAreaHandle;
import org.rappsilber.peaklist.MgfStyleTitleParser;
import org.rappsilber.utils.UpdateableInteger;
import org.rappsilber.utils.XiFDRUtils;
import rappsilber.config.RunConfigFile;

/**
 *
 * @author lfischer
 */
public class FDRGUI extends javax.swing.JFrame {

    private boolean stopMaximizing = false;

    private JoinedThreadedTextOuput m_status = new JoinedThreadedTextOuput();

    private JTextAreaHandle loggingOutput;

    private FDRResult m_result;

    private FDRSettingsPanel fdrSettings;

    private OfflineFDR m_fdr;
//    /** denotes required fields in the CSV */
//    public String missingColumn = "!! !MISSING! !!";
//    public String optionalColumn = "  OPTIONAL   ";
//    public String[] csvColumns = new String[]{missingColumn};
//    public String[] csvColumnsOptional = new String[]{optionalColumn};

    private XiFDRUtils.FDRCSVOUT csvout;
    /**
     * Creates new form
     */
    public FDRGUI() {
        initComponents();
        ckDBSize.setVisible(false);
//        cmbPeakListFormat.setVisible(false);
//        lblPeaklistExtension.setVisible(false);
//        cmbMzMLScan2ID.setVisible(false);
//        lblMzMLScan2ID.setVisible(false);
        
        

        // setup the logging to the text-field
        loggingOutput = new JTextAreaHandle(txtLog);
        loggingOutput.setFilter(new Filter() {

            public boolean isLoggable(LogRecord record) {
                return true;
            }
        });
        loggingOutput.setLevel(Level.ALL);

        Handler tabRiser = riseLoggignTabOnError();
        
        // make sure we display one of the possible fdrsettings-panel
        changeFDRSettings(null);

        //Logger.getLogger("rappsilber").addHandler(loggingOutput);
        Logger.getLogger("org.rappsilber").setLevel(Level.ALL);
        Logger.getLogger("org.rappsilber").addHandler(loggingOutput);
        Logger.getLogger("org.rappsilber").addHandler(tabRiser);

        this.setTitle("");
        txtXiFDRVersion.setText(OfflineFDR.getXiFDRVersion().toString());
        //tpInput.removeTabAt(0);
        tblPepLength.getModel().addTableModelListener(new AutoAddTableModelListener());

        // setup a context-menu with copy past-functions for all test-fields
        GenericTextPopUpMenu gtpm = new GenericTextPopUpMenu();
        gtpm.installContextMenu(this);
        gtpm.installContextMenu(fdrSettingsComplete);

        this.toFront();

        // register, where we want to have status messages printed
        m_status.addLoggerOutput(Logger.getLogger(this.getClass().getName()));
        m_status.addTextOutput(txtStatus);

//        fdrSettings=fdrSettingsSimple;
        // make spinners select all test on gaining the focus
        // it's a bit of a dirty hack - but the best I could find
        ArrayList<Container> con = new ArrayList<Container>();
        HashSet<Container> spinSet = new HashSet<Container>();
        ArrayList<Container> spins = new ArrayList<Container>();
        con.add(this);
        con.add(fdrSettingsComplete);
        con.add(fdrSettingsMedium);
        con.add(fdrSettingsSimple);
        makeSpinnersSelectTextOnEnter(con, spinSet, spins);

        // hide the DB-sizes panel
        ckDBSizeActionPerformed(null);
        // hide the FDR-groups panel
        ckDefineGroupsActionPerformed(null);


        ActionListener stopMaximizingAl = new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                stopMaximizing = true;
                getFdr().stopMaximizing();
            }
        };

        fdrSettingsComplete.addStopMaximizingListener(stopMaximizingAl);
        fdrSettingsMedium.addStopMaximizingListener(stopMaximizingAl);
        fdrSettingsSimple.addStopMaximizingListener(stopMaximizingAl);

        ActionListener calcListener = new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                if (e.getActionCommand().contentEquals("boost")) {
                    getFdrSettingsComplete().btnStopBoost.setEnabled(true);
                    fdrSettingsSimple.btnStopBoost.setEnabled(true);
                    fdrSettingsMedium.btnStopBoost.setEnabled(true);
                }
                basicCalc();
            }
        };

        fdrSettingsSimple.addCalcListener(calcListener);
        fdrSettingsMedium.addCalcListener(calcListener);
        fdrSettingsComplete.addCalcListener(calcListener);

        spDistanceGroup.setSpecialValueText("No distance group");
        spDistanceGroup.setValue((Integer) 0);

        setMZIdentMLOwner();

        csvSelect.addAddListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                csvSelectAddActionPerformed(evt);
            }
        });

        fbFolderAddActionListener();
        
        cbLevel.setSelectedItem(Level.INFO);
        cbLevelActionPerformed(new ActionEvent(cbLevel, 0, ""));

        getDBFDR.setFdrgui(this);
        writeDB.setFdrgui(this);
        
        parseChageLog();
        try {
            if (this.getDBFDR.getSearch.getConnection() == null)
                this.tpInput.remove(0);
        } catch (Exception e) {
                this.tpInput.remove(0);
        }
        
        
    }

    public Handler riseLoggignTabOnError() {
        // add a hanlder that shows the log - tab when evver a warning or a an error happens
        java.util.logging.Handler tabRiser = new java.util.logging.Handler() {
            
            {
                this.setFilter(new Filter() {
                    
                    public boolean isLoggable(LogRecord record) {
                        return true;
                    }
                });
            }
            
            @Override
            public void publish(LogRecord record) {
                if (record.getLevel().intValue() >= Level.WARNING.intValue()) {
                    // find the log tab
                    for (int i = 0; i < jTabbedPane1.getTabCount(); i++) {
                        if (jTabbedPane1.getTitleAt(i).contentEquals("Log")) {
                            final int tabID = i;
                            // warnings are yellow
                            Color pc = Color.YELLOW;
                            // severe is red
                            if (record.getLevel().intValue() >= Level.SEVERE.intValue()) {
                                pc = Color.RED;
                            }
                            final Color c = pc;
                            SwingUtilities.invokeLater(new Runnable() {
                                public void run() {
                                    jTabbedPane1.setBackgroundAt(tabID, c);
                                    // if severe raise the log tab
                                    if (c == Color.RED) {
                                        jTabbedPane1.setSelectedIndex(tabID);
                                    }
                                    jTabbedPane1.addChangeListener(new ChangeListener() {
                                        @Override
                                        public void stateChanged(ChangeEvent e) {
                                            jTabbedPane1.removeChangeListener(this);
                                            jTabbedPane1.setBackgroundAt(tabID, null);
                                        }
                                    });
                                }
                            });
                            break;
                        }

                    }
                    final String message = record.getMessage();
                    final int level = record.getLevel().intValue();
                    SwingUtilities.invokeLater(new Runnable() {
                        @Override
                        public void run() {
                            if (level >= Level.SEVERE.intValue()) {
                                JOptionPane.showMessageDialog(FDRGUI.this, message, "Error", JOptionPane.ERROR_MESSAGE);
                            } else {
                                JOptionPane.showMessageDialog(FDRGUI.this, message, "Warning", JOptionPane.WARNING_MESSAGE);
                            }
                        }
                    });
                    
                }
            }
            
            @Override
            public void flush() {
                
            }
            
            @Override
            public void close() throws SecurityException {
                
            }
        };
        return tabRiser;
    }

    public void makeSpinnersSelectTextOnEnter(ArrayList<Container> con, HashSet<Container> spinSet, ArrayList<Container> spins) {
        // get all the spinners in all containers recursivly
        for (int i = 0; i < con.size(); i++) {
            for (Component c : con.get(i).getComponents()) {
                if (c instanceof JSpinner) {
                    if (!spinSet.contains(c)) {
                        spins.add((Container) c);
                        spinSet.add((Container) c);
                    }
                } else if (c instanceof Container) {
                    con.add((Container) c);
                }
            }
        }
        // now install the a handler on all text-fields within all spinners that selects all
        for (int i = 0; i < spins.size(); i++) {
            Container spin = spins.get(i);
            for (Component c : spin.getComponents()) {
                if (c instanceof JTextField) {
                    ((JTextField) c).addFocusListener(new java.awt.event.FocusAdapter() {
                        public void focusGained(final java.awt.event.FocusEvent evt) {
                            javax.swing.SwingUtilities.invokeLater(new Runnable() {
                                public void run() {
                                    // on focus gain select all text
                                    ((JTextField) evt.getSource()).selectAll();
                                }
                            });
                        }
                    });
                } else if (c instanceof Container) {
                    if (!spinSet.contains(c)) {
                        spins.add((Container) c);
                        spinSet.add((Container) c);
                    }
                }
            }
        }
    }

    public void setMZIdentMLOwner() {
        txtmzIdentOwnerFirst.setText(LocalProperties.getProperty("mzIdenMLOwnerFirst", txtmzIdentOwnerFirst.getText()));
        txtmzIdentOwnerLast.setText(LocalProperties.getProperty("mzIdenMLOwnerLast", txtmzIdentOwnerLast.getText()));
        txtmzIdentOwnerEmail.setText(LocalProperties.getProperty("mzIdenMLOwnerEmail", txtmzIdentOwnerEmail.getText()));
        txtmzIdentAdress.setText(LocalProperties.getProperty("mzIdenMLOwnerAddress", txtmzIdentAdress.getText()));
        txtmzIdentOwnerOrg.setText(LocalProperties.getProperty("mzIdenMLOwnerOrg", txtmzIdentOwnerOrg.getText()));
    }

    public void fbFolderAddActionListener() {
        fbFolder.addActionListener(new ActionListener() {
            public boolean setting = false;

            @Override
            public void actionPerformed(ActionEvent e) {
                if (setting) {
                    return;
                }
                setting = true;
                String filename = fbFolder.getText();
                setCsvOut(XiFDRUtils.splitFilename(filename));
                if (getCsvOut().tsv)
                    rbTSV.setSelected(true);
                if (getCsvOut().tsv)
                    rbCSV.setSelected(true);

                fbFolder.setFile(getCsvOut().folder + File.separator + getCsvOut().basename + getCsvOut().extension);
                setting = false;

            }
        });
        fbFolder.setLocalPropertyKey("XiFDR_LAST_CSV_OUT_FOLDER");
    }
    
    public void parseChageLog() {
        final Properties properties = new Properties();
        String propertyFile ="xifdrproject.properties";
        try {
            properties.load(this.getClass().getResourceAsStream(propertyFile));
        } catch (Exception e) {
            try {
                properties.load(this.getClass().getClassLoader().getResourceAsStream(propertyFile));                
            }catch (Exception ex) {
                Logger.getLogger(this.getClass().getName()).log(Level.WARNING,"Could not parse changelog",ex);
            }
        }
        String v = properties.getProperty("xifdr.changelog");
        txtchangelog.setText(v);
        
    }

    public void setTitle(String title) {
        super.setTitle("xiFDR (" + OfflineFDR.getXiFDRVersion() + ")" + title);

    }

    public void prepareFDRCalculation() throws NumberFormatException {
        ProteinGroupLink.MIN_DISTANCE_FOR_LONG = ((Number) spDistanceGroup.getValue()).intValue();
        getFdr().setSettings(fdrSettings);
        getFdr().setMinimumPeptideLength(fdrSettings.getMinPeptideLength());
        getFdr().setGroupPPIByHasInternal(cbGroupPPIBySelf.isSelected());
        getFdr().setGroupLinksByHasInternal(cbGroupLinksBySelf.isSelected());
        getFdr().setGroupPepPairByHasInternal(cbGroupPepsBySelf.isSelected());
        getFdr().setGroupPSMsByRun(cbGroupRun.isSelected());
        
        if (ckDBSize.isSelected()) {
            getFdr().setTargetDBSize((Double) spTargetDB.getValue());
            getFdr().setDecoyDBSize((Double) spDecoyDB.getValue());
            getFdr().setTargetProtDBSize((Double) spTargetDBProt.getValue());
            getFdr().setDecoyProtDBSize((Double) spDecoyDBProt.getValue());
            getFdr().setTargetProtDBSize((Double) spTargetDBProt.getValue());
            getFdr().setDecoyProtDBSize((Double) spDecoyDBProt.getValue());
        } else {
            getFdr().setTargetDBSize(999999999);
            getFdr().setDecoyDBSize(999999999);
            getFdr().setTargetProtDBSize(999999999);
            getFdr().setDecoyProtDBSize(999999999);
        }

        getFdr().setPsm_directional(fdrSettings.isPSMDirectional());
        getFdr().setPeptides_directional(fdrSettings.isPeptidePairDirectional());
        getFdr().setLinks_directional(fdrSettings.isLinkDirectional());
        getFdr().setPpi_directional(fdrSettings.isPPIDirectional());

        getFdr().setMinPepPerProteinGroup(fdrSettings.getMinProteinPepCount());
        getFdr().setMinPepPerProteinGroupLink(fdrSettings.getMinLinkPepCount());
        getFdr().setMinPepPerProteinGroupPair(fdrSettings.getMinPPIPepCount());

        getFdr().setMinDecoys(fdrSettings.getMinTD());

        int lengthSteps = tblPepLength.getRowCount();
        int[] pepLength = new int[lengthSteps];
        for (int i = 0; i < lengthSteps - 1; i++) {
            pepLength[i] = Integer.parseInt(tblPepLength.getValueAt(i, 0).toString());
        }

        getFdr().setLengthGroups(pepLength);
    }

    protected String fdrLevelSummary(FDRResultLevel fdrl) {
        StringBuilder sbSummary = new StringBuilder();
        Set<String> ids = fdrl.getGroupIDs();
        for (String fg : ids) {
            SubGroupFdrInfo sg = (SubGroupFdrInfo) fdrl.getGroup(fg);
            String gn = sg.fdrGroup;
            sbSummary.append(gn);
            sbSummary.append(" : ");
            sbSummary.append(String.format("%6d", sg.results.size()));
            sbSummary.append("[");
            sbSummary.append(String.format("%6.2f", sg.higherFDR * 100));
            sbSummary.append(",");
            sbSummary.append(String.format("%6.2f", sg.lowerFDR * 100));
            sbSummary.append("]");
            sbSummary.append(sg.filteredResult.size());
            sbSummary.append(" \n");
        }
        return sbSummary.toString();
    }

    /**
     * @return the m_fdr
     */
    public OfflineFDR getFdr() {
        return m_fdr;
    }

    /**
     * @param m_fdr the m_fdr to set
     */
    public void setFdr(OfflineFDR m_fdr) {
        this.m_fdr = m_fdr;
        ckGroupByCrossLinkerStubs.setEnabled(m_fdr.stubsFound());
    }

    /**
     * @return the m_result
     */
    public FDRResult getResult() {
        return m_result;
    }

    /**
     * @param m_result the m_result to set
     */
    public void setResult(FDRResult m_result) {
        this.m_result = m_result;
    }

    private class NeededOptionalComboBoxCellEditor extends DefaultCellEditor {

        EditorDelegate neededDelegate;
        EditorDelegate optionalDelegate;
        JComponent neededEditorComponent;
        JComponent optionalEditorComponent;

        public NeededOptionalComboBoxCellEditor(final JComboBox needed, final JComboBox optional) {
            super(needed);
            this.neededDelegate = delegate;
            this.neededEditorComponent = needed;

            optionalDelegate = new EditorDelegate() {
                public void setValue(Object value) {
                    optional.setSelectedItem(value);
                }

                public Object getCellEditorValue() {
                    return optional.getSelectedItem();
                }

                public boolean shouldSelectCell(EventObject anEvent) {
                    if (anEvent instanceof MouseEvent) {
                        MouseEvent e = (MouseEvent) anEvent;
                        return e.getID() != MouseEvent.MOUSE_DRAGGED;
                    }
                    return true;
                }

                public boolean stopCellEditing() {
                    if (optional.isEditable()) {
                        // Commit edited value.
                        optional.actionPerformed(new ActionEvent(NeededOptionalComboBoxCellEditor.this, 0, ""));
                    }
                    return super.stopCellEditing();
                }
            };
            optional.addActionListener(optionalDelegate);
            optionalEditorComponent = optional;

        }

        @Override
        public Component getTableCellEditorComponent(JTable table, Object value,
                boolean isSelected,
                int row, int column) {
            TableModel tm = table.getModel();
            if (Boolean.TRUE.equals(tm.getValueAt(row, 1))) {
                delegate = optionalDelegate;
                editorComponent = optionalEditorComponent;
            } else {
                delegate = neededDelegate;
                editorComponent = neededEditorComponent;
            }
            return super.getTableCellEditorComponent(table, value, isSelected, row, column);
        }
    }

    protected File getMzIdentMLOutput() {
        return fbMzIdentMLOut.getFile();
    }

    protected String getMzIdentMLOwnerLast() {
        return txtmzIdentOwnerLast.getText();
    }

    protected String getMzIdentMLOwnerFirst() {
        return txtmzIdentOwnerFirst.getText();
    }

    protected String getMzIdentMLOwnerEmail() {
        return txtmzIdentOwnerEmail.getText();
    }

    protected String getMzIdentMLOwnerOrg() {
        return txtmzIdentOwnerOrg.getText();
    }

    protected String getMzIdentMLOwnerAdress() {
        return txtmzIdentAdress.getText();
    }

    protected void exportMZIdentML() {
        if (!java.nio.charset.Charset.defaultCharset().displayName().contentEquals("UTF-8")) {
            JOptionPane.showMessageDialog(this, "Probably wrong file encoding: start xiFDRDB with: -Dfile.encoding=UTF-8 "); 
        }
                
        try {
            final File f = getMzIdentMLOutput();
            if (f == null) {
                setStatus("No output file selected");
                return;
            }
            final OfflineFDR fdr = getFdr();
            final FDRResult res = getResult();
            final MZIdentMLOwner mzo = new MZIdentMLOwner(getMzIdentMLOwnerFirst(), getMzIdentMLOwnerLast(), getMzIdentMLOwnerEmail(), getMzIdentMLOwnerOrg(), getMzIdentMLOwnerAdress());
            final String extension = cmbPeakListFormat.getSelectedItem().toString();
            final String template  = cmbMzMLScan2ID.getSelectedItem().toString();
            // if (f.canWrite())
            if ( fdr instanceof MZIdentXLFDR) {
                setStatus("Writing mzIdenML");
                writeMZIdentML();
                setStatus("done writing mzIdentML");
            } else if (getFdr()instanceof XiInFDR) {
                setEnableWrite(false);
                Runnable runnable = new Runnable() {
                    public void run() {
                        
                        try {
                            setStatus("Exporting mzIdenML");
                            MZIdentMLExport export = new MZIdentMLExport(fdr, res, mzo);
                            export.setForceExtension(extension);
                            export.setMZMLTemplate(template);
                            export.convertFile(f.getAbsolutePath());
                            
                            setStatus("done exporting mzIdentML");
                        } catch (Exception ex) {
                            setStatus("Error exporting mzIdentML");
                            Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, "Error writing mzIdentML", ex);
                        } finally {
                            setEnableWrite(true);
                        }
                    }
                };
                Thread t = new Thread(runnable);
                t.setName(t.getName() + " Writing mzIdentML");
                t.start();
            }
        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this, "Error while writing the xml-file:" + ex + "\n" + RArrayUtils.toString(ex.getStackTrace(), "\n"));
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Error while writing the xml-file", ex);
            setStatus("error:" + ex);
        }
    }    
    
    protected void writeMZIdentML() {
        try {
            final File f = getMzIdentMLOutput();
            // if (f.canWrite())
            if (getFdr() instanceof MZIdentXLFDR) {
                ((MZIdentXLFDR) getFdr()).writeMZIdentML(f.getAbsolutePath(), getResult());
            } else {
                JOptionPane.showMessageDialog(this, "mzIdentML export currently only supported with mzIdentML source files");
            }
        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this, "Error while writing the xml-file:" + ex + "\n" + RArrayUtils.toString(ex.getStackTrace(), "\n"));
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Error while writing the xml-file", ex);
            setStatus("error:" + ex);

        }
    }

    private boolean mzScoreDirectionChanged = false;

    protected void writeResults() {
        final OfflineFDR ofdr = getFdr();
        ofdr.setOutputLocale(lpCsvOutLocal.getSelectLocale());

        setEnableRead(false);
        setEnableCalc(false);
        setEnableWrite(false);
        final boolean prepostaa = ckPrePostAA.isSelected();
        
        final String sep = getSeparator(ofdr);
        

        Runnable runnable = new Runnable() {
            public void run() {
                try {
                    
                    setStatus("start writing");
                    Logger.getLogger(this.getClass().getName()).log(Level.INFO, "start writing");
                    //                    ofdr.writeFiles(txtFolder.getText(), txtBaseName.getText(), 0.05, 0.05, 0.05, 0.05, new int[]{4, 9});
                    String folder = getCsvOut().folder;
                    String basename = getCsvOut().basename;
                    if (getCsvOut().extension != null && !csvout.extension.isEmpty()) {
                        ofdr.writeFiles(folder, basename,getCsvOut().extension, sep, getResult());
                    } else {
                        ofdr.writeFiles(folder, basename, sep, getResult());
                    }
                    setStatus("finished writing: " + ofdr.summaryString(getResult()));

                } catch (FileNotFoundException ex) {
                    Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
                    setStatus("error: " + ex);
                } catch (Exception e) {
                    Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, e);
                    setStatus("error: " + e);
                }
                setEnableRead(true);
                setEnableCalc(true);
                setEnableWrite(true);
            }

        };
        Thread t = new Thread(runnable);
        t.setName("Write CSV");
        t.start();
    }

    protected String getSeparator(final OfflineFDR ofdr) {
        final String sep = rbTSV.isSelected() ? "\t" : (((DecimalFormat) ofdr.getNumberFormat()).getDecimalFormatSymbols().getDecimalSeparator()==','? ";" : ",");
        return sep;
    }

    
    protected void parseSuppliedMGF() throws IOException {
        UpdateableInteger ambiguousPeakListEntry = new UpdateableInteger(0);
        int noPeakListEntry = 0;
        ArrayList<MgfStyleTitleParser> peakLookup = new ArrayList<>();
        if (flPeakLists.getFiles().length > 0) {
            setStatus("Parsing peaklist");
            for (File f : flPeakLists.getFiles()) {
                MgfStyleTitleParser p = new MgfStyleTitleParser();
                p.parseFile(f);
                peakLookup.add(p);
            }
        }

        // try to recover peaklist informations (name and index)
        if (!peakLookup.isEmpty()) {
            setStatus("updating PSMs");
            int current = 0;
            int oldCurrent = -100;
            int all = getResult().psmFDR.size();
            int oldPerc = -1;
            for (PSM p : getResult().psmFDR) {
                int percent = (int) (current++ *1000.0 /all);
                if (percent>oldPerc&& current-oldCurrent>100) {
                    setStatus("updating PSMs ("+percent/10.0+"%)" );
                    oldPerc=percent;
                    oldCurrent = current;
                }
                if (p.getPeakListName() == null || p.getPeakListName().trim().isEmpty()) {
                    MgfStyleTitleParser.ParseEnrtry  i = null;
                    for (MgfStyleTitleParser pl : peakLookup) {
                        if (pl.getParsedFile().getName().contains(p.getRun())) {
                            i = pl.findScanIndex(p.getRun(), p.getScan());
                            if (i == null)
                                continue;
                            if (i.index == -1) {
                                ambiguousPeakListEntry.value++;
                                break;
                            } else if ( i.index > 0) {
                                setExperimentalPSMValues(p, pl, i);
                                break;
                            }
                        }
                    }
                    if (i == null) {
                        for (MgfStyleTitleParser pl : peakLookup) {
                            if (!pl.getParsedFile().getName().contains(p.getRun())) {
                                i = pl.findScanIndex(p.getRun(), p.getScan());
                                if (i == null)
                                    continue;
                                if (i.index == -1) {
                                    ambiguousPeakListEntry.value++;
                                    break;
                                } else if ( i.index > 0) {
                                    setExperimentalPSMValues(p, pl, i);
                                    break;
                                }
                            }
                        }
                    }
                    if (i == null) {
                        noPeakListEntry++;
                    }
                } else if (p.getFileScanIndex()== null) {
                    MgfStyleTitleParser.ParseEnrtry i = null;
                    for (MgfStyleTitleParser pl : peakLookup) {
                        String name = pl.getParsedFile().getName();
                        if (name.contentEquals(p.getPeakListName())) {
                            i = pl.findScanIndex(p.getRun(), p.getScan());
                            if (i.index == -1) {
                                ambiguousPeakListEntry.value++;
                                break;
                            } else if ( i.index > 0) {
                                setExperimentalPSMValues(p, pl, i);
                            }
                        }
                    }
                    if (i == null) {
                        noPeakListEntry++;
                    }
                }
            }
            String error = "";
            if (ambiguousPeakListEntry.value >0) {
                error="Scans found with ambiguous sources: " + ambiguousPeakListEntry +"\n";
            }
            if (noPeakListEntry >0) {
                error="Scans without detected sources: " + noPeakListEntry +"\n";
            }
            if (error.length()>0) {
                setStatus(error);
                JOptionPane.showMessageDialog(rootPane, error, "Not usable for XiView", JOptionPane.INFORMATION_MESSAGE);
            } else {
                setStatus("PSMs updated");
            }

        }
    }

    protected void setExperimentalPSMValues(PSM p, MgfStyleTitleParser pl, MgfStyleTitleParser.ParseEnrtry i) {
        p.setPeakListName(pl.getParsedFile().getName());
        p.setFileScanIndex(i.index);
        if (p.getPeakListName() == null || p.getPeakListName().isEmpty())
            p.setPeakListName(pl.fileName);
        if (p.getExperimentalMZ() == 0.0 || Double.isNaN(p.getExperimentalMZ()))
            p.setExperimentalMZ(i.expMZ);
        if (p.getExpCharge() == 0 )
            p.setExpCharge(i.expCharge);
    }
    
//    public void readCSV() {
//        try {
//
//            if (csvSelect.getFile() == null) {
//                JOptionPane.showMessageDialog(this, "No file selected", "no File Selected", JOptionPane.ERROR_MESSAGE);
//                return;
//            }
//
//            setStatus("Start");
//            final CSVinFDR ofdr = new CSVinFDR();
//
//
//
//            setFdr(ofdr);
//
//            addCSV(ofdr,null);
//        } catch (Exception ex) {
//            Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
//            setStatus("error:" + ex);
//            setEnableRead(true);
//        }
//
//    }
    
    protected CSVinFDR innerBatchReadCSV(File config, File fasta,CsvCondition filter,OfflineFDR.Normalisation normalisation) throws IOException, ParseException {
        Iterator<CsvParser> csvs = csvSelect.iterator();
        CsvParser csv = csvs.next();
        String forwardPattern = csvSelect.getForwardPattern();

        CSVinFDR ofdr = null;
        if (config != null && fasta != null) {
            ofdr = new XiCSVinFDR();
            ((XiCSVinFDR)ofdr).setConfig(new RunConfigFile(csvSelect.fbConfigIn.getFile()));
            ArrayList<String> fastas = new ArrayList<>(1);
            fastas.add(csvSelect.fbFastaIn.getFile().getAbsolutePath());
            ((XiCSVinFDR)ofdr).setFastas(fastas);
        }else {
            ofdr = new CSVinFDR();
        }
        if (forwardPattern != null && !forwardPattern.isEmpty())
            ofdr.setForwardPattern(forwardPattern);
        if (addCSV(ofdr, null, csv,filter)) {

            while (csvs.hasNext()) {
                csv = csvs.next();
                CSVinFDR nextfdr = null;

                if (config != null && fasta != null) {
                    nextfdr = new XiCSVinFDR();
                    ((XiCSVinFDR)nextfdr).setConfig(new RunConfigFile(csvSelect.fbConfigIn.getFile()));
                    ArrayList<String> fastas = new ArrayList<>(1);
                    fastas.add(csvSelect.fbFastaIn.getFile().getAbsolutePath());
                    ((XiCSVinFDR)nextfdr).setFastas(fastas);
                }else {
                    nextfdr = new CSVinFDR();
                }
                if (forwardPattern != null && !forwardPattern.isEmpty())
                    nextfdr.setForwardPattern(forwardPattern);


                if (!addCSV(nextfdr, ofdr, csv,filter)) {
                    return null;
                }
            }
            if (normalisation != OfflineFDR.Normalisation.None) {
                ofdr.normalizePSMs(normalisation);
            }
            return ofdr;
        }
        return null;
    }
    
    public void readAllCSV() {


        if (csvSelect.getFile() == null) {
            JOptionPane.showMessageDialog(this, "No file selected", "no File Selected", JOptionPane.ERROR_MESSAGE);
            return;
        }
        final CsvCondition filter = csvSelect.getFilter();
        final File config=csvSelect.fbConfigIn.getFile();
        final File fasta=csvSelect.fbFastaIn.getFile();
        final OfflineFDR.Normalisation normalisation = csvSelect.doNormalize();
        setEnableRead(false);
        setEnableCalc(false);
        setEnableWrite(false);

        Runnable runnable = new Runnable() {
            public void run() {
                try {

                    setStatus("Start");
                    setFdr(innerBatchReadCSV(config, fasta, filter, normalisation));
                    setEnableCalc(true);
                    setStatus("finished reading");
                } catch (Exception ex) {
                    Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
                    setStatus("error:" + ex);
                }
                setEnableRead(true);
                //  setEnableWrite(true);
            }

        };
        Thread t = new Thread(runnable);
        t.setName("Reading From CSV");
        t.start();

    }

    public void readAdditionalCSV() {
        if (csvSelect.getFile() == null) {
            JOptionPane.showMessageDialog(this, "No file selected", "no File Selected", JOptionPane.ERROR_MESSAGE);
            return;
        }

        setEnableRead(false);
        setEnableCalc(false);
        setEnableWrite(false);
        final CsvCondition filter = csvSelect.getFilter();
        final File config=csvSelect.fbConfigIn.getFile();
        final File fasta=csvSelect.fbFastaIn.getFile();
        final OfflineFDR.Normalisation normalisation = csvSelect.doNormalize();

        Runnable runnable = new Runnable() {
            public void run() {
                try {

                    setStatus("Start");

                    CSVinFDR nextfdr = innerBatchReadCSV(config, fasta, filter, normalisation);
                    
                    getFdr().normaliseAndAddPsmList(nextfdr, normalisation);
                    
                    setEnableCalc(true);
                    setStatus("finished reading");

                } catch (Exception ex) {
                    Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
                    setStatus("error:" + ex);
                    setEnableRead(true);
                }
                setEnableRead(true);
                //  setEnableWrite(true);
            }

        };
        Thread t = new Thread(runnable);
        t.setName("Reading From CSV");
        t.start();


    }

    protected boolean addCSV(final CSVinFDR ofdr, final CSVinFDR addto, final CsvParser csv, CsvCondition filter) throws IOException, ParseException {
        ofdr.setPSMScoreHighBetter(csvSelect.higherIsBetter());

        setStatus("Read from " + csv.getInputFile().getName());
        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "read from " + csv.getInputFile().getAbsolutePath());
        if (addto == null) {
            PSM.resetAdditionalColumnNames();
        }
        if (!ofdr.readCSV(csv, filter)) {
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Could not read from " + csv.getInputFile().getAbsolutePath());
            setStatus("Error reading from " + csv.getInputFile().getName());            
            return false;
        }
        if (addto != null) {
            addto.add(ofdr);
            addto.addSource(ofdr.getSources(), ofdr.getFilter());
        }
        return true;
    }

    public void readMZIdentML() {
        try {
            setEnableRead(false);
            setEnableCalc(false);
            setEnableWrite(false);

            final File mzIdentIn = fbMZIdentMLIn.getFile();
            if (mzIdentIn == null) {
                setStatus("no file selected");
                setEnableRead(true);
                //btnRead.setEnabled(true);
                return;
            } else if (!mzIdentIn.canRead()) {
                setStatus("Can't read selected file");
                setEnableRead(true);
                return;
            }

            setStatus("Start");
            final MZIdentXLFDR ofdr = new MZIdentXLFDR();

            setFdr(ofdr);

            getFdr().setPSMScoreHighBetter(rbMZHighBetter.isSelected());
//            ofdr.setCrosslinkedDonorModAcc(txtCrossLinkedDonorModCvParam.getText());
//            ofdr.setCrosslinkedReceptorModAcc(txtCrossLinkedReceptorModCvParam.getText());
//            ofdr.setCrosslinkedSIIAcc(txtCrossLinkedPepCvParam.getText());
            ofdr.setPSMScore(cbMZMatchScoreName.getSelectedItem().toString());

            Runnable runnable = new Runnable() {
                public void run() {
                    try {
                        setStatus("Read from " + mzIdentIn.getName());
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "read from " + mzIdentIn.getAbsolutePath());
                        ofdr.readMzIdentML(mzIdentIn);
                        setEnableCalc(true);
                        setStatus("finished reading");
                    } catch (Exception ex) {
                        Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
                        setStatus("error:" + ex);
                    }
                    setEnableRead(true);
                    //  setEnableWrite(true);
                }
            };
            new Thread(runnable).start();
        } catch (Exception ex) {
            Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
            setStatus("error:" + ex);
            setEnableRead(true);
        }

    }

    protected void calculateFDR() {
        final OfflineFDR ofdr = getFdr();
        setEnableCalc(false);
        setEnableWrite(false);
        setEnableRead(false);
        final Double saftyfactor = (Double) fdrSettings.getReportFactor();
        ofdr.setSettings(fdrSettings);
        prepareFDRCalculation();

        Runnable runnable = new Runnable() {
            final Double psmfdr = (Double) getFdrSettings().getPSMFDR();
            final Double pepfdr = (Double) getFdrSettings().getPeptidePairFDR();
            final Double protfdr = (Double) getFdrSettings().getProteinGroupFDR();
            final Double linkfdr = (Double) getFdrSettings().getProteinGroupLinkFDR();
            final Double ppifdr = (Double) getFdrSettings().getProteinGroupPairFDR();
            final boolean filterToUniquePSM = getFdrSettings().filterToUniquePSM();

            public void run() {
                try {
                    innerFDRCalculation(psmfdr, pepfdr, protfdr, linkfdr, ppifdr, ofdr, saftyfactor, filterToUniquePSM);

                } catch (Exception e) {
                    Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, e);
                    setStatus("error:" + e);
                }
                setEnableCalc(true);
                setEnableRead(true);
                //setEnableWrite(true);
            }
        };
        Thread r = new Thread(runnable);
        r.setName(r.getName() + " - " + "FDR");
        r.start();
    }

    protected void innerFDRCalculation(Double psmfdr, Double pepfdr, Double protfdr, Double linkfdr, Double ppifdr, final OfflineFDR ofdr, Double saftyfactor, boolean filterToUniquePSM) {
        setStatus("Start");
        setStatus("Calculating fdr");
        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Calculating fdr");
        
        double fm = Runtime.getRuntime().freeMemory();
        double mm = Runtime.getRuntime().maxMemory();
        double tm = Runtime.getRuntime().totalMemory();
        double um = tm-fm;
        if (um/tm>0.5) {
            // not enough memory - will have to give up on the old result already now
            if (getResult()!=null)
                getResult().clear();
            //setResult(null);
        }
            
//        final FDRResult result = ofdr.calculateFDR(psmfdr, pepfdr, protfdr, linkfdr, ppifdr, saftyfactor, ckIgnoreGroups1.isSelected(), true, filterToUniquePSM);
        ofdr.setIgnoreGroupsSetting(ckIgnoreGroups1.isSelected());
        final FDRResult result = ofdr.calculateFDR(fdrSettings, true);
        setResult(result);
        EventQueue.invokeLater(new Runnable() {
            public void run() {
                reportResultToSumaryTab(ofdr, result);
            }
        });
    }

    protected void reportResultToSumaryTab(OfflineFDR ofdr, FDRResult result) {
        setStatus("finished");
        setStatus(ofdr.summaryString(getResult()));
        setEnableWrite(true);
//        HashMap<String, Double> summary = ofdr.summaryString(m_result);
        txtSumInput.setText(Integer.toString(ofdr.getInputPSMs().size()));

        int sumXLPSM = 0;
        int sumLinearPSM = 0;

        int targetXLPSM = 0;
        int targetLinearPSM = 0;
        int decoyXLPSM = 0;
        int decoyLinearPSM = 0;
        for (PSM psm : getResult().psmFDR.filteredResults()) {
            if (psm.isLinear()) {
                sumLinearPSM++;
                if (!psm.isDecoy()) {
                    targetLinearPSM++;
                }

            } else {
                sumXLPSM++;
                if (!psm.isDecoy()) {
                    targetXLPSM++;
                }
            }
        }
        int sumPSM = sumXLPSM + sumLinearPSM;

        int sumXLPepPairs = 0;
        int sumLinearPepPairs = 0;

        int targetXLPepPairs = 0;
        int targetLinearPepPairs = 0;
        int decoyXLPepPairs = 0;
        int decoyLinearPepPairs = 0;
        for (PeptidePair pp : result.peptidePairFDR.filteredResults()) {
            if (pp.isLinear()) {
                sumLinearPepPairs++;
                if (!pp.isDecoy()) {
                    targetLinearPepPairs++;
                }

            } else {
                sumXLPepPairs++;
                if (!pp.isDecoy()) {
                    targetXLPepPairs++;
                }
            }
        }
        int sumPepPairs = sumXLPepPairs + sumLinearPepPairs;

        int sumLinksBetweenDecoy = 0;
        int sumLinksInternalDecoy = 0;
        int sumLinksBetweenDD = 0;
        int sumLinksInternalDD = 0;
        int sumLinksInternalTarget = 0;
        int sumLinksBetweenTarget = 0;
        int sumLinks = getResult().proteinGroupLinkFDR.getResultCount();
        for (ProteinGroupLink l : getResult().proteinGroupLinkFDR.filteredResults()) {
            if (l.isDecoy()) {
                if (l.isInternal) {
                    sumLinksInternalDecoy++;
                    if (l.isDD()) {
                        sumLinksInternalDD++;
                    }
                } else {
                    sumLinksBetweenDecoy++;
                    if (l.isDD()) {
                        sumLinksBetweenDD++;
                    }
                }
            } else if (l.isInternal) {
                sumLinksInternalTarget++;
            } else {
                sumLinksBetweenTarget++;
            }

        }

        int sumLinksBetweenDecoyUF = 0;
        int sumLinksInternalDecoyUF = 0;
        int sumLinksBetweenDDUF = 0;
        int sumLinksInternalDDUF = 0;
        int sumLinksInternalTargetUF = 0;
        int sumLinksBetweenTargetUF = 0;
        for (ProteinGroupLink l : getResult().proteinGroupLinkFDR) {
            if (l.isDecoy()) {
                if (l.isInternal) {
                    sumLinksInternalDecoyUF++;
                    if (l.isDD()) {
                        sumLinksInternalDDUF++;
                    }
                } else {
                    sumLinksBetweenDecoyUF++;
                    if (l.isDD()) {
                        sumLinksBetweenDDUF++;
                    }
                }
            } else if (l.isInternal) {
                sumLinksInternalTargetUF++;
            } else {
                sumLinksBetweenTargetUF++;
            }

        }

        
        
        int sumProteinGroupPairs = getResult().proteinGroupPairFDR.getResultCount();
        int sumProteinGroupPairsBetweenDecoy = 0;
        int sumProteinGroupPairsInternalDecoy = 0;
        int sumProteinGroupPairsBetweenDD = 0;
        int sumProteinGroupPairsInternalDD = 0;
        int sumProteinGroupPairsInternalTarget = 0;
        int sumProteinGroupPairsBetweenTarget = 0;
        for (ProteinGroupPair pgl : getResult().proteinGroupPairFDR.filteredResults()) {
            if (pgl.isDecoy()) {
                if (pgl.isInternal()) {
                    sumProteinGroupPairsInternalDecoy++;
                    if (pgl.isDD())
                        sumProteinGroupPairsInternalDD++;
                } else {
                    sumProteinGroupPairsBetweenDecoy++;
                    if (pgl.isDD())
                        sumProteinGroupPairsBetweenDD++;
                }
            } else if (pgl.isInternal()) {
                sumProteinGroupPairsInternalTarget++;
            } else {
                sumProteinGroupPairsBetweenTarget++;
            }

        }
//        Integer sumLinksProtGroups =  ofdr.getFDRProteinGroups().size();

        String[] nice = MiscUtils.arrayToStringWithDifferenceOrientedFormat(new double[]{result.psmFDR.getHigherFDR() * 100, result.psmFDR.getLowerFDR() * 100}, 1);
        txtSumPSM.setText(sumPSM + " [" + nice[0] + "%," + nice[1] + "%]");
        txtSumPSM.setToolTipText(fdrLevelSummary(result.psmFDR));
        txtSumPSMXL.setText(sumXLPSM + " (" + (targetXLPSM) + " Target)");
        txtSumPSMLinear.setText(sumLinearPSM + " (" + (targetLinearPSM) + " Target)");

        nice = MiscUtils.arrayToStringWithDifferenceOrientedFormat(new double[]{result.peptidePairFDR.getHigherFDR() * 100, result.peptidePairFDR.getLowerFDR() * 100}, 1);
        txtSumPepPairs.setText(sumPepPairs + " [" + nice[0] + "%," + nice[1] + "%]");
        txtSumPepPairs.setToolTipText(fdrLevelSummary(result.peptidePairFDR));
        txtSumPepPairsXL.setText(sumXLPepPairs + " (" + targetXLPepPairs + " Target)");
        txtSumPepPairsLinear.setText(sumLinearPepPairs + " (" + targetLinearPepPairs + " Target)");

        nice = MiscUtils.arrayToStringWithDifferenceOrientedFormat(new double[]{result.proteinGroupLinkFDR.getHigherFDR() * 100, result.proteinGroupLinkFDR.getLowerFDR() * 100}, 1);
        txtSumLinks.setText(sumLinks + " [" + nice[0] + "%," + nice[1] + "%]");
        txtSumLinks.setToolTipText(fdrLevelSummary(result.proteinGroupLinkFDR));
        nice = MiscUtils.arrayToStringWithDifferenceOrientedFormat(new double[]{(sumLinksBetweenDecoyUF-2*sumLinksBetweenDDUF)/(double)sumLinksBetweenTargetUF * 100, (sumLinksBetweenDecoyUF-2*sumLinksBetweenDDUF+1)/(double)sumLinksBetweenTargetUF * 100}, 1);
        txtSumLinksBetween.setText((sumLinksBetweenTarget + sumLinksBetweenDecoy) + " (" + sumLinksBetweenTarget + " TT) [" + nice[0] + "%," + nice[1] + "%]");
        nice = MiscUtils.arrayToStringWithDifferenceOrientedFormat(new double[]{(sumLinksInternalDecoyUF-2*sumLinksInternalDDUF)/(double)sumLinksInternalTargetUF * 100, (sumLinksInternalDecoyUF-2*sumLinksInternalDDUF+1)/(double)sumLinksInternalTargetUF * 100}, 1);
        txtSumLinksInternal.setText((sumLinksInternalDecoy + sumLinksInternalTarget) + " (" + sumLinksInternalTarget + " TT)  [" + nice[0] + "%," + nice[1] + "%]");

        nice = MiscUtils.arrayToStringWithDifferenceOrientedFormat(new double[]{result.proteinGroupFDR.getHigherFDR() * 100, result.proteinGroupFDR.getLowerFDR() * 100}, 1);
        txtSumProtGroups.setText(Integer.toString(getResult().proteinGroupFDR.getResultCount()) + " [" + nice[0] + "%," + nice[1] + "%]");
        txtSumProtGroups.setToolTipText(fdrLevelSummary(result.proteinGroupFDR));

        nice = MiscUtils.arrayToStringWithDifferenceOrientedFormat(new double[]{result.proteinGroupPairFDR.getHigherFDR() * 100, result.proteinGroupPairFDR.getLowerFDR() * 100}, 1);
        txtSumProtGroupPairs.setText(sumProteinGroupPairs + " [" + nice[0] + "%," + nice[1] + "%]");
        txtSumProtGroupPairs.setToolTipText(fdrLevelSummary(result.proteinGroupPairFDR));
        nice = MiscUtils.arrayToStringWithDifferenceOrientedFormat(new double[]{(sumProteinGroupPairsBetweenDecoy-2*sumProteinGroupPairsBetweenDD)/(double)sumProteinGroupPairsBetweenTarget * 100, (sumProteinGroupPairsBetweenDecoy-2*sumProteinGroupPairsBetweenDD+1)/(double)sumProteinGroupPairsBetweenTarget * 100}, 1);
        txtSumProtGroupPairsBetween.setText((sumProteinGroupPairsBetweenTarget + sumProteinGroupPairsBetweenDecoy) + " (" + sumProteinGroupPairsBetweenTarget + " TT) [" + nice[0] + "%," + nice[1] + "%]");
        nice = MiscUtils.arrayToStringWithDifferenceOrientedFormat(new double[]{(sumProteinGroupPairsInternalDecoy-2*sumProteinGroupPairsInternalDD)/(double)sumProteinGroupPairsInternalTarget * 100, (sumProteinGroupPairsInternalDecoy-2*sumProteinGroupPairsInternalDD+1)/(double)sumProteinGroupPairsInternalTarget * 100}, 1);
        txtSumProtGroupPairsInternal.setText((sumProteinGroupPairsInternalDecoy + sumProteinGroupPairsInternalTarget) + " (" + sumProteinGroupPairsInternalTarget + " TT) [" + nice[0] + "%," + nice[1] + "%]");

// Logger.getLogger(this.getClass().getName()).log(Level.INFO, "finished writing");
    }

    private void basicCalc() {
        if (getFdr() != null) {
            FDRSettingsImpl.transferSettings(fdrSettings, fdrSettingsComplete);
            FDRSettingsImpl.transferSettings(fdrSettings, fdrSettingsMedium);
            FDRSettingsImpl.transferSettings(fdrSettings, fdrSettingsSimple);
            final FDRSettingsImpl settings = new FDRSettingsImpl(fdrSettings);
            
            getFdr().setMinDecoys(settings.getMinTD());
            
            settings.setGroupByCrosslinkerStubs(ckGroupByCrossLinkerStubs.isEnabled() && ckGroupByCrossLinkerStubs.isSelected());
                    
            if (fdrSettings.doOptimize() != null) {
                final OfflineFDR.FDRLevel l = settings.doOptimize();
                if (l == OfflineFDR.FDRLevel.PSM && !(settings.boostMinFragments() || settings.boostDeltaScore() || settings.boostPepCoverage()))
                    JOptionPane.showMessageDialog(this, "Boosting of that Level currently not supported!");
                else { 
                    Thread ml = new Thread() {
                        public void run() {
                            maximise(l, settings);
                            getFdrSettingsComplete().btnStopBoost.setEnabled(false);
                            fdrSettingsMedium.btnStopBoost.setEnabled(false);
                            fdrSettingsSimple.btnStopBoost.setEnabled(false);
                            
                        }
                    };
                    ml.start();
                }
                            
            } else {
                calculateFDR();
            }
        }
    }

    public void setStatus(final String status) {
        m_status.write(status);
//        m_statusmessages.put(Thread.currentThread(),status);
//        Runnable setModel = new Runnable() {
//            public void run() {
//                String message = "";
//                for (Thread t : m_statusmessages.keySet())  {
//                    message+=" - " + m_statusmessages.get(t);
//                    if (!t.isAlive()) {
//                        m_statusmessages.remove(t);
//                    }
//                }
//                txtStatus.setText(message.substring(3));
//                set
//                
//            }
//        };
//        javax.swing.SwingUtilities.invokeLater(setModel);
    }

    public void setEnableRead(final boolean enable) {
        Runnable setModel = new Runnable() {
            public void run() {
                csvSelect.setEnabled(enable);
                csvSelect.setEnableAdd(enable && m_fdr != null);
                btnReadMZIdent.setEnabled(enable && ((JPanel) getFdrSettings()).isEnabled());
                getDBFDR.setEnableRead(enable);
                getDBFDR.setEnableAdd(enable && getFdr() != null);
            }
        };
        javax.swing.SwingUtilities.invokeLater(setModel);
    }

    public void setEnableCalc(final boolean enable) {
        Runnable setModel = new Runnable() {
            public void run() {
                ((JPanel) getFdrSettings()).setEnabled(enable);
//                btnMaxLink.setEnabled(enable);
//                btnMaxPPI.setEnabled(enable);
                getFdrSettingsComplete().setEnabled(enable);
                fdrSettingsMedium.setEnabled(enable);
                fdrSettingsSimple.setEnabled(enable);
                btnCalcRanges.setEnabled(enable);
            }
        };
        javax.swing.SwingUtilities.invokeLater(setModel);
    }

    public void setEnableWrite(final boolean enable) {
        final boolean isMzIdent = getFdr() instanceof MZIdentXLFDR;
        Runnable setModel = new Runnable() {
            public void run() {
                btnWrite.setEnabled(enable);
                writeDB.setEnableWrite(getFdr() instanceof DBinFDR && enable);
                ckPrePostAA.setEnabled(getFdr() instanceof DBinFDR);

                btnWriteMzIdentML.setEnabled(enable && ((isMzIdent) || getFdr() instanceof XiInFDR));
            }
        };
        javax.swing.SwingUtilities.invokeLater(setModel);
    }

    public void clearResults() {
        txtSumLinks.setText("");
        txtSumLinksBetween.setText("");
        txtSumLinksInternal.setText("");
        txtSumPSM.setText("");
        txtSumPSMLinear.setText("");
        txtSumPSMXL.setText("");
        txtSumPepPairs.setText("");
        txtSumPepPairsLinear.setText("");
        txtSumPepPairsXL.setText("");
        txtSumProtGroupPairs.setText("");
        txtSumProtGroupPairsBetween.setText("");
        txtSumProtGroupPairsInternal.setText("");
        txtSumProtGroups.setText("");
    }

    public void maximise(OfflineFDR.FDRLevel level, FDRSettings settings) {
        clearResults();
        prepareFDRCalculation();
       
        setEnableRead(false);
        setEnableCalc(false);
        setEnableWrite(false);
        m_fdr.setIgnoreGroupsSetting(ckIgnoreGroups1.isSelected());
        
        final MaximisingStatus result = m_fdr.maximise(settings, level, settings.getBoostBetween(), new MaximizingUpdate() {
            @Override
            public void setStatus(final MaximisingStatus state) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        getFdrSettingsComplete().setMinPeptideStubFilter(state.showMinStubs);
                        getFdrSettingsComplete().setMinPeptideDoubletFilter(state.showMinDoublets);
                        getFdrSettingsComplete().setMinPeptideFragmentsFilter(state.showMinFrags);
                        getFdrSettingsComplete().setMinDeltaScoreFilter(state.showDelta);
                        getFdrSettingsComplete().setMinPeptideCoverageFilter(state.showPepCoverage);
                        getFdrSettingsComplete().setPeptidePairFDR(state.showPepFDR);
                        getFdrSettingsComplete().setPSMFDR(state.showPSMFDR);
                        getFdrSettingsComplete().setPeptidePairFDR(state.showPepFDR);
                        getFdrSettingsComplete().setProteinGroupFDR(state.showProtFDR);
                        getFdrSettingsComplete().setProteinGroupLinkFDR(state.showLinkFDR);

                        fdrSettingsMedium.setPeptidePairFDR(state.showPepFDR);
                        fdrSettingsMedium.setPSMFDR(state.showPSMFDR);
                        fdrSettingsMedium.setPeptidePairFDR(state.showPepFDR);
                        fdrSettingsMedium.setProteinGroupFDR(state.showProtFDR);
                        fdrSettingsMedium.setProteinGroupLinkFDR(state.showLinkFDR);
                        

                        txtSumPSM.setText(state.showPSMCount);
                        txtSumPepPairs.setText(state.showPepCount);
                        txtSumProtGroups.setText(state.showProtCount);
                        txtSumLinks.setText(state.showLinkCount);
                        txtSumProtGroupPairs.setText(state.showPPICount);

                        txtSumLinksBetween.setText(state.showLinkCountBetween);

                        txtSumProtGroupPairsBetween.setText(state.showPPICountBetween);
                    }
                });
            }

            @Override
            public void setStatusText(String text) {
                FDRGUI.this.setStatus(text);
            }

            @Override
            public void reportError(final String text, final Exception ex) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        JOptionPane.showMessageDialog(FDRGUI.this, text + "\n" + ex, "Error Maximizing", JOptionPane.ERROR_MESSAGE);
                    }
                });
            }
        });
        setEnableRead(true);
        setEnableCalc(true);
        if (result != null) {
            setEnableWrite(true);
        
            setResult(result.result);
            javax.swing.SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                    reportResultToSumaryTab(m_fdr, result.result);
                    JOptionPane.showMessageDialog(rootPane, "found " + result.resultCount + "("+ result.resultCountBetween +" between) matches for following settings \n"
                            + "\nPSM fdr:    " + result.showPSMFDR
                            + "\nPeptide fdr:" + result.showPepFDR
                            + "\nProtein FDR:" + result.showProtFDR
                            + "\nLink FDR:" + result.showLinkFDR, "best parameters found for max protein group pairs", JOptionPane.INFORMATION_MESSAGE);
                }
            });        
        }
    }
    

    private void changeFDRSettings(java.awt.event.ActionEvent evt) {
        if (rbFDRSimple.isSelected()) {
            if (fdrSettings != null) {
                fdrSettingsSimple.setAll(fdrSettings);
            }
            if (spFDRSettingsWrapper.getComponentCount() > 3) {
                spFDRSettingsWrapper.remove(spFDRSettingsWrapper.getComponent(3));
            }
            spFDRSettingsWrapper.setViewportView(fdrSettingsSimple);
            fdrSettings = fdrSettingsSimple;
        } else if (rbFDRMedium.isSelected()) {
            if (fdrSettings != null) {
                fdrSettingsMedium.setAll(fdrSettings);
            }
            if (spFDRSettingsWrapper.getComponentCount() > 3) {
                spFDRSettingsWrapper.remove(spFDRSettingsWrapper.getComponent(3));
            }
            spFDRSettingsWrapper.setViewportView(fdrSettingsMedium);
            fdrSettings = fdrSettingsMedium;
        } else if (rbFDRComplete.isSelected()) {
            if (fdrSettings != null) {
                fdrSettingsComplete.setAll(fdrSettings);
            }
            if (spFDRSettingsWrapper.getComponentCount() > 3) {
                spFDRSettingsWrapper.remove(spFDRSettingsWrapper.getComponent(3));
            }
            spFDRSettingsWrapper.setViewportView(fdrSettingsComplete);
            fdrSettings = fdrSettingsComplete;
        }
        // TODO add your handling code here:
    }

    public void setFDRSettings(final FDRSettings settings) {
        EventQueue.invokeLater(new Runnable() {
            @Override
            public void run() {
                rbFDRMedium.setSelected(true);
                rbFDRComplete.setSelected(false);
                rbFDRSimple.setSelected(false);
                rbFDRCompleteActionPerformed(null);

                fdrSettingsSimple.setAll(settings);
                fdrSettingsMedium.setAll(settings);
                getFdrSettingsComplete().setAll(settings);
            }
        });
    }
    
    public void setInput(final String fileName) {
        EventQueue.invokeLater(new Runnable() {
            @Override
            public void run() {
                csvSelect.setInputFile(fileName);
                tpInput.setSelectedComponent(pCSVInput);
            }
        });
    }

    public void setXiConfig(String fileName) {
        csvSelect.setXiConfig(fileName);
    }

    public void setFasta(final String fileName) {
        EventQueue.invokeLater(new Runnable() {
            @Override
            public void run() {
                csvSelect.setFasta(fileName);
            }
        });
    }

    public void setCsvFlagModifications(final boolean flag) {
        EventQueue.invokeLater(new Runnable() {
            @Override
            public void run() {
                csvSelect.flagModifications(flag);
                }
        });
    }
    
    public void setCSVOutFile() {
        setCSVOutFile(fbFolder.getFile().getAbsolutePath());
        setCsvOut(XiFDRUtils.splitFilename(fbFolder.getFile().getAbsolutePath()));

    }
    
    public void setCSVOutFile(String path) {
        fbFolder.setFile(path);
        
    }
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        bgSeparator = new javax.swing.ButtonGroup();
        bgScoreDirectionMzIdentML = new javax.swing.ButtonGroup();
        cbCSVHeaders = new javax.swing.JComboBox();
        cbCSVHeaderOptional = new javax.swing.JComboBox();
        bgFDRSettingType = new javax.swing.ButtonGroup();
        fdrSettingsComplete = new org.rappsilber.fdr.gui.components.FDRSettingsComplete();
        fdrSettingsSimple = new org.rappsilber.fdr.gui.components.FDRSettingsSimple();
        fdrSettingsMedium = new org.rappsilber.fdr.gui.components.FDRSettingsMedium();
        rbFDRSimple = new javax.swing.JRadioButton();
        jScrollPane6 = new javax.swing.JScrollPane();
        jTabbedPane1 = new javax.swing.JTabbedPane();
        jPanel2 = new javax.swing.JPanel();
        tpInput = new javax.swing.JTabbedPane();
        pDatabase = new javax.swing.JPanel();
        getDBFDR = new org.rappsilber.fdr.gui.components.GetDBFDR();
        pCSVInput = new javax.swing.JPanel();
        csvSelect = new org.rappsilber.fdr.gui.components.CSVSelection();
        pMZIdentMLInput = new javax.swing.JPanel();
        fbMZIdentMLIn = new org.rappsilber.gui.components.FileBrowser();
        jLabel19 = new javax.swing.JLabel();
        btnReadMZIdent = new javax.swing.JButton();
        jLabel22 = new javax.swing.JLabel();
        rbMZHighBetter = new javax.swing.JRadioButton();
        rbMZLowBetter = new javax.swing.JRadioButton();
        cbMZMatchScoreName = new javax.swing.JComboBox();
        spExtraLong = new javax.swing.JPanel();
        pDatabseSize = new javax.swing.JPanel();
        spDecoyDBProt = new javax.swing.JSpinner();
        spTargetDBProt = new javax.swing.JSpinner();
        lblProtein = new javax.swing.JLabel();
        lblPeptide = new javax.swing.JLabel();
        spTargetDB = new javax.swing.JSpinner();
        spDecoyDB = new javax.swing.JSpinner();
        lblDecoyDB = new javax.swing.JLabel();
        lblTargetDB = new javax.swing.JLabel();
        spDecoyDBLinks = new javax.swing.JSpinner();
        spTargetDBLinks = new javax.swing.JSpinner();
        lblLinkDB = new javax.swing.JLabel();
        ckDBSize = new javax.swing.JCheckBox();
        rbFDRMedium = new javax.swing.JRadioButton();
        rbFDRComplete = new javax.swing.JRadioButton();
        spFDRSettingsWrapper = new javax.swing.JScrollPane();
        pFDRGroups = new javax.swing.JPanel();
        spPepLength = new javax.swing.JScrollPane();
        tblPepLength = new javax.swing.JTable();
        ckIgnoreGroups1 = new javax.swing.JCheckBox();
        spDistanceGroup = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        jLabel1 = new javax.swing.JLabel();
        cbGroupPepsBySelf = new javax.swing.JCheckBox();
        cbGroupLinksBySelf = new javax.swing.JCheckBox();
        cbGroupPPIBySelf = new javax.swing.JCheckBox();
        cbGroupRun = new javax.swing.JCheckBox();
        ckGroupByCrossLinkerStubs = new javax.swing.JCheckBox();
        ckDefineGroups = new javax.swing.JCheckBox();
        pResult = new javax.swing.JPanel();
        tpResult = new javax.swing.JTabbedPane();
        jPanel8 = new javax.swing.JPanel();
        jLabel12 = new javax.swing.JLabel();
        txtSumInput = new javax.swing.JTextField();
        jPanel9 = new javax.swing.JPanel();
        btnPSMInfo = new javax.swing.JButton();
        jPanel18 = new javax.swing.JPanel();
        jLabel34 = new javax.swing.JLabel();
        jLabel35 = new javax.swing.JLabel();
        lblSumBetween1 = new javax.swing.JLabel();
        lblSumInternal1 = new javax.swing.JLabel();
        jPanel15 = new javax.swing.JPanel();
        jLabel13 = new javax.swing.JLabel();
        txtSumPSM = new javax.swing.JTextField();
        txtSumPSMXL = new javax.swing.JTextField();
        txtSumPSMLinear = new javax.swing.JTextField();
        jPanel16 = new javax.swing.JPanel();
        jLabel14 = new javax.swing.JLabel();
        txtSumPepPairs = new javax.swing.JTextField();
        txtSumPepPairsXL = new javax.swing.JTextField();
        txtSumPepPairsLinear = new javax.swing.JTextField();
        btnPepInfo = new javax.swing.JButton();
        jPanel20 = new javax.swing.JPanel();
        jLabel15 = new javax.swing.JLabel();
        txtSumProtGroups = new javax.swing.JTextField();
        jLabel30 = new javax.swing.JLabel();
        jLabel31 = new javax.swing.JLabel();
        btnProtInfo = new javax.swing.JButton();
        jPanel12 = new javax.swing.JPanel();
        jLabel33 = new javax.swing.JLabel();
        jLabel32 = new javax.swing.JLabel();
        lblSumBetween = new javax.swing.JLabel();
        lblSumInternal = new javax.swing.JLabel();
        jPanel21 = new javax.swing.JPanel();
        jLabel16 = new javax.swing.JLabel();
        txtSumLinks = new javax.swing.JTextField();
        txtSumLinksBetween = new javax.swing.JTextField();
        txtSumLinksInternal = new javax.swing.JTextField();
        btnLinkInfo = new javax.swing.JButton();
        jPanel22 = new javax.swing.JPanel();
        jLabel17 = new javax.swing.JLabel();
        txtSumProtGroupPairs = new javax.swing.JTextField();
        txtSumProtGroupPairsBetween = new javax.swing.JTextField();
        txtSumProtGroupPairsInternal = new javax.swing.JTextField();
        btnPPIInfo = new javax.swing.JButton();
        jPanel10 = new javax.swing.JPanel();
        rbTSV = new javax.swing.JRadioButton();
        rbCSV = new javax.swing.JRadioButton();
        btnWrite = new javax.swing.JButton();
        fbFolder = new org.rappsilber.gui.components.FileBrowser();
        jLabel9 = new javax.swing.JLabel();
        ckPrePostAA = new javax.swing.JCheckBox();
        lpCsvOutLocal = new org.rappsilber.gui.components.LocalPicker();
        jLabel28 = new javax.swing.JLabel();
        btnCalcRanges = new javax.swing.JButton();
        jPanel14 = new javax.swing.JPanel();
        fbMzIdentMLOut = new org.rappsilber.gui.components.FileBrowser();
        btnWriteMzIdentML = new javax.swing.JButton();
        txtmzIdentOwnerFirst = new javax.swing.JTextField();
        jLabel2 = new javax.swing.JLabel();
        txtmzIdentOwnerLast = new javax.swing.JTextField();
        jLabel3 = new javax.swing.JLabel();
        txtmzIdentOwnerEmail = new javax.swing.JTextField();
        jScrollPane1 = new javax.swing.JScrollPane();
        txtmzIdentAdress = new javax.swing.JTextArea();
        jLabel4 = new javax.swing.JLabel();
        jLabel5 = new javax.swing.JLabel();
        txtmzIdentOwnerOrg = new javax.swing.JTextField();
        lblPeaklistExtension = new javax.swing.JLabel();
        cmbPeakListFormat = new javax.swing.JComboBox();
        cmbMzMLScan2ID = new javax.swing.JComboBox<>();
        lblMzMLScan2ID = new javax.swing.JLabel();
        pPeakListLookup = new javax.swing.JPanel();
        flPeakLists = new org.rappsilber.gui.components.FileList();
        btnPeakListParse = new javax.swing.JButton();
        writeDB = new org.rappsilber.fdr.gui.components.WriteToDB();
        pLog = new javax.swing.JPanel();
        spLog = new javax.swing.JScrollPane();
        txtLog = new javax.swing.JTextArea();
        memory2 = new org.rappsilber.gui.components.memory.Memory();
        cbLevel = new javax.swing.JComboBox<>();
        pAbout = new javax.swing.JPanel();
        jTabbedPane4 = new javax.swing.JTabbedPane();
        pVersion = new javax.swing.JPanel();
        jLabel29 = new javax.swing.JLabel();
        txtXiFDRVersion = new javax.swing.JTextField();
        jScrollPane2 = new javax.swing.JScrollPane();
        txtchangelog = new javax.swing.JTextArea();
        jScrollPane3 = new javax.swing.JScrollPane();
        editAbout = new javax.swing.JEditorPane();
        jScrollPane5 = new javax.swing.JScrollPane();
        editAboutCSV = new javax.swing.JEditorPane();
        jScrollPane4 = new javax.swing.JScrollPane();
        editAboutMzIdentML = new javax.swing.JEditorPane();
        jSplitPane1 = new javax.swing.JSplitPane();
        txtStatus = new javax.swing.JTextField();
        memory3 = new org.rappsilber.gui.components.memory.Memory();

        cbCSVHeaders.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Item 1", "Item 2", "Item 3", "Item 4" }));

        cbCSVHeaderOptional.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Item 1", "Item 2", "Item 3", "Item 4" }));

        fdrSettingsSimple.setBorder(javax.swing.BorderFactory.createEtchedBorder());

        bgFDRSettingType.add(rbFDRSimple);
        rbFDRSimple.setText("Simple FDR");
        rbFDRSimple.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                changFDRSettingsInterface(evt);
            }
        });

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("Cross-Link FDR");
        addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                formFocusGained(evt);
            }
        });
        addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                formMouseClicked(evt);
            }
        });

        jTabbedPane1.setPreferredSize(new java.awt.Dimension(700, 400));

        javax.swing.GroupLayout pDatabaseLayout = new javax.swing.GroupLayout(pDatabase);
        pDatabase.setLayout(pDatabaseLayout);
        pDatabaseLayout.setHorizontalGroup(
            pDatabaseLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pDatabaseLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(getDBFDR, javax.swing.GroupLayout.DEFAULT_SIZE, 853, Short.MAX_VALUE)
                .addContainerGap())
        );
        pDatabaseLayout.setVerticalGroup(
            pDatabaseLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pDatabaseLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(getDBFDR, javax.swing.GroupLayout.DEFAULT_SIZE, 361, Short.MAX_VALUE)
                .addContainerGap())
        );

        tpInput.addTab("Database", pDatabase);

        csvSelect.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                csvSelectActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout pCSVInputLayout = new javax.swing.GroupLayout(pCSVInput);
        pCSVInput.setLayout(pCSVInputLayout);
        pCSVInputLayout.setHorizontalGroup(
            pCSVInputLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pCSVInputLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(csvSelect, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        pCSVInputLayout.setVerticalGroup(
            pCSVInputLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pCSVInputLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(csvSelect, javax.swing.GroupLayout.DEFAULT_SIZE, 361, Short.MAX_VALUE)
                .addContainerGap())
        );

        tpInput.addTab("CSV", pCSVInput);

        fbMZIdentMLIn.setDescription("MZIdentML-Files");
        fbMZIdentMLIn.setExtensions(new String[] {"mzid"});

        jLabel19.setText("mzid-File");

        btnReadMZIdent.setText("Read");
        btnReadMZIdent.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnReadMZIdentActionPerformed(evt);
            }
        });

        jLabel22.setText("PSM-Score");

        bgScoreDirectionMzIdentML.add(rbMZHighBetter);
        rbMZHighBetter.setSelected(true);
        rbMZHighBetter.setText("High better");
        rbMZHighBetter.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbMZHighBetterActionPerformed(evt);
            }
        });

        bgScoreDirectionMzIdentML.add(rbMZLowBetter);
        rbMZLowBetter.setText("Low better");
        rbMZLowBetter.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbMZLowBetterActionPerformed(evt);
            }
        });

        cbMZMatchScoreName.setEditable(true);
        cbMZMatchScoreName.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Score", "pvalue", "FDRScore", "xi:score", "MS:1001143", "search engine specific score for PSMs" }));
        cbMZMatchScoreName.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cbMZMatchScoreNameActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout pMZIdentMLInputLayout = new javax.swing.GroupLayout(pMZIdentMLInput);
        pMZIdentMLInput.setLayout(pMZIdentMLInputLayout);
        pMZIdentMLInputLayout.setHorizontalGroup(
            pMZIdentMLInputLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pMZIdentMLInputLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(pMZIdentMLInputLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(pMZIdentMLInputLayout.createSequentialGroup()
                        .addComponent(jLabel19)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fbMZIdentMLIn, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(btnReadMZIdent))
                    .addGroup(pMZIdentMLInputLayout.createSequentialGroup()
                        .addComponent(jLabel22)
                        .addGap(239, 239, 239)
                        .addGroup(pMZIdentMLInputLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(pMZIdentMLInputLayout.createSequentialGroup()
                                .addGroup(pMZIdentMLInputLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(rbMZLowBetter)
                                    .addComponent(rbMZHighBetter))
                                .addGap(88, 433, Short.MAX_VALUE))
                            .addComponent(cbMZMatchScoreName, javax.swing.GroupLayout.Alignment.TRAILING, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
                .addContainerGap())
        );
        pMZIdentMLInputLayout.setVerticalGroup(
            pMZIdentMLInputLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pMZIdentMLInputLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(pMZIdentMLInputLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel22)
                    .addComponent(cbMZMatchScoreName, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(1, 1, 1)
                .addComponent(rbMZHighBetter)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(rbMZLowBetter)
                .addGap(18, 18, 18)
                .addGroup(pMZIdentMLInputLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(pMZIdentMLInputLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(btnReadMZIdent)
                        .addComponent(jLabel19))
                    .addComponent(fbMZIdentMLIn, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(259, Short.MAX_VALUE))
        );

        tpInput.addTab("mzIdentML", pMZIdentMLInput);

        javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
        jPanel2.setLayout(jPanel2Layout);
        jPanel2Layout.setHorizontalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addComponent(tpInput)
                .addContainerGap())
        );
        jPanel2Layout.setVerticalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(tpInput)
                .addContainerGap())
        );

        jTabbedPane1.addTab("Input", jPanel2);

        pDatabseSize.setBorder(javax.swing.BorderFactory.createTitledBorder("Database Sizes"));
        pDatabseSize.setEnabled(false);

        spDecoyDBProt.setModel(new javax.swing.SpinnerNumberModel(9.99999999E8d, 0.001d, null, 1.0d));
        spDecoyDBProt.setEnabled(false);

        spTargetDBProt.setModel(new javax.swing.SpinnerNumberModel(9.99999999E8d, 0.001d, null, 1.0d));
        spTargetDBProt.setEnabled(false);

        lblProtein.setText("Protein");
        lblProtein.setEnabled(false);

        lblPeptide.setText("Peptide");
        lblPeptide.setEnabled(false);

        spTargetDB.setModel(new javax.swing.SpinnerNumberModel(9.99999999E8d, 0.001d, null, 1.0d));
        spTargetDB.setEnabled(false);

        spDecoyDB.setModel(new javax.swing.SpinnerNumberModel(9.99999999E8d, 0.001d, null, 1.0d));
        spDecoyDB.setEnabled(false);

        lblDecoyDB.setText("Decoy");
        lblDecoyDB.setEnabled(false);

        lblTargetDB.setText("Target");
        lblTargetDB.setEnabled(false);

        spDecoyDBLinks.setModel(new javax.swing.SpinnerNumberModel(9.99999999E8d, 0.001d, null, 1.0d));
        spDecoyDBLinks.setEnabled(false);

        spTargetDBLinks.setModel(new javax.swing.SpinnerNumberModel(9.99999999E8d, 0.001d, null, 1.0d));
        spTargetDBLinks.setEnabled(false);

        lblLinkDB.setText("Links");
        lblLinkDB.setEnabled(false);

        javax.swing.GroupLayout pDatabseSizeLayout = new javax.swing.GroupLayout(pDatabseSize);
        pDatabseSize.setLayout(pDatabseSizeLayout);
        pDatabseSizeLayout.setHorizontalGroup(
            pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pDatabseSizeLayout.createSequentialGroup()
                .addGap(78, 78, 78)
                .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(lblDecoyDB, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(lblTargetDB, javax.swing.GroupLayout.PREFERRED_SIZE, 56, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(spDecoyDB, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(lblPeptide, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(spTargetDB, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(spDecoyDBLinks, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(spTargetDBLinks, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(lblLinkDB, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(spDecoyDBProt, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(spTargetDBProt, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(lblProtein, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(384, Short.MAX_VALUE))
        );
        pDatabseSizeLayout.setVerticalGroup(
            pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pDatabseSizeLayout.createSequentialGroup()
                .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(lblPeptide)
                        .addComponent(lblProtein))
                    .addComponent(lblLinkDB))
                .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(spTargetDBLinks, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(pDatabseSizeLayout.createSequentialGroup()
                        .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(spTargetDB, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(lblTargetDB)
                            .addComponent(spTargetDBProt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(lblDecoyDB)
                            .addComponent(spDecoyDB, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addGroup(pDatabseSizeLayout.createSequentialGroup()
                        .addGap(26, 26, 26)
                        .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spDecoyDBLinks, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(spDecoyDBProt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))))
        );

        ckDBSize.setText("Define databse size");
        ckDBSize.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckDBSizeActionPerformed(evt);
            }
        });

        bgFDRSettingType.add(rbFDRMedium);
        rbFDRMedium.setSelected(true);
        rbFDRMedium.setText("Reduced");
        rbFDRMedium.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbFDRMediumActionPerformed(evt);
            }
        });

        bgFDRSettingType.add(rbFDRComplete);
        rbFDRComplete.setText("Complete FDR");
        rbFDRComplete.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbFDRCompleteActionPerformed(evt);
            }
        });

        spPepLength.setMinimumSize(new java.awt.Dimension(50, 50));
        spPepLength.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                spPepLengthFocusLost(evt);
            }
        });

        tblPepLength.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {
                { new Integer(0)},
                {null}
            },
            new String [] {
                "Peptide Length"
            }
        ) {
            Class[] types = new Class [] {
                java.lang.Integer.class
            };

            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }
        });
        spPepLength.setViewportView(tblPepLength);

        ckIgnoreGroups1.setText("Ignore Groups");
        ckIgnoreGroups1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckIgnoreGroups1ActionPerformed(evt);
            }
        });

        spDistanceGroup.setSpecialValueText("No distance group");

        jLabel1.setText("Sequnce Distance");

        cbGroupPepsBySelf.setText("Peptides by self-links support");

        cbGroupLinksBySelf.setText("Links by self-links support");

        cbGroupPPIBySelf.setText("PPIs by self-links support");

        cbGroupRun.setText("PSMs by Run");

        ckGroupByCrossLinkerStubs.setText("MS Cleavable XL Grouping");
        ckGroupByCrossLinkerStubs.setToolTipText("If the present PSMs will be groupped by presents of peptides with crosslinker stub fragments");

        javax.swing.GroupLayout pFDRGroupsLayout = new javax.swing.GroupLayout(pFDRGroups);
        pFDRGroups.setLayout(pFDRGroupsLayout);
        pFDRGroupsLayout.setHorizontalGroup(
            pFDRGroupsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pFDRGroupsLayout.createSequentialGroup()
                .addGroup(pFDRGroupsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(cbGroupPPIBySelf)
                    .addComponent(cbGroupLinksBySelf)
                    .addComponent(cbGroupRun)
                    .addComponent(spPepLength, javax.swing.GroupLayout.PREFERRED_SIZE, 226, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(pFDRGroupsLayout.createSequentialGroup()
                        .addComponent(jLabel1)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(spDistanceGroup, javax.swing.GroupLayout.PREFERRED_SIZE, 86, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(ckIgnoreGroups1, javax.swing.GroupLayout.PREFERRED_SIZE, 226, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(cbGroupPepsBySelf)
                    .addComponent(ckGroupByCrossLinkerStubs))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        pFDRGroupsLayout.setVerticalGroup(
            pFDRGroupsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pFDRGroupsLayout.createSequentialGroup()
                .addComponent(ckIgnoreGroups1)
                .addGap(18, 18, 18)
                .addComponent(spPepLength, javax.swing.GroupLayout.PREFERRED_SIZE, 69, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(cbGroupRun)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(cbGroupPepsBySelf)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(cbGroupLinksBySelf)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(cbGroupPPIBySelf)
                .addGap(4, 4, 4)
                .addComponent(ckGroupByCrossLinkerStubs)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(pFDRGroupsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(spDistanceGroup, javax.swing.GroupLayout.PREFERRED_SIZE, 27, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, Short.MAX_VALUE))
        );

        ckDefineGroups.setText("Define Groups");
        ckDefineGroups.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckDefineGroupsActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout spExtraLongLayout = new javax.swing.GroupLayout(spExtraLong);
        spExtraLong.setLayout(spExtraLongLayout);
        spExtraLongLayout.setHorizontalGroup(
            spExtraLongLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(spExtraLongLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(spExtraLongLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(pDatabseSize, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(spExtraLongLayout.createSequentialGroup()
                        .addComponent(ckDBSize)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(spExtraLongLayout.createSequentialGroup()
                        .addComponent(rbFDRMedium)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(rbFDRComplete)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(ckDefineGroups, javax.swing.GroupLayout.PREFERRED_SIZE, 134, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, spExtraLongLayout.createSequentialGroup()
                        .addComponent(spFDRSettingsWrapper)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(pFDRGroups, javax.swing.GroupLayout.PREFERRED_SIZE, 227, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap())
        );
        spExtraLongLayout.setVerticalGroup(
            spExtraLongLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(spExtraLongLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(spExtraLongLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(rbFDRComplete)
                    .addComponent(ckDefineGroups)
                    .addComponent(rbFDRMedium))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(spExtraLongLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(spFDRSettingsWrapper)
                    .addComponent(pFDRGroups, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(ckDBSize)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(pDatabseSize, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(20, 20, 20))
        );

        jTabbedPane1.addTab("FDR Settings", spExtraLong);

        jLabel12.setText("PSM input");

        jPanel9.setBorder(javax.swing.BorderFactory.createTitledBorder("AfterFDR"));

        btnPSMInfo.setText("+");
        btnPSMInfo.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnPSMInfoActionPerformed(evt);
            }
        });

        jPanel18.setLayout(new java.awt.GridLayout(1, 4, 15, 0));
        jPanel18.add(jLabel34);
        jPanel18.add(jLabel35);

        lblSumBetween1.setText("Cross-Linked");
        lblSumBetween1.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jPanel18.add(lblSumBetween1);

        lblSumInternal1.setText("Linear");
        jPanel18.add(lblSumInternal1);

        jPanel15.setLayout(new java.awt.GridLayout(1, 0, 10, 0));

        jLabel13.setText("PSMs");
        jPanel15.add(jLabel13);

        txtSumPSM.setMaximumSize(null);
        txtSumPSM.setMinimumSize(null);
        txtSumPSM.setPreferredSize(null);
        txtSumPSM.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtSumPSMActionPerformed(evt);
            }
        });
        jPanel15.add(txtSumPSM);

        txtSumPSMXL.setMaximumSize(null);
        txtSumPSMXL.setMinimumSize(null);
        txtSumPSMXL.setPreferredSize(null);
        txtSumPSMXL.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtSumPSMXLActionPerformed(evt);
            }
        });
        jPanel15.add(txtSumPSMXL);

        txtSumPSMLinear.setMaximumSize(null);
        txtSumPSMLinear.setMinimumSize(null);
        txtSumPSMLinear.setPreferredSize(null);
        txtSumPSMLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtSumPSMLinearActionPerformed(evt);
            }
        });
        jPanel15.add(txtSumPSMLinear);

        jPanel16.setLayout(new java.awt.GridLayout(1, 0, 10, 0));

        jLabel14.setText("Peptides");
        jPanel16.add(jLabel14);

        txtSumPepPairs.setMaximumSize(null);
        txtSumPepPairs.setMinimumSize(null);
        txtSumPepPairs.setPreferredSize(null);
        jPanel16.add(txtSumPepPairs);

        txtSumPepPairsXL.setMaximumSize(null);
        txtSumPepPairsXL.setMinimumSize(null);
        txtSumPepPairsXL.setPreferredSize(null);
        jPanel16.add(txtSumPepPairsXL);

        txtSumPepPairsLinear.setMaximumSize(null);
        txtSumPepPairsLinear.setMinimumSize(null);
        txtSumPepPairsLinear.setPreferredSize(null);
        jPanel16.add(txtSumPepPairsLinear);

        btnPepInfo.setText("+");
        btnPepInfo.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnPepInfoActionPerformed(evt);
            }
        });

        jPanel20.setLayout(new java.awt.GridLayout(1, 0, 10, 0));

        jLabel15.setText("Protein Groups");
        jPanel20.add(jLabel15);

        txtSumProtGroups.setMaximumSize(null);
        txtSumProtGroups.setMinimumSize(null);
        txtSumProtGroups.setPreferredSize(null);
        jPanel20.add(txtSumProtGroups);
        jPanel20.add(jLabel30);
        jPanel20.add(jLabel31);

        btnProtInfo.setText("+");
        btnProtInfo.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnProtInfoActionPerformed(evt);
            }
        });

        jPanel12.setLayout(new java.awt.GridLayout(1, 4, 15, 0));
        jPanel12.add(jLabel33);
        jPanel12.add(jLabel32);

        lblSumBetween.setText("Between");
        lblSumBetween.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jPanel12.add(lblSumBetween);

        lblSumInternal.setText("Within");
        jPanel12.add(lblSumInternal);

        jPanel21.setLayout(new java.awt.GridLayout(1, 0, 10, 0));

        jLabel16.setText("Residue Pairs");
        jPanel21.add(jLabel16);

        txtSumLinks.setMaximumSize(null);
        txtSumLinks.setMinimumSize(null);
        txtSumLinks.setPreferredSize(null);
        jPanel21.add(txtSumLinks);

        txtSumLinksBetween.setMaximumSize(null);
        txtSumLinksBetween.setMinimumSize(null);
        txtSumLinksBetween.setPreferredSize(null);
        jPanel21.add(txtSumLinksBetween);

        txtSumLinksInternal.setMaximumSize(null);
        txtSumLinksInternal.setMinimumSize(null);
        txtSumLinksInternal.setPreferredSize(null);
        jPanel21.add(txtSumLinksInternal);

        btnLinkInfo.setText("+");
        btnLinkInfo.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnLinkInfoActionPerformed(evt);
            }
        });

        jPanel22.setLayout(new java.awt.GridLayout(1, 0, 10, 0));

        jLabel17.setText("Protein Group Pairs");
        jPanel22.add(jLabel17);

        txtSumProtGroupPairs.setMaximumSize(null);
        txtSumProtGroupPairs.setMinimumSize(null);
        txtSumProtGroupPairs.setPreferredSize(null);
        jPanel22.add(txtSumProtGroupPairs);

        txtSumProtGroupPairsBetween.setMaximumSize(null);
        txtSumProtGroupPairsBetween.setMinimumSize(null);
        txtSumProtGroupPairsBetween.setPreferredSize(null);
        jPanel22.add(txtSumProtGroupPairsBetween);

        txtSumProtGroupPairsInternal.setMaximumSize(null);
        txtSumProtGroupPairsInternal.setMinimumSize(null);
        txtSumProtGroupPairsInternal.setPreferredSize(null);
        jPanel22.add(txtSumProtGroupPairsInternal);

        btnPPIInfo.setText("+");
        btnPPIInfo.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnPPIInfoActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel9Layout = new javax.swing.GroupLayout(jPanel9);
        jPanel9.setLayout(jPanel9Layout);
        jPanel9Layout.setHorizontalGroup(
            jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel9Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel9Layout.createSequentialGroup()
                        .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(jPanel22, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE)
                            .addComponent(jPanel21, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jPanel12, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(btnLinkInfo)
                            .addComponent(btnPPIInfo)))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel9Layout.createSequentialGroup()
                        .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(jPanel18, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jPanel20, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, 698, Short.MAX_VALUE)
                            .addComponent(jPanel16, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jPanel15, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(btnPSMInfo)
                            .addComponent(btnPepInfo)
                            .addComponent(btnProtInfo))))
                .addContainerGap())
        );
        jPanel9Layout.setVerticalGroup(
            jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel9Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jPanel18, javax.swing.GroupLayout.PREFERRED_SIZE, 20, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(btnProtInfo)
                    .addGroup(jPanel9Layout.createSequentialGroup()
                        .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(btnPSMInfo, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jPanel15, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(jPanel16, javax.swing.GroupLayout.PREFERRED_SIZE, 25, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(btnPepInfo))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jPanel20, javax.swing.GroupLayout.PREFERRED_SIZE, 25, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jPanel12, javax.swing.GroupLayout.PREFERRED_SIZE, 18, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                    .addComponent(btnLinkInfo, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jPanel21, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED, 32, Short.MAX_VALUE)
                .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(btnPPIInfo, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jPanel22, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );

        javax.swing.GroupLayout jPanel8Layout = new javax.swing.GroupLayout(jPanel8);
        jPanel8.setLayout(jPanel8Layout);
        jPanel8Layout.setHorizontalGroup(
            jPanel8Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel8Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel8Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(jLabel12, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(txtSumInput))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jPanel9, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        jPanel8Layout.setVerticalGroup(
            jPanel8Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel8Layout.createSequentialGroup()
                .addGroup(jPanel8Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel8Layout.createSequentialGroup()
                        .addComponent(jLabel12)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(txtSumInput, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(jPanel9, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        tpResult.addTab("Summary", jPanel8);

        bgSeparator.add(rbTSV);
        rbTSV.setText("Tab Separated");
        rbTSV.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbTSVActionPerformed(evt);
            }
        });

        bgSeparator.add(rbCSV);
        rbCSV.setSelected(true);
        rbCSV.setText("Comma Separated");
        rbCSV.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbCSVActionPerformed(evt);
            }
        });

        btnWrite.setText("Write");
        btnWrite.setEnabled(false);
        btnWrite.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnWriteActionPerformed(evt);
            }
        });

        fbFolder.setAutoAddDefaultExtension("csv");
        fbFolder.setExtensions(new String[] {"tsv", "txt", "csv"});
        fbFolder.setLoad(false);

        jLabel9.setText("Output");

        ckPrePostAA.setText("Pre and post amino-acids");
        ckPrePostAA.setEnabled(false);

        lpCsvOutLocal.setDefaultLocal(java.util.Locale.ENGLISH);
        lpCsvOutLocal.setMinimumSize(new java.awt.Dimension(80, 27));
        lpCsvOutLocal.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                lpCsvOutLocalActionPerformed(evt);
            }
        });

        jLabel28.setText("Language");

        btnCalcRanges.setText("Calculate Ranges");
        btnCalcRanges.setEnabled(false);
        btnCalcRanges.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnCalcRangesActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel10Layout = new javax.swing.GroupLayout(jPanel10);
        jPanel10.setLayout(jPanel10Layout);
        jPanel10Layout.setHorizontalGroup(
            jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel10Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel9)
                    .addComponent(jLabel28))
                .addGap(26, 26, 26)
                .addGroup(jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(fbFolder, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(jPanel10Layout.createSequentialGroup()
                        .addComponent(rbTSV)
                        .addGap(95, 95, 95)
                        .addComponent(ckPrePostAA)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(jPanel10Layout.createSequentialGroup()
                        .addComponent(rbCSV)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(btnCalcRanges)
                        .addGap(18, 18, 18)
                        .addComponent(btnWrite))
                    .addComponent(lpCsvOutLocal, javax.swing.GroupLayout.DEFAULT_SIZE, 757, Short.MAX_VALUE))
                .addGap(23, 23, 23))
        );
        jPanel10Layout.setVerticalGroup(
            jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel10Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(fbFolder, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel9))
                .addGap(18, 18, 18)
                .addGroup(jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(lpCsvOutLocal, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel28))
                .addGap(10, 10, 10)
                .addGroup(jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(rbTSV)
                    .addComponent(ckPrePostAA))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(rbCSV)
                    .addComponent(btnWrite)
                    .addComponent(btnCalcRanges))
                .addContainerGap(243, Short.MAX_VALUE))
        );

        tpResult.addTab("CSV/TSV", jPanel10);

        fbMzIdentMLOut.setDescription("mzIdentML files");
        fbMzIdentMLOut.setExtensions(new String[] {"mzid"});
        fbMzIdentMLOut.setLoad(false);
        fbMzIdentMLOut.setLocalPropertyKey("LastAccessedmzIdentMLFolder");

        btnWriteMzIdentML.setText("Write");
        btnWriteMzIdentML.setEnabled(false);
        btnWriteMzIdentML.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnWriteMzIdentMLActionPerformed(evt);
            }
        });

        txtmzIdentOwnerFirst.setText("First");
        txtmzIdentOwnerFirst.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtmzIdentOwnerFirstActionPerformed(evt);
            }
        });

        jLabel2.setText("Document Owner");

        txtmzIdentOwnerLast.setText("Last");
        txtmzIdentOwnerLast.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtmzIdentOwnerLastActionPerformed(evt);
            }
        });

        jLabel3.setText("E-Mail");

        txtmzIdentOwnerEmail.setText("reseacher@organisation.org");
        txtmzIdentOwnerEmail.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtmzIdentOwnerEmailActionPerformed(evt);
            }
        });

        txtmzIdentAdress.setColumns(20);
        txtmzIdentAdress.setRows(5);
        jScrollPane1.setViewportView(txtmzIdentAdress);

        jLabel4.setText("Address");

        jLabel5.setText("Organisation");

        txtmzIdentOwnerOrg.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtmzIdentOwnerOrgActionPerformed(evt);
            }
        });

        lblPeaklistExtension.setText("Peak-list file extension");

        cmbPeakListFormat.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "mgf", "mzML", "raw", " " }));

        cmbMzMLScan2ID.setEditable(true);
        cmbMzMLScan2ID.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "controllerType=0 controllerNumber=1 scan=%s%", "controllerType=0 controllerNumber=1 scan=%s-1%", "scan=%s%", "scan=%s-1%", "%s%" }));
        cmbMzMLScan2ID.setToolTipText("%s% will be replaced by the actuall scan-number\n%s-1% will be replaced by scannumber-1\n%s+1% will be replaced by scannumber+1\n");

        lblMzMLScan2ID.setText("mzML spectrumID template");

        javax.swing.GroupLayout jPanel14Layout = new javax.swing.GroupLayout(jPanel14);
        jPanel14.setLayout(jPanel14Layout);
        jPanel14Layout.setHorizontalGroup(
            jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel14Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(fbMzIdentMLOut, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(jPanel14Layout.createSequentialGroup()
                        .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel2)
                            .addComponent(jLabel3)
                            .addComponent(jLabel4)
                            .addComponent(jLabel5)
                            .addComponent(lblPeaklistExtension)
                            .addComponent(lblMzMLScan2ID))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(jPanel14Layout.createSequentialGroup()
                                .addComponent(cmbMzMLScan2ID, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addGap(30, 30, 30))
                            .addGroup(jPanel14Layout.createSequentialGroup()
                                .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                    .addComponent(jScrollPane1)
                                    .addGroup(jPanel14Layout.createSequentialGroup()
                                        .addComponent(txtmzIdentOwnerFirst, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addGap(18, 18, 18)
                                        .addComponent(txtmzIdentOwnerLast, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE))
                                    .addComponent(txtmzIdentOwnerEmail)
                                    .addComponent(txtmzIdentOwnerOrg)
                                    .addComponent(cmbPeakListFormat, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED, 350, Short.MAX_VALUE)
                                .addComponent(btnWriteMzIdentML)))))
                .addContainerGap())
        );
        jPanel14Layout.setVerticalGroup(
            jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel14Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(fbMzIdentMLOut, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(btnWriteMzIdentML)
                    .addComponent(jLabel2)
                    .addComponent(txtmzIdentOwnerFirst, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(txtmzIdentOwnerLast, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel3)
                    .addComponent(txtmzIdentOwnerEmail, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(txtmzIdentOwnerOrg, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel5))
                .addGap(9, 9, 9)
                .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel4))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(lblPeaklistExtension)
                    .addComponent(cmbPeakListFormat, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(cmbMzMLScan2ID, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(lblMzMLScan2ID)))
        );

        tpResult.addTab("mzIdentML", jPanel14);

        flPeakLists.setToolTipText("if peak list file and peak list name are not readable from the input these can be suplemented here");
        flPeakLists.setDescription("PeakLists");
        flPeakLists.setExtensions(new String[] {".apl", ".mgf"});

        btnPeakListParse.setText("Parse");
        btnPeakListParse.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnPeakListParseActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout pPeakListLookupLayout = new javax.swing.GroupLayout(pPeakListLookup);
        pPeakListLookup.setLayout(pPeakListLookupLayout);
        pPeakListLookupLayout.setHorizontalGroup(
            pPeakListLookupLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, pPeakListLookupLayout.createSequentialGroup()
                .addContainerGap(805, Short.MAX_VALUE)
                .addComponent(btnPeakListParse)
                .addContainerGap())
            .addGroup(pPeakListLookupLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(pPeakListLookupLayout.createSequentialGroup()
                    .addContainerGap()
                    .addComponent(flPeakLists, javax.swing.GroupLayout.DEFAULT_SIZE, 865, Short.MAX_VALUE)
                    .addContainerGap()))
        );
        pPeakListLookupLayout.setVerticalGroup(
            pPeakListLookupLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, pPeakListLookupLayout.createSequentialGroup()
                .addContainerGap(348, Short.MAX_VALUE)
                .addComponent(btnPeakListParse)
                .addContainerGap())
            .addGroup(pPeakListLookupLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(pPeakListLookupLayout.createSequentialGroup()
                    .addContainerGap()
                    .addComponent(flPeakLists, javax.swing.GroupLayout.PREFERRED_SIZE, 284, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addContainerGap(89, Short.MAX_VALUE)))
        );

        tpResult.addTab("PeakListLookup", pPeakListLookup);
        tpResult.addTab("To Database", writeDB);

        javax.swing.GroupLayout pResultLayout = new javax.swing.GroupLayout(pResult);
        pResult.setLayout(pResultLayout);
        pResultLayout.setHorizontalGroup(
            pResultLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(tpResult, javax.swing.GroupLayout.Alignment.TRAILING)
        );
        pResultLayout.setVerticalGroup(
            pResultLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pResultLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(tpResult)
                .addContainerGap())
        );

        jTabbedPane1.addTab("Results", pResult);

        txtLog.setColumns(20);
        txtLog.setRows(5);
        spLog.setViewportView(txtLog);

        cbLevel.setModel(new javax.swing.DefaultComboBoxModel<Level>(new Level[] { Level.ALL,Level.FINEST,Level.FINER,Level.FINE,Level.INFO,Level.WARNING,Level.SEVERE,Level.OFF}));
        cbLevel.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cbLevelActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout pLogLayout = new javax.swing.GroupLayout(pLog);
        pLog.setLayout(pLogLayout);
        pLogLayout.setHorizontalGroup(
            pLogLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pLogLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(pLogLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(spLog)
                    .addComponent(memory2, javax.swing.GroupLayout.DEFAULT_SIZE, 870, Short.MAX_VALUE)
                    .addComponent(cbLevel, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );
        pLogLayout.setVerticalGroup(
            pLogLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, pLogLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(memory2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(9, 9, 9)
                .addComponent(cbLevel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(spLog, javax.swing.GroupLayout.DEFAULT_SIZE, 348, Short.MAX_VALUE)
                .addContainerGap())
        );

        jTabbedPane1.addTab("Log", pLog);

        jLabel29.setText("Version");

        txtXiFDRVersion.setEditable(false);

        txtchangelog.setColumns(20);
        txtchangelog.setRows(5);
        jScrollPane2.setViewportView(txtchangelog);

        javax.swing.GroupLayout pVersionLayout = new javax.swing.GroupLayout(pVersion);
        pVersion.setLayout(pVersionLayout);
        pVersionLayout.setHorizontalGroup(
            pVersionLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pVersionLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(pVersionLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jScrollPane2, javax.swing.GroupLayout.DEFAULT_SIZE, 865, Short.MAX_VALUE)
                    .addGroup(pVersionLayout.createSequentialGroup()
                        .addComponent(jLabel29)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(txtXiFDRVersion, javax.swing.GroupLayout.PREFERRED_SIZE, 119, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap())
        );
        pVersionLayout.setVerticalGroup(
            pVersionLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pVersionLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(pVersionLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel29)
                    .addComponent(txtXiFDRVersion, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jScrollPane2, javax.swing.GroupLayout.DEFAULT_SIZE, 336, Short.MAX_VALUE)
                .addContainerGap())
        );

        jTabbedPane4.addTab("Version", pVersion);

        editAbout.setEditable(false);
        jScrollPane3.setViewportView(editAbout);

        jTabbedPane4.addTab("General", jScrollPane3);

        jScrollPane5.setViewportView(editAboutCSV);

        jTabbedPane4.addTab("CSV-Import", jScrollPane5);

        jScrollPane4.setViewportView(editAboutMzIdentML);

        jTabbedPane4.addTab("mzIdentML-import", jScrollPane4);

        javax.swing.GroupLayout pAboutLayout = new javax.swing.GroupLayout(pAbout);
        pAbout.setLayout(pAboutLayout);
        pAboutLayout.setHorizontalGroup(
            pAboutLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jTabbedPane4)
        );
        pAboutLayout.setVerticalGroup(
            pAboutLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pAboutLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jTabbedPane4)
                .addContainerGap())
        );

        jTabbedPane1.addTab("About", pAbout);

        jScrollPane6.setViewportView(jTabbedPane1);

        jSplitPane1.setDividerLocation(600);
        jSplitPane1.setDividerSize(5);

        txtStatus.setEditable(false);
        txtStatus.setText("status");
        txtStatus.setPreferredSize(new java.awt.Dimension(200, 19));
        jSplitPane1.setLeftComponent(txtStatus);

        memory3.setShowGCButton(false);
        jSplitPane1.setRightComponent(memory3);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane6, javax.swing.GroupLayout.DEFAULT_SIZE, 902, Short.MAX_VALUE)
            .addComponent(jSplitPane1)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane6, javax.swing.GroupLayout.DEFAULT_SIZE, 466, Short.MAX_VALUE)
                .addGap(18, 18, 18)
                .addComponent(jSplitPane1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void btnWriteActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnWriteActionPerformed
        writeResults();
    }//GEN-LAST:event_btnWriteActionPerformed

    private void txtSumPSMActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtSumPSMActionPerformed
    }//GEN-LAST:event_txtSumPSMActionPerformed

    private void btnReadMZIdentActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnReadMZIdentActionPerformed
        readMZIdentML();
        ckPrePostAA.setVisible(false);

    }//GEN-LAST:event_btnReadMZIdentActionPerformed

    private void cbMZMatchScoreNameActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cbMZMatchScoreNameActionPerformed

        if (!mzScoreDirectionChanged) {
            String name = cbMZMatchScoreName.getModel().getSelectedItem().toString();
            if (name.toLowerCase().contains("pvalue") || name.toLowerCase().contains("fdr")) {
                rbMZLowBetter.setSelected(true);
            }
        }

    }//GEN-LAST:event_cbMZMatchScoreNameActionPerformed

    private void rbMZHighBetterActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbMZHighBetterActionPerformed
        mzScoreDirectionChanged = true;
    }//GEN-LAST:event_rbMZHighBetterActionPerformed

    private void rbMZLowBetterActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbMZLowBetterActionPerformed
        mzScoreDirectionChanged = true;
    }//GEN-LAST:event_rbMZLowBetterActionPerformed

    private void btnWriteMzIdentMLActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnWriteMzIdentMLActionPerformed
        exportMZIdentML();
        LocalProperties.setProperty(MZIdentMLOwner.propertyFirst, txtmzIdentOwnerFirst.getText());
        LocalProperties.setProperty(MZIdentMLOwner.propertyLast, txtmzIdentOwnerLast.getText());
        LocalProperties.setProperty(MZIdentMLOwner.propertyEMail, txtmzIdentOwnerEmail.getText());
        LocalProperties.setProperty(MZIdentMLOwner.propertyAddress, txtmzIdentAdress.getText());
        LocalProperties.setProperty(MZIdentMLOwner.propertyOrg, txtmzIdentOwnerOrg.getText());

    }//GEN-LAST:event_btnWriteMzIdentMLActionPerformed

    private void txtSumPSMXLActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtSumPSMXLActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_txtSumPSMXLActionPerformed

    private void txtSumPSMLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtSumPSMLinearActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_txtSumPSMLinearActionPerformed

    private void spPepLengthFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_spPepLengthFocusLost
    }//GEN-LAST:event_spPepLengthFocusLost

    private void ckDefineGroupsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckDefineGroupsActionPerformed

        pFDRGroups.setVisible(ckDefineGroups.isSelected());
        tblPepLength.setEnabled(!ckIgnoreGroups1.isSelected());
        spPepLength.setEnabled(!ckIgnoreGroups1.isSelected());

    }//GEN-LAST:event_ckDefineGroupsActionPerformed

    private void formFocusGained(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_formFocusGained


    }//GEN-LAST:event_formFocusGained

    private void formMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_formMouseClicked
        // TODO add your handling code here:
        this.toFront();
    }//GEN-LAST:event_formMouseClicked

    private void btnPSMInfoActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnPSMInfoActionPerformed
        if (getResult() != null) {
            new FDRLevelInformations(getResult().psmFDR, "PSM FDR").setVisible(true);
        }
    }//GEN-LAST:event_btnPSMInfoActionPerformed

    private void btnPepInfoActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnPepInfoActionPerformed
        if (getResult() != null) {
            new FDRLevelInformations(getResult().peptidePairFDR, "Peptide-Pair FDR").setVisible(true);
        }
    }//GEN-LAST:event_btnPepInfoActionPerformed

    private void btnProtInfoActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnProtInfoActionPerformed
        if (getResult() != null) {
            new FDRLevelInformations(getResult().proteinGroupFDR, "Protein Group FDR").setVisible(true);
        }
    }//GEN-LAST:event_btnProtInfoActionPerformed

    private void btnLinkInfoActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnLinkInfoActionPerformed
        if (getResult() != null) {
            new FDRLevelInformations(getResult().proteinGroupLinkFDR, "Protein Group Link FDR").setVisible(true);
        }
    }//GEN-LAST:event_btnLinkInfoActionPerformed

    private void btnPPIInfoActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnPPIInfoActionPerformed
        if (getResult() != null) {
            new FDRLevelInformations(getResult().proteinGroupPairFDR, "Protein Group Pairs FDR").setVisible(true);
        }
    }//GEN-LAST:event_btnPPIInfoActionPerformed

    private void ckDBSizeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckDBSizeActionPerformed

        pDatabseSize.setVisible(ckDBSize.isSelected());

        spTargetDB.setEnabled(ckDBSize.isSelected());
        spDecoyDB.setEnabled(ckDBSize.isSelected());
        spTargetDBLinks.setEnabled(ckDBSize.isSelected());
        spDecoyDBLinks.setEnabled(ckDBSize.isSelected());
        spTargetDBProt.setEnabled(ckDBSize.isSelected());
        spDecoyDBProt.setEnabled(ckDBSize.isSelected());
        lblDecoyDB.setEnabled(ckDBSize.isSelected());
        lblTargetDB.setEnabled(ckDBSize.isSelected());
        lblLinkDB.setEnabled(ckDBSize.isSelected());
        lblPeptide.setEnabled(ckDBSize.isSelected());
        lblProtein.setEnabled(ckDBSize.isSelected());

        spTargetDB.setVisible(ckDBSize.isSelected());
        spDecoyDB.setVisible(ckDBSize.isSelected());
        spTargetDBLinks.setVisible(ckDBSize.isSelected());
        spDecoyDBLinks.setVisible(ckDBSize.isSelected());
        spTargetDBProt.setVisible(ckDBSize.isSelected());
        spDecoyDBProt.setVisible(ckDBSize.isSelected());
        lblDecoyDB.setVisible(ckDBSize.isSelected());
        lblTargetDB.setVisible(ckDBSize.isSelected());
        lblLinkDB.setVisible(ckDBSize.isSelected());
        lblPeptide.setVisible(ckDBSize.isSelected());
        lblProtein.setVisible(ckDBSize.isSelected());

    }//GEN-LAST:event_ckDBSizeActionPerformed

    private void changFDRSettingsInterface(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_changFDRSettingsInterface

        changeFDRSettings(evt);

    }//GEN-LAST:event_changFDRSettingsInterface

    private void rbFDRCompleteActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbFDRCompleteActionPerformed
        changeFDRSettings(evt);
    }//GEN-LAST:event_rbFDRCompleteActionPerformed

    private void ckIgnoreGroups1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckIgnoreGroups1ActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_ckIgnoreGroups1ActionPerformed

    private void csvSelectAddActionPerformed(java.awt.event.ActionEvent evt) {
        readAdditionalCSV();
    }

    private void csvSelectActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_csvSelectActionPerformed
        readAllCSV();
    }//GEN-LAST:event_csvSelectActionPerformed

    private void txtmzIdentOwnerFirstActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtmzIdentOwnerFirstActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_txtmzIdentOwnerFirstActionPerformed

    private void txtmzIdentOwnerLastActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtmzIdentOwnerLastActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_txtmzIdentOwnerLastActionPerformed

    private void txtmzIdentOwnerEmailActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtmzIdentOwnerEmailActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_txtmzIdentOwnerEmailActionPerformed

    private void txtmzIdentOwnerOrgActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtmzIdentOwnerOrgActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_txtmzIdentOwnerOrgActionPerformed

    private void rbTSVActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbTSVActionPerformed
        if (rbTSV.isSelected()) {
            fbFolder.setAutoAddDefaultExtension("txt");
        }
    }//GEN-LAST:event_rbTSVActionPerformed

    private void rbCSVActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbCSVActionPerformed
        if (rbCSV.isSelected()) {
            fbFolder.setAutoAddDefaultExtension("csv");
        }
    }//GEN-LAST:event_rbCSVActionPerformed

    private void cbLevelActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cbLevelActionPerformed
        Runnable runnable = new Runnable() {
            public void run() {
                Handler[] handlers
                        = Logger.getLogger("").getHandlers();
                for (int index = 0; index < handlers.length; index++) {
                    handlers[index].setLevel((Level) cbLevel.getSelectedItem());
                }
                loggingOutput.setLevel((Level) cbLevel.getSelectedItem());
            }
        };
        Thread t = new Thread(runnable);
        t.setName("Thread-" + t.getId()+ " UpdateLogLevel");
        t.setDaemon(true);
        t.start();
    }//GEN-LAST:event_cbLevelActionPerformed

    private void lpCsvOutLocalActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_lpCsvOutLocalActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_lpCsvOutLocalActionPerformed

    private void btnPeakListParseActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnPeakListParseActionPerformed
        Runnable runnable = new Runnable() {
            public void run() {
                try {
                    parseSuppliedMGF();
                } catch (IOException ex) {
                    Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, "Error parsing peaklist", ex);
                    setStatus("Error parsing PeakList");
                } finally {
                    SwingUtilities.invokeLater(new Runnable() {
                        @Override
                        public void run() {
                            btnPeakListParse.setEnabled(true);
                        }
                    });
                }
            }
        };
        btnPeakListParse.setEnabled(false);
        Thread peaklistparser = new Thread(runnable);
        peaklistparser.setName("PeakListParser" + peaklistparser.getName());
        peaklistparser.start();
    }//GEN-LAST:event_btnPeakListParseActionPerformed

    private void rbFDRMediumActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbFDRMediumActionPerformed
        changeFDRSettings(evt);
    }//GEN-LAST:event_rbFDRMediumActionPerformed

    private void btnCalcRangesActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnCalcRangesActionPerformed
        if (fbFolder.getFile() != null) {
            CalculateRanges cr = new CalculateRanges(this);
            cr.setVisible(true);
        } else {
            setStatus("No output file defined");
        }
    }//GEN-LAST:event_btnCalcRangesActionPerformed

//    private void fdrSpinnerMaximumCheck(JSpinner sp, double max) {                                   
//        SpinnerModel sm = sp.getModel();
//        Double fdr = (Double) sm.getValue();
//        if (fdr > max)
//            sm.setValue(max);
//    }                                  
    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        
        startGUI();
        
    }
    
    public static FDRGUI startGUI() {

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
            java.util.logging.Logger.getLogger(FDRGUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(FDRGUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(FDRGUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(FDRGUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        //</editor-fold>
        
        final FDRGUI gui =  new FDRGUI();
        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                gui.setVisible(true);
            }
        });
        return gui;
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.ButtonGroup bgFDRSettingType;
    private javax.swing.ButtonGroup bgScoreDirectionMzIdentML;
    private javax.swing.ButtonGroup bgSeparator;
    private javax.swing.JButton btnCalcRanges;
    private javax.swing.JButton btnLinkInfo;
    private javax.swing.JButton btnPPIInfo;
    private javax.swing.JButton btnPSMInfo;
    private javax.swing.JButton btnPeakListParse;
    private javax.swing.JButton btnPepInfo;
    private javax.swing.JButton btnProtInfo;
    private javax.swing.JButton btnReadMZIdent;
    protected javax.swing.JButton btnWrite;
    protected javax.swing.JButton btnWriteMzIdentML;
    private javax.swing.JComboBox cbCSVHeaderOptional;
    private javax.swing.JComboBox cbCSVHeaders;
    private javax.swing.JCheckBox cbGroupLinksBySelf;
    private javax.swing.JCheckBox cbGroupPPIBySelf;
    private javax.swing.JCheckBox cbGroupPepsBySelf;
    private javax.swing.JCheckBox cbGroupRun;
    private javax.swing.JComboBox<Level> cbLevel;
    private javax.swing.JComboBox cbMZMatchScoreName;
    private javax.swing.JCheckBox ckDBSize;
    private javax.swing.JCheckBox ckDefineGroups;
    private javax.swing.JCheckBox ckGroupByCrossLinkerStubs;
    private javax.swing.JCheckBox ckIgnoreGroups1;
    public javax.swing.JCheckBox ckPrePostAA;
    protected javax.swing.JComboBox<String> cmbMzMLScan2ID;
    public javax.swing.JComboBox cmbPeakListFormat;
    protected org.rappsilber.fdr.gui.components.CSVSelection csvSelect;
    private javax.swing.JEditorPane editAbout;
    private javax.swing.JEditorPane editAboutCSV;
    private javax.swing.JEditorPane editAboutMzIdentML;
    private org.rappsilber.gui.components.FileBrowser fbFolder;
    private org.rappsilber.gui.components.FileBrowser fbMZIdentMLIn;
    private org.rappsilber.gui.components.FileBrowser fbMzIdentMLOut;
    private org.rappsilber.fdr.gui.components.FDRSettingsComplete fdrSettingsComplete;
    private org.rappsilber.fdr.gui.components.FDRSettingsMedium fdrSettingsMedium;
    private org.rappsilber.fdr.gui.components.FDRSettingsSimple fdrSettingsSimple;
    private org.rappsilber.gui.components.FileList flPeakLists;
    private org.rappsilber.fdr.gui.components.GetDBFDR getDBFDR;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel12;
    private javax.swing.JLabel jLabel13;
    private javax.swing.JLabel jLabel14;
    private javax.swing.JLabel jLabel15;
    private javax.swing.JLabel jLabel16;
    private javax.swing.JLabel jLabel17;
    private javax.swing.JLabel jLabel19;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel22;
    private javax.swing.JLabel jLabel28;
    private javax.swing.JLabel jLabel29;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel30;
    private javax.swing.JLabel jLabel31;
    private javax.swing.JLabel jLabel32;
    private javax.swing.JLabel jLabel33;
    private javax.swing.JLabel jLabel34;
    private javax.swing.JLabel jLabel35;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel9;
    private javax.swing.JPanel jPanel10;
    private javax.swing.JPanel jPanel12;
    private javax.swing.JPanel jPanel14;
    private javax.swing.JPanel jPanel15;
    private javax.swing.JPanel jPanel16;
    private javax.swing.JPanel jPanel18;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel20;
    private javax.swing.JPanel jPanel21;
    private javax.swing.JPanel jPanel22;
    private javax.swing.JPanel jPanel8;
    private javax.swing.JPanel jPanel9;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JScrollPane jScrollPane3;
    private javax.swing.JScrollPane jScrollPane4;
    private javax.swing.JScrollPane jScrollPane5;
    private javax.swing.JScrollPane jScrollPane6;
    private javax.swing.JSplitPane jSplitPane1;
    private javax.swing.JTabbedPane jTabbedPane1;
    private javax.swing.JTabbedPane jTabbedPane4;
    private javax.swing.JLabel lblDecoyDB;
    private javax.swing.JLabel lblLinkDB;
    protected javax.swing.JLabel lblMzMLScan2ID;
    public javax.swing.JLabel lblPeaklistExtension;
    private javax.swing.JLabel lblPeptide;
    private javax.swing.JLabel lblProtein;
    private javax.swing.JLabel lblSumBetween;
    private javax.swing.JLabel lblSumBetween1;
    private javax.swing.JLabel lblSumInternal;
    private javax.swing.JLabel lblSumInternal1;
    private javax.swing.JLabel lblTargetDB;
    private org.rappsilber.gui.components.LocalPicker lpCsvOutLocal;
    private org.rappsilber.gui.components.memory.Memory memory2;
    private org.rappsilber.gui.components.memory.Memory memory3;
    private javax.swing.JPanel pAbout;
    private javax.swing.JPanel pCSVInput;
    private javax.swing.JPanel pDatabase;
    private javax.swing.JPanel pDatabseSize;
    private javax.swing.JPanel pFDRGroups;
    private javax.swing.JPanel pLog;
    private javax.swing.JPanel pMZIdentMLInput;
    private javax.swing.JPanel pPeakListLookup;
    private javax.swing.JPanel pResult;
    private javax.swing.JPanel pVersion;
    private javax.swing.JRadioButton rbCSV;
    private javax.swing.JRadioButton rbFDRComplete;
    private javax.swing.JRadioButton rbFDRMedium;
    private javax.swing.JRadioButton rbFDRSimple;
    private javax.swing.JRadioButton rbMZHighBetter;
    private javax.swing.JRadioButton rbMZLowBetter;
    private javax.swing.JRadioButton rbTSV;
    private javax.swing.JSpinner spDecoyDB;
    private javax.swing.JSpinner spDecoyDBLinks;
    private javax.swing.JSpinner spDecoyDBProt;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner spDistanceGroup;
    private javax.swing.JPanel spExtraLong;
    private javax.swing.JScrollPane spFDRSettingsWrapper;
    private javax.swing.JScrollPane spLog;
    private javax.swing.JScrollPane spPepLength;
    private javax.swing.JSpinner spTargetDB;
    private javax.swing.JSpinner spTargetDBLinks;
    private javax.swing.JSpinner spTargetDBProt;
    private javax.swing.JTable tblPepLength;
    protected javax.swing.JTabbedPane tpInput;
    protected javax.swing.JTabbedPane tpResult;
    private javax.swing.JTextArea txtLog;
    private javax.swing.JTextField txtStatus;
    private javax.swing.JTextField txtSumInput;
    public javax.swing.JTextField txtSumLinks;
    public javax.swing.JTextField txtSumLinksBetween;
    public javax.swing.JTextField txtSumLinksInternal;
    public javax.swing.JTextField txtSumPSM;
    public javax.swing.JTextField txtSumPSMLinear;
    public javax.swing.JTextField txtSumPSMXL;
    public javax.swing.JTextField txtSumPepPairs;
    public javax.swing.JTextField txtSumPepPairsLinear;
    public javax.swing.JTextField txtSumPepPairsXL;
    public javax.swing.JTextField txtSumProtGroupPairs;
    public javax.swing.JTextField txtSumProtGroupPairsBetween;
    public javax.swing.JTextField txtSumProtGroupPairsInternal;
    public javax.swing.JTextField txtSumProtGroups;
    private javax.swing.JTextField txtXiFDRVersion;
    private javax.swing.JTextArea txtchangelog;
    private javax.swing.JTextArea txtmzIdentAdress;
    private javax.swing.JTextField txtmzIdentOwnerEmail;
    private javax.swing.JTextField txtmzIdentOwnerFirst;
    private javax.swing.JTextField txtmzIdentOwnerLast;
    private javax.swing.JTextField txtmzIdentOwnerOrg;
    private org.rappsilber.fdr.gui.components.WriteToDB writeDB;
    // End of variables declaration//GEN-END:variables

    /**
     * @return the fdrSettings
     */
    public FDRSettingsPanel getFdrSettings() {
        return fdrSettings;
    }

    /**
     * @param fdrSettings the fdrSettings to set
     */
    public void setFdrSettings(FDRSettingsPanel fdrSettings) {
        this.fdrSettings = fdrSettings;
    }

    /**
     * @return the csvout
     */
    public XiFDRUtils.FDRCSVOUT getCsvOut() {
        return csvout;
    }

    /**
     * @param csvout the csvout to set
     */
    public void setCsvOut(XiFDRUtils.FDRCSVOUT csvout) {
        this.csvout = csvout;
    }

    /**
     * @return the fdrSettingsComplete
     */
    public org.rappsilber.fdr.gui.components.FDRSettingsComplete getFdrSettingsComplete() {
        return fdrSettingsComplete;
    }

    /**
     * @return the fdrSettingsComplete
     */
    public org.rappsilber.fdr.gui.components.FDRSettingsMedium getFdrSettingsMedium() {
        return fdrSettingsMedium;
    }

}
