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
package org.rappsilber.fdr.gui.components;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.HeadlessException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.EventObject;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Locale;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.DefaultCellEditor;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JTable;
import javax.swing.UIManager;
import javax.swing.table.TableColumn;
import javax.swing.table.TableModel;
import org.rappsilber.data.csv.CsvParser;
import org.rappsilber.gui.components.JoinedThreadedTextOuput;
import org.rappsilber.data.csv.ColumnAlternatives;
import org.rappsilber.data.csv.condition.CsvCondition;
import org.rappsilber.data.csv.gui.filter.ConditionList;
import org.rappsilber.fdr.CSVinFDR;
import org.rappsilber.fdr.OfflineFDR;
import org.rappsilber.utils.Version;

/**
 *
 * @author lfischer
 */
public class CSVSelection extends javax.swing.JPanel implements Iterable<CsvParser>{

    private ArrayList<java.awt.event.ActionListener> m_actionlisteners = new ArrayList<ActionListener>();

    private ArrayList<java.awt.event.ActionListener> m_addlisteners = new ArrayList<ActionListener>();
    
    public String missingColumn = "!! !MISSING! !!";
    public String optionalColumn = "  OPTIONAL   ";
    public String[] csvColumns = new String[]{missingColumn};
    public String[] csvColumnsOptional = new String[]{optionalColumn};


    private JoinedThreadedTextOuput m_status = null;

    public void flagModifications(boolean flag) {
        this.ckCSVMarkModifications.setSelected(flag);
    }
    

    
    private class NeededOptionalComboBoxCellEditor extends DefaultCellEditor {

        DefaultCellEditor.EditorDelegate neededDelegate;
        DefaultCellEditor.EditorDelegate optionalDelegate;
        JComponent neededEditorComponent;
        JComponent optionalEditorComponent;

        public NeededOptionalComboBoxCellEditor(final JComboBox needed, final JComboBox optional) {
            super(needed);
            this.neededDelegate = delegate;
            this.neededEditorComponent = needed;

            optionalDelegate = new DefaultCellEditor.EditorDelegate() {
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
                        optional.actionPerformed(new ActionEvent(this, 0, ""));
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
    
    public OfflineFDR.Normalisation doNormalize() {
        return cbNormalize.getItemAt(cbNormalize.getSelectedIndex());
    }
    /**
     * Creates new form CSVSelection
     */
    public CSVSelection() {
        initComponents();
        fbCsvIn.setMultipleFiles(true);
        fbCsvIn.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                testInputFile();
                ckFilter.setEnabled(true);
            }
        });       
        
        fbCsvIn.setLocalPropertyKey("xiFDR_CSV_IN_LAST_FOLDER");

        ckCSVHasHeader.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                testInputFile();
            }
        });
                
        
        cbCSVHeaders.setModel(new DefaultComboBoxModel(csvColumns));
        cbCSVHeaderOptional.setModel(new DefaultComboBoxModel(csvColumnsOptional));

        //spAdditional.setVisible(false);
        
        tblCSVColumns.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {
                {"run",  Boolean.TRUE, "raw file name that the spectrum derived from", ""},
                {"scan",  Boolean.TRUE, "scan number within that run", null},
                {"psmid",  Boolean.TRUE, "a unique ID for the PSM - if not given will be defined based on run and scan", null},
                {"peak list file",  Boolean.TRUE, "the name of the actual peak file containing the spectrum (needed for mzML export)", null},
                {"peak list index",  Boolean.TRUE, "an index for the spectrum within the file(needed for mzML export)", null},
                {"rank",  Boolean.TRUE, "The rank of the PSM (e.g. 1=best match for a spectrum ;2 second best ...)", null},
                {"peptide1", null, "sequence of the first peptide", null},
                {"peptide2", null, "sequence of the second peptide", null},
                {"peptide length 1",  Boolean.TRUE, "length (in amino acids) of the first peptide", null},
                {"peptide length 2",  Boolean.TRUE, "length (in amino acids) of the second peptide", null},
                {"peptide link 1", null, "which residue of the first peptide does the cross-linker attach to", null},
                {"peptide link 2", null, "which residue of the second peptide does the cross-linker attach to", null},
                {"is decoy 1", null, "is the first peptide from the decoy database", null},
                {"is decoy 2", null, "is the second peptide from the decoy database", null},
                {"precursor charge", null, "charge state of the precursor ion", null},
                {"score",  Boolean.TRUE, "score of the spectrum match", null},
                {"score ratio",  Boolean.TRUE, "if there is a joined score given for the match how to separate the score for cases where each peptide individually has to be considered (only affects protein fdr)", null},
                {"peptide1 score",  Boolean.TRUE, "a score for the first peptide", null},
                {"peptide2 score",  Boolean.TRUE, "a score for the second peptide", null},
                {"accession1", null, "protein accession number(s) for the source of the first peptide", null},
                {"accession2", null, "protein accession number(s) for the source of the second peptide", null},
                {"name1",  Boolean.TRUE, "name of the first protein", null},
                {"name2",  Boolean.TRUE, "name of the second protein", null},
                {"description1",  Boolean.TRUE, "description of the first protein", null},
                {"description2",  Boolean.TRUE, "description of the second protein", null},
                {"peptide position 1", null, "position(s) of the first peptide in the protein(s)", null},
                {"peptide position 2", null, "position(s) of the second peptide in the protein(s)", null},
                {"crosslinker",  Boolean.TRUE, "name of the cross-linker involved in this PSM", null},
                {"crossLinkerModMass",  Boolean.TRUE, "mass difference between the sum of the non-cross-linked peptides and teh cross-linked peptides", null},
                {"experimental mz",  Boolean.TRUE, "experimental precursor M/Z", null},
                {"calculated mass",  Boolean.TRUE, "calculated mass of the precursor", null},
                {"negative grouping",  Boolean.TRUE, "if some matches have an inherently higher chance to be false positive then they can be flaged here", null},
                {"positive grouping",  Boolean.TRUE, "if some matches have an inherently lower chance to be false positive then they can be flaged here", null},
                {"info",  Boolean.TRUE, "arbitrary info field to be forwarded to the results table", null},
                {"delta score",  Boolean.TRUE, "the delta score of the match", null},
                {"peptide coverage1",  Boolean.TRUE, "how well is peptide explained", null},
                {"peptide coverage2",  Boolean.TRUE, "how well is peptide 2 explained", null},
                {"minimum peptide coverage",  Boolean.TRUE, "how well is the wors explained peptide explained", null},
                {"peptide1 fragments",  Boolean.TRUE, "How many fragments where assigned to peptide 1", null},
                {"peptide2 fragments",  Boolean.TRUE, "How many fragments where assigned to peptide 2", null},
                {"peptides with stubs",  Boolean.TRUE, "For how many peptides was an cross-linker stub identified", null},
                {"peptides with doublets",  Boolean.TRUE, "For how many peptides was an doublet of cross-linker stubs identified", null}
            },
            new String [] {
                "Column", "Optional", "Description", "Name in CSV"
            }
        ) {
            Class[] types = new Class [] {
                java.lang.String.class, java.lang.Boolean.class, java.lang.String.class, java.lang.String.class
            };
            boolean[] canEdit = new boolean [] {
                false, false, false, true
            };

            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }

            public boolean isCellEditable(int rowIndex, int columnIndex) {
                return canEdit [columnIndex];
            }
        });        
        TableColumn columnNamesColumn = tblCSVColumns.getColumnModel().getColumn(3);
        columnNamesColumn.setCellEditor(new NeededOptionalComboBoxCellEditor(cbCSVHeaders, cbCSVHeaderOptional));
        resetColumnMappings();

        
    }

    public void resetColumnMappings() {
        TableModel tm = tblCSVColumns.getModel();
        for (int r = 0; r < tblCSVColumns.getRowCount(); r++) {
            if (Boolean.TRUE.equals(tm.getValueAt(r, 1))) {
                tblCSVColumns.getModel().setValueAt(optionalColumn, r, 3);
            } else {
                tblCSVColumns.getModel().setValueAt(missingColumn, r, 3);
            }
        }
    }
    
    public void guessXiVersion() {
        Pattern p = Pattern.compile(".*(?:xi|xisearch)[\\-_]?([0-9]\\.[0-9]+.*).(?:csv|tsv|txt)(?:.gz)?$", Pattern.CASE_INSENSITIVE);
        File[] files  =  fbCsvIn.getFiles();
        Matcher m = p.matcher(files[0].getName());
        Version v = null;
        if (m.matches()) {
            try {
                v = new Version(m.group(1));
                if (txtXiVersion.getText().trim().length()>0 && new Version(txtXiVersion.getText().trim()).compareTo(v) != 0) {
                    txtXiVersion.setText("");
                    return;
                }
            } catch (Exception e) {
                txtXiVersion.setText("");
                return;
            }
        }
        if (files.length>1 && v != null) {
            for (int f=1; f<files.length; f++) {
                m = p.matcher(files[f].getName());
                if (m.matches()) {
                    try {
                        Version v2 = new Version(m.group(1));
                        if (v2.compareTo(v) != 0) {
                            v = null;
                            break;
                        }
                    } catch (Exception e) {
                        v = null;
                        break;
                    }
                } else {
                    v = null;
                    break;
                }
                
            }
        }
        if (v!=null)
            txtXiVersion.setText(v.toString());
    }
    
    public Version getXiVersion(){
        if (txtXiVersion.getText().trim().length() > 0)
            return new Version(txtXiVersion.getText());
        return null;
    }
    
    public void setXiVersion(String xiVersion) {
        if (xiVersion != null && xiVersion.trim().length()>0) {
            txtXiVersion.setText(xiVersion);
        }
    }

    public void testInputFile() {
        File f = fbCsvIn.getFile();
        if (f != null && f.canRead()) {
            CsvParser csv;
            try {
                csv = CsvParser.guessCsv(f, 50);
                String delimiter = csv.getDelimiter();
                if (delimiter.contentEquals(",")) {
                    delimiter = "Comma(,)";
                } else if (delimiter.contentEquals(";")) {
                    delimiter = "Semicolon(;)";
                } else if (delimiter.contentEquals("\t")) {
                    delimiter = "Tab";
                } else if (delimiter.contentEquals(" ")) {
                    delimiter = "Space";
                }
                cbCSVDElimiter.setSelectedItem(delimiter);
                cbCSVQuote.setSelectedItem(csv.getQuote());
                csv.close();
                readColumns();
                guessXiVersion();

            } catch (FileNotFoundException ex) {
                Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
                setStatus(ex.toString());
            } catch (IOException ex) {
                Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
                setStatus(ex.toString());
            }
        }

    }

    protected void readColumns() throws IOException {
        File f = fbCsvIn.getFile();
        if (f != null && f.canRead()) {
            String delimiter = cbCSVDElimiter.getSelectedItem().toString();
            String quote = cbCSVQuote.getSelectedItem().toString();
            if (delimiter.contentEquals("Comma(,)")) {
                delimiter = ",";
            } else if (delimiter.contentEquals("Tab")) {
                delimiter = "\t";
            } else if (delimiter.contentEquals("Semicolon(;)")) {
                delimiter = ";";
            } else if (delimiter.contentEquals("Space")) {
                delimiter = " ";
            }
            CsvParser csv = new CsvParser(delimiter.charAt(0), quote.charAt(0));
            ColumnAlternatives.setupAlternatives(csv,CSVinFDR.DEFAULT_COLUMN_MAPPING);
            csv.openFile(f);
            csv.next();
            csvColumns = new String[csv.getMaxColumns() + 1];
            csvColumnsOptional = new String[csv.getMaxColumns() + 1];
            csvColumns[0] = missingColumn;
            csvColumnsOptional[0] = optionalColumn;

            
            if (ckCSVHasHeader.isSelected()) {
                csv.setCurrentLineAsHeader();
                HashMap<String,String[]> excludeMappings = new HashMap<>(2);
                if (ckSmartMatch.isSelected()) {
                    excludeMappings.put("psmid",new String[]{"spectra to matched","matchrank","spectrum quality score"});
                    ColumnAlternatives.levenshteinMatchHeadersAlternatives(csv, excludeMappings,0.7);
                }
                for (int c = 0; c < csv.getMaxColumns(); c++) {
                    csvColumns[c + 1] = csv.getHeader(c);
                    csvColumnsOptional[c + 1] = csv.getHeader(c);
                }
                TableModel tm = tblCSVColumns.getModel();
                for (int r = 0; r < tm.getRowCount(); r++) {
                    Integer rc = csv.getColumn(tm.getValueAt(r, 0).toString());
                    if (rc != null) {
                        tm.setValueAt(csv.getHeader(rc), r, 3);
                    } else {
                        if (tm.getValueAt(r, 1)!= null && ((Boolean)tm.getValueAt(r, 1))) {
                            tm.setValueAt(optionalColumn, r, 3);
                        } else
                            tm.setValueAt(missingColumn, r, 3);
                    }
                }

            } else {
                for (int c = 0; c < csv.getMaxColumns(); c++) {
                    csvColumns[c + 1] = Integer.toString(c);
                }
            }

            cbCSVHeaders.setModel(new DefaultComboBoxModel(csvColumns));
            cbCSVHeaderOptional.setModel(new DefaultComboBoxModel(csvColumnsOptional));
            if (filter.setCsvParser(csv) >0) {
                btnFilter.setIcon(UIManager.getIcon("OptionPane.errorIcon"));
                btnFilter.setToolTipText("Filter no longer match to file");
            }
            csv.close();
        }
    }
 
    public Locale getLocale() {
        return localPicker1.getSelectLocale();
    }
    private void setStatus(final String status) {
        if (m_status != null)
            m_status.write(status);
    }
    
    public void setStatusInterface(final JoinedThreadedTextOuput status ) {
        m_status = status;
    }
    
    public File getFile() {
        return fbCsvIn.getFile();
    }

    public File[] getFiles() {
        return fbCsvIn.getFiles();
    }
    
    public String getForwardPattern() {
        if (!txtForwardColumns.getText().trim().isEmpty()) {
            return txtForwardColumns.getText().trim();
        } else {
            return null;
        }
    }

    public void doActionPerformed() {
        ActionEvent ae = new ActionEvent(this, ActionEvent.ACTION_PERFORMED, "file selected",Calendar.getInstance().getTimeInMillis(), 0);
        for (java.awt.event.ActionListener al : m_actionlisteners) 
            al.actionPerformed(null);
    }

    public void addActionListener(java.awt.event.ActionListener al) {
        this.m_actionlisteners.add(al);
    }
    

    public void removeActionListener(java.awt.event.ActionListener al) {
        this.m_actionlisteners.remove(al);
    }


    protected void doAddPerformed() {
        ActionEvent ae = new ActionEvent(this, ActionEvent.ACTION_PERFORMED, "file selected",Calendar.getInstance().getTimeInMillis(), 0);
        for (java.awt.event.ActionListener al : m_addlisteners) 
            al.actionPerformed(null);
    }

    
    public void addAddListener(java.awt.event.ActionListener al) {
        this.m_addlisteners.add(al);
    }
    

    public void removeAddListener(java.awt.event.ActionListener al) {
        this.m_addlisteners.remove(al);
    }
    
    
    public String getDelimiter() {
        String delim = cbCSVDElimiter.getSelectedItem().toString();
        if (delim.startsWith("Comma")) {
            delim = ",";
        } else if (delim.startsWith("Semicolon")) {
            delim = ";";
        } else if (delim.startsWith("Tab")) {
            delim = "\t";
        } else if (delim.startsWith("Space")) {
            delim = " ";
        }
        return delim;
    }

    public char getQuote() {
        return cbCSVQuote.getSelectedItem().toString().charAt(0);
    }
    
    public boolean higherIsBetter() {
        return rbCSVHighBetter.isSelected();
    }
    
    public boolean hasHeader() {
        return ckCSVHasHeader.isSelected();
    }
    
    
    public CsvParser getCsv() throws IOException {
        final TableModel tm = tblCSVColumns.getModel();
        HashMap<String,String> file2Expected = new HashMap();
        CsvParser csv = new CsvParser();
        for (int r = 0; r < tm.getRowCount(); r++) {
            if (tm.getValueAt(r, 3) != missingColumn) {
                file2Expected.put(tm.getValueAt(r, 3).toString(), tm.getValueAt(r, 0).toString());
                csv.setAlternative(tm.getValueAt(r, 3).toString(), tm.getValueAt(r, 0).toString());
            }
        }
        csv.renameHeader(file2Expected);
        String delimiter = getDelimiter();
        csv.setDelimiter(delimiter.charAt(0));
        csv.setQuote(getQuote());
        csv.openFile(getFile(), hasHeader());
        csv.setLocale(localPicker1.getSelectLocale());
        return csv;
    }
    
    public Iterator<CsvParser> iterator() {
        final File[] inputs = getFiles();
        final TableModel tm = tblCSVColumns.getModel();
        String preDelimiter = getDelimiter();
        final char delimiter = preDelimiter.charAt(0);
        final char quote = getQuote();
        final boolean hasHeader =hasHeader();
        return new Iterator<CsvParser>() {
            int n=0;
            
            @Override
            public boolean hasNext() {
                return n<inputs.length;
            }

            @Override
            public CsvParser next()  {
                File nf = inputs[n++];
                CsvParser csv = new CsvParser();
                for (int r = 0; r < tm.getRowCount(); r++) {
                    if (tm.getValueAt(r, 3) != missingColumn && tm.getValueAt(r, 3) != optionalColumn) {
                        csv.setAlternative(tm.getValueAt(r, 3).toString(), tm.getValueAt(r, 0).toString());
                    }
                }
                csv.setDelimiter(delimiter);
                csv.setQuote(quote);
                csv.setLocale(localPicker1.getSelectLocale());
                try {                
                    csv.openFile(nf, hasHeader);
                } catch (IOException ex) {
                    Logger.getLogger(CSVSelection.class.getName()).log(Level.SEVERE, null, ex);
                    throw new RuntimeException(ex);
                }
                return csv;
            }
        };
        
    }

    
    JFrame windowFilter;
    
    public void showFilter() {
        if (windowFilter == null){
            windowFilter = new JFrame();
            windowFilter.setName("XiFDR - CSV Filter");
            windowFilter.setLayout(new BorderLayout());
            ConditionList csv = filter;
            windowFilter.setPreferredSize(csv.getPreferredSize());
            windowFilter.add(csv, BorderLayout.CENTER);
            JButton btnok = new JButton("OK");
            btnok.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    windowFilter.setVisible(false);
                }
            });
            windowFilter.add(csv, BorderLayout.CENTER);
            windowFilter.add(btnok, BorderLayout.PAGE_END);
            windowFilter.pack();
            Dimension s = windowFilter.getSize();
            s.height=s.height*2;
            windowFilter.setSize(s);
            windowFilter.setVisible(true);
            
        } else {
            windowFilter.setVisible(true);
        }
    }
    
    public CsvCondition getFilter() throws ParseException{
        if (ckFilter.isSelected()) {
            return filter.getCondition();
        }
        return null;
    }
    
    public void setEnabled(boolean enabled) {
        super.setEnabled(enabled);
        btnReadCsv.setEnabled(enabled);
    }

    public void setEnableAdd(boolean enabled) {
        btnAddCSV.setEnabled(enabled);
    }

    
    
    public void addAdditionalOption(Component c) {
        spAdditional.setVisible(true);
        pAdditional.add(c);
                
    }
    
    public void setInputFile(String path) {
        fbCsvIn.setFile(path);
        testInputFile();
        ckFilter.setEnabled(true);
        
    }

    public void setXiConfig(String path) {
        fbConfigIn.setFile(path);
    }
    
    public void setFasta(String path) {
        fbFastaIn.setFile(path);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        cbCSVHeaders = new javax.swing.JComboBox();
        cbCSVHeaderOptional = new javax.swing.JComboBox();
        bgScoreDirectionCSV = new javax.swing.ButtonGroup();
        filter = new org.rappsilber.data.csv.gui.filter.ConditionList();
        pAdditional = new javax.swing.JPanel();
        jLabel26 = new javax.swing.JLabel();
        cbCSVQuote = new javax.swing.JComboBox();
        cbCSVDElimiter = new javax.swing.JComboBox();
        rbCSVHighBetter = new javax.swing.JRadioButton();
        jScrollPane2 = new javax.swing.JScrollPane();
        tblCSVColumns = new javax.swing.JTable();
        ckCSVHasHeader = new javax.swing.JCheckBox();
        fbCsvIn = new org.rappsilber.gui.components.FileBrowser();
        rbCSVLowBetter = new javax.swing.JRadioButton();
        jLabel11 = new javax.swing.JLabel();
        btnReadCsv = new javax.swing.JButton();
        jLabel27 = new javax.swing.JLabel();
        btnAddCSV = new javax.swing.JButton();
        ckSmartMatch = new javax.swing.JCheckBox();
        spAdditional = new javax.swing.JScrollPane();
        pXiConfig = new javax.swing.JPanel();
        fbConfigIn = new org.rappsilber.gui.components.FileBrowser();
        jLabel12 = new javax.swing.JLabel();
        jLabel13 = new javax.swing.JLabel();
        fbFastaIn = new org.rappsilber.gui.components.FileBrowser();
        ckCSVMarkModifications = new javax.swing.JCheckBox();
        txtForwardColumns = new javax.swing.JTextField();
        jLabel1 = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        cbNormalize = new javax.swing.JComboBox<>();
        txtXiVersion = new javax.swing.JTextField();
        jLabel4 = new javax.swing.JLabel();
        ckFilter = new javax.swing.JCheckBox();
        btnFilter = new javax.swing.JButton();
        localPicker1 = new org.rappsilber.gui.components.LocalPicker();
        jLabel28 = new javax.swing.JLabel();
        txtDecoyPrefix = new javax.swing.JTextField();
        jLabel3 = new javax.swing.JLabel();

        cbCSVHeaders.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Item 1", "Item 2", "Item 3", "Item 4" }));

        cbCSVHeaderOptional.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Item 1", "Item 2", "Item 3", "Item 4" }));

        pAdditional.setPreferredSize(new java.awt.Dimension(300, 90));
        pAdditional.setLayout(new javax.swing.BoxLayout(pAdditional, javax.swing.BoxLayout.Y_AXIS));

        jLabel26.setText("Delimiter");

        cbCSVQuote.setEditable(true);
        cbCSVQuote.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "\"", "'" }));
        cbCSVQuote.setPreferredSize(new java.awt.Dimension(100, 24));
        cbCSVQuote.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cbCSVQuoteActionPerformed(evt);
            }
        });

        cbCSVDElimiter.setEditable(true);
        cbCSVDElimiter.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Comma(,)", "Semicolon(;)", "Tab", "Space", "|", " " }));
        cbCSVDElimiter.setPreferredSize(new java.awt.Dimension(100, 24));
        cbCSVDElimiter.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cbCSVDElimiterActionPerformed(evt);
            }
        });

        bgScoreDirectionCSV.add(rbCSVHighBetter);
        rbCSVHighBetter.setSelected(true);
        rbCSVHighBetter.setText("High Score better");

        tblCSVColumns.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {
                {"run",  new Boolean(true), "raw file name that the spectrum derived from", ""},
                {"scan",  new Boolean(true), "scan number within that run", null},
                {"psmid",  new Boolean(true), "a unique ID for the PSM - if not given will be defined based on run and scan", null},
                {"peak list file",  new Boolean(true), "the name of the actual peak file containing the spectrum (needed for mzML export)", null},
                {"peak list index",  new Boolean(true), "an index for the spectrum within the file(needed for mzML export)", null},
                {"rank",  new Boolean(true), "The rank of the PSM (e.g. 1=best match for a spectrum ;2 second best ...)", null},
                {"peptide1", null, "sequence of the first peptide", null},
                {"peptide2", null, "sequence of the second peptide", null},
                {"peptide length 1",  new Boolean(true), "length (in amino acids) of the first peptide", null},
                {"peptide length 2",  new Boolean(true), "length (in amino acids) of the second peptide", null},
                {"peptide link 1", null, "which residue of the first peptide does the cross-linker attach to", null},
                {"peptide link 2", null, "which residue of the second peptide does the cross-linker attach to", null},
                {"is decoy 1", null, "is the first peptide from the decoy database", null},
                {"is decoy 2", null, "is the second peptide from the decoy database", null},
                {"precursor charge", null, "charge state of the precursor ion", null},
                {"score",  new Boolean(true), "score of the spectrum match", null},
                {"score ratio",  new Boolean(true), "if there is a joined score given for the match how to separate the score for cases where each peptide individually has to be considered (only affects protein fdr)", null},
                {"peptide1 score",  new Boolean(true), "a score for the first peptide", null},
                {"peptide2 score",  new Boolean(true), "a score for the second peptide", null},
                {"accession1", null, "protein accession number(s) for the source of the first peptide", null},
                {"accession2", null, "protein accession number(s) for the source of the second peptide", null},
                {"description1",  new Boolean(true), "description of the first protein", null},
                {"description2",  new Boolean(true), "description of the second protein", null},
                {"peptide position 1", null, "position(s) of the first peptide in the protein(s)", null},
                {"peptide position 2", null, "position(s) of the second peptide in the protein(s)", null},
                {"crosslinker",  new Boolean(true), "name of the cross-linker involved in this PSM", null},
                {"crossLinkerModMass",  new Boolean(true), "mass difference between the sum of the non-cross-linked peptides and teh cross-linked peptides", null},
                {"experimental mz",  new Boolean(true), "experimental precursor M/Z", null},
                {"calculated mass",  new Boolean(true), "calculated mass of the precursor", null},
                {"negative grouping",  new Boolean(true), "if some matches have an inherently higher chance to be false positive then they can be flaged here", null},
                {"positive grouping",  new Boolean(true), "if some matches have an inherently lower chance to be false positive then they can be flaged here", null},
                {"info",  new Boolean(true), "arbitrary info field to be forwarded to the results table", null},
                {"delta score",  new Boolean(true), "the delta score of the match", null},
                {"peptide coverage1",  new Boolean(true), "how well is peptide explained", null},
                {"peptide coverage2",  new Boolean(true), "how well is peptide 2 explained", null},
                {"minimum peptide coverage",  new Boolean(true), "how well is the wors explained peptide explained", null},
                {"peptide1 fragments",  new Boolean(true), "How many fragments where assigned to peptide 1", null},
                {"peptide2 fragments",  new Boolean(true), "How many fragments where assigned to peptide 2", null}
            },
            new String [] {
                "Column", "Optional", "Description", "Name in CSV"
            }
        ) {
            Class[] types = new Class [] {
                java.lang.String.class, java.lang.Boolean.class, java.lang.String.class, java.lang.String.class
            };
            boolean[] canEdit = new boolean [] {
                false, false, false, true
            };

            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }

            public boolean isCellEditable(int rowIndex, int columnIndex) {
                return canEdit [columnIndex];
            }
        });
        tblCSVColumns.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_LAST_COLUMN);
        jScrollPane2.setViewportView(tblCSVColumns);

        ckCSVHasHeader.setSelected(true);
        ckCSVHasHeader.setText("hasHeader");
        ckCSVHasHeader.setToolTipText("first row of input is considered header");
        ckCSVHasHeader.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckCSVHasHeaderActionPerformed(evt);
            }
        });

        fbCsvIn.setDescription("Text Files");
        fbCsvIn.setExtensions(new String[] {"csv", "txt", "tsv"});

        bgScoreDirectionCSV.add(rbCSVLowBetter);
        rbCSVLowBetter.setText("Lower score better");

        jLabel11.setText("CSV-File");

        btnReadCsv.setText("Read");
        btnReadCsv.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnReadCsvActionPerformed(evt);
            }
        });

        jLabel27.setText("Quote");

        btnAddCSV.setText("Add");
        btnAddCSV.setEnabled(false);
        btnAddCSV.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnAddCSVActionPerformed(evt);
            }
        });

        ckSmartMatch.setText("tolerant column matching");
        ckSmartMatch.setToolTipText("Needs Manual control:tries to select the best header by string distance");
        ckSmartMatch.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckSmartMatchActionPerformed(evt);
            }
        });

        pXiConfig.setPreferredSize(new java.awt.Dimension(300, 107));

        fbConfigIn.setDescription("Xi-Config");
        fbConfigIn.setExtensions(new String[] {"config", "conf", "json"});

        jLabel12.setText("xi Config");

        jLabel13.setText("FASTA");

        fbFastaIn.setDescription("Text Files");
        fbFastaIn.setExtensions(new String[] {"fasta", "txt"});

        ckCSVMarkModifications.setText("Flag Modifications");

        txtForwardColumns.setToolTipText("A Regular Expression defining what addtional columns should be forwarded to the PSM-file");

        jLabel1.setText("Forward Columns");

        jLabel2.setText("Normalize");

        cbNormalize.setModel(new DefaultComboBoxModel<OfflineFDR.Normalisation>(OfflineFDR.Normalisation.values()));
        cbNormalize.setToolTipText("How to normalize the input data");

        txtXiVersion.setToolTipText("xi version to be written into the mzIdentML file");

        jLabel4.setText("xi Version");

        javax.swing.GroupLayout pXiConfigLayout = new javax.swing.GroupLayout(pXiConfig);
        pXiConfig.setLayout(pXiConfigLayout);
        pXiConfigLayout.setHorizontalGroup(
            pXiConfigLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pXiConfigLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(pXiConfigLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel12)
                    .addComponent(jLabel13))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(pXiConfigLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(pXiConfigLayout.createSequentialGroup()
                        .addComponent(ckCSVMarkModifications, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jLabel2)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(cbNormalize, 0, 209, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jLabel1)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(txtForwardColumns, javax.swing.GroupLayout.DEFAULT_SIZE, 128, Short.MAX_VALUE))
                    .addComponent(fbFastaIn, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(pXiConfigLayout.createSequentialGroup()
                        .addComponent(fbConfigIn, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jLabel4)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(txtXiVersion, javax.swing.GroupLayout.PREFERRED_SIZE, 64, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap())
        );
        pXiConfigLayout.setVerticalGroup(
            pXiConfigLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pXiConfigLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(pXiConfigLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.CENTER)
                    .addComponent(fbConfigIn, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel12)
                    .addComponent(txtXiVersion, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel4))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(pXiConfigLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.CENTER)
                    .addComponent(fbFastaIn, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel13))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(pXiConfigLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ckCSVMarkModifications)
                    .addComponent(txtForwardColumns, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel1)
                    .addComponent(cbNormalize, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel2))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        spAdditional.setViewportView(pXiConfig);

        ckFilter.setEnabled(false);
        ckFilter.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckFilterActionPerformed(evt);
            }
        });

        btnFilter.setText("Filter");
        btnFilter.setEnabled(false);
        btnFilter.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnFilterActionPerformed(evt);
            }
        });

        localPicker1.setDefaultLocal(java.util.Locale.ENGLISH);
        localPicker1.setMinimumSize(new java.awt.Dimension(80, 27));
        localPicker1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                localPicker1ActionPerformed(evt);
            }
        });

        jLabel28.setText("Language");

        txtDecoyPrefix.setText("auto");
        txtDecoyPrefix.setToolTipText("<html><body><p>In order to accurately match decoy proteins to their corresponding target proteins, <br/>it's important to remove any prefix included in the accession numbers for the decoys.<br/> \nFailing to do so can prevent proper assignment of target-decoy matches as self-links, resulting in false FDR estimates.\u0000<br/></p><p>\nBy default REV_, REVERSE_, DECOY_, DECOY:, RANDOM_, and RAN_ are recognised.</p></body></html>");

        jLabel3.setText("decoy prefix:");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(13, 13, 13)
                        .addComponent(jLabel11))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(jLabel26))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(ckFilter)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(btnFilter)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(ckCSVHasHeader)
                            .addComponent(cbCSVDElimiter, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(ckSmartMatch)
                                .addGap(18, 18, 18)
                                .addComponent(jLabel3)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(txtDecoyPrefix, javax.swing.GroupLayout.PREFERRED_SIZE, 59, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGap(0, 0, Short.MAX_VALUE))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jLabel27)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(cbCSVQuote, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jLabel28)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(localPicker1, javax.swing.GroupLayout.PREFERRED_SIZE, 1, Short.MAX_VALUE)))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(rbCSVHighBetter)
                            .addComponent(rbCSVLowBetter))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(btnReadCsv, javax.swing.GroupLayout.PREFERRED_SIZE, 76, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(btnAddCSV)))
                    .addComponent(fbCsvIn, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jScrollPane2)
                    .addComponent(spAdditional, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE)))
        );

        layout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {btnAddCSV, btnReadCsv});

        layout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {cbCSVDElimiter, cbCSVQuote});

        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(18, 18, 18)
                        .addComponent(jLabel11))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(fbCsvIn, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(btnAddCSV)
                    .addComponent(rbCSVHighBetter)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(jLabel26)
                        .addComponent(cbCSVDElimiter, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jLabel27)
                        .addComponent(cbCSVQuote, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(localPicker1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jLabel28)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.CENTER)
                    .addComponent(btnReadCsv)
                    .addComponent(ckCSVHasHeader)
                    .addComponent(btnFilter)
                    .addComponent(ckFilter)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(rbCSVLowBetter)
                        .addComponent(ckSmartMatch)
                        .addComponent(txtDecoyPrefix, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jLabel3)))
                .addGap(18, 18, 18)
                .addComponent(spAdditional, javax.swing.GroupLayout.DEFAULT_SIZE, 105, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jScrollPane2, javax.swing.GroupLayout.DEFAULT_SIZE, 201, Short.MAX_VALUE))
        );
    }// </editor-fold>//GEN-END:initComponents

    private void cbCSVQuoteActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cbCSVQuoteActionPerformed
        String quote = cbCSVQuote.getSelectedItem().toString();
        if (quote.length() > 1) {
            JOptionPane.showMessageDialog(this, "Quote \"" + quote + "\" is not supported. \n\nThe quote character has to be a single character", "Unsupoorted quote", JOptionPane.ERROR_MESSAGE);
            cbCSVQuote.requestFocusInWindow();
        } else if (quote.isEmpty()) {
            JOptionPane.showMessageDialog(this, "Sorry we need a quote character", "Unsupoorted quote", JOptionPane.ERROR_MESSAGE);
            cbCSVQuote.setSelectedItem("\"");
            cbCSVQuote.requestFocusInWindow();
        }
        try {
            readColumns();
        } catch (IOException ex) {
            Logger.getLogger(CSVSelection.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }//GEN-LAST:event_cbCSVQuoteActionPerformed

    private void cbCSVDElimiterActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cbCSVDElimiterActionPerformed

        String delimiter = cbCSVDElimiter.getSelectedItem().toString();
        if (delimiter.length() > 1) {
            if (!delimiter.matches("(Tab|Comma|Comma\\(,\\)|Space|Semicolon|Semicolon\\(;\\))")) {
                JOptionPane.showMessageDialog(this, "Delimiter \"" + delimiter + "\" is not supported. \n\nOnly single character delimiter are supported", "Unsupoorted delimiter", JOptionPane.ERROR_MESSAGE);
                cbCSVDElimiter.requestFocusInWindow();
            }
        } else if (delimiter.isEmpty()) {
            JOptionPane.showMessageDialog(this, "Sorry we need a delimiter", "Unsupoorted delimiter", JOptionPane.ERROR_MESSAGE);
            cbCSVDElimiter.setSelectedItem("Comma(,)");
            cbCSVDElimiter.requestFocusInWindow();
        } else {
            if (delimiter.contentEquals(",")) {
                cbCSVDElimiter.setSelectedItem("Comma(,)");
            } else if (delimiter.contentEquals(";")) {
                cbCSVDElimiter.setSelectedItem("Semicolon(;)");
            } else if (delimiter.contentEquals(" ")) {
                cbCSVDElimiter.setSelectedItem("Space");
            } else if (delimiter.contentEquals("\t")) {
                cbCSVDElimiter.setSelectedItem("Tab");
            }
        }
        try {
            readColumns();
        } catch (IOException ex) {
            Logger.getLogger(CSVSelection.class.getName()).log(Level.SEVERE, null, ex);
        }
    }//GEN-LAST:event_cbCSVDElimiterActionPerformed

    private void btnReadCsvActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnReadCsvActionPerformed

            if (!testInput()) 
                return;
        
        doActionPerformed();
    }//GEN-LAST:event_btnReadCsvActionPerformed

    protected boolean testInput() throws HeadlessException {
        boolean scoreSelected = false;
        boolean peptide1ScoreSelected = false;
        boolean peptide2ScoreSelected = false;
        boolean psmidSelected = false;
        boolean runSelected = false;
        boolean scanSelected = false;
        final TableModel tm = tblCSVColumns.getModel();
        // check, whether we have all the neded columns
        for (int r = 0; r < tm.getRowCount(); r++) {
            if ((!Boolean.TRUE.equals(tm.getValueAt(r, 1))) && tm.getValueAt(r, 3) == missingColumn) {
                JOptionPane.showMessageDialog(this, "No column for " + tm.getValueAt(r, 0) + " selected", "Missing Column", JOptionPane.WARNING_MESSAGE);
                return false;
            }
            if (tm.getValueAt(r, 0).toString().contentEquals("Peptide1 Score") && tm.getValueAt(r, 3) != optionalColumn) {
                peptide1ScoreSelected = true;
            }
            if (tm.getValueAt(r, 0).toString().contentEquals("Peptide2 Score") && tm.getValueAt(r, 3) != optionalColumn) {
                peptide2ScoreSelected = true;
            }
            if (tm.getValueAt(r, 0).toString().contentEquals("score") && tm.getValueAt(r, 3) != optionalColumn) {
                scoreSelected = true;
            }
            if (tm.getValueAt(r, 0).toString().contentEquals("run") && tm.getValueAt(r, 3) != optionalColumn) {
                runSelected = true;
            }
            if (tm.getValueAt(r, 0).toString().contentEquals("scan") && tm.getValueAt(r, 3) != optionalColumn) {
                scanSelected = true;
            }
            if (tm.getValueAt(r, 0).toString().contentEquals("psmid") && tm.getValueAt(r, 3) != optionalColumn) {
                psmidSelected = true;
            }
        }
        if (!(scoreSelected || (peptide1ScoreSelected && peptide2ScoreSelected))) {
            JOptionPane.showMessageDialog(this, "At least a PSM-level-score (Score-column) \n"
                    + "or two peptide scores for a PSM is required.", "Missing Column", JOptionPane.WARNING_MESSAGE);
            return false;
        }
        return true;
    }
    
    public String getDecoyPrefix() {
        String prefix = txtDecoyPrefix.getText();
        if (prefix != null && prefix.trim().toLowerCase().contentEquals("auto")) {
            prefix = null;
        }
        return prefix;
    }

    private void btnAddCSVActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnAddCSVActionPerformed
            if (!testInput()) return;

            doAddPerformed();
    }//GEN-LAST:event_btnAddCSVActionPerformed

    private void ckSmartMatchActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckSmartMatchActionPerformed
        resetColumnMappings();
        
        testInputFile();
    }//GEN-LAST:event_ckSmartMatchActionPerformed

    private void ckCSVHasHeaderActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckCSVHasHeaderActionPerformed
        ckSmartMatch.setEnabled(ckCSVHasHeader.isSelected());
    }//GEN-LAST:event_ckCSVHasHeaderActionPerformed

    private void ckFilterActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckFilterActionPerformed
        btnFilter.setEnabled(ckFilter.isSelected());
        if (ckFilter.isSelected()) {
            showFilter();
        }
    }//GEN-LAST:event_ckFilterActionPerformed

    private void btnFilterActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnFilterActionPerformed
        showFilter();
        btnFilter.setIcon(null);
        btnFilter.setToolTipText("open filter dialog");
        
    }//GEN-LAST:event_btnFilterActionPerformed

    private void localPicker1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_localPicker1ActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_localPicker1ActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.ButtonGroup bgScoreDirectionCSV;
    private javax.swing.JButton btnAddCSV;
    private javax.swing.JButton btnFilter;
    private javax.swing.JButton btnReadCsv;
    private javax.swing.JComboBox cbCSVDElimiter;
    private javax.swing.JComboBox cbCSVHeaderOptional;
    private javax.swing.JComboBox cbCSVHeaders;
    private javax.swing.JComboBox cbCSVQuote;
    private javax.swing.JComboBox<OfflineFDR.Normalisation> cbNormalize;
    private javax.swing.JCheckBox ckCSVHasHeader;
    public javax.swing.JCheckBox ckCSVMarkModifications;
    private javax.swing.JCheckBox ckFilter;
    private javax.swing.JCheckBox ckSmartMatch;
    public org.rappsilber.gui.components.FileBrowser fbConfigIn;
    private org.rappsilber.gui.components.FileBrowser fbCsvIn;
    public org.rappsilber.gui.components.FileBrowser fbFastaIn;
    private org.rappsilber.data.csv.gui.filter.ConditionList filter;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel11;
    private javax.swing.JLabel jLabel12;
    private javax.swing.JLabel jLabel13;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel26;
    private javax.swing.JLabel jLabel27;
    private javax.swing.JLabel jLabel28;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JScrollPane jScrollPane2;
    private org.rappsilber.gui.components.LocalPicker localPicker1;
    private javax.swing.JPanel pAdditional;
    private javax.swing.JPanel pXiConfig;
    private javax.swing.JRadioButton rbCSVHighBetter;
    private javax.swing.JRadioButton rbCSVLowBetter;
    public javax.swing.JScrollPane spAdditional;
    public javax.swing.JTable tblCSVColumns;
    private javax.swing.JTextField txtDecoyPrefix;
    private javax.swing.JTextField txtForwardColumns;
    private javax.swing.JTextField txtXiVersion;
    // End of variables declaration//GEN-END:variables
}
