/*
 * Copyright 2021 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
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
package org.rappsilber.fdr.dataimport;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.json.simple.parser.ParseException;
import org.json.simple.parser.JSONParser;
import org.rappsilber.utils.ms.Composition;
import rappsilber.config.AbstractRunConfig;
import rappsilber.config.ConfigurationParserException;
import rappsilber.config.RunConfig;
import rappsilber.ms.ToleranceUnit;
import rappsilber.ms.crosslinker.AsymetricSingleAminoAcidRestrictedCrossLinker;
import rappsilber.ms.crosslinker.CrossLinker;
import rappsilber.ms.crosslinker.SymetricSingleAminoAcidRestrictedCrossLinker;
import rappsilber.ms.sequence.AminoAcid;
import rappsilber.ms.sequence.AminoModification;
import rappsilber.ms.sequence.digest.AAConstrainedDigestion;


/**
 * Parses the minimal needed information from a xi2config json
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class Xi2Xi1Config extends AbstractRunConfig{
    HashMap<String, Double> default_xl_masses = new HashMap<>();
    HashMap<String, CrossLinker> default_xl_xi1 = new HashMap<>();
    HashMap<String, Xi2Crosslinker> default_xl_xi2 = new HashMap<>();
    public boolean isModX = false;
    {   
        AminoAcid[]  KSTY = new AminoAcid[]{AminoAcid.K, AminoAcid.S, AminoAcid.T, AminoAcid.Y};
        AminoAcid[]  DE = new AminoAcid[]{AminoAcid.D, AminoAcid.E};
        AminoAcid[]  X = new AminoAcid[]{AminoAcid.X};
        SymetricSingleAminoAcidRestrictedCrossLinker BS3 = new SymetricSingleAminoAcidRestrictedCrossLinker("BS3", 138.06807961, 138.06807961, KSTY);
        BS3.setlinksNTerm(true);
        SymetricSingleAminoAcidRestrictedCrossLinker DSSO = new SymetricSingleAminoAcidRestrictedCrossLinker("DSSO", 158.0038, 158.0038, KSTY);
        BS3.setlinksNTerm(true);
        AsymetricSingleAminoAcidRestrictedCrossLinker EDC = new AsymetricSingleAminoAcidRestrictedCrossLinker(
                "EDC", -18.01056027, -18.01056027, KSTY, DE);
        EDC.setlinksNTerm(true);
        EDC.setLinksCTermSecondary(true);
        AsymetricSingleAminoAcidRestrictedCrossLinker SDA = new AsymetricSingleAminoAcidRestrictedCrossLinker(
                "SDA", 82.04186484, 82.04186484, KSTY, X);
        SDA.setlinksNTerm(true);
                
        default_xl_xi1.put("BS3", BS3);
        default_xl_xi1.put("EDC", EDC);
        default_xl_xi1.put("SDA", SDA);
        default_xl_xi1.put("DSSO", DSSO);

        String[]  KSTY2 = new String[]{"K", "S", "T", "Y", "nterm"};
        String[]  DE2 = new String[]{"D", "E", "cterm"};
        String[]  X2 = new String[]{"X"};

        default_xl_xi2.put("BS3", 
                new Xi2Crosslinker(BS3.getName(), BS3.getCrossLinkedMass(), new String[][] {KSTY2,KSTY2}));
        default_xl_xi2.put("EDC", 
                new Xi2Crosslinker(EDC.getName(), EDC.getCrossLinkedMass(), new String[][] {KSTY2,DE2}));
        default_xl_xi2.put("SDA", 
                new Xi2Crosslinker(SDA.getName(), SDA.getCrossLinkedMass(), new String[][] {KSTY2,X2}));
        default_xl_xi2.put("DSSO", 
                new Xi2Crosslinker(DSSO.getName(), DSSO.getCrossLinkedMass(), new String[][] {KSTY2,KSTY2}));

        for (CrossLinker xl: default_xl_xi1.values()){
            default_xl_masses.put(xl.getName(), xl.getCrossLinkedMass());        
        }
        
    }
    
    
    /**
     * Java representation of the xiSEARCH2 crosslinker definition. 
     * Only the here needed information are retained
     */
    public class Xi2Crosslinker {
        public String name;
        public Double mass ;
        String[][] specificity = new String[2][];
        
        /**
         * initialise just with name and get the rest from default definitions.
         * @param name Name of the crosslinker
         */
        public Xi2Crosslinker(String name) {
            this.name = name;
            if (default_xl_masses.containsKey(name.toUpperCase())) {
                this.mass = default_xl_masses.get(name.toUpperCase());
            } else {
                this.mass = Double.NaN;
            }
        }
        
        /**
         * Initialise by name and mass.
         * @param name Name of the crosslinker
         * @param mass mass of the reacted crosslinker
         */
        public Xi2Crosslinker(String name, double mass) {
            this.name = name;
            this.mass = mass;
        }

        /**
         * Initialise by name, mass, and specificity.
         * @param name Name of the crosslinker
         * @param mass mass of the reacted crosslinker
         * @param specificity where can the crosslinker react
         */
        public Xi2Crosslinker(String name, double mass, String[][] specificity) {
            this.name = name;
            this.mass = mass;
            this.specificity = specificity;
        }

        /**
         * initialise based on a json crosslinker defintion derived map
         * @param m 
         */
        public Xi2Crosslinker(Map m) {
            this.name = m.get("name").toString();
            this.mass = (Double) m.get("mass");
            Object s = m.get("specificity");
            if (s instanceof String) {
                this.specificity[0] = new String[] {s.toString()};
                this.specificity[1] = this.specificity[0];
            } else {
                List sa = (List) s;
                if (sa.get(0) instanceof String) {
                    this.specificity[0] = new String[sa.size()];
                    for (int i =0; i<sa.size();i++)
                        this.specificity[0][i]=sa.get(i).toString();
                    this.specificity[1] = this.specificity[0];
                } else {
                    List sa0 = (List)sa.get(0);
                    this.specificity[0] = new String[sa0.size()];
                    for (int i =0; i<sa0.size();i++)
                        this.specificity[0][i]=sa0.get(i).toString();
                    if (sa.size() == 1)
                        this.specificity[1] = this.specificity[0];
                    else {
                        List sa1 = (List)sa.get(1);
                        this.specificity[1] = new String[sa1.size()];
                        for (int i =0; i<sa1.size();i++)
                            this.specificity[1][i]=sa1.get(i).toString();
                    }
                }
            }
        }
        
        /**
         * Convert to a xiSEARCH1 crosslinker definition.
         * @return
         * @throws java.text.ParseException
         * @throws ConfigurationParserException 
         */
        public String toXi1Crosslinker() throws java.text.ParseException, ConfigurationParserException {
            StringBuilder sb = new StringBuilder("crosslinker:AsymetricSingleAminoAcidRestrictedCrossLinker:NAME:");
            sb.append(this.name);
            sb.append(";FIRSTLINKEDAMINOACIDS:");
            boolean first = true;
            for (String aa: this.specificity[0]) {
                if (first) {
                    first = false;
                    sb.append(aa);
                } else {
                    sb.append(",").append(aa);
                }
            }
            sb.append(";SECONDLINKEDAMINOACIDS:");
            first = true;
            for (String aa: this.specificity[0]) {
                if (first) {
                    first = false;
                    sb.append(aa);
                } else {
                    sb.append(",").append(aa);
                }
            }
            sb.append(";MASS:").append(this.mass);
            return sb.toString();
            //return AsymetricSingleAminoAcidRestrictedCrossLinker.parseArgs(sb.toString(), DUMMYCONFIG);
        }
        
    }
    public class Xi2Modification {
        String symbol;
        String type;
        double mass;
        String[] specificity;

        public Xi2Modification(String symbol, String modtype) {
            this.symbol = symbol;
            this.type = modtype;
        }
        public Xi2Modification(String symbol, String modtype, double mass) {
            this.symbol = symbol;
            this.type = modtype;
            this.mass = mass;
        }

        public Xi2Modification(Map m) {
            this.mass = 0;
            this.symbol = m.get("name").toString();
            this.type = m.get("type").toString();
            if (m.get("specificity") instanceof List) {
                List specList = (List) m.get("specificity");
                this.specificity  =new String[specList.size()];
                int sid = 0;
                for (Object s : specList) {
                    specificity[sid++]=s.toString();
                }
            } else {
                this.specificity  =new String[1];
                this.specificity[0] = m.get("specificity").toString();
            }
            if (m.containsKey("mass"))
                mass = (Double) m.get("mass");
            else
                mass = Composition.formula2mass(m.get("composition").toString());
            
        }
        
        public String toxi1Mod() throws java.text.ParseException {
            StringBuilder sb = new StringBuilder("modification:variable::SYMBOLEXT:");
            sb.append(symbol);
            sb.append(";MODIFIED:");
            boolean first  =true;
            for (String aa: this.specificity) {
                if (first) {
                    first = false;
                    sb.append(aa);
                } else {
                    sb.append(",").append(aa);
                }
            }
                    
            sb.append(";DELTAMASS:").append(this.mass);
            return sb.toString();
            //return AminoModification.parseArgs(sb.toString(), DUMMYCONFIG);
        }
        
        
    }
    public ArrayList<Xi2Crosslinker> xi2crosslinker = new ArrayList<>();
    public ArrayList<Xi2Modification> xi2modifications = new ArrayList<>();
    
    public String textConfig;

    public Xi2Xi1Config() {
    }
    
    public Xi2Xi1Config(String config) {
        try {
            loadFromString(config);
        } catch (ParseException |java.text.ParseException|ConfigurationParserException ex) {
            Logger.getLogger(Xi2Xi1Config.class.getName()).log(Level.SEVERE, "Error converting xi2 config to xi1 config", ex);
        }
        
    }

    private void loadFromString(String config) throws java.text.ParseException, ConfigurationParserException, ParseException {
        this.xi2crosslinker = new ArrayList<>(1);
        this.xi2modifications = new ArrayList<>();
        this.textConfig = config;
        JSONParser prsr = new JSONParser();
        Map json = (Map)prsr.parse(config);
        Object xls = parseListSetting(json.get("crosslinker"));
        Object mod_peptide_syntax = json.get("mod_peptide_syntax");
        if (mod_peptide_syntax != null && mod_peptide_syntax.toString().contentEquals("modX")) {
            this.isModX = true;
        }
        List crosslinker;
        if (xls instanceof List)
            crosslinker = (List) xls;
        else {
            crosslinker = new ArrayList(1);
            crosslinker.add(xls);
        }

        for (Object xlo : crosslinker) {
            String name;
            Xi2Crosslinker xl;
            if (xlo instanceof Map) {
                xl = new Xi2Crosslinker((Map)xlo);
            } else {
                name = xlo.toString();
                xl = default_xl_xi2.get(name);
            }
            this.xi2crosslinker.add(xl);                
            this.evaluateConfigLine(xl.toXi1Crosslinker());
            //this.addCrossLinker( new AsymetricSingleAminoAcidRestrictedCrossLinker(xl.name,xl.mass,xl.mass,PrimaryLinkableAminoAcids, SecondaryLinkableAminoAcids));
        }

        List modifications = parseListSetting(((Map)json.get("modification")).get("modifications"));
        for (Object mo : modifications) {
            Map m = (Map)mo;                
            Xi2Modification mod = new Xi2Modification(m);
            this.xi2modifications.add(mod);
            this.evaluateConfigLine(mod.toxi1Mod());
        }
        // tolerances
        this.setPrecoursorTolerance(new ToleranceUnit(json.get("ms1_tol").toString()));
        this.setFragmentTolerance(new ToleranceUnit(json.get("ms2_tol").toString()));
        
        // enzyme info
        Map digestConfjson = (Map)json.get("digestion");
        List enzymes = parseListSetting(digestConfjson.get("enzymes"));
        Object e0 = enzymes.get(0);
        AAConstrainedDigestion xi1digest = null;
        if (enzymes.get(0) instanceof Map) {
            Map enzyme = (Map) e0;
            List nspecificity = parseListSetting(enzyme.get("nterminal_of"));
            List cspecificity = parseListSetting(enzyme.get("cterminal_of"));
            List rspecificity = parseListSetting(enzyme.get("restraining"));
            AminoAcid[] nterm = 
                    nspecificity == null ? new AminoAcid[0] : new AminoAcid[nspecificity.size()];
            AminoAcid[] cterm = 
                    cspecificity == null ? new AminoAcid[0] : new AminoAcid[cspecificity.size()];
            AminoAcid[] rest = 
                    rspecificity == null ? new AminoAcid[0] : new AminoAcid[rspecificity.size()];
            int sid = 0;
            if (nspecificity != null) for (Object a : nspecificity) {
                nterm[sid++] = this.getAminoAcid(a.toString());
            }
            sid = 0;
            if (cspecificity != null) for (Object a : cspecificity) {
                cterm[sid++] = this.getAminoAcid(a.toString());
            }
            sid = 0;
            if (rspecificity != null) for (Object a : rspecificity) {
                rest[sid++] = this.getAminoAcid(a.toString());
            }

            String name = enzyme.get("name").toString();
            xi1digest = new AAConstrainedDigestion(nterm, cterm, rest,rest, this);
            xi1digest.setName(name);
        } else {
            String se = e0.toString();
            if (se.toLowerCase().contentEquals("trypsin")) {
                AminoAcid[] nterm = new AminoAcid[]{AminoAcid.K, AminoAcid.R};
                AminoAcid[] cterm = new AminoAcid[0];
                AminoAcid[] rest = new AminoAcid[]{AminoAcid.P};
                String name = "trypsin";
                xi1digest = new AAConstrainedDigestion(nterm, cterm, rest,rest, this);
                xi1digest.setName(name);
            } else if (se.toLowerCase().contentEquals("asp_n")) {
                AminoAcid[] nterm = new AminoAcid[]{AminoAcid.D};
                AminoAcid[] cterm = new AminoAcid[0];
                AminoAcid[] rest = new AminoAcid[0];
                String name = "asp-n";
                xi1digest = new AAConstrainedDigestion(nterm, cterm, rest,rest, this);
                xi1digest.setName(name);
            } 
        }
        xi1digest.setMaxMissCleavages(Integer.parseInt(digestConfjson.getOrDefault("missed_cleavages", 0).toString()));
        this.setDigestion(xi1digest);
        
        // basic ions
        Map fragmentation  = (Map)json.get("fragmentation");
        List nterm = parseListSetting(fragmentation.getOrDefault("nterm_ions","b"));
        List cterm = parseListSetting(fragmentation.getOrDefault("cterm_ions","y"));
        ArrayList allIons = new ArrayList(nterm);
        allIons.addAll(cterm);
        for (Object i : allIons) {
            String ion = i.toString();
            try {
                if (ion.contentEquals("a")) {
                    this.getFragmentMethods().add(rappsilber.ms.sequence.ions.AIon.class.getMethod("fragment", rappsilber.ms.sequence.Peptide.class));
                } else if (ion.contentEquals("b")) {
                    this.getFragmentMethods().add(rappsilber.ms.sequence.ions.BIon.class.getMethod("fragment", rappsilber.ms.sequence.Peptide.class));
                } else if (ion.contentEquals("c")) {
                    this.getFragmentMethods().add(rappsilber.ms.sequence.ions.CIon.class.getMethod("fragment", rappsilber.ms.sequence.Peptide.class));
                } else if (ion.contentEquals("x")) {
                    this.getFragmentMethods().add(rappsilber.ms.sequence.ions.XIon.class.getMethod("fragment", rappsilber.ms.sequence.Peptide.class));
                } else if (ion.contentEquals("y")) {
                    this.getFragmentMethods().add(rappsilber.ms.sequence.ions.YIon.class.getMethod("fragment", rappsilber.ms.sequence.Peptide.class));
                } else if (ion.contentEquals("z")) {
                    this.getFragmentMethods().add(rappsilber.ms.sequence.ions.ZIon.class.getMethod("fragment", rappsilber.ms.sequence.Peptide.class));
                } 
            } catch (NoSuchMethodException ex) {
                Logger.getLogger(Xi2Xi1Config.class.getName()).log(Level.SEVERE, null, ex);
            } catch (SecurityException ex) {
                Logger.getLogger(Xi2Xi1Config.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    public Xi2Xi1Config(File f) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(f));
        StringBuilder sb = new StringBuilder();
        String line;
        while ((line = br.readLine())!= null)
            sb.append(line);
        try{
            loadFromString(sb.toString());
        } catch (ParseException |java.text.ParseException|ConfigurationParserException ex) {
            throw new IOException("Error converting xi2 config to xi1 config",ex);
        }
}
    
    public Xi2Xi1Config(RunConfig config) {
        if (config instanceof Xi2Xi1Config) {
            this.xi2crosslinker = new ArrayList<>(((Xi2Xi1Config)config).xi2crosslinker);
            this.xi2modifications = new ArrayList<>(((Xi2Xi1Config)config).xi2modifications);
        }
        for (AminoModification am :  config.getAllModifications()) {
            this.addKnownModification(am);
        }
        for (CrossLinker xl : config.getCrossLinker())
            this.addCrossLinker(xl);
        this.setDigestion(config.getDigestion_method());
        this.setFragmentTolerance(config.getFragmentTolerance());
        this.setPrecoursorTolerance(config.getPrecousorTolerance());
    }

    
    public Xi2Xi1Config add(Xi2Xi1Config o) {
        return add(o,true);
    }

    public Xi2Xi1Config add(Xi2Xi1Config o, boolean self) {
        Xi2Xi1Config ret;
        if (self) {
            ret = this;
        } else {
            ret = new Xi2Xi1Config(this);
        }
        ArrayList<Xi2Modification> addNewMod = new ArrayList<>();
        oloop: for (Xi2Modification mo : o.xi2modifications) {
            for (Xi2Modification mt : ret.xi2modifications) {
                if (mo.symbol.contentEquals(mt.symbol) && mo.type.contentEquals(mt.type))
                    continue oloop;
            }
            addNewMod.add(mo);
            try {
                this.evaluateConfigLine(mo.toxi1Mod());
            } catch (java.text.ParseException ex) {
                Logger.getLogger(Xi2Xi1Config.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        ret.xi2modifications.addAll(addNewMod);
        
        ArrayList<Xi2Crosslinker> addNewXL = new ArrayList<>();
        oloop: for (Xi2Crosslinker co : o.xi2crosslinker) {
            for (Xi2Crosslinker ct : ret.xi2crosslinker) {
                if (co.name.contentEquals(ct.name) && Math.round(co.mass * 1000) == Math.round(ct.mass * 1000))
                    continue oloop;
            }
            addNewXL.add(co);
            try {
                ret.evaluateConfigLine(co.toXi1Crosslinker());
                //ret.addCrossLinker(co.toXi1Crosslinker());
            } catch (java.text.ParseException ex) {
                Logger.getLogger(Xi2Xi1Config.class.getName()).log(Level.SEVERE, null, ex);
            } catch (ConfigurationParserException ex) {
                Logger.getLogger(Xi2Xi1Config.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        ret.xi2crosslinker.addAll(addNewXL);
        
        return ret;
    }
    
    public List parseListSetting(Object json) {
        if (json == null)
            return null;
        if (json instanceof List)
            return (List)json;
        List ret = new ArrayList();
        ret.add(json);
        return ret;
    }
    
    public String peptPeptide(String basepeptide, int[] mods) {
        StringBuffer ret = new StringBuffer();
        if (mods[0] != 0)
            ret.append(this.xi2modifications.get(mods[0]).symbol);
        for (int i = 0; i<basepeptide.length(); i++) {
            ret.append(basepeptide.substring(i,i+1));
            if (mods[2+i] != 0)
                ret.append(this.xi2modifications.get(mods[0]).symbol);
        }
        if (mods[1] != 0)
            ret.append(this.xi2modifications.get(mods[1]).symbol);
        return ret.toString();
        
    }
    

    public static void main(String[] args) {
        String conf = "{\n" +
                "  \"reporting_requirements\": {\n" +
                "    \"report_top_ranking_only\": true\n" +
                "  },\n" +
                "  \"digestion\": {\n" +
                "    \"enzymes\": [\"trypsin\"],\n" +
                "    \"missed_cleavages\": 2,\n" +
                "    \"min_peptide_length\": 2\n" +
                "  },\n" +
                "  \"isotope_error_ximpa\": 2,\n" +
                "  \"ms1_tol\": \"3ppm\",\n" +
                "  \"ms2_tol\": \"5ppm\",\n" +
                "  \"top_n_alpha_scores\": 10,\n" +
                "  \"top_n_alpha_beta_scores\": 10,\n" +
                "  \"crosslinker\": \"BS3\",\n" +
                "  \"conservative_n_multi_loss\": 3,\n" +
                "  \"denoise_alpha\": {\n" +
                "    \"top_n\": 10,\n" +
                "    \"bin_size\": 100\n" +
                "  },\n" +
                "  \"denoise_alpha_beta\": {\n" +
                "    \"top_n\": 20,\n" +
                "    \"bin_size\": 100\n" +
                "  },\n" +
                "  \"fragmentation\": {\n" +
                "    \"nterm_ions\": [\"b\"],\n" +
                "    \"cterm_ions\": [\"y\"],\n" +
                "    \"add_precursor\": true,\n" +
                "    \"max_nloss\": 4,\n" +
                "    \"match_missing_monoisotopic\": true,\n" +
                "    \"losses\": [\n" +
                "      {\n" +
                "        \"name\": \"H2O\",\n" +
                "        \"specificity\": [\"S\", \"T\", \"D\", \"E\", \"cterm\"],\n" +
                "        \"composition\": \"\",\n" +
                "        \"mass\": 18.01056027\n" +
                "      },\n" +
                "      {\n" +
                "        \"name\": \"NH3\",\n" +
                "        \"specificity\": [\"R\", \"K\", \"N\", \"Q\", \"nterm\"],\n" +
                "        \"composition\": \"\",\n" +
                "        \"mass\": 17.02654493\n" +
                "      },\n" +
                "      {\n" +
                "        \"name\": \"CH3SOH\",\n" +
                "        \"specificity\": [\"oxM\"],\n" +
                "        \"composition\": \"\",\n" +
                "        \"mass\": 63.99828547\n" +
                "      }\n" +
                "    ]\n" +
                "  },\n" +
                "  \"modification\": {\n" +
                "    \"max_var_protein_mods\": 3,\n" +
                "    \"max_modified_peps\": 20,\n" +
                "    \"modifications\": [\n" +
                "      {\n" +
                "        \"name\": \"cm\",\n" +
                "        \"specificity\": [\"C\"],\n" +
                "        \"type\": \"fixed\",\n" +
                "        \"composition\": \"C2H3N1O1\"\n" +
                "      },\n" +
                "      {\n" +
                "        \"name\": \"ox\",\n" +
                "        \"specificity\": [\"M\"],\n" +
                "        \"type\": \"variable\",\n" +
                "        \"composition\": \"O1\"\n" +
                "      }\n" +
                "    ]\n" +
                "  }\n" +
                "}";
        Xi2Xi1Config xiconf = new Xi2Xi1Config(conf);
        for (Xi2Xi1Config.Xi2Crosslinker cl : xiconf.xi2crosslinker) {
            System.out.println("found crosslinker " + cl.name +"(" + cl.mass +")");
        }
        for (Xi2Xi1Config.Xi2Modification mod : xiconf.xi2modifications) {
            System.out.println("found modification " + mod.symbol +"(" + mod.mass +")");
        }
        
    }
}
