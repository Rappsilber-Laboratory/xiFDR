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
package org.rappsilber.fdr;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.sql.SQLException;
import org.rappsilber.fdr.gui.FDRGUI;
import java.lang.management.ManagementFactory;
import java.lang.management.RuntimeMXBean;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Proxy class that looks at the given arguments and then forwards the whole 
 * command-line to the appropriate class.
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class CommandLine {
    public static List<String> getVMArguments() {
        RuntimeMXBean runtimeMxBean = ManagementFactory.getRuntimeMXBean();
        return runtimeMxBean.getInputArguments();        
    }
    
    public static void restart(String args[]) {
        List<String> vmArgs = getVMArguments();
        boolean addopenlang = true;
        boolean addopenutil = true;
        boolean forceutf8 = true;
        ProcessBuilder builder = new ProcessBuilder();
        LinkedList<String> newargs = new LinkedList<>();
        for (String a : vmArgs) {
            if (a.contains("java.base/java.lang=ALL-UNNAMED"))
                addopenlang = false;
            if (a.contains("java.base/java.util=ALL-UNNAMED"))
                addopenutil = false;                
            if (a.contains("java.base/java.util=ALL-UNNAMED"))
                addopenutil = false;   
            if (a.contains("file.encoding=UTF-8")) {
                forceutf8 = false;
            }
        }
        int java_ver = org.rappsilber.utils.JavaUtils.getJavaMajorVersion();

        newargs.add(org.rappsilber.utils.JavaUtils.findJava());
        // -jar xiFDR.jar
        newargs.add("-Dfile.encoding=UTF-8");
        if (java_ver >= 11 && addopenlang) {
            newargs.add("--add-opens");
            newargs.add("java.base/java.lang=ALL-UNNAMED");
        }
        if (java_ver >= 16 && addopenutil) {
            newargs.add("--add-opens");
            newargs.add("java.base/java.util=ALL-UNNAMED");
        }
        if (java_ver >= 16 && addopenutil) {
            newargs.add("--add-opens");
            newargs.add("java.base/java.util=ALL-UNNAMED");
        }
        if (addopenutil && java_ver>= 9) {
            newargs.add("--add-opens");
            newargs.add("java.base/sun.reflect.annotation=ALL-UNNAMED");
        }
        newargs.add("-cp");
        newargs.add(System.getProperty("java.class.path"));
        newargs.addAll(vmArgs);
        newargs.add(CommandLine.class.getCanonicalName());
        newargs.add("--got-restarted");
        newargs.addAll(Arrays.asList(args));
        Process p = null;
        if (forceutf8 || (addopenlang && java_ver >= 11) || (addopenutil && java_ver >= 16)) {
            Logger.getLogger(CommandLine.class.getName()).log(Level.INFO, "some suggested args are missing - trying to restart");
            try {
                Logger.getLogger(CommandLine.class.getName()).log(Level.INFO, "Arguments: " + newargs.toString());
                builder.command(newargs);
                p = builder.start();
                final InputStream es = p.getErrorStream();
                final InputStream is = p.getInputStream();
                Runnable osTransfer = new Runnable() {
                    @Override
                    public void run() {
                        byte[] buf = new byte[8192];
                        int length;
                        try {
                            while ((length = is.read(buf)) != -1) {
                                System.out.write(buf, 0, length);
                            }
                        } catch (IOException ex) {
                            Logger.getLogger(CommandLine.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }
                };
                Runnable esTransfer = new Runnable() {
                    @Override
                    public void run() {
                        try {
                            byte[] buf = new byte[8192];
                            int length;
                            while ((length = es.read(buf)) != -1) {
                                System.err.write(buf, 0, length);
                            }
                        } catch (IOException ex) {
                            Logger.getLogger(CommandLine.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }
                };
                Thread stdForward = new Thread(osTransfer, "STD Forward");
                stdForward.start();
                Thread errForward = new Thread(esTransfer, "ERR Forward");
                errForward.start();
                
                    
            } catch (IOException ex) {
                Logger.getLogger(CommandLine.class.getName()).log(Level.SEVERE, "Could not restart with suggested options - will try without", ex);

            }
            try {        
                int ret  = p.waitFor();
                System.exit(ret);
            } catch (InterruptedException ex) {
                Logger.getLogger(CommandLine.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
                
    }
    
    public static void main(String args[]) throws SQLException, FileNotFoundException {
        boolean xicsv  = false;
        boolean csv = false;
        boolean gui = true;
        boolean db = false;
        boolean dorestart=false;

        for (String a : args) {
            if (a.startsWith("--xiconfig=") || a.contentEquals("--gui")) {
                xicsv=true;
                gui = false;
            } else if (a.startsWith("--dbids=")) {
                db = true;
                gui = false;
            } else if (!a.startsWith("-")) {
                csv = true;
                gui = false;
            } else if (a.contentEquals("--got-restarted")) {
                dorestart = false;
            } else {
                gui = false;
            }
        }
        if (dorestart) {
            restart(args);
        }
        if (db && !csv) {
            DBinFDR.main(args);
            return;
        }
        if (xicsv && !db) {
            XiCSVinFDR.main(args);
            return;
        }
        if (csv && (!db) && (!xicsv)) {
            CSVinFDR.main(args);
            return;
        }
        if (gui) {
            FDRGUI.main(args);
            return;
        }
        new XiCSVinFDR().printUsage();
    }
}
