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

package org.rappsilber.config;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Properties;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Just small wrapper around java properties to persist some settings - predominantly last used folders.
 * @author Lutz Fischer &lt;lfischer at staffmail.ed.ac.uk&gt;
 */
public  class LocalProperties {
    private static final long serialVersionUID = -2872695155848651883L;
    
    /**
     * folder under which to store persistent settings
     */
    static String homeDir = System.getProperty("user.home");
    /**
     * where to remember the last opened folders
     */
    static String userPropertiesFile = homeDir + "/.Xlink.properties";

    static Properties localProperties = new Properties();
    /**
     * some default names for properties 
     */
    public static final String LAST_MSM_FOLDER = "LastMSMFolder";
    /**
     * some default names for properties 
     */
    public static final String LAST_SEQUNECE_FOLDER = "LastSequenceFolder";
    
    /**
     * load the default property file
     */
    static {
        try {
            localProperties.load(new FileInputStream(userPropertiesFile));
        } catch (IOException ex) {
        }
    }

    /**
     * The last folder that contained the last selected peak-list
     * @return the folder that was last accessed. If no previous folder is stored the user-home-folder is returned
     */
    public static File getLastMSMFolder() {
        return getFolder(LAST_MSM_FOLDER);
    }

    /**
     * The last folder that contained the last selected peak-list
     * @param path set the last used folder
     * @return the previously stored value
     */
    public static Object setLastMSMFolder(String path) {
        return setFolder(LAST_MSM_FOLDER, path);
    }
    
    /**
     * The last folder that contained the last selected peak-list
     * @param path set the last used folder
     * @return the previously stored value
     */
    public static Object setLastMSMFolder(File path) {
        return setFolder(LAST_MSM_FOLDER, path);
    }


    /**
     * Returns the last accessed folder for the given key
     * 
     * @param key a identifier to for the type of folder/action
     * @return the folder that was last accessed under that key. If no previous folder is stored the user-home-folder is returned
     */
    public static File getFolder(String key) {
        if (key != null)
            return new File(localProperties.getProperty(key,homeDir));
        else
            return new File(homeDir);
    }

    /**
     * Returns the last accessed folder for the given key
     * 
     * @param key a identifier to for the type of folder/action
     * @param path the path to store under the given key
     * @return the previously stored value
     */
    public static Object setFolder(String key, File path) {
        return setProperty(key, path.getAbsolutePath());
    }

    /**
     * Returns the last accessed folder for the given key
     * 
     * @param key a identifier to for the type of folder/action
     * @param path the path to store under the given key
     * @return the previously stored value
     */
    public static Object setFolder(String key, String path) {
        return setProperty(key, path);
    }


    /**
     * Writes an arbitrary key value pair into the default-properties file
     * @param key the key under what to write the value
     * @param value the (surprise) value to store
     * @return the previously stored value
     */
    public static synchronized Object setProperty(String key, String value)  {
        String old = localProperties.getProperty(key);
        if ((old == null && value != null) || old.contentEquals(value)) {
            Object ret = localProperties.setProperty(key, value);
            try {
                localProperties.store(new FileOutputStream(userPropertiesFile), "XLink local properties file");
            } catch (FileNotFoundException ex) {
                Logger.getLogger(LocalProperties.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                Logger.getLogger(LocalProperties.class.getName()).log(Level.SEVERE, null, ex);
            }
            return ret;
        } else
            if (value == null) {
                return localProperties.remove(key);
            } else 
                return localProperties.setProperty(key, value);
    }


    /**
     * returns the value of the given property
     * @param key name/key of the property
     * @return value of the property
     */
    public static String getProperty(String key) {
        return localProperties.getProperty(key);
    }

    /**
     * returns the value of the given property
     * @param key name/key of the property
     * @param defaultValue default value to be returned if nothing was stored
     * @return the stored value for this key or the default value if nothing was stored
     */
    public static String getProperty(String key, String defaultValue) {
        return localProperties.getProperty(key, defaultValue);
    }

}
