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
package org.rappsilber.fdr.utils;

import java.io.IOException;
import java.io.Writer;

/**
 * A writer, that replaces every occurrence of a string in the the data 
 * to be written with another string.
 * To do so a minimum of data is temporarily stored in a buffer and
 * replacement is done on that buffer.
 * @author lfischer@staffmail.ed.ac.uk
 * 
 */
public class StreamReplaceWriter extends Writer{
    /**
     * The actual writer that where the data are forwarded to.
     */
    private Writer innerWriter;
    /**
     * the search-term
     */
    private String search;
    /**
     * String that should replace every occurrence of the search term
     */
    private String replace;
    /**
     * Buffer that keeps enough of the data to make sure that we can replace every 
     * occurrence of search with replace
     */
    private StringBuffer sb = new StringBuffer();

    /**
     * 
     * @param innerWriter where are the - possibly modified - data to be written 
     * to.
     * @param search term to be searched
     * @param replace what replaces search
     */
    public StreamReplaceWriter(Writer innerWriter, String search, String replace) {
        this.innerWriter = innerWriter;
        this.search = search;
        this.replace = replace;
    }
    
    @Override
    public void write(char[] cbuf, int off, int len) throws IOException {
        
        char[] out = new char[len];
        System.arraycopy(cbuf, off, out, 0, len);
        sb.append(out);
        int pos = 0;
        int writeOut =-1;
        //replace every occurence of search with replace
        while ((pos = sb.indexOf(search)) >=0) {
            sb.replace(pos, pos+search.length(), replace);
            pos = pos+replace.length()-1;
            writeOut = pos;
        }
        // we can write everything out that was replaced or is more then search-length away from the end
        writeOut = Math.max(writeOut, sb.length()-search.length());
        
        // write out everything that we can
        if (writeOut >0) {
            String wos = sb.substring(0,writeOut);
            sb.delete(0, writeOut);
            innerWriter.write(wos);
        }
        
        
    }

    /**
     * Forwards the flush to the inner writer
     * We cant really write out everything in our buffer - as we don't know yet
     * whether anything in it has to be replaced.<br/>
     * So it is not really a proper flush of this writer.
     * @throws IOException 
     */
    @Override
    public void flush() throws IOException {
        innerWriter.flush();
    }

    @Override
    public void close() throws IOException {
        innerWriter.write(sb.toString());
        sb.setLength(0);
        innerWriter.close();
    }

    
    
    
    
}
