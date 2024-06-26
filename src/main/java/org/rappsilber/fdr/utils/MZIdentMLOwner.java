/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr.utils;

import java.util.Objects;
import org.rappsilber.config.LocalProperties;

/**
 *
 * @author lfischer
 */
public class MZIdentMLOwner {
    public String first;
    public String last;
    public String email;
    public String org;
    public String address;

    public final static String propertyFirst = "mzIdenMLOwnerFirst";
    public final static String propertyLast = "mzIdenMLOwnerLast";
    public final static String propertyEMail = "mzIdenMLOwnerEmail";
    public final static String propertyAddress = "mzIdenMLOwnerAddress";
    public final static String propertyOrg = "mzIdenMLOwnerOrg";
    
    public MZIdentMLOwner() {
        
    }
    public MZIdentMLOwner(String first, String last, String email, String org, String adress) {
        this.first = first;
        this.last = last;
        this.email = email;
        this.org = org;
        this.address = adress;
    }
    
    public void readLast() {
        first = LocalProperties.getProperty(propertyFirst, "");
        last = LocalProperties.getProperty(propertyLast, "");
        email = LocalProperties.getProperty(propertyEMail, "");
        address = LocalProperties.getProperty(propertyAddress, "");
        org = LocalProperties.getProperty(MZIdentMLOwner.propertyOrg, "");

    }
    
    public boolean isEmpty() {
        return (this.first == null  || this.first.isEmpty()) &&
                (this.last == null  || this.last.isEmpty())  &&
                (this.email == null || this.email.isEmpty()) &&
                (this.org == null   || this.email.isEmpty()) &&
                (this.address == null|| this.address.isEmpty());
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof MZIdentMLOwner) {
            MZIdentMLOwner o = (MZIdentMLOwner) obj;
            return ((this.first == null && o.first == null) || (this.first != null && this.first.contentEquals(o.first))) &&
                    ((this.last == null && o.last == null) || (this.last != null && this.last.contentEquals(o.last))) &&
                    ((this.email == null && o.email == null) || (this.email != null && this.email.contentEquals(o.email))) &&
                    ((this.org == null && o.org == null) || (this.org != null && this.org.contentEquals(o.org))) &&
                    ((this.address == null && o.address == null) || (this.address != null && this.address.contentEquals(o.address)));
        }
        return false;
    }

    @Override
    public int hashCode() {
        // somewhat autogenerated
        int hash = 5;
        hash = 97 * hash + this.first.hashCode();
        hash = 97 * hash + this.last.hashCode();
        hash = 97 * hash + this.email.hashCode();
        hash = 97 * hash + this.org.hashCode();
        hash = 97 * hash + this.address.hashCode();
        return hash;
    }

    
}
