/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package comparebacngssvcalls;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import org.apache.commons.lang3.Range;

/**
 *
 * @author wim
 */
public class Deletion {
    
    Range<Integer> location;
    String chromosome;
    
    Boolean common;
    
    
    TreeMap<String, Double> ratStrainsGenotype;
    

    public Deletion(Range<Integer> location, String chromosome) {
        this.location = location;
        this.chromosome = chromosome;
        
        ratStrainsGenotype = new TreeMap<String, Double>();
    }
    
    public void updateLocationEnd(Integer locationEnd)
    {
        location = Range.between(location.getMinimum(), locationEnd);    
    }  
    

    public Range<Integer> getLocation() {
        return location;
    }

    public void setLocation(Range<Integer> location) {
        this.location = location;
    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }
    
    public Integer getSize()
    {
        return location.getMaximum() - location.getMinimum();
    }
    
    public void addRatGenotype(String ratName, Double genotype)
    {
        ratStrainsGenotype.put(ratName, genotype);
    }

    public TreeMap<String, Double> getRatStrainsGenotype() {
        return ratStrainsGenotype;
    }

    public Boolean getCommon() {
        return common;
    }

    public void setCommon(Boolean common) {
        this.common = common;
    }
    
    
    
    

    
    
    
    
    
    
    
    
}
