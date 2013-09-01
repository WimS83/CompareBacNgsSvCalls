/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package comparebacngssvcalls;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.Map;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.apache.commons.lang3.Range;

/**
 *
 * @author wim
 */
public class BamAlignmentRegionsParser {
    
    
    File bamFile;

    public BamAlignmentRegionsParser(File bamFile) {
        this.bamFile = bamFile;
    }
    
    
    
    
    
    
    public EnumMap<RatChromosomes, ArrayList<Range<Integer>>> readBacRegions() throws IOException {

        BufferedReader bacRegionsBR = new BufferedReader(new FileReader(bamFile));
        
        
        EnumMap<RatChromosomes, ArrayList<Range<Integer>>> alignmentsListPerChromosome = new EnumMap<RatChromosomes, ArrayList<Range<Integer>>>(RatChromosomes.class);
        for(RatChromosomes ratChromosomes: RatChromosomes.values())
        {
            ArrayList alignment = new ArrayList<Range<Integer>>();
            alignmentsListPerChromosome.put(ratChromosomes, alignment);
        }
        
        
        SAMFileReader bamreader = new SAMFileReader(bamFile);
        

       for(SAMRecord samRecord : bamreader)
        {
            //skip alignment with low mapping qual
            if(samRecord.getMappingQuality() < 60 ){continue;}
            
            //get the start of the alignment
            String currentChrom = samRecord.getReferenceName();
            if(!currentChrom.contains("chr")){currentChrom = "chr"+currentChrom;}
            
            RatChromosomes ratChromosome;
            try{                
                ratChromosome = RatChromosomes.valueOf(currentChrom);
            }
            //continue to next line if unknown chromosome
            catch(IllegalArgumentException  ex)
            {
                continue;
            }                       
            
            
            Integer begin = samRecord.getAlignmentStart();  
            Integer end = samRecord.getAlignmentEnd();

            Range<Integer> range = Range.between(begin, end);
            alignmentsListPerChromosome.get(ratChromosome).add(range);
        }

        return alignmentsListPerChromosome;

    }
    
    
    
}
