/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package comparebacngssvcalls;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.apache.commons.lang3.Range;

/**
 *
 * @author wim
 */
public class PiclParser {

    File piclFile;
    BufferedReader piclBR;
    
    
    
    public PiclParser(File piclFile) throws FileNotFoundException {
        this.piclFile = piclFile;
        
        piclBR = new BufferedReader(new FileReader(piclFile));    
        
    }
    
    public  EnumMap<RatChromosomes, ArrayList<Deletion>> readPiclSVCalls() throws IOException {
       
       // Map<String, ArrayList<Deletion>> chromPiclSVCallsMap = new HashMap<String, ArrayList<Deletion>>();
        
        //initialize map to store list of deletions per ratchromosome
        EnumMap<RatChromosomes, ArrayList<Deletion>> deletionListPerChromosome = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);
        for(RatChromosomes ratChromosomes: RatChromosomes.values())
        {
            ArrayList deletions = new ArrayList<Deletion>();
            deletionListPerChromosome.put(ratChromosomes, deletions);
        }
        
        

        String line = piclBR.readLine();

        while (line != null) {
            String[] splitLine = line.split("\t");
            Integer begin = new Integer(splitLine[1]);
            Integer end = new Integer(splitLine[5]);
            String type = splitLine[9];            
            if(!type.equalsIgnoreCase("deletion")){
                line = piclBR.readLine();
                continue;
            }               
            
            Integer readPairSupport = new Integer(splitLine[7]);
            if(readPairSupport < 4 )
            {
                 line = piclBR.readLine();
                continue;
                
            }
            
            
            String chrom = splitLine[0];
            
            if(!chrom.contains("chr")){chrom = "chr"+chrom;}
            
            RatChromosomes ratChromosome;
            try{                
                ratChromosome = RatChromosomes.valueOf(chrom);
            }
            //continue to next line if unknown chromosome
            catch(IllegalArgumentException  ex)
            {
                line = piclBR.readLine();
                continue;
            }                       
            
            Deletion deletion = new Deletion(Range.between(begin, end), chrom);            
            deletionListPerChromosome.get(ratChromosome).add(deletion);      

            line = piclBR.readLine();
            String blaat = "blaat";
        }
        
       
       


        return deletionListPerChromosome;
        
    }
    
   
    
    
    
    
    
}
