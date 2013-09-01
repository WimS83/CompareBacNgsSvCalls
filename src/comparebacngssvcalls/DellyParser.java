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
import org.apache.commons.lang3.Range;

/**
 *
 * @author wim
 */
public class DellyParser {
    
    File dellyFile;
    BufferedReader dellyBR;

    public DellyParser(File dellyFile) throws FileNotFoundException {
        this.dellyFile = dellyFile;
         dellyBR = new BufferedReader(new FileReader(dellyFile));    
    }
    
      public  EnumMap<RatChromosomes, ArrayList<Deletion>> readDellySVCalls() throws IOException {
       
       // Map<String, ArrayList<Deletion>> chromPiclSVCallsMap = new HashMap<String, ArrayList<Deletion>>();
        
        //initialize map to store list of deletions per ratchromosome
        EnumMap<RatChromosomes, ArrayList<Deletion>> deletionListPerChromosome = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);
        for(RatChromosomes ratChromosomes: RatChromosomes.values())
        {
            ArrayList deletions = new ArrayList<Deletion>();
            deletionListPerChromosome.put(ratChromosomes, deletions);
        }
        
        

        String line = dellyBR.readLine();

        while (line != null) {
            String[] splitLine = line.split("\t");
            Integer begin = new Integer(splitLine[1]);
            Integer end = new Integer(splitLine[2]);
//            String type = splitLine[9];            
//            if(!type.equalsIgnoreCase("deletion")){continue;}       
            
            Integer supportingPairs = new Integer(splitLine[4]);
            if(supportingPairs < 4 ){
                line = dellyBR.readLine();
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
                line = dellyBR.readLine();
                continue;
            }                       
            
            Deletion deletion = new Deletion(Range.between(begin, end), chrom);            
            deletionListPerChromosome.get(ratChromosome).add(deletion);      

            line = dellyBR.readLine();
            String blaat = "blaat";
        }
        
       
       


        return deletionListPerChromosome;
        
    }
    
    
    
}
