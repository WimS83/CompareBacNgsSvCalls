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
public class CNVnatorParser {
    
    File cnvNatorFile;
    BufferedReader cnvBR;
        
    public CNVnatorParser(File cnvNatorFile) throws FileNotFoundException {
        this.cnvNatorFile = cnvNatorFile;
        
        cnvBR = new BufferedReader(new FileReader(cnvNatorFile));    
        
    }
    
    public  EnumMap<RatChromosomes, ArrayList<Deletion>> readCNVSVCalls() throws IOException {
       
       // Map<String, ArrayList<Deletion>> chromPiclSVCallsMap = new HashMap<String, ArrayList<Deletion>>();
        
        //initialize map to store list of deletions per ratchromosome
        EnumMap<RatChromosomes, ArrayList<Deletion>> deletionListPerChromosome = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);
        for(RatChromosomes ratChromosomes: RatChromosomes.values())
        {
            ArrayList deletions = new ArrayList<Deletion>();
            deletionListPerChromosome.put(ratChromosomes, deletions);
        }
        
        

        String line = cnvBR.readLine();

        while (line != null) {
            String[] splitLine = line.split("\t");
            String location = splitLine[1];
            String[] locationSplit = location.split(":");
            String locationBeginEnd = locationSplit[1];
            String[] locationBeginEndSplit = locationBeginEnd.split("-");
            
            Integer begin = new Integer(locationBeginEndSplit[0]);
            Integer end = new Integer(locationBeginEndSplit[1]);
            String type = splitLine[0];            
            if(!type.equalsIgnoreCase("deletion")){
                line = cnvBR.readLine();
                continue;
            }   
            
//            Double mapq0Fraction  = new Double(splitLine[8]);
//            if(mapq0Fraction.compareTo(new Double(1)) ==0  )
//            {
//                line = cnvBR.readLine();
//                continue;
//            }
            
            
            
//            Integer readPairSupport = new Integer(splitLine[7]);
//            if(readPairSupport < 4 )
//            {
//                 line = cnvBR.readLine();
//                continue;
//                
//            }
            
            
            String chrom = locationSplit[0];
            
            if(!chrom.contains("chr")){chrom = "chr"+chrom;}
            
            RatChromosomes ratChromosome;
            try{                
                ratChromosome = RatChromosomes.valueOf(chrom);
            }
            //continue to next line if unknown chromosome
            catch(IllegalArgumentException  ex)
            {
                line = cnvBR.readLine();
                continue;
            }                       
            
            Deletion deletion = new Deletion(Range.between(begin, end), chrom);            
            deletionListPerChromosome.get(ratChromosome).add(deletion);      

            line = cnvBR.readLine();
            String blaat = "blaat";
        }
        
       
       


        return deletionListPerChromosome;
        
    }
    
    
    
    
}
