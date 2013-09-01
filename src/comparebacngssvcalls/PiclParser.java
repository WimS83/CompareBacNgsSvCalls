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
            if(!type.equalsIgnoreCase("deletion")){continue;}               
            
            String chrom = splitLine[0];
            
            if(!chrom.contains("chr")){chrom = "chr"+chrom;}
            
            RatChromosomes ratChromosome;
            try{                
                ratChromosome = RatChromosomes.valueOf(chrom);
            }
            //continue to next line if unknown chromosome
            catch(IllegalArgumentException  ex)
            {
                continue;
            }                       
            
            Deletion deletion = new Deletion(Range.between(begin, end), chrom);            
            deletionListPerChromosome.get(ratChromosome).add(deletion);      

            line = piclBR.readLine();
            String blaat = "blaat";
        }
        
        TreeMap<Integer, Integer> deletionsCountPer100bpWindow=  storeDeletionsSizes(deletionListPerChromosome, 100, 1000000);
        printDeletionCount(deletionsCountPer100bpWindow);
       


        return deletionListPerChromosome;
        
    }
    
    private void printDeletions(EnumMap<RatChromosomes, ArrayList<Deletion>> deletionListPerChromosome) {
        
        for(RatChromosomes ratChromosomes : deletionListPerChromosome.keySet())
        {
            for(Deletion deletion : deletionListPerChromosome.get(ratChromosomes))
            {
                System.out.println(deletion.getChromosome()+ "\t"+deletion.getLocation().getMinimum()+"\t"+deletion.getLocation().getMaximum()+"\t"+"deletion"+"\t"+deletion.getSize());         
            }        
        }
        
    }
     
    private TreeMap<Integer, Integer> storeDeletionsSizes( EnumMap<RatChromosomes, ArrayList<Deletion>> deletionListPerChromosome, int binSize, int max) {
        
        TreeMap<Integer, Integer> deletionsSizesCountMap = new TreeMap<Integer,Integer>();
        
        //initialize the countMap with zero for each bin up to the max
        Integer i = 0;        
        while(i < max)
        {
            deletionsSizesCountMap.put(i, 0);
            i = i +binSize;
        }
        
        //loop over the deletions and increase the counters for the bin in which the deletions size falls for every deletion
        for(RatChromosomes ratChromosome : deletionListPerChromosome.keySet())
        {
            for(Deletion deletion : deletionListPerChromosome.get(ratChromosome))
            {
                if(deletion.getSize() > max) {continue;}

                Integer modulus = deletion.getSize() % binSize;
                Integer oldCount = deletionsSizesCountMap.get(deletion.getSize()-modulus);
                deletionsSizesCountMap.put(deletion.getSize()-modulus, oldCount+1);   
            }       
        }
        
        return deletionsSizesCountMap;       
    }
    
     private void printDeletionCount(TreeMap<Integer, Integer> deletionsSizesCountMap) {
        
        for(Integer deletionSize : deletionsSizesCountMap.keySet())
        {
            Integer count = deletionsSizesCountMap.get(deletionSize);
            System.out.println(deletionSize+"\t"+count);
        
        }
    }
    
    
    
    
    
    
}
