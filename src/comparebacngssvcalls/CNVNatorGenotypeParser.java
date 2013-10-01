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
public class CNVNatorGenotypeParser {
    
    File cnvNatorGenotypeFile;
    BufferedReader cnvGenotypeBR;
    
    String ratStrain;
        
    public CNVNatorGenotypeParser(File cnvNatorGenotypeFile, String ratStrain) throws FileNotFoundException {
        this.cnvNatorGenotypeFile = cnvNatorGenotypeFile;
        this.ratStrain = ratStrain;
        
        cnvGenotypeBR = new BufferedReader(new FileReader(cnvNatorGenotypeFile));    
        
    }
    
    public  EnumMap<RatChromosomes, ArrayList<Deletion>> readCNVSGenotypes() throws IOException {
       
       // Map<String, ArrayList<Deletion>> chromPiclSVCallsMap = new HashMap<String, ArrayList<Deletion>>();
        
        //initialize map to store list of deletions per ratchromosome
        EnumMap<RatChromosomes, ArrayList<Deletion>> deletionListPerChromosome = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);
        for(RatChromosomes ratChromosomes: RatChromosomes.values())
        {
            ArrayList deletions = new ArrayList<Deletion>();
            deletionListPerChromosome.put(ratChromosomes, deletions);
        }
        
        

        String line = cnvGenotypeBR.readLine();

        while (line != null) {
            
            
            if(line.contains("Assuming"))
            {
                line = cnvGenotypeBR.readLine();
                continue;
            }
            
            
            
            
            String[] splitLine = line.split(" ");
            String location = splitLine[1];
            String[] locationSplit = location.split(":");
            String locationBeginEnd = locationSplit[1];
            String[] locationBeginEndSplit = locationBeginEnd.split("-");
            
            Integer begin = new Integer(locationBeginEndSplit[0]);
            Integer end = new Integer(locationBeginEndSplit[1]);
          
            Double genotype = new Double(splitLine[3]);
            
            
            
//            Double mapq0Fraction  = new Double(splitLine[8]);
//            if(mapq0Fraction.compareTo(new Double(1)) ==0  )
//            {
//                line = cnvGenotypeBR.readLine();
//                continue;
//            }
            
            
            
//            Integer readPairSupport = new Integer(splitLine[7]);
//            if(readPairSupport < 4 )
//            {
//                 line = cnvGenotypeBR.readLine();
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
                line = cnvGenotypeBR.readLine();
                continue;
            }                       
            
            Deletion deletion = new Deletion(Range.between(begin, end), chrom);
            deletion.addRatGenotype(ratStrain, genotype);
            deletionListPerChromosome.get(ratChromosome).add(deletion);      

            line = cnvGenotypeBR.readLine();
            String blaat = "blaat";
        }
        
       
       


        return deletionListPerChromosome;
        
    }
    
    
}
