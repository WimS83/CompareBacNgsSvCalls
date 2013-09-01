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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.lang3.Range;

/**
 *
 * @author wim
 */
public class CompareBacNgsSvCalls {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException, IOException {
        // TODO code application logic here

        //File PiclFilteredOutput = new File("/home/sge_share_fedor13/wim/ratClusterWorkSpace/LE/LE_c4_deletion.csv");
        File PiclFilteredOutput = new File("/home/wim/Analysis/ratfounder/NGS/Picl_SV_calls/LE/LE_c4_deletion.csv");
        
        //File PiclFilteredOutput = new File("/home/wim/Analysis/ratfounder/NGS/Picl_SV_calls/LE/Illumina/LE_ILLUMINA/LE_Illumina_c4.csv");   
        

        File inBacSVCallsFile = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/13BACs_vs_10SolidRatNGS_GenotypeConcordance/13BACS_k400_w5000/13BACS_OnlyHighQualSorted_deletionsMergedWindow1000.bed");

       // File betweenBacSVCallsFile = new File("/home/wim/Analysis/ratfounder/bac/rnor5/betweenBacCalls/betweenBacCalls.txt");

        File bacRegionsFile = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/13BACs_vs_10SolidRatNGS_GenotypeConcordance/13BACS_k400_w5000/bed/13BACS_OnlyHighQualSorted_bamToBedMerged.bed");

        //File BacContigsBamFile = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/13BACs_vs_10SolidRatNGS_GenotypeConcordance/13BACS_k400_w5000/13BACS_OnlyHighQualSorted.bam"); 
        File BacContigsBamFile = new File("/home/wim/Analysis/ratfounder/bac/rnor5/13BACS_OnlyHighQualSorted.bam");
        
        
        
        //parse the deletions from the bac alignment
        BamIndelParser bamIndelParser = new BamIndelParser(BacContigsBamFile);        
        EnumMap<RatChromosomes, ArrayList<Deletion>> bacBasedDeletionCalls = bamIndelParser.parseIndelInAlignments();
        
        //parse the deletions from the Picl calls
        PiclParser piclParser = new PiclParser(PiclFilteredOutput );
        EnumMap<RatChromosomes, ArrayList<Deletion>> piclBasedDeletionsCalls = piclParser.readPiclSVCalls();
        
        //parse the bacRegions
        BamAlignmentRegionsParser bamAlignmentRegionsParser = new BamAlignmentRegionsParser(BacContigsBamFile);
        EnumMap<RatChromosomes, ArrayList<Range<Integer>>> bacAlignments = bamAlignmentRegionsParser.readBacRegions();
       
        EnumMap<RatChromosomes, ArrayList<Deletion>> piclBasedDeletionsCallsInBacRegions =  getPiclCallsInBacRegions(bacAlignments, piclBasedDeletionsCalls);
        
        
        calculateOverlap("PiclCalls", piclBasedDeletionsCallsInBacRegions,"bacCalls", bacBasedDeletionCalls);
        
        calculateOverlap("bacCalls", bacBasedDeletionCalls, "PiclCalls", piclBasedDeletionsCallsInBacRegions);
        
        
//       
//       calculatePrecision(inBacCalls, piclCallsInBacRegions);
//       calculateRecall(inBacCalls, piclCallsInBacRegions);      

    }
    

    

    private static EnumMap<RatChromosomes, ArrayList<Deletion>> getPiclCallsInBacRegions(EnumMap<RatChromosomes, ArrayList<Range<Integer>>> bacAlignments , EnumMap<RatChromosomes, ArrayList<Deletion>> piclCallsMap) {
        
       EnumMap<RatChromosomes, ArrayList<Deletion>> deletionInBacAligmentListPerChromosome = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);
        for(RatChromosomes ratChromosomes: RatChromosomes.values())
        {
            ArrayList deletions = new ArrayList<Deletion>();
            deletionInBacAligmentListPerChromosome.put(ratChromosomes, deletions);
        }
        
        for(RatChromosomes ratChromosome: RatChromosomes.values())
        {
            List<Deletion> piclDeletions = piclCallsMap.get(ratChromosome);
            List<Range<Integer>> bacAlignmentsList = bacAlignments.get(ratChromosome);
            
            for(Deletion deletion: piclDeletions)
            {
                Boolean containedInBacAlignment = false;
                
                for(Range<Integer> bacAlignment: bacAlignmentsList)
                {
                    if(bacAlignment.containsRange(deletion.getLocation()))
                    {
                        containedInBacAlignment = true;
                    }
                
                }
                
                if(containedInBacAlignment)
                {
                    deletionInBacAligmentListPerChromosome.get(ratChromosome).add(deletion);
                }               
            }           
        }   
        
       return deletionInBacAligmentListPerChromosome;
       
    }
    
    private static void calculateOverlap(String nameSetA, EnumMap<RatChromosomes, ArrayList<Deletion>> deletionSetA, String nameSetB, EnumMap<RatChromosomes, ArrayList<Deletion>> deletionSetB)
    {
        
       //initialize the maps to store the concordant and discordant call
       EnumMap<RatChromosomes, ArrayList<Deletion>> concordantWithSetB = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);       
       for(RatChromosomes ratChromosomes: RatChromosomes.values())
       {
           ArrayList<Deletion> deletions = new ArrayList<Deletion>();
           concordantWithSetB.put(ratChromosomes, deletions);
       } 
       
       EnumMap<RatChromosomes, ArrayList<Deletion>> discordantWithSetB = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);
       for(RatChromosomes ratChromosomes: RatChromosomes.values())
       {
           ArrayList<Deletion> deletions = new ArrayList<Deletion>();
           discordantWithSetB.put(ratChromosomes, deletions);
       }          
       
       //initialize the concordant and discordant counters
       Integer concordantWithBCounter = 0;
       Integer discordantWithBCounter = 0;
       
       for(RatChromosomes ratChromosome : RatChromosomes.values())
       {
           List<Deletion> deletionListSetA = deletionSetA.get(ratChromosome);
           List<Deletion> deletionListSetB = deletionSetB.get(ratChromosome);
           
           for(Deletion deletionA : deletionListSetA)
           {
               Boolean overLap = false;
               for(Deletion deletionB: deletionListSetB )
               {
                   if(deletionA.getLocation().isOverlappedBy(deletionB.getLocation()))
                   {
                       overLap = true;
                   }
                   
               }
               if(overLap == true)
               {
                   concordantWithBCounter++;
                   concordantWithSetB.get(ratChromosome).add(deletionA);
               
               }
               else
               {
                   discordantWithBCounter++;
                   discordantWithSetB.get(ratChromosome).add(deletionA);
               }
           
           }   
       }
       
       System.out.println(concordantWithBCounter +" from set "+nameSetA +" concordant with "+nameSetB);
       System.out.println(discordantWithBCounter +" from set "+nameSetA +" disconcordant with "+nameSetB);
       
       Double overlapPercentage = new Double(concordantWithBCounter) / ( new Double(concordantWithBCounter) + new Double(discordantWithBCounter)) * new Double(100);
       
       System.out.println(overlapPercentage +" % overlap from set "+nameSetA +"  with "+nameSetB);
       
        
        
        
        
    
    }       
    
    
    

}
