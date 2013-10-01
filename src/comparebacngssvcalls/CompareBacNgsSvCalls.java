/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package comparebacngssvcalls;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.TreeMap;
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

        File PiclFilteredOutput = new File("/home/sge_share_fedor13/wim/ratClusterWorkSpace/LE/LE_c4_deletion.csv");
        //File PiclFilteredOutput = new File("/home/wim/Analysis/ratfounder/NGS/Picl_SV_calls/LE/LE_deletion.csv");

        //File PiclFilteredOutput = new File("/home/wim/Analysis/ratfounder/NGS/Picl_SV_calls/LE/Illumina/LE_ILLUMINA/LE_Illumina_c4.csv");   


        File inBacSVCallsFile = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/13BACs_vs_10SolidRatNGS_GenotypeConcordance/13BACS_k400_w5000/13BACS_OnlyHighQualSorted_deletionsMergedWindow1000.bed");

        // File betweenBacSVCallsFile = new File("/home/wim/Analysis/ratfounder/bac/rnor5/betweenBacCalls/betweenBacCalls.txt");

        File bacRegionsFile = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/13BACs_vs_10SolidRatNGS_GenotypeConcordance/13BACS_k400_w5000/bed/13BACS_OnlyHighQualSorted_bamToBedMerged.bed");

        File BacContigsBamFile = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/13BACs_vs_10SolidRatNGS_GenotypeConcordance/13BACS_k400_w5000/13BACS_OnlyHighQualSorted.bam");
        //File BacContigsBamFile = new File("/home/wim/Analysis/ratfounder/bac/rnor5/13BACS_OnlyHighQualSorted.bam");

        File dellyDeletionFile = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/dellyCalls/LE/minQ30/del_onlyCalls.txt");

        File breakDancerFile = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/breakdancer_calls/LE_try2/LE_calls.csv");

        File cnvNatorFile = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/LE/LE_calls.out");




        //parse the bacRegions
        BamAlignmentRegionsParser bamAlignmentRegionsParser = new BamAlignmentRegionsParser(BacContigsBamFile);
        EnumMap<RatChromosomes, ArrayList<Range<Integer>>> bacAlignments = bamAlignmentRegionsParser.readBacRegions();


        //parse the deletions from the bac alignment
        BamIndelParser bamIndelParser = new BamIndelParser(BacContigsBamFile);
        EnumMap<RatChromosomes, ArrayList<Deletion>> bacBasedDeletionCalls = bamIndelParser.parseIndelInAlignments(100);

        //parse the deletions from the Picl calls
        PiclParser piclParser = new PiclParser(PiclFilteredOutput);
        EnumMap<RatChromosomes, ArrayList<Deletion>> piclBasedDeletionsCalls = piclParser.readPiclSVCalls();
        EnumMap<RatChromosomes, ArrayList<Deletion>> piclBasedDeletionsCallsInBacRegions = getCallsInBacRegions(bacAlignments, piclBasedDeletionsCalls);

        //parse the deletions from the Delly calls
        DellyParser dellyParser = new DellyParser(dellyDeletionFile);
        EnumMap<RatChromosomes, ArrayList<Deletion>> dellyBasedDeletionsCalls = dellyParser.readDellySVCalls();
        EnumMap<RatChromosomes, ArrayList<Deletion>> dellyBasedDeletionsCallsInBacRegions = getCallsInBacRegions(bacAlignments, dellyBasedDeletionsCalls);

        //parse the deletions from the cnvnator calls
        CNVnatorParser cNVnatorParser = new CNVnatorParser(cnvNatorFile);
        EnumMap<RatChromosomes, ArrayList<Deletion>> cnvNatorBasedDeletionsCalls = cNVnatorParser.readCNVSVCalls();
        EnumMap<RatChromosomes, ArrayList<Deletion>> cnvNatorBasedDeletionsCallsInBacRegions = getCallsInBacRegions(bacAlignments, cnvNatorBasedDeletionsCalls);

        //parse the deletions from the Delly calls
        BreakdancerParser breakDancerParser = new BreakdancerParser(breakDancerFile);
        EnumMap<RatChromosomes, ArrayList<Deletion>> breakdancerBasedDeletionsCalls = breakDancerParser.readBreakdancerSVCalls();
        EnumMap<RatChromosomes, ArrayList<Deletion>> breakdancerBasedDeletionsCallsInBacRegions = getCallsInBacRegions(bacAlignments, breakdancerBasedDeletionsCalls);


//        calculateOverlap("breakdancerCalls", breakdancerBasedDeletionsCallsInBacRegions,"bacCalls", bacBasedDeletionCalls);
//        
//        calculateOverlap("bacCalls", bacBasedDeletionCalls, "breakdancerCalls", breakdancerBasedDeletionsCallsInBacRegions);


       // calculateOverlap("PiclCalls", piclBasedDeletionsCallsInBacRegions,"bacCalls", bacBasedDeletionCalls);
        
//       calculateOverlap("bacCalls", bacBasedDeletionCalls, "PiclCalls", piclBasedDeletionsCallsInBacRegions);
//        

//        calculateOverlap("DellyCalls", dellyBasedDeletionsCallsInBacRegions,"bacCalls", bacBasedDeletionCalls);
//        
//        calculateOverlap("bacCalls", bacBasedDeletionCalls, "DellyCalls", dellyBasedDeletionsCallsInBacRegions);

//          calculateOverlap("breakdancerCalls", breakdancerBasedDeletionsCallsInBacRegions,"piclCalls", piclBasedDeletionsCallsInBacRegions);
//        
//          calculateOverlap("piclCalls", piclBasedDeletionsCallsInBacRegions, "breakdancerCalls", breakdancerBasedDeletionsCallsInBacRegions);

        calculateOverlap("breakdancerCalls", breakdancerBasedDeletionsCalls, "piclCalls", piclBasedDeletionsCalls);

        calculateOverlap("piclCalls", piclBasedDeletionsCalls, "breakdancerCalls", breakdancerBasedDeletionsCalls);
        
//        calculateOverlapDifference3Sets("bacBasedDeletionCalls", bacBasedDeletionCalls, "piclCalls", piclBasedDeletionsCallsInBacRegions, "breakDancerCalls", breakdancerBasedDeletionsCallsInBacRegions);




//        //print the bacBased deletionCount in 100bp windows
//        TreeMap<Integer, Integer> bacDeletionCount100bpWindows = storeDeletionsSizes(bacBasedDeletionCalls, 100, 6000);        
//        printDeletionCount(bacDeletionCount100bpWindows);       
//        
//        //print the picl based deletionCount in 100bp windows 
//        TreeMap<Integer, Integer> piclDeletionsCountPer100bpWindow=  storeDeletionsSizes(piclBasedDeletionsCalls, 100, 1000000);
//        printDeletionCount(piclDeletionsCountPer100bpWindow);



//       
//       calculatePrecision(inBacCalls, piclCallsInBacRegions);
//       calculateRecall(inBacCalls, piclCallsInBacRegions);      

    }

    private static EnumMap<RatChromosomes, ArrayList<Deletion>> getCallsInBacRegions(EnumMap<RatChromosomes, ArrayList<Range<Integer>>> bacAlignments, EnumMap<RatChromosomes, ArrayList<Deletion>> piclCallsMap) {

        EnumMap<RatChromosomes, ArrayList<Deletion>> deletionInBacAligmentListPerChromosome = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);
        for (RatChromosomes ratChromosomes : RatChromosomes.values()) {
            ArrayList deletions = new ArrayList<Deletion>();
            deletionInBacAligmentListPerChromosome.put(ratChromosomes, deletions);
        }

        for (RatChromosomes ratChromosome : RatChromosomes.values()) {
            List<Deletion> piclDeletions = piclCallsMap.get(ratChromosome);
            List<Range<Integer>> bacAlignmentsList = bacAlignments.get(ratChromosome);

            for (Deletion deletion : piclDeletions) {
                Boolean containedInBacAlignment = false;

                for (Range<Integer> bacAlignment : bacAlignmentsList) {
                    if (bacAlignment.containsRange(deletion.getLocation())) {
                        containedInBacAlignment = true;
                    }

                }

                if (containedInBacAlignment) {
                    deletionInBacAligmentListPerChromosome.get(ratChromosome).add(deletion);
                }
            }
        }

        return deletionInBacAligmentListPerChromosome;

    }

    private static void calculateOverlap(String nameSetA, EnumMap<RatChromosomes, ArrayList<Deletion>> deletionSetA, String nameSetB, EnumMap<RatChromosomes, ArrayList<Deletion>> deletionSetB) {

        //initialize the maps to store the concordant and discordant call
        EnumMap<RatChromosomes, ArrayList<Deletion>> concordantWithSetB = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);
        for (RatChromosomes ratChromosomes : RatChromosomes.values()) {
            ArrayList<Deletion> deletions = new ArrayList<Deletion>();
            concordantWithSetB.put(ratChromosomes, deletions);
        }

        EnumMap<RatChromosomes, ArrayList<Deletion>> discordantWithSetB = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);
        for (RatChromosomes ratChromosomes : RatChromosomes.values()) {
            ArrayList<Deletion> deletions = new ArrayList<Deletion>();
            discordantWithSetB.put(ratChromosomes, deletions);
        }

        //initialize the concordant and discordant counters
        Integer concordantWithBCounter = 0;
        Integer discordantWithBCounter = 0;

        for (RatChromosomes ratChromosome : RatChromosomes.values()) {
            List<Deletion> deletionListSetA = deletionSetA.get(ratChromosome);
            List<Deletion> deletionListSetB = deletionSetB.get(ratChromosome);

            for (Deletion deletionA : deletionListSetA) {
                Boolean overLap = false;
                for (Deletion deletionB : deletionListSetB) {
//                    if (deletionA.getLocation().isOverlappedBy(deletionB.getLocation())) {
//                        overLap = true;
//                    }
                    
                    if (deletionA.getLocation().isOverlappedBy(deletionB.getLocation())) 
                    {
                        Range<Integer> intersection = deletionA.getLocation().intersectionWith(deletionB.getLocation());
                        Double overlapLenght = new Double(intersection.getMaximum() - intersection.getMinimum());
                        Double deletionALenght = new Double(deletionA.getLocation().getMaximum() - deletionA.getLocation().getMinimum());
                        Double deletionBLenght = new Double(deletionB.getLocation().getMaximum() - deletionB.getLocation().getMinimum());
                        Double overlapFractionOfA = overlapLenght / deletionALenght;
                        Double overlapFractionOfB = overlapLenght / deletionBLenght;
                        Double cutOff = new Double("0.5");
                        
                        if(overlapFractionOfA > cutOff && overlapFractionOfB > cutOff)
                        //if(overlapFractionOfA > cutOff)
                        {
                            overLap = true;
                        }
                    }
                    

                }
                if (overLap == true) {
                    concordantWithBCounter++;
                    concordantWithSetB.get(ratChromosome).add(deletionA);

                } else {
                    discordantWithBCounter++;
                    discordantWithSetB.get(ratChromosome).add(deletionA);
                }

            }
        }

        System.out.println(concordantWithBCounter + " from set " + nameSetA + " concordant with " + nameSetB);
        System.out.println(discordantWithBCounter + " from set " + nameSetA + " disconcordant with " + nameSetB);

        Double overlapPercentage = new Double(concordantWithBCounter) / (new Double(concordantWithBCounter) + new Double(discordantWithBCounter)) * new Double(100);

        System.out.println(overlapPercentage + " % overlap from set " + nameSetA + "  with " + nameSetB);

     //  printDeletions(concordantWithSetB, "concordantWith "+nameSetB);
       //printDeletions(discordantWithSetB, "disconcordantWith "+nameSetB);




    }

    private static void calculateOverlapDifference3Sets(String nameSetA, EnumMap<RatChromosomes, ArrayList<Deletion>> deletionSetA, String nameSetB, EnumMap<RatChromosomes, ArrayList<Deletion>> deletionSetB, String nameSetC, EnumMap<RatChromosomes, ArrayList<Deletion>> deletionSetC) {

        //initialize the maps to store the concordant and discordant call
        EnumMap<RatChromosomes, ArrayList<Deletion>> calledInBothBandCandCordandantWithA = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);
        for (RatChromosomes ratChromosomes : RatChromosomes.values()) {
            ArrayList<Deletion> deletions = new ArrayList<Deletion>();
            calledInBothBandCandCordandantWithA.put(ratChromosomes, deletions);
        }

        EnumMap<RatChromosomes, ArrayList<Deletion>> calledInBorCandConcordantWithA = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);
        for (RatChromosomes ratChromosomes : RatChromosomes.values()) {
            ArrayList<Deletion> deletions = new ArrayList<Deletion>();
            calledInBorCandConcordantWithA.put(ratChromosomes, deletions);
        }



        //initialize the concordant and discordant counters
        Integer bothConconcordantCounter = 0;
        Integer oneConcordantCoutner = 0;
        Integer noneConcordantCoutner = 0;

        for (RatChromosomes ratChromosome : RatChromosomes.values()) {
            List<Deletion> deletionListSetA = deletionSetA.get(ratChromosome);
            List<Deletion> deletionListSetB = deletionSetB.get(ratChromosome);
            List<Deletion> deletionListSetC = deletionSetC.get(ratChromosome);

           





            for (Deletion deletionA : deletionListSetA) 
            {
                 Boolean overLapBwithA = false;
                 Boolean overLapCwithA = false;
                
                for (Deletion deletionB : deletionListSetB) 
                {
                    if (deletionA.getLocation().isOverlappedBy(deletionB.getLocation())) 
                    {
                        overLapBwithA = true;
                    }

                }

                for (Deletion deletionC : deletionListSetC) 
                {
                    if (deletionC.getLocation().isOverlappedBy(deletionA.getLocation())) 
                    {
                        overLapCwithA = true;
                    }

                }

                if (overLapBwithA && overLapCwithA) 
                {
                    bothConconcordantCounter++;
                    calledInBothBandCandCordandantWithA.get(ratChromosome).add(deletionA);

                } else {
                    if (overLapBwithA || overLapCwithA) 
                    {
                        oneConcordantCoutner++;
                        calledInBorCandConcordantWithA.get(ratChromosome).add(deletionA);
                    }
                    else
                    {
                        noneConcordantCoutner++;
                    
                    }                   
                }
            }
            

        }


        System.out.println(bothConconcordantCounter  + " from set " + nameSetB + " and set "+ nameSetC + "concordant with " + nameSetA);
        System.out.println(oneConcordantCoutner  + " from set " + nameSetB + " or set "+ nameSetC + "concordant with " + nameSetA);
        System.out.println(noneConcordantCoutner  + " from set " + nameSetB + " and set "+ nameSetC + "disconcordant with " + nameSetA);
        
        //printDeletions(calledInBothBandCandCordandantWithA, "Picl and BreakDancer concordant with "+nameSetA);
        printDeletions(calledInBorCandConcordantWithA, "Picl or BreakDancer concordant with "+nameSetA);
        
        

//        Double overlapPercentage = new Double(concordantWithBCounter) / (new Double(concordantWithBCounter) + new Double(discordantWithBCounter)) * new Double(100);
//
//        System.out.println(overlapPercentage
//                + " % overlap from set " + nameSetA + "  with " + nameSetB);
//
//        printDeletions(calledInBothBandCandCordandantWithA,
//                "concordantWith " + nameSetB);
//        printDeletions(calledInBorCandConcordantWithA,
//                "disconcordantWith " + nameSetB);




    }

    private static void printDeletions(EnumMap<RatChromosomes, ArrayList<Deletion>> deletionListPerChromosome, String setName) {

        for (RatChromosomes ratChromosomes : deletionListPerChromosome.keySet()) {
            for (Deletion deletion : deletionListPerChromosome.get(ratChromosomes)) {
                System.out.println(deletion.getChromosome() + "\t" + deletion.getLocation().getMinimum() + "\t" + deletion.getLocation().getMaximum() + "\t" + "deletion" + "\t" + deletion.getSize() + "\t" + setName);
            }
        }

    }

    private static TreeMap<Integer, Integer> storeDeletionsSizes(EnumMap<RatChromosomes, ArrayList<Deletion>> deletionListPerChromosome, int binSize, int max) {

        TreeMap<Integer, Integer> deletionsSizesCountMap = new TreeMap<Integer, Integer>();

        //initialize the countMap with zero for each bin up to the max
        Integer i = 0;
        while (i < max) {
            deletionsSizesCountMap.put(i, 0);
            i = i + binSize;
        }

        //loop over the deletions and increase the counters for the bin in which the deletions size falls for every deletion
        for (RatChromosomes ratChromosome : deletionListPerChromosome.keySet()) {
            for (Deletion deletion : deletionListPerChromosome.get(ratChromosome)) {
                if (deletion.getSize() > max) {
                    continue;
                }

                Integer modulus = deletion.getSize() % binSize;
                Integer oldCount = deletionsSizesCountMap.get(deletion.getSize() - modulus);
                deletionsSizesCountMap.put(deletion.getSize() - modulus, oldCount + 1);
            }
        }

        return deletionsSizesCountMap;
    }

    private static void printDeletionCount(TreeMap<Integer, Integer> deletionsSizesCountMap) {

        for (Integer deletionSize : deletionsSizesCountMap.keySet()) {
            Integer count = deletionsSizesCountMap.get(deletionSize);
            System.out.println(deletionSize + "\t" + count);

        }
    }
}
