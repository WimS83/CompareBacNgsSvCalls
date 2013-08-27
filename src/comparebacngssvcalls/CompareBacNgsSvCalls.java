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
        File PiclFilteredOutput = new File("/home/sge_share_fedor13/wim/ratClusterWorkSpace/LE_ILLUMINA/LE_Illumina_c4.csv");
        

        File inBacSVCallsFile = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/13BACs_vs_10SolidRatNGS_GenotypeConcordance/13BACS_k400_w5000/13BACS_OnlyHighQualSorted_deletionsMergedWindow1000.bed");

       // File betweenBacSVCallsFile = new File("/home/wim/Analysis/ratfounder/bac/rnor5/betweenBacCalls/betweenBacCalls.txt");

        File bacRegionsFile = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/13BACs_vs_10SolidRatNGS_GenotypeConcordance/13BACS_k400_w5000/bed/13BACS_OnlyHighQualSorted_bamToBedMerged.bed");


        BufferedReader piclBR = new BufferedReader(new FileReader(PiclFilteredOutput));
        BufferedReader inBacBR = new BufferedReader(new FileReader(inBacSVCallsFile));
    //    BufferedReader betweenBacBR = new BufferedReader(new FileReader(betweenBacSVCallsFile));
        BufferedReader bacRegionsBR = new BufferedReader(new FileReader(bacRegionsFile));

        

       Map<String, ArrayList<Range<Integer>>> bacRegions =   readBacRegions(bacRegionsBR);
       Map<String, ArrayList<Range<Integer>>> piclCalls =  readPiclSVCalls(piclBR);
       Map<String, ArrayList<Range<Integer>>> inBacCalls = readInBacSVCalls(inBacBR);
       
       Map<String, ArrayList<Range<Integer>>> piclCallsInBacRegions =  getPiclCallsInBacRegions(bacRegions, piclCalls);
       
       calculatePrecision(inBacCalls, piclCallsInBacRegions);
       calculateRecall(inBacCalls, piclCallsInBacRegions);
       
       


    }

    private static  Map<String, ArrayList<Range<Integer>>> readInBacSVCalls(BufferedReader inBacBR) throws IOException {

        Map<String, ArrayList<Range<Integer>>> inBacSVCalls = new HashMap<String, ArrayList<Range<Integer>>>();

        String line = inBacBR.readLine();

        while (line != null) {
            String[] splitLine = line.split("\t");
            String chrom = splitLine[0];
            Integer begin = new Integer(splitLine[1]);
            Integer end = new Integer(splitLine[2]);
            String type = splitLine[3];     
            
            

            Range<Integer> range = Range.between(begin, end);

            if (!inBacSVCalls.containsKey(chrom)) {
                ArrayList<Range<Integer>> chromSVCallList = new ArrayList<Range<Integer>>();
                inBacSVCalls.put(chrom, chromSVCallList);
            }
            inBacSVCalls.get(chrom).add(range);


            String blaat = "blaat";
         


            line = inBacBR.readLine();

        }

        return inBacSVCalls;
    }

    private static  Map<String, ArrayList<Range<Integer>>> readBacRegions(BufferedReader bacRegionsBR) throws IOException {

        Map<String, ArrayList<Range<Integer>>> chromBacRegionsMap = new HashMap<String, ArrayList<Range<Integer>>>();

        String line = bacRegionsBR.readLine();

        while (line != null) {
            String[] splitLine = line.split("\t");
            String chrom = splitLine[0];
            Integer begin = new Integer(splitLine[1]);
            Integer end = new Integer(splitLine[2]);

            Range<Integer> range = Range.between(begin, end);

            if (!chromBacRegionsMap.containsKey(chrom)) {
                ArrayList<Range<Integer>> chromBacRegionsList = new ArrayList<Range<Integer>>();
                chromBacRegionsMap.put(chrom, chromBacRegionsList);
            }
            chromBacRegionsMap.get(chrom).add(range);

             line = bacRegionsBR.readLine();
            String blaat = "blaat";
        }


        return chromBacRegionsMap;

    }

    private static  Map<String, ArrayList<Range<Integer>>> readPiclSVCalls(BufferedReader piclBR) throws IOException {
       
        Map<String, ArrayList<Range<Integer>>> chromPiclSVCallsMap = new HashMap<String, ArrayList<Range<Integer>>>();

        String line = piclBR.readLine();

        while (line != null) {
            String[] splitLine = line.split("\t");
            String chrom = splitLine[0];
            Integer begin = new Integer(splitLine[1]);
            Integer end = new Integer(splitLine[5]);

            Range<Integer> range = Range.between(begin, end);

            if (!chromPiclSVCallsMap.containsKey(chrom)) {
                ArrayList<Range<Integer>> chromBacRegionsList = new ArrayList<Range<Integer>>();
                chromPiclSVCallsMap.put(chrom, chromBacRegionsList);
            }
            chromPiclSVCallsMap.get(chrom).add(range);

             line = piclBR.readLine();
            String blaat = "blaat";
        }


        return chromPiclSVCallsMap;
        
    }

    private static Map<String, ArrayList<Range<Integer>>> getPiclCallsInBacRegions(Map<String, ArrayList<Range<Integer>>> bacRegionsMap, Map<String, ArrayList<Range<Integer>>> piclCallsMap) {
        
        Map<String, ArrayList<Range<Integer>>> piclCallsInBacRegionsMap = new HashMap<String, ArrayList<Range<Integer>>>();
        
        for(String chrom : piclCallsMap.keySet())
        {
            //skip the picl cals if they are on a chromosome where no bac is mapped\
            if(!bacRegionsMap.containsKey(chrom)){continue;}
            
            List<Range<Integer>> bacRegions = bacRegionsMap.get(chrom);
            List<Range<Integer>> piclCalls = piclCallsMap.get(chrom);
            
            for(Range<Integer> piclCall : piclCalls)
            {
                Boolean containedInBacRegions = false;
                
                for(Range<Integer> bacRegion : bacRegions)
                {
                    if(bacRegion.containsRange(piclCall))
                    {
                        containedInBacRegions = true;
                    }
                }
                
                if(containedInBacRegions == true)
                {
                    if(!piclCallsInBacRegionsMap.containsKey(chrom))
                    {
                        ArrayList<Range<Integer>> chromPiclInBacCallsList = new ArrayList<Range<Integer>>();
                        piclCallsInBacRegionsMap.put(chrom, chromPiclInBacCallsList);
                    }
                    piclCallsInBacRegionsMap.get(chrom).add(piclCall);       
                    
                    String blaat = "blaat";
                }            
            
            }
            String blaat = "blaat";
            
        
        
        }
        
       return piclCallsInBacRegionsMap;            
    }

    private static void calculatePrecision(Map<String, ArrayList<Range<Integer>>> inBacCalls, Map<String, ArrayList<Range<Integer>>> piclCallsInBacRegions) {
       
        
       Map<String, ArrayList<Range<Integer>>> piclCallsConcordantWithBacCalls = new HashMap<String, ArrayList<Range<Integer>>>();
       Map<String, ArrayList<Range<Integer>>> piclCallsDisConcordantWithBacCalls = new HashMap<String, ArrayList<Range<Integer>>>();
       
       Integer piclSVWasMatchedCounter = 0;
       Integer piclSVWasMissedCounter = 0;
        
        for(String chrom : piclCallsInBacRegions.keySet())
        {
            //skip the picl cals if they are on a chromosome where no bac is mapped\
            if(!inBacCalls.containsKey(chrom))
            {
                piclSVWasMissedCounter = piclSVWasMissedCounter + piclCallsInBacRegions.get(chrom).size();
                continue;
            }
            
            List<Range<Integer>> bacSVCalls = inBacCalls.get(chrom);
            List<Range<Integer>> piclCalls = piclCallsInBacRegions.get(chrom);
            
            for(Range<Integer> piclCall : piclCalls)
            {
                Boolean concordantWithABacSVCall = false;
                
                for(Range<Integer> bacSVCall : bacSVCalls)
                {
                    if(bacSVCall.isOverlappedBy(piclCall))
                    {
                        concordantWithABacSVCall = true;
                    }
                }
                
                if(concordantWithABacSVCall == true)
                {
                    piclSVWasMatchedCounter++;
                    if(!piclCallsConcordantWithBacCalls.containsKey(chrom))
                    {
                        ArrayList<Range<Integer>> chromPiclConcordantList = new ArrayList<Range<Integer>>();
                        piclCallsConcordantWithBacCalls.put(chrom, chromPiclConcordantList);
                    }
                    piclCallsConcordantWithBacCalls.get(chrom).add(piclCall);   
                    
                    System.out.println(chrom+"\t"+piclCall.getMinimum()+"\t"+piclCall.getMaximum()+"\t"+"piclCall and bacBacCall !");
                    
                    String blaat = "blaat";
                }
                else
                {
                    piclSVWasMissedCounter++;
                    if(!piclCallsDisConcordantWithBacCalls.containsKey(chrom))
                    {
                        ArrayList<Range<Integer>> chromPiclDiscordantList = new ArrayList<Range<Integer>>();
                        piclCallsDisConcordantWithBacCalls.put(chrom, chromPiclDiscordantList);
                    }
                    piclCallsDisConcordantWithBacCalls.get(chrom).add(piclCall);       
                    
                    System.out.println(chrom+"\t"+piclCall.getMinimum()+"\t"+piclCall.getMaximum()+"\t"+"piclCall_noBacCall");
                    String blaat = "blaat";
                
                }
            
            }
            String blaat = "blaat";
            
        
        
        }
        
        Double precision = new Double(piclSVWasMatchedCounter) / ( new Double(piclSVWasMatchedCounter) + new Double(piclSVWasMissedCounter)) * new Double(100);
        
        
        System.out.println("Picl SVs that overlap Bac SV: "+ piclSVWasMatchedCounter);
        System.out.println("Picl SVs that do not overlap Bac SV: "+ piclSVWasMissedCounter);
        System.out.println("Precision is : "+ precision);
        
        
        
        String blaat = "blaat";        
        
    }
    
    
    private static void calculateRecall(Map<String, ArrayList<Range<Integer>>> inBacCalls, Map<String, ArrayList<Range<Integer>>> piclCallsInBacRegions) {
       
        
       Map<String, ArrayList<Range<Integer>>> BAcCallsConcordantWithPiclalls = new HashMap<String, ArrayList<Range<Integer>>>();
       Map<String, ArrayList<Range<Integer>>> baclCallsDisConcordantWithPiclCalls = new HashMap<String, ArrayList<Range<Integer>>>();
       
       Integer bacSVWasMatchedCounter = 0;
       Integer bacSVWasMissedCounter = 0;
        
        for(String chrom : inBacCalls.keySet())
        {
            //skip the picl cals if they are on a chromosome where no bac is mapped\
            if(!piclCallsInBacRegions.containsKey(chrom))
            {
                baclCallsDisConcordantWithPiclCalls.put(chrom, inBacCalls.get(chrom));
                bacSVWasMissedCounter = bacSVWasMissedCounter + inBacCalls.get(chrom).size();
                
                System.out.println("no picl call within a bac on chrom "+ chrom );
                
                continue;
            }
            
            List<Range<Integer>> bacSVCalls = inBacCalls.get(chrom);
            List<Range<Integer>> piclCalls = piclCallsInBacRegions.get(chrom);
            
            for(Range<Integer> bacSVCall : bacSVCalls)
            {
                Boolean concordantWithAPiclSVCall = false;
                
                for(Range<Integer> piclSVCall : piclCalls)
                {
                    if(bacSVCall.isOverlappedBy(piclSVCall))
                    {
                        concordantWithAPiclSVCall = true;
                    }
                }
                
                if(concordantWithAPiclSVCall == true)
                {
                    bacSVWasMatchedCounter++;
                    if(!BAcCallsConcordantWithPiclalls.containsKey(chrom))
                    {
                        ArrayList<Range<Integer>> chromBacConcordantList = new ArrayList<Range<Integer>>();
                        BAcCallsConcordantWithPiclalls.put(chrom, chromBacConcordantList);
                    }
                    BAcCallsConcordantWithPiclalls.get(chrom).add(bacSVCall);       
                    
                    String blaat = "blaat";
                }
                else
                {
                    bacSVWasMissedCounter++;
                    if(!baclCallsDisConcordantWithPiclCalls.containsKey(chrom))
                    {
                        ArrayList<Range<Integer>> chromBacDiscordantList = new ArrayList<Range<Integer>>();
                        baclCallsDisConcordantWithPiclCalls.put(chrom, chromBacDiscordantList);
                    }
                    baclCallsDisConcordantWithPiclCalls.get(chrom).add(bacSVCall);       
                    
                    //System.out.println(chrom+"\t"+bacSVCall.getMinimum()+"\t"+bacSVCall.getMaximum()+"\t"+"bacCall_noNGSCall");
                    String blaat = "blaat";
                
                }
            
            }
            String blaat = "blaat";
            
        
        
        }
        
        
        
        for(String chrom : baclCallsDisConcordantWithPiclCalls.keySet())
        {
            for(Range<Integer> bacSVCall : baclCallsDisConcordantWithPiclCalls.get(chrom))
            {
                 System.out.println(chrom+"\t"+bacSVCall.getMinimum()+"\t"+bacSVCall.getMaximum()+"\t"+"bacCall_noNGSCall");
            }
        
        }
        
        String blaat = "blaat";
        
        Double recall = new Double(bacSVWasMatchedCounter) / ( new Double(bacSVWasMatchedCounter) + new Double(bacSVWasMissedCounter)) * new Double(100);
        
        
        System.out.println("Bac SVs that overlap Picl SV: "+ bacSVWasMatchedCounter);
        System.out.println("Bac SVs that do not overlap Picl SV: "+ bacSVWasMissedCounter);
        System.out.println("Recall is : "+ recall);
       
        
        
        
    }
    
   
}
