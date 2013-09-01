/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package comparebacngssvcalls;

import java.io.File;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.apache.commons.lang3.Range;

/**
 *
 * @author wim
 */
public class BamIndelParser {
    
    
    File bamFile;   

    public BamIndelParser(File bamFile) {
        this.bamFile = bamFile;
        
        
        
        String blaat = "";
        
        
        
    }       
     
    
    public EnumMap<RatChromosomes, ArrayList<Deletion>> parseIndelInAlignments() {
        SAMFileReader bamreader = new SAMFileReader(bamFile);
        
        //initialize map to store list of deletions per ratchromosome
        EnumMap<RatChromosomes, ArrayList<Deletion>> deletionListPerChromosome = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);
        for(RatChromosomes ratChromosomes: RatChromosomes.values())
        {
            ArrayList deletions = new ArrayList<Deletion>();
            deletionListPerChromosome.put(ratChromosomes, deletions);
        }
        
        
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
            
            
            Integer currentReferencePos = samRecord.getAlignmentStart();
            
            for (CigarElement cigarElement :samRecord.getCigar().getCigarElements())
            {
                if(cigarElement.getOperator() == CigarOperator.D)
                {
                    Integer cigarElementStart = currentReferencePos;
                    Integer cigarElementLenght = cigarElement.getLength();
                    Integer cigarElementEnd = currentReferencePos+cigarElementLenght;
                    
                    if(cigarElementLenght > 100 )
                    {
                        //System.out.println("deletion found from " + currentChrom +":"  + cigarElementStart + "-" + cigarElementEnd );
                        //System.out.println(currentChrom+ "\t"+cigarElementStart+"\t"+cigarElementEnd+"\t"+"deletion");
                        //if there are previous deletions check is this one falls within 100 bp of the previous one. If so update the of the last one to be the end of this one 
                        ArrayList<Deletion> deletionList = deletionListPerChromosome.get(ratChromosome);
                        
                        if(deletionList.size() > 0)
                        {
                            Deletion previousDeletion = deletionList.get(deletionList.size()-1 );
                            Integer previousDeletionEnd = previousDeletion.getLocation().getMaximum();
                            
                            
                            if(previousDeletion.getChromosome().equals(currentChrom) &&  cigarElementStart - previousDeletionEnd  < 1000   )
                            {
                                previousDeletion.updateLocationEnd(cigarElementEnd);
                            }
                            else
                            {
                                 Deletion deletion = new Deletion(Range.between(cigarElementStart, cigarElementEnd), currentChrom);
                                deletionList.add(deletion);
                            }
                        }
                        else
                        {
                             Deletion deletion = new Deletion(Range.between(cigarElementStart, cigarElementEnd), currentChrom);
                             deletionList.add(deletion);
                        }                        
                        
                    }  
                }
                if(cigarElement.getOperator() == CigarOperator.I)
                {
                    Integer cigarElementStart = currentReferencePos;
                    Integer cigarElementLenght = cigarElement.getLength();
                   
                    
                    if(cigarElementLenght > 100 )
                    {
                        //System.out.println("insertion found at " + currentChrom +":"  + cigarElementStart );
                        //System.out.println(currentChrom+ "\t"+(cigarElementStart-1)+"\t"+cigarElementStart+"\t"+"insertion");
                    }  
                }  
                
                
                //add the lenght of this cigarElement to the current pos to go to the start of the next element
                if(cigarElement.getOperator().consumesReferenceBases())
                {
                     currentReferencePos = currentReferencePos + cigarElement.getLength();
                } 
            }            
        }   
        
        
        printDeletions(deletionListPerChromosome);
        TreeMap<Integer, Integer> deletionCount100bpWindows = storeDeletionsSizes(deletionListPerChromosome, 100, 6000);        
        printDeletionCount(deletionCount100bpWindows);       
       
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
     
     
     
     
    

    private static void parseIndelsBetweenAlignments(File inputBam) {
        
        
        SAMFileReader bamreader = new SAMFileReader(inputBam);
        
        //put the alignments in a map entry per bac
        Map<String, ArrayList<SAMRecord>> mapContigAlignments = new HashMap<String, ArrayList<SAMRecord>>();        
        for(SAMRecord samRecord : bamreader)
        {
            String readName  = samRecord.getReadName();
            
            if(!mapContigAlignments.containsKey(readName))
            {
                 ArrayList<SAMRecord> alignmentList = new ArrayList<SAMRecord>();
                 mapContigAlignments.put(readName, alignmentList);
            }
            
            mapContigAlignments.get(readName).add(samRecord);          
            
            
            String blaat = "blaat";
            
        
        }
        
        
        for(String bac : mapContigAlignments.keySet())
        {
            ArrayList<SAMRecord> records = mapContigAlignments.get(bac);
            Integer numberOfAlignments = records.size();
            System.out.println();
            System.out.println(bac +"  has " + numberOfAlignments + " alignments" );
            
            for(SAMRecord sAMRecord : records)
            {
                String chromName = sAMRecord.getReferenceName();
                Integer start = sAMRecord.getAlignmentStart();
                Integer stop = sAMRecord.getAlignmentEnd();
                
                System.out.println(chromName+":"+start+"-"+stop);  
            }
            
            
            TreeMap<Integer, SAMRecord> mapStartAlignments = new TreeMap<Integer, SAMRecord>();   
            
            for(SAMRecord sAMRecord : mapContigAlignments.get(bac))
            {
                mapStartAlignments.put(sAMRecord.getAlignmentStart(), sAMRecord);                
            }
            
            SAMRecord previousRecord = null;
            
            for (Integer key : mapStartAlignments.keySet()) { 
                 
                 SAMRecord currentRecord = mapStartAlignments.get(key);
                 
                 if(previousRecord != null)
                 {
                     Integer previousRecordEnd = previousRecord.getAlignmentEnd();
                     Integer currentRecordStart = currentRecord.getAlignmentStart();
                     System.out.println("found deletion on chrom "+ previousRecord.getReferenceName()+ ":" + previousRecordEnd+"-"+currentRecordStart);
                 
                 }
                 
                 previousRecord = currentRecord;
                 
                 
                 
                 
            // do something
            }
            
            
        }
        
        
        
        
        
        String blaat = "blaat";
    }

   

    

   
    
    
}
