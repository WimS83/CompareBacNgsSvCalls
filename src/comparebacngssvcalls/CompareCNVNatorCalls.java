/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package comparebacngssvcalls;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;

/**
 *
 * @author wim
 */
public class CompareCNVNatorCalls {

    public CompareCNVNatorCalls() throws FileNotFoundException, IOException {        
        
       
        File cnvNatorLE = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/LE/LE_calls.out");
        File cnvNatorACI = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/ACI/ACI_calls.out");
        File cnvNatorBNLX = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/BNLX/BNLX_calls.out");
        File cnvNatorBNSSN = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/BNSSN/BNSSN_calls.out");
        File cnvNatorBUF = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/BUF/BUF_calls.out");
        File cnvNatorF344 = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/F334/F334_calls.out");
        File cnvNatorM520 = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/M520/M520_calls.out");
        File cnvNatorMR = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/MR/MR_calls.out");
        File cnvNatorSHR = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/SHR/SHR_calls.out");
        File cnvNatorWKY = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/WKY/WKY_calls.out");
        File cnvNatorWN = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/WN/WN_calls.out");
        
        File nonOverlappingCalls = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/nonOverlappingCalls/nonOverlappingCalls.txt");
        
        FileWriter fw = new FileWriter(nonOverlappingCalls);
	BufferedWriter bw = new BufferedWriter(fw);        
        
        //map to store all non overlapping calls
        EnumMap<RatChromosomes, ArrayList<Deletion>> totalNonOverlappingDeletions = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);
        for(RatChromosomes ratChromosomes: RatChromosomes.values())
        {
            ArrayList deletions = new ArrayList<Deletion>();
            totalNonOverlappingDeletions.put(ratChromosomes, deletions);
        }        
        
        
        List<File> cnvNatorFiles = new ArrayList<File>();
        
        cnvNatorFiles.add(cnvNatorLE);
        cnvNatorFiles.add(cnvNatorACI);
        cnvNatorFiles.add(cnvNatorBNLX);
        cnvNatorFiles.add(cnvNatorBNSSN);
        cnvNatorFiles.add(cnvNatorBUF);
        cnvNatorFiles.add(cnvNatorF344);
        cnvNatorFiles.add(cnvNatorM520);
        cnvNatorFiles.add(cnvNatorMR);
        cnvNatorFiles.add(cnvNatorSHR);
        cnvNatorFiles.add(cnvNatorWKY);
        cnvNatorFiles.add(cnvNatorWN);
        
        
        for(File cnvNatorFile : cnvNatorFiles)
        {        
            CNVnatorParser cNVnatorParser = new CNVnatorParser(cnvNatorFile);
            EnumMap<RatChromosomes, ArrayList<Deletion>> cnvDeletionCalls = cNVnatorParser.readCNVSVCalls();        
            
            addNonOverlappingDeletions(totalNonOverlappingDeletions, cnvDeletionCalls );
        }
        
        
        for(RatChromosomes ratChromosomes: RatChromosomes.values())
        {
            for( Deletion nonOverlappingDeletion :    totalNonOverlappingDeletions.get(ratChromosomes))
            {   
                System.out.println(nonOverlappingDeletion.chromosome+":"+nonOverlappingDeletion.getLocation().getMinimum()+"-"+nonOverlappingDeletion.getLocation().getMaximum());
                
                bw.write(nonOverlappingDeletion.chromosome+":"+nonOverlappingDeletion.getLocation().getMinimum()+"-"+nonOverlappingDeletion.getLocation().getMaximum());
                bw.write("\n");
                        
            }
        }
        
        bw.close();
        
        
   
        
        String blaat = "blaat";
        
    }
    
     private void addNonOverlappingDeletions(EnumMap<RatChromosomes, ArrayList<Deletion>> totalNonOverlappingDeletions, EnumMap<RatChromosomes, ArrayList<Deletion>> toAddCalls) {
       
        Integer addedDeletions = 0; 
       
        for(RatChromosomes ratChromosome: RatChromosomes.values())
        {
            ArrayList<Deletion> knownDeletions = totalNonOverlappingDeletions.get(ratChromosome);
            ArrayList<Deletion> toAddDeletions = toAddCalls.get(ratChromosome);
            
            for(Deletion toAdDeletion : toAddDeletions)
            {
                Boolean newDeletion = true;
                
                for(Deletion knownDeletion : knownDeletions)
                {
                    if(knownDeletion.getLocation().isOverlappedBy(toAdDeletion.getLocation()))
                    {
                        newDeletion = false;
                        
//                        System.out.println("to add deletion "+toAdDeletion.getChromosome()+ ":"+ toAdDeletion.getLocation().getMinimum() + "-" + toAdDeletion.getLocation().getMaximum());
//                        System.out.println("overlaps known deletion "+knownDeletion.getChromosome()+ ":"+ knownDeletion.getLocation().getMinimum() + "-" + knownDeletion.getLocation().getMaximum());
//                        System.out.println("\n");
                    }
                }
                
                if(newDeletion)
                {
                    knownDeletions.add(toAdDeletion);   
                    addedDeletions++;
                }            
            }
            
            totalNonOverlappingDeletions.put(ratChromosome, knownDeletions);
        
        }
        
         System.out.println("added "+addedDeletions + " new deletions");
         
    }
    
    
    
    
    
    public static void main(String[] args) throws FileNotFoundException, IOException 
    {
        new CompareCNVNatorCalls();
    }

   
    
}
