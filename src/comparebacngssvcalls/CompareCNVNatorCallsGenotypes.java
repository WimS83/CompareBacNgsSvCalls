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

/**
 *
 * @author wim
 */
public class CompareCNVNatorCallsGenotypes {
    
    
    public static void main(String[] args) throws FileNotFoundException, IOException 
    {
            
        
        File ACI_genotypes = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/ACI/ACI_genotyes.out");              
        File BNLX_genotypes = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/BNLX/BNLX_genotyes.out");      
        File LE_genotypes = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/LE/LE_genotyes.out");      
        File BNSSN_genotypes = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/BNSSN/BNSSN_genotyes.out");      
        File BUF_genotypes = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/BUF/BUF_genotyes.out");      
        File F344_genotypes = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/F334/F334_genotyes.out");      
        File MR_genotypes = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/MR/MR_genotyes.out");      
        File SHR_genotypes = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/SHR/SHR_genotyes.out");      
        File WKY_genotypes = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/WKY/WKY_genotyes.out");      
        File WN_genotypes = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/WN/WN_genotyes.out");      
        File M520_genotypes = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/M520/M520_genotyes.out");     
        
        List<File> genotypeFiles = new ArrayList<File>();
        genotypeFiles.add(ACI_genotypes);
        genotypeFiles.add(BNLX_genotypes);
        genotypeFiles.add(LE_genotypes);
        genotypeFiles.add(BNSSN_genotypes);
        genotypeFiles.add(BUF_genotypes);
        genotypeFiles.add(F344_genotypes);
        genotypeFiles.add(MR_genotypes);
        genotypeFiles.add(SHR_genotypes);
        genotypeFiles.add(WKY_genotypes);
        genotypeFiles.add(WN_genotypes);        
        genotypeFiles.add(M520_genotypes);
        
        
        EnumMap<RatChromosomes, ArrayList<Deletion>> deletionListPerChromosome = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);
        for(RatChromosomes ratChromosomes: RatChromosomes.values())
        {
            ArrayList deletions = new ArrayList<Deletion>();
            deletionListPerChromosome.put(ratChromosomes, deletions);
        }
        
        
        
        
        for(File genotypeFile: genotypeFiles)
        {
            String fileName = genotypeFile.getName();
            String[] fileNameSplit = fileName.split("_");
            String strain = fileNameSplit[0];
            
            
            CNVNatorGenotypeParser cNVNatorGenotypeParser = new CNVNatorGenotypeParser(genotypeFile, strain);
            cNVNatorGenotypeParser.readCNVSGenotypes();
        
        }
        
        
        
        
        
        
        
        
        
    
    
    
    
    
    
    }
    
}
