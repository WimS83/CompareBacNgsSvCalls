/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package comparebacngssvcalls;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.TreeMap;
import org.apache.commons.lang3.Range;
import org.apache.commons.lang3.StringUtils;


/**
 *
 * @author wim
 */
public class CompareCNVNatorCallsgenotypes2 {
    
    
    public static void main(String[] args) throws FileNotFoundException, IOException 
    {
        
        //input file with all the cnv genotype scores for the rats
        File genotypesFile = new File("/home/sge_share_fedor12/wim/Analysis/ratfounders/raw_data/cnvnator_calls/100bpWindows/genotypes.txt");   
        
        //output file for the normal heatmap
        File heatMapCSVFile = new File("/home/wim/heatMap.csv");  
        FileWriter fw = new FileWriter(heatMapCSVFile.getAbsoluteFile());        
        BufferedWriter bw = new BufferedWriter(fw);
        
        //output file for the category heatmap
        File heatMapCategoryCSVFile = new File("/home/wim/heatMapCategory.csv");  
        FileWriter fwCategory = new FileWriter(heatMapCategoryCSVFile.getAbsoluteFile());        
        BufferedWriter bwCategory = new BufferedWriter(fwCategory);
        
        
        //map of all the rat strains
        TreeMap<String, Double> genotypesMapDummy = new TreeMap<String, Double>();
        genotypesMapDummy.put("ACI",null);
        genotypesMapDummy.put("BNLX",null);
        genotypesMapDummy.put("LE",null);
        genotypesMapDummy.put("BNSSN",null);
        genotypesMapDummy.put("BUF",null);
        genotypesMapDummy.put("F344", null);
        genotypesMapDummy.put("MR",null);
        genotypesMapDummy.put("SHR",null);
        genotypesMapDummy.put("WKY", null);
        genotypesMapDummy.put("WN",null);
        genotypesMapDummy.put("M520",null);
        
        //write the header line
        List<String> ratStrains = new ArrayList<>();
        for(String ratStrain : genotypesMapDummy.keySet())
        {
            ratStrains.add(ratStrain);                      
        }        
        bw.write(StringUtils.join(ratStrains, "\t"));
        bw.write("\n");
        bwCategory.write(StringUtils.join(ratStrains, "\t"));
        bwCategory.write("\n");
        
        
        //list to store all the deletions
        EnumMap<RatChromosomes, ArrayList<Deletion>> deletionListPerChromosome = new EnumMap<RatChromosomes, ArrayList<Deletion>>(RatChromosomes.class);
        for(RatChromosomes ratChromosomes: RatChromosomes.values())
        {
            ArrayList deletions = new ArrayList<Deletion>();
            deletionListPerChromosome.put(ratChromosomes, deletions);
        }
        

        BufferedReader cnvGenotypeBR = new BufferedReader(new FileReader(genotypesFile));      
        String line = cnvGenotypeBR.readLine();
        
        Integer deletionCounter = 0;
        Integer deletionCounterDifferentInStrains = 0;       

        //loop over aal the lines in the input file 
        while (line != null) {              
            
            deletionCounter++;
            
            Boolean common = true;            
            
            String[] splitLine = line.split("\\s+");
            String location = splitLine[1];
            String[] locationSplit = location.split(":");
            String locationBeginEnd = locationSplit[1];
            String[] locationBeginEndSplit = locationBeginEnd.split("-");
            
            Integer begin = new Integer(locationBeginEndSplit[0]);
            Integer end = new Integer(locationBeginEndSplit[1]);
          
            Double ACI_genotype = new Double(splitLine[3]);
            Double BNLX_genotype = new Double(splitLine[8]);
            Double LE_genotype = new Double(splitLine[13]);
            Double BNSSN_genotype = new Double(splitLine[18]);
            Double BUF_genotype = new Double(splitLine[23]);
            Double F344_genotype = new Double(splitLine[28]);
            Double MR_genotype = new Double(splitLine[33]);
            Double SHR_genotype = new Double(splitLine[38]);
            Double WKY_genotype = new Double(splitLine[43]);
            Double WN_genotype = new Double(splitLine[48]);
            Double M520_genotype = new Double(splitLine[53]);
            
            TreeMap<String, Double> genotypesMap = new TreeMap<String, Double>();
            genotypesMap.put("ACI",ACI_genotype);
            genotypesMap.put("BNLX",BNLX_genotype);
            genotypesMap.put("LE",LE_genotype);
            genotypesMap.put("BNSSN",BNSSN_genotype);
            genotypesMap.put("BUF",BUF_genotype);
            genotypesMap.put("F344", F344_genotype);
            genotypesMap.put("MR",MR_genotype);
            genotypesMap.put("SHR",SHR_genotype);
            genotypesMap.put("WKY", WKY_genotype);
            genotypesMap.put("WN",WN_genotype);
            genotypesMap.put("M520",M520_genotype);
            
            Boolean categoryDeletion = false;
            Boolean categoryNormal = false;
            Boolean categoryDuplication = false;
            
            StringBuilder sbHeatMap = new StringBuilder();    
            StringBuilder sbHeatMapCategory = new StringBuilder();    
            
            List<String> genotypeValues = new ArrayList<String>();
            List<String> genotypeValuesInCategorie = new ArrayList<String>();
            
            
            
            
            //loop over all the rat strains genotypes for this cnv
            for(String strain : genotypesMap.keySet())
            {
                Double genotype = genotypesMap.get(strain);
                Integer category = null;
                
                genotypeValues.add(genotype.toString());
                
                if(genotype <= new Double("1.5")){categoryDeletion = true;category = 0;}
                if(genotype > new Double("1.5") && genotype < new Double("2.5") ){categoryNormal = true;category = 1;}
                if(genotype >= new Double("2.5")){categoryDuplication = true;category = 2;}
                
//                if(category == null)
//                {
//                    String blaat = "blaat";
//                }
                
                genotypeValuesInCategorie.add(category.toString());
                
            }
            
            
            sbHeatMap.append(StringUtils.join(genotypeValues, "\t"));
            sbHeatMapCategory.append(StringUtils.join(genotypeValuesInCategorie, "\t"));            
           
            //if al least 2 categories are represented in the genotypes            
            if((categoryDeletion && categoryNormal) ||  (categoryDeletion && categoryDuplication) || (categoryNormal && categoryDuplication) )
            {    
               common = false;
               
               sbHeatMap.append("\n");      
               sbHeatMapCategory.append("\n");
               bw.write(sbHeatMap.toString());    
               bwCategory.write(sbHeatMapCategory.toString());
               
               
               deletionCounterDifferentInStrains++;
               
            }                 
         
            
            
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
            deletion.setCommon(common);
            
           // deletion.addRatGenotype(ratStrain, genotype);
            deletionListPerChromosome.get(ratChromosome).add(deletion);      

            line = cnvGenotypeBR.readLine();
            String blaat = "blaat";
        }
        
        bw.close();
        bwCategory.close();
        
        
        
//        double[][] cnvvGenotypesArray = new double[deletionCounterDifferentInStrains][11];
//        
//        Integer xAxisCounter = 0;
//        Integer yAxisCounter = 0;
//        
//        
//        for( RatChromosomes ratChromosomes :   deletionListPerChromosome.keySet())
//        {
//            for(Deletion deletion : deletionListPerChromosome.get(ratChromosomes))
//            {
//                if(deletion.getCommon()){continue;}
//                
//                TreeMap<String, Double> ratGenotypes = deletion.getRatStrainsGenotype();
//                
//                for(String ratStrain : ratGenotypes.keySet())
//                {
//                    Double genotype = ratGenotypes.get(ratStrain);
//                    
//                    cnvvGenotypesArray[xAxisCounter][yAxisCounter] = genotype;
//                    
//                    
//                    
//                    yAxisCounter++;
//                
//                }
//                
//                yAxisCounter = 0;
//                xAxisCounter++;            
//            }
//        }
//        
        
        
        
       
        
        
        
        
         System.out.println("Deletion count "+ deletionCounter);
         System.out.println("Deletion count different in strains "+ deletionCounterDifferentInStrains);
        
    }
     
    
     
     
    
}
