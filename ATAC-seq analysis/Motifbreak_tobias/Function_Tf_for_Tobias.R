

#### Function to get the TF_list in a TOBIAS friendly format

# we need 4 things:
# - The ortology between your Species and Human, taken from Biomart (From your organism:
#    Atributes are : ENSEMBL geneID, Human.gene.stable.ID, Human.gene.name)
# - The TF list of your organism. Taken from here: http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/download
# - The output folder of tobias TFBS: for example TOBIAS_snakemake/results_tobias/TFBS/
# - The desired name of the output file (its TF_list_for_TOBIAS.txt)

# here it goes:

TF_list_for_TOBIAS <- function(Orthology,TF_list, output_file = "TF_list_for_TOBIAS.txt", folder_of_motifs){
  
  TFs <- read.delim(TF_list, stringsAsFactors = F) # data taken from here: http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Astyanax_mexicanus_TF
  
  TFs <- TFs[,c("Symbol","Ensembl")] # we select only those columns that interest us
  
  orthologs <- read.delim(Orthology, stringsAsFactors = F) # we load the human orthologs of these TF
  
  #colnames(TFs)
  colnames(orthologs)[1] <- "Ensembl" # we tweak a little the colnames so we can do a proper merge of the 2 tables
  
  TFs_with_orth <- merge(TFs, orthologs)
  
  TF_list_jaspar <- list.dirs(path = folder_of_motifs, full.names = F, recursive = F) # this is the total list of jaspar TFBS
  
  ## We will use this list of jaspar motifs to generate the table that we need.
  
  TF_final_table <- TFs_with_orth
  
  for (i in TF_list_jaspar) {
    
    symbol <- gsub("_.*$", "", i)
    
    finding <- grep(pattern = paste0("^",symbol,"$"), x = TF_final_table$Human.gene.name, ignore.case = T )
    
    if (!is.null(finding)) { 
      
      TF_final_table[finding,"Symbol"] <- i
      
    }
    
  }
  
  # finally we will keep only those entries that have the proper format:
  
  TF_final_table_filt <- TF_final_table[grepl("^.*_MA.*$", TF_final_table$Symbol), c(1,2)] 
  
  # TF_final_table_filt$Symbol <- gsub("_MA.*$", "", TF_final_table_filt$Symbol)
  
  write.table(TF_final_table_filt, file = output_file, sep = "\t",
              quote = F, col.names = T, row.names = F)
  
}
TF_list_for_TOBIAS(Orthology = "data/Astyanax_mexicanus_TF_human_orthologs.txt", TF_list = "data/Astyanax_mexicanus_TF.txt",
                   folder_of_motifs = "TOBIAS_snakemake/results_tobias/TFBS/", output_file = "foo")
