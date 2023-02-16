generate_ME_plot <- function(Directories, top_n_motifs = 10, Conditions, let.size = 14){
  
  if(!length(Directories)==length(Conditions)){
    stop("Conditions shoould be the same length than Directories, they will be their names in the plot!")
  }
  # loading the libraries
  
  require(ggplot2)
  require(marge)
  require(dplyr)
  require(viridis)
  
  for (folder in Directories) {
    
    if (folder==Directories[1]) {
      
      known_motifs <- marge::read_known_results(folder)
      
      known_motifs_combined <- known_motifs %>% arrange(desc(log_p_value)) %>%  head(top_n_motifs)
      
    } else {
      known_motifs_i <- marge::read_known_results(folder)
      
      if (nrow(known_motifs_i) < top_n_motifs) {
        
        stop(paste("The file", folder, "has less than", top_n_motifs, "rows. Consider reducing top_n_motifs."))
        
      }
      
      known_motifs_i <- known_motifs_i %>% arrange(desc(log_p_value)) %>%  head(n=top_n_motifs)
      
      known_motifs_combined <- rbind(known_motifs_combined, known_motifs_i)
    }
  }
  
  known_motifs_combined$Condition <- factor(rep(Conditions,each = top_n_motifs), levels = rev(Conditions))
  
  known_motifs_combined <- known_motifs_combined[order(known_motifs_combined$Condition, known_motifs_combined$motif_name),]
  
  known_motifs_combined$motif_name <- factor(known_motifs_combined$motif_name, levels = unique(known_motifs_combined$motif_name))
  gg <- known_motifs_combined %>% 
    ggplot(aes(x=motif_name, y = Condition, size = tgt_pct*100, color = log_p_value)) + geom_point() +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = let.size ))+ 
    theme(axis.text.y.left = element_text(hjust = 1, size = let.size ))+ 
    scale_size_continuous(name = "%Peaks with \n motif",range = c(3,10)) + scale_color_viridis("-log10(P value)") +
    xlab("")+ylab("")
  
  return(gg)
  
}