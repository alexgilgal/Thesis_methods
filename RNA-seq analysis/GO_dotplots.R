### GO plot function 

generate_GO_plot <- function(Files, top_n_GOs = 10, Conditions, let.size = 14 ){
  
  if(!length(Files)==length(Conditions)){
    stop("Conditions shoould be the same length than Files, they will be their names in the plot!")
  }
  # loading the libraries
  
  require(ggplot2)
  require(dplyr)
  require(viridis)
  require(vroom)
  
  # selecttion GO terms. Taking into account common terms in all the possible combinations
  
  used_terms <- NULL
  for (i in 1:length(Files)) {
    
    i_table <- vroom(Files[i])
    
    i_table <- i_table %>% filter(Annotated < 3000)
    
    used_terms <- unique(c(used_terms, i_table$Term[1:top_n_GOs]))
  }
  
  for (i in 1:length(Files)) {
    
    if (i==1) {
      
      GO_table <- read.delim(Files[i], sep=" ")
      
      GO_table$Condition <- Conditions[i]
      
      GO_table_combined <- GO_table  %>% filter(Annotated < 3000) %>% arrange(elim) %>%  filter(Term %in% used_terms)
      
    } else {
      GO_table_i <- read.delim(Files[i], sep=" ")
      
      if (nrow(GO_table_i) < top_n_GOs) {
        
        warning(paste("The file", Files[i], "has less than", top_n_GOs, "rows. Consider reducing top_n_GOs."))
        
      }
      
      # we keep the used GO terms in order to keep all those terms that are apearing, even though they are not at the current file top10
      

      
      GO_table_i <- GO_table_i %>% filter(Annotated < 3000)  %>% arrange(elim)%>% filter(Term %in% used_terms)
      
      GO_table_i$Condition <- Conditions[i]
      
      GO_table_combined <- rbind(GO_table_combined, GO_table_i)
    }
  }
  
  GO_table_combined$elim <- -log10(as.numeric(GO_table_combined$elim))
  
  GO_table_combined$elim[is.na(GO_table_combined$elim)] <- max(GO_table_combined$elim, na.rm = T)
  
  GO_table_combined <- GO_table_combined[order(GO_table_combined$Condition, GO_table_combined$Term),]
  
  GO_table_combined$Term <- factor(GO_table_combined$Term, levels = unique(GO_table_combined$Term))
  
  gg <- GO_table_combined  %>% 
    ggplot(aes(x=Condition, y = Term, size = Significant, color = elim)) + geom_point() +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = let.size ))+ 
    theme(axis.text.y.left = element_text(hjust = 1, size = let.size ),
          legend.title=element_text(size=let.size),
          legend.text=element_text(size=let.size))+ 
    scale_size_continuous(name = "Genes affected",range = c(3,10)) + scale_color_viridis("-log10(P-value)") +
    xlab("")+ylab("")
  
  return(gg)
  
}
