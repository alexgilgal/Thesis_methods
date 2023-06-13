# function that looks into orthology families to retrieve the orthologous zebrafish genes of amphioxus

# v is a character vector of Bralan gene IDS
# fams are the gene families 
amphi2zebra <- function(v, fams) {
  result <- character(0)
  for(a in v) {
    
     # for(x in 1:nrow(fams)) {
      genes <- unlist(strsplit(as.character(fams[grep(a, fams[,2]),2]), ":"))
      # genes <- unlist(strsplit(as.character(fams[x,2]), ":"))
      if(a %in% genes) {
        genesZ <- unlist(strsplit(as.character(fams[grep(a, fams[,2]),3]), ":"))
        result <- c(result, paste(genesZ, sep = ", "))
        notfound <- 0
      } else {
        notfound <- 1
      }
    # }
    if(notfound) result <- c(result, "NO")
  }
  return(result)
}
