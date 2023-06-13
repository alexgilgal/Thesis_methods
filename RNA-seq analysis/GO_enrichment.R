############## TopGO Gene Ontology ##############

library(topGO)

# Loading files with ENSEMBL gene ID

ensm_file <- dir(path = 'analisys/results_ens94', pattern = 'ENS_ID', full.names = T)

# Loading translation from ENS_id to Gene_symbol

eqiv <- read.delim('~/cavefish/Cavefish_Trans_Gene_name_ENS99.txt', header = T, stringsAsFactors = F, na.strings = '')

# Loading Ens_id to GO file and preparing it to the correct format (list)

id_go <-  read.delim('~/cavefish/Ens99_geneID_GOTerm.txt', header = T, stringsAsFactors = F, na.strings = '')

# sum(is.na(id_go$GO.term.accession))/nrow(id_go) # 30% of ENS_ids dont have assoc GO

id_go <- id_go[!is.na(id_go$GO.term.accession),]


gene.to.GO <- split(id_go$GO.term.accession, id_go$Gene.stable.ID) 
gene.to.GO <- lapply(gene.to.GO, unique) # to remove duplicates

# str(gene.to.GO)

library(topGO)
library(xtable)

### For loop that will go file by file making GO using top GO

All_cave_genes <- unique(read.delim('~/cavefish/Cavefish_Trans_Gene_name_ENS99.txt',
                             header = T, stringsAsFactors = F, na.strings = '')[,1])

for (i in ensm_file) {
  
message(i)
  
# Loading all the "background genes"

universe <- All_cave_genes

names(universe) <- universe

# Genes of interest are from each file evaluated by for loop

genelist <- read.delim(i, header = F, stringsAsFactors = F)[,1]

# Universe need to be a named factor, with 2 factors (selected or not , 0,1) and the name beeing the ENS_id

universe <- as.factor(as.integer(universe %in% genelist))

names(universe) <- All_cave_genes

# generation of topGO object. We take only BP ontology. Node size is referenced to how many genes must 
# be annotated in a GO node to be considered

sampleGOdata <- new("topGOdata",
                    description = "cavefish_ens_94", ontology = "BP",
                    allGenes = universe, nodeSize = 10,
                    annot =annFUN.gene2GO, gene2GO = gene.to.GO
)

# we use fisher test and the elim parameter (which yields more specific nodes)

resultFisher <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")


## GenTable retrieves a table with significants GO terms, but no with inside genes.

GOTable <- GenTable(object = sampleGOdata, elim = resultFisher, topNodes = sum(resultFisher@score < 0.01),
                    numChar = 5000)

# genesInTerm gives us all genes inside all GOTable ids, in the form of named list

go_w_terms <- genesInTerm(sampleGOdata, whichGO=GOTable$GO.ID)

# we transform this to a named character vectors which converts ENs_id to Gene_symbol and separe
# everything by comma

GO_with_genesym <- sapply(go_w_terms, function(x){
  
  # we obtain the genes that have annotation from Ens to Symbol 
  # then we deduplicate results, and sort them alphabetically.
  # Since genesInTerm retrieves ALL genes inside a node, we MUST take only significant genes.
  # Finally we paste all of them using collapse
  
  filt_genes <- x[x %in% genelist]
  
  symbol <- sort(unique(eqiv[eqiv$Gene.stable.ID %in% filt_genes, "Gene.name"]))

  
  paste(symbol, collapse = ', ')
  
})

### Modify GO tables and write them

GOTable$Genes <- GO_with_genesym

GOTable <- GOTable[,c(1,2,7,3:6)]

print(xtable(GOTable, auto = T), type = 'html',
      file = paste('analisys/results_ens94/Gene_ontology_2022_09/BP_', gsub('analisys/results_ens94/', '', gsub(pattern = '_ENS_ID.txt', replacement = '', i)),
                   '.html', sep = ''))

write.table(GOTable, file =  paste('analisys/results_ens94/Gene_ontology_2022_09/BP_TXT_', gsub('analisys/results_ens94/', '', gsub(pattern = '_ENS_ID.txt', replacement = '', i)),
                                 '.txt', sep = ''))


### Finally, we will convert All significant Ens_id to Gene_symbol, improving my sanity

gene_symbols <- eqiv[eqiv$Gene.stable.ID %in% genelist, "Gene.name"]

gene_symbols <- gene_symbols[!is.na(gene_symbols)]

write.table(sort(unique(gene_symbols)), 
            file =  paste(gsub(pattern = '_ENS_ID.txt', replacement = '_Gene_sym.txt', i), sep = ''),
            sep = "\t", quote = F, col.names = F, row.names = F)

}
