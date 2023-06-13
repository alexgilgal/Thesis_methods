library(dplyr)
library(motifbreakR)
library(BSgenome)
library(BSgenome.Amexicanus.Ensembl.Astmex2)
library(MotifDb)
library(BiocParallel)
library(ggplot2)
library(vroom)



## getting the genome of Astyanax

#available.genomes() # of course, our genome is not available

# Making a 2bit file for the genome (run only once)

# genome_fa <- Biostrings::readDNAStringSet("~/genomes/AstMex2.0/AstMex_v2_fixed_names.fa")
# 
# rtracklayer::export.2bit(genome_fa, "data/AstMex2.0_fixed_names.2bit")
# 
# forgeBSgenomeDataPkg(x="data/bsgenomesrc/AstMex2.0_fixed_names-seed", seqs_srcdir = "data/bsgenomesrc/",
#                      destdir = "data/Astmex2_bsgenome/",verbose = T)

### reading the vcf that we produced using CRISP

# the file results/Genotyped_variant_calls_CRISP.vcf.gz was generated using CRISP and bgziped. Then we used tabix -p vcf file
#in order to generate an index file

#  variants <- snps.from.file(file = "results/Genotyped_variant_calls_CRISP_filtered_v2.vcf.gz", format = "vcf",
#                             search.genome = BSgenome.Amexicanus.Ensembl.Astmex2)
# 
#  #we need to trim elements out of boundaries (?)
# 
#  variants_mod <- GenomicRanges::trim(variants,use.names = T)
# 
# # for some reason this does nothing, we will remove the values manually
# 
#  variants_mod <- variants_mod[!(seqnames(variants_mod) %in% c("Scaffold2973", "Scaffold3656", "Scaffold4631")),]
#  
#  motifs <- query(MotifDb, "jaspar")
#  
#  results <- motifbreakR(snpList = variants_mod, filterp = T,
#                         pwmList = motifs, threshold = 1e-4,
#                         method = "ic",
#                         bkg=c(A=0.25, C=0.25, G=0.25, T=0.25),BPPARAM = BiocParallel::MulticoreParam(32))
#  
#  
#  saveRDS(results, file = paste0("results/motifbreakR_results", "_2022_10_03",".rds"))

#### We can load the rds instead of computing everything #########

results <- readRDS("results/motifbreakR_results_2022_10_03.rds")

# results_pval <- calculatePvalue(results)

#plot some examples

results <- results[!(is.na(results$geneSymbol)),]

klf4 <- results[results$geneSymbol == "KLF4" & results$effect == "strong", ]

ctcf <- results[results$geneSymbol == "CTCF", ]

options(ucscChromosomeNames=FALSE)

plotMB(klf4, rsid = "chr3:16890413_A/G") #klf4 improved
plotMB(results, rsid = "chr1:1729503_G/A") #crx improved
plotMB(ctcf,rsid="chr24:13347898_G/C")

## We will see the top 20 most mutated elements:

TF_ids <- as.data.frame(results[!is.na(results$geneSymbol) & results$effect == "strong", ]$geneSymbol)

colnames(TF_ids) <- "geneSymbol"

Tf_counts <- TF_ids %>% count(geneSymbol)

Tf_counts <- Tf_counts[base::order(Tf_counts$n, decreasing = T), ]


top_20 <- Tf_counts[1:20,]

top_20$geneSymbol <- factor(top_20$geneSymbol, levels = top_20$geneSymbol)

Tf_counts$geneSymbol <- factor(Tf_counts$geneSymbol, levels = Tf_counts$geneSymbol)

Tf_counts %>% ggplot(aes(x=geneSymbol, y = n)) + geom_col() + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

top_20 %>% ggplot(aes(x=geneSymbol, y = n)) + geom_col() + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### Lets try to do GSEA with this list of ordered genes

library(fgsea)
library(reactome.db)
library(org.Hs.eg.db)

# We need to transform the data to ENTREZ. since we are using Human data with homer, it shouldnt be a problem

hs <- org.Hs.eg.db
my.symbols <- Tf_counts$geneSymbol
Tf_Entrez <- AnnotationDbi::select(hs, 
       keys = as.character(my.symbols),
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

# merging tables

colnames(Tf_Entrez) <- c("geneSymbol", "ENTREZID")

Tf_counts_merged <- merge(Tf_counts, Tf_Entrez)

Tf_counts_merged <- Tf_counts_merged[!is.na(Tf_counts_merged$ENTREZID),]

Tf_counts_merged$rank <- base::order(Tf_counts_merged$n, decreasing = T)

pathways <- reactomePathways(Tf_counts_merged$ENTREZID)


named_ranks <- Tf_counts_merged$rank

names(named_ranks) <- Tf_counts_merged$ENTREZID

fgsea_result <- fgsea(pathways, named_ranks, eps=0.0, nproc = 12, scoreType = "pos", minSize = 5)

fgsea_result_filt <- fgsea_result[fgsea_result$pval < 0.05,]

plotEnrichment(pathways[["Circadian Clock"]], named_ranks) + labs(title = "Circadian rithm")

plotEnrichment(pathways[["Signal Transduction"]], named_ranks) + labs(title = "Signal Transduction")

plotEnrichment(pathways[["Transcriptional regulation of white adipocyte differentiation"]], named_ranks) + labs(title = "Adipocyte differentiation")


## now we will use GO terms not reactome pathways

Bp <- gmtPathways("data/c5.go.bp.v7.2.entrez.xls")

fgsea_bp <- fgsea(Bp, named_ranks, eps=0.0, nproc = 12, minSize = 5)

fgsea_bp_filt <- fgsea_bp[fgsea_bp$pval<0.05,]

fgsea_bp_filt <- fgsea_bp_filt[, leadingEdge := mapIdsList(x=org.Hs.eg.db,
                                                           keys =leadingEdge ,
                                                           keytype = "ENTREZID",
                                                           column = "SYMBOL")]

for (j in 1:nrow(fgsea_bp_filt)) {
  
  # print(paste('Row', j, 'of',nrow(table)))
  
  Entrez.ids <- unlist( fgsea_bp_filt[j,leadingEdge])

  fgsea_bp_filt[j,"leadingEdge"] <- base::paste(Entrez.ids,collapse = ", ")
  
  
}

fgsea_bp_filt <- fgsea_bp_filt[order(fgsea_bp_filt$pval), ]

plotEnrichment(Bp[["GO_CIRCADIAN_RHYTHM"]], named_ranks) + labs(title = "Circadian rithm")
GO_NEUROGENESIS

plotEnrichment(Bp[["GO_NEUROGENESIS"]], named_ranks) + labs(title = "Neurogenesis")
plotEnrichment(Bp[["GO_RESPONSE_TO_INSULIN"]], named_ranks) + labs(title = "Response to Insulin")


###### Crossing the info with TOBIAS #########

results <- readRDS("results/motifbreakR_results.rds")

# results_pval <- calculatePvalue(results)

#plot some examples

results <- results[!(is.na(results$geneSymbol)),]

# first wee ned to take all the footprint info that tobias has for ALL TFs.

# we have computed it using "cat */*overview.txt > All_TFBS_overviews.txt" in TOBIAS_snakemake/results_tobias/TFBS/

# Then we have use intersect bed

tobias_motifbreakr <- vroom("results/results_tobias_motifbreaksr_merge.bed", col_names = F)

# we can take only those results that overlaps with QTL regions

# tobias_motifbreakr <- vroom("results/results_tobias_motifbreaksr_merge_QTL_overlaping.bed", col_names = F)


results_tobias <- vroom("TOBIAS_snakemake/results_tobias/TFBS/AhrArnt_MA0006.1_overview.txt", col_names = T)

colnames(tobias_motifbreakr) <- c(colnames(results_tobias),
                                 c("chr", "start", "end", "score", "strand",colnames(mcols(results))), "foo")

# we can start to filter things

tobias_motifbreakr_filt <- as.data.frame(tobias_motifbreakr[tobias_motifbreakr$effect == "strong",])

hist(tobias_motifbreakr_filt$Cave_Surface_log2fc, breaks = 200)

tobias_motifbreakr_filt <- tobias_motifbreakr_filt[abs(tobias_motifbreakr_filt$Cave_Surface_log2fc) > log2(1.5),]

tobias_motifbreakr_filt_Not_binding <- tobias_motifbreakr_filt[!(tobias_motifbreakr_filt$Cave_bound & tobias_motifbreakr_filt$Surface_bound),]

# we filter by motif bounding in both populations

# tobias_motifbreakr_filt <- tobias_motifbreakr_filt[tobias_motifbreakr_filt$Cave_bound &
                                                     # tobias_motifbreakr_filt$Surface_bound,]


# Now we want to use this cross information in order to know broken regulation loops between TFs

# we have used the GO:0006355 term ()regulation of transcription from our most recent GO database

GO_ens <- vroom("~/cavefish/ENSEMBL_GO_2022_90.txt")

astyanax_TFs <- GO_ens %>% filter(grepl("GO:0006355", `GO term accession`))

tobias_motifbreakr_filt_TF_targets <- tobias_motifbreakr_filt[tobias_motifbreakr_filt$gene_id %in% astyanax_TFs$`Gene stable ID`,]

# # Now lets try to filter and keep only TFBS affected in Tobias that are Affected also in motifbreakR

keep <- NULL

for (i in 1:nrow(tobias_motifbreakr_filt_TF_targets)) {

  # (tobias_motifbreakr_filt$geneSymbol[i] %in% tobias_motifbreakr_filt$TFBS_name[i])
  
  # tobias_motifbreakr_filt$genename is the name of the gene that is near a TFBS altered
  
  # tobias_motifbreakr_filt$geneSymbol is the name of the TF that binds to that TFBS

  keep <- c(keep, grepl(pattern =tobias_motifbreakr_filt_TF_targets$geneSymbol[i], x =  tobias_motifbreakr_filt_TF_targets$TFBS_name[i],
        ignore.case = T))

}

tobias_motifbreakr_filt_TF_2_TF <- tobias_motifbreakr_filt_TF_targets[keep,]

## We can also explore TF-nonTF loops 

keep <- NULL

for (i in 1:nrow(tobias_motifbreakr_filt)) {
  
  # (tobias_motifbreakr_filt$geneSymbol[i] %in% tobias_motifbreakr_filt$TFBS_name[i])
  
  # tobias_motifbreakr_filt$genename is the name of the gene that is near a TFBS altered
  
  # tobias_motifbreakr_filt$geneSymbol is the name of the TF that binds to that TFBS
  
  keep <- c(keep, grepl(pattern =tobias_motifbreakr_filt$geneSymbol[i], x =  tobias_motifbreakr_filt$TFBS_name[i],
                        ignore.case = T))
  
}

tobias_motifbreakr_filt_TF_nonTF <- tobias_motifbreakr_filt[keep & !(tobias_motifbreakr_filt$gene_id %in% astyanax_TFs$`Gene stable ID`) ,]


### Lets keep those regulation loops in which the TF changes and has a motifbreakR instance:

tobias_motifbreakr_filt_v2 <- tobias_motifbreakr_filt[keep,]


tobias_motifbreakr_filt_v2 <- tobias_motifbreakr_filt_v2[!is.na(tobias_motifbreakr_filt_v2$gene_id),]

write.table(tobias_motifbreakr_filt_v2,
            file = "results/results_tobias_motifbreakr/results_tobias_motifbreakr_filtered.bed",
            quote = F, row.names = F, sep = "\t", col.names = F)

# We have generated the next file using: 
# intersectBed -wa -a results_tobias_motifbreakr_filtered.bed /
# -b ~/cavefish/ATAC/new_genome/analysis_cave_vs_surface/results_idr/All_ATAC_IDR_diff_peaks.bed > /
# results_tobias_motifbreakr_Within_diff_peaks.bed

tobias_motifbreakr_filt_v2_diff_peaks <- vroom("results/results_tobias_motifbreakr/results_tobias_motifbreakr_Within_diff_peaks.bed", col_names = F)
colnames(tobias_motifbreakr_filt_v2_diff_peaks) <- colnames(tobias_motifbreakr_filt_v2)

# Lets see which are the TFs more affected by these changes:

top_tobias_motifbreakr_diff_TFs <- tobias_motifbreakr_filt_v2_diff_peaks %>% count(geneSymbol) %>% dplyr::arrange(desc(n)) %>% filter(n > 1)

top_tobias_motifbreakr_diff_TFs$geneSymbol <- factor(top_tobias_motifbreakr_diff_TFs$geneSymbol,
                                                     levels = top_tobias_motifbreakr_diff_TFs$geneSymbol[order(top_tobias_motifbreakr_diff_TFs$n, decreasing = T)])

top_tobias_motifbreakr_diff_TFs %>% 
ggplot(aes(y=geneSymbol, x = n)) +
  geom_col(fill= "#123f5d") + 
  theme_classic() +
  theme(axis.text.x = element_text( vjust = 0.5, hjust=1, size = 24, face = "bold"),
        axis.text.y = element_text(size=24, face = "bold")) + xlab("") + ylab("")

ggsave("results_idr/Top_Tf_diff_mutated.pdf", height = 12, width = 10 )

"#123f5d"


### Network visualization ####

### Lets try generate a cool network graphic using this tutorial as example: https://kateto.net/network-visualization
library(igraph)
library(dplyr)
# We will select those relationships between TFs

tobias_motifbreakr_filt_v2_network <- tobias_motifbreakr_filt_v2_diff_peaks %>% 
  dplyr::select(geneSymbol, gene_name, Cave_Surface_log2fc, alleleDiff) %>% filter(!is.na(gene_name) & !is.na(geneSymbol)) %>% 
  mutate(alleleDiff = ifelse(alleleDiff >0, "Improved", "Disrupted"),
         geneSymbol = tolower(geneSymbol),
         gene_name = tolower(gene_name))



g <- graph_from_data_frame(d=tobias_motifbreakr_filt_v2_network)

E(g)$color <- ifelse(E(g)$alleleDiff == "Improved", "blue", "red")

E(g)$width <- abs(E(g)$Cave_Surface_log2fc)


l <- norm_coords(layout_with_graphopt(g), ymin=-1, ymax=1, xmin=-1, xmax=1)


plot(g, edge.arrow.size=0.2, layaout = l,
     vertex.size = 5, edge.length = 10,
     vertex.label.cex = 1)


### This is a bit messy, lets try to find genes that have lateral line function

GO_ens <- vroom("~/cavefish/ENSEMBL_GO_2022_90.txt")

Lat_line_genes <- GO_ens %>%
  filter(grepl("sensory organ", `GO term name`)) %>%
  filter(!duplicated(`Gene stable ID`)) #%>% 
  # select(`Gene stable ID`,`Gene name`)


## now lets cross the info that we have

tobias_motifbreakr_filt_v2_lat_line <- tobias_motifbreakr_filt_v2[tobias_motifbreakr_filt_v2$gene_id %in% Lat_line_genes$`Gene stable ID`,]

# not a single gene of lateral line is in our dataset :_( 


#### topGO enrichment ####


library(topGO)
library(xtable)
library(knitr)
library(kableExtra)

# Loading translation from ENS_id to Gene_symbol

eqiv <- read.delim('~/cavefish/Cavefish_Trans_Gene_name_ENS99.txt', header = T, stringsAsFactors = F, na.strings = '')

# Loading Ens_id to GO file and preparing it to the correct format (list)

id_go <-  read.delim('~/cavefish/Ens99_geneID_GOTerm.txt', header = T, stringsAsFactors = F, na.strings = '')

# sum(is.na(id_go$GO.term.accession))/nrow(id_go) # 30% of ENS_ids dont have assoc GO

id_go <- id_go[!is.na(id_go$GO.term.accession),]


gene.to.GO <- split(id_go$GO.term.accession, id_go$Gene.stable.ID) 
gene.to.GO <- lapply(gene.to.GO, unique) # to remove duplicates

# str(gene.to.GO)

### For loop that will go file by file making GO using top GO

All_cave_genes <- read.delim('~/cavefish/Cavefish_Trans_Gene_name_ENS99.txt',
                             header = T, stringsAsFactors = F, na.strings = '')[,1]

# Loading all the "background genes"

universe <- All_cave_genes

names(universe) <- universe

# Genes of interest are extracted from the tobias_motifbreakr_filt_v2


genelist <- unique(tobias_motifbreakr_filt_v2$gene_id)

genelist <- genelist[!is.na(genelist)]
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

resultFisher <- runTest(sampleGOdata,algorithm = "elim", statistic = "fisher")


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


write.table(GOTable, file = "results/results_tobias_motifbreakr/BP_TXT_Genes_affected_tobias_motifbreakr.txt")

kable(GOTable) %>%
  kable_styling(bootstrap_options = c( "hover", "bordered", "responsive"), font_size = 18) %>%
  save_kable(file ="results/results_tobias_motifbreakr/BP_Genes_affected_tobias_motifbreakr.html")



