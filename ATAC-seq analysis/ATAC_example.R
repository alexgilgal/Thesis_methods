## ATACseq script

### Library loading

library(DESeq2)
library(BiocParallel)
library(ggplot2)

### Sample Table generation

directory <- "my_counts_folder"

n_threads <- 4 # put a number of cores to use

# Use your own data and follow this guide: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


sampleTable <- read.csv("sampleTable.csv") 

# (Optional) You can remove low counts peaks. A popular filter is to ensure at least X samples with a count of 10 or more

keep <- rowSums(counts(dds) >= 10) >= X
dds <- dds[keep,]

# Load the data from HTseq like files, which are stored in my_counts_folder directory

dds<- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~condition)

dds <- DESeq(dds_allrep, parallel = T,  BPPARAM = MulticoreParam(n_threads))


### PCA plot generation

rlog_data <- rlog(dds)

g <- plotPCA(rlog_data, intgroup= c('condition'), ntop = 1000, returnData = T)

ggplot(g) + geom_point(aes(PC1, PC2, color = stage, shape = condition), size = 3)+
  ggtitle("PCA of ATAC-seq ")+ 
  theme_bw()+
  theme(plot.title = element_text( size=16, hjust=0.5))

ggsave("results/PCA_ATAC.pdf", height = 8, width = 8)


#### Extraction of DE information of each condition

## Definition of statistical thresholds

logfc_thres <- 1

sig_thres <- 0.05

## For loop which will iterate all the stages

  result_deseq <- results(dds, contrast = c("condition", "treatment", "control")) # What you want to compare and that reference condition is the LAST
  
  
  plot(result_deseq$log2FoldChange, -log10(result_deseq$padj),
       pch = 19, cex = 0.3, col = 'grey', xlab = 'Log2FC', ylab = '-log10(P Adj)',
       main = "ATAC-seq volcano plot")
  
  sig_genes <- subset(result_deseq, abs(log2FoldChange) > logfc_thres )
  
  sig_genes <- subset(sig_genes, padj < sig_thres)
  
  points(sig_genes$log2FoldChange, -log10(sig_genes$padj),
         col = 'lightgreen', pch = 19, cex = 0.31)
  
  print(paste('Number of DEs:',nrow(sig_genes)))
  
  print(paste('Upregulated:',nrow(sig_genes[sig_genes$log2FoldChange>0,])))
  print(paste('Downregulated:',nrow(sig_genes[sig_genes$log2FoldChange<0,])))
  
  print(paste('NAs:',sum(is.na(result_deseq$padj))))
  
  print(paste('total features:',length(result_deseq$padj)))
  
  
  #### Writing the information
  
  # Write table, be carefull with paste statement
  
  write.table(sig_genes[sig_genes$log2FoldChange>0,], file = "results/ATAC_upregulated_peaks.txt", quote = F, sep = "\t")
  
  write.table(sig_genes[sig_genes$log2FoldChange<0,], file = "results/ATAC_downregulated_peaks.txt", quote = F, sep = "\t")
  
  write.table(sig_genes, file = "results_idr/ATAC_All_differential.txt", quote = F, sep = "\t")
  
  ### Obtaining position information about those peaks
  
  bed_counts <- read.delim(file = "control_rep1.bed", header = F, # read a bed file of the peaks
                           stringsAsFactors = F)
  
  upreg <- rownames(sig_genes[sig_genes$log2FoldChange>0,])
  
  downreg <- rownames(sig_genes[sig_genes$log2FoldChange<0,])
  
  upreg_bed <- bed_counts[bed_counts$V4 %in% upreg, ]
  
  downreg_bed <-  bed_counts[bed_counts$V4 %in% downreg, ]
  
  write.table(upreg_bed, file = "results/ATAC_upregulated_peaks.bed",
              row.names = F, quote = F, col.names = F, sep = '\t')
  
  write.table(downreg_bed, file = "results/ATAC_downregulated_peaks.bed",
              row.names = F, quote = F, col.names = F, sep = '\t')
  
  write.table(rbind(upreg_bed, downreg_bed), file = "results/ATAC_All_differential_peaks.bed",
              row.names = F, quote = F, col.names = F, sep = '\t')
  
