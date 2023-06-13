

### Selecting those TF that are in peaks with mutations inside. We will perform GO of the filtered sets.

setwd("~/cavefish/ATAC/new_genome/analysis_cave_vs_surface/results/motifs/")

Up_24 <- read.delim('SNP_ATAC_Cave_24hpf_sobreexp/knownResults.txt', header = T, stringsAsFactors = F)[1:20,1]

Up_80 <- read.delim('SNP_ATAC_Cave_80epi_sobreexp/knownResults.txt', header = T, stringsAsFactors = F)[1:20,1]

Up_5ss <- read.delim('SNP_ATAC_Cave_5ss_sobreexp/knownResults.txt', header = T, stringsAsFactors = F)[1:20,1]
  
Up_48 <- read.delim('SNP_ATAC_Cave_48hpf_sobreexp/knownResults.txt', header = T, stringsAsFactors = F)[1:20,1]


Down_24 <- read.delim('SNP_ATAC_Cave_24hpf_subexp/knownResults.txt', header = T, stringsAsFactors = F)[1:20,1]

Down_80 <- read.delim('SNP_ATAC_Cave_80epi_subexp/knownResults.txt', header = T, stringsAsFactors = F)[1:20,1]

Down_5ss <- read.delim('SNP_ATAC_Cave_5ss_subexp/knownResults.txt', header = T, stringsAsFactors = F)[1:20,1]

Down_48 <- read.delim('SNP_ATAC_Cave_48hpf_subexp/knownResults.txt', header = T, stringsAsFactors = F)[1:20,1]


### We make a pull of TF depending on their status, UP, or down. We will eliminate those that 
## are both up and down, by stage.

Down_24_u <- Down_24[!(Down_24 %in% Up_24)]

Down_48_u <- Down_48[!(Down_48 %in% Up_48)]

Down_5ss_u <- Down_5ss[!(Down_5ss %in% Up_5ss)]

Down_80_u <- Down_80[!(Down_80 %in% Up_80)]

down_u_jl <- unique(c(Down_24_u, Down_48_u, Down_5ss_u, Down_80_u))

down_u_jl <- gsub(pattern = '\\(\\d*\\D*\\S*',replacement =  '', x = down_u_jl)

write.table(down_u_jl, file = 'TFs_in_mutated_down_peaks.txt', quote = F, col.names = F, row.names = F)


Up_24_u <- Up_24[!(Up_24 %in% Down_24)]

Up_48_u <- Up_48[!(Up_48 %in% Down_48)]

Up_5ss_u <- Up_5ss[!(Up_5ss %in% Down_5ss)]

Up_80_u <- Up_80[!(Up_80 %in% Down_80)]

Up_u_jl <- unique(c(Up_24_u, Up_48_u, Up_5ss_u, Up_80_u))

Up_u_jl <- gsub(pattern = '\\(\\d*\\D*\\S*',replacement =  '', x = Up_u_jl)

write.table(Up_u_jl, file = 'TFs_in_mutated_Up_peaks.txt', quote = F, col.names = F, row.names = F)



#### Now we will perform the same for TFs that are around mutations (hopefully affected by those mutations)

snp_up_24 <- read.delim('Variants_inSNP_ATAC_Cave_24hpf_sobreexp/knownResults.txt', header = T, stringsAsFactors = F)[1:5,1]

snp_up_80 <- read.delim('Variants_inSNP_ATAC_Cave_80epi_sobreexp/knownResults.txt', header = T, stringsAsFactors = F)[1:5,1]

snp_up_5ss <- read.delim('Variants_inSNP_ATAC_Cave_5ss_sobreexp/knownResults.txt', header = T, stringsAsFactors = F)[1:5,1]

snp_up_48 <- read.delim('Variants_inSNP_ATAC_Cave_48hpf_sobreexp/knownResults.txt', header = T, stringsAsFactors = F)[1:5,1]

snp_down_24 <- read.delim('Variants_inSNP_ATAC_Cave_24hpf_subexp/knownResults.txt', header = T, stringsAsFactors = F)[1:5,1]

snp_down_80 <- read.delim('Variants_inSNP_ATAC_Cave_80epi_subexp/knownResults.txt', header = T, stringsAsFactors = F)[1:5,1]

snp_down_5ss <- read.delim('Variants_inSNP_ATAC_Cave_5ss_subexp/knownResults.txt', header = T, stringsAsFactors = F)[1:5,1]

snp_down_48 <- read.delim('Variants_inSNP_ATAC_Cave_48hpf_subexp/knownResults.txt', header = T, stringsAsFactors = F)[1:5,1]

snp_up <- unique(c(snp_up_24, snp_up_80, snp_up_48, snp_up_5ss))

snp_down <- unique(c(snp_down_24, snp_down_80, snp_down_48, snp_down_5ss))

snp_up <-  gsub(pattern = '\\(\\d*\\D*\\S*',replacement =  '', x = snp_up)

snp_down <-  gsub(pattern = '\\(\\d*\\D*\\S*',replacement =  '', x = snp_down)

### David txt table to html


library(xtable)



# txt_table <- read.delim('/data/Zebra_drugs/bio_problem/Upregulated/BP_Bio_all_sobrexp.txt')
# 
# txt_table
# 
# out_table_x <- xtable(txt_table, auto = T)
# 
# print(out_table_x, type = 'html', file = 'example.html')

tables <- dir(pattern = 'BP_')

for (i in tables) {
  print(i)
  
  raw <- read.delim(i)
  
  print(xtable(raw, auto = T), type = 'html',
        file = paste(gsub(pattern = '.txt', replacement = '', i), '.html', sep = ''))
  
}

