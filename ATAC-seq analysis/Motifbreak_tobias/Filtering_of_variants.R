####### VCF filtering

setwd('cavefish/ATAC/new_genome/analysis_cave_vs_surface/results/')

### loading vcf

library(data.table)

library("VariantAnnotation")

data <- readVcf('Genotyped_variant_calls_CRISP_filtered_v2.vcf' )

head(data)

genodata <- geno(data)$GT

colnames(genodata)

tail(genodata)

mod <- genodata == '.'

dim(mod)

genodata[mod] <- '0'

dim(genodata)

genodata_mod <- matrix(NA, nrow = 45968, ncol = 17)

for (i in seq(nrow(genodata))) {
  
  row <- genodata[i,]
  
  new_row <- NULL
  
  for (j in seq(length(row))) {
    
    suma <- sum(as.integer(unlist(strsplit(row[j], '/'))))
    
    new_row <- c(new_row, suma)
  }
  
  genodata_mod[i,] <- new_row
}


 colnames(genodata_mod) <- colnames(genodata)

 rownames(genodata_mod) <- row.names(genodata) 
 
keep <- NULL

sum_cave <- NULL

sum_surface <- NULL

for (i in seq(nrow(genodata_mod))) {
  
  row <- genodata_mod[i,]
  
  may_i_keep <-  sum(row[1:8]) >= 50 & sum(row[9:17]) < 66
  
  keep <- c(keep, may_i_keep)
  
  sum_cave <- c(sum_cave, sum(row[1:8]))
  
  sum_surface <- c(sum_surface, sum(row[9:17]))
}

sum(keep)

hist(sum_cave, freq = F, col=rgb(1,0,0,0.5))

hist(sum_surface, freq = F, add = T, col=rgb(0,0,1,0.5))

filtered_variants <- genodata[keep,]

filter_positions <- rownames(filtered_variants)

filter_positions <- gsub('_[A,T,C,G]*/[A,T,C,G]*', '', filter_positions)

chr <- NULL

pos <- NULL

for (i in seq(length(filter_positions))) {

  row <- filter_positions[i]
  
  chr <- c(chr, strsplit(row, ':')[[1]][1])
  
  pos <- c(pos, as.numeric(strsplit(row, ':')[[1]][2]))
  
    
}

bed <- data.frame(chr = chr, start = pos, end = pos+1)

write.table(bed, 'filtered_variants.bed', quote = F, row.names = F, col.names = F, sep = '\t')

### Enrichment of mutations in peaks that change

## We have a total of 15943 peaks that change

## We have a total of 8207 peaks that have a mutation inside

## We have a total of 792 peaks that change AND have a mutation

## we have a total of 556282 peaks

binom.test(792, 15943, p = (8207/556282), alternative = 'greater')

library(ggplot2)

data <- data.frame( Group = c('All peaks','Differentially regulated peaks') ,Percentage = c(1.475331,4.967697))

ggplot(data, aes(Group, Percentage, fill=Group )) + geom_bar(stat = 'Identity', width = 0.4) + theme_minimal()

ggsave('Mutations_barplot.pdf', scale = 2)
  
g + geom_bar( stat = 'Percentage')

###### Analysis of coverage #######

## File genrated using 
# samtools depth -b ~/cavefish/ATAC/new_genome//analysis_cave_vs_surface/results/filtered_variants.bed 
#  ~/cavefish/ATAC/new_genome/raw_files/*/*bam > Coverage_of_filtered_variants.txt

coverage <- read.delim('~/cavefish/ATAC/new_genome/raw_files/Coverage_of_filtered_variants.txt',
                       header = F, stringsAsFactors = F)

coverage_PA <- apply(coverage[,c(3:10)], 1, mean)


coverage_SF <- apply(coverage[,c(11:19)], 1, mean)

positions <- paste(coverage[,1], coverage[,2], sep = '-')

predata <- c(coverage_PA, coverage_SF)

data <- data.frame(Coverage = predata, pos = rep(1:length(coverage_PA),2), 
                   type =  as.factor(rep(c("Cave", "Surface"), each = length(coverage_PA))))

library(ggplot2)

ggplot(data = data, aes(x=pos, y = Coverage, col = type)) + geom_smooth()

ggplot(data = data, aes(x=type, y = Coverage, fill = type)) + geom_boxplot()

ggplot(data = data) + geom_histogram(aes(x=Coverage, fill = type), bins = 200) + 
  geom_vline(xintercept = 40, lty = 2,lwd = 1, alpha = 0.5) + theme_bw() + 
  geom_vline(xintercept = mean(coverage_PA), lty = 2,lwd = 1, alpha = 0.5, colour = 'red') + 
  geom_vline(xintercept = mean(coverage_SF), lty = 2,lwd = 1, alpha = 0.5, colour = 'blue')

qqnorm(coverage_PA)
qqline(coverage_PA,col='red')

qqnorm(coverage_SF)
qqline(coverage_SF,col='red')

wilcox.test(coverage_PA, coverage_SF, alternative = 'less')

######### Distribution of the mutations inside ALL the ATAC peaks #######

library(ggplot2)

setwd('~/cavefish/ATAC/new_genome/analysis_cave_vs_surface/')

## loading the data

data <- read.delim('data/Mut_in_all_ATAC_peaks.bed', header = F, stringsAsFactors = F)

data_trim <- data[,c(2:3,5:6)]


# computing distance of the mutation to the center of the peak, in percentage

distance <- apply(data_trim, 1, function(x){
  # where is the center of the atac peak?
  
  center_peak <- round((x[3] + x[4])/2)
  
  # distance to the center of the peak of the mutation, in percentage
  
  distance2center <- ((center_peak - x[1]) / ((x[4] - x[3])/2))*100
  
})

ggdata <- data.frame(Distance = distance)

ggplot(ggdata, aes(x = ggdata$Distance ,y = (..count..)*100/sum(..count..))) +
  geom_histogram(bins = 200, fill = 'royalblue1') + theme_minimal() +
  ggtitle('Relative Position of Mutations') + theme(plot.title = element_text(hjust = 0.5)) +
  xlab('Relative distance to center') + ylab('Frequency') +
  geom_vline(aes(xintercept=0), lwd = 1, lty = 2, col = 'grey30')+ 
  theme(plot.title = element_text(size = 22, face = 'bold'),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=16, face = 'bold'), 
        legend.text=element_text(size=14))

ggsave('results/Position_of_mutations_in_all_ATAC.pdf')

######### Distribution of the mutations inside ALL the ATAC peaks #######

library(ggplot2)

setwd('~/cavefish/ATAC/new_genome/analysis_cave_vs_surface/')

## loading the data

data <- read.delim('results/Mut_in_diff_peaks.bed', header = F, stringsAsFactors = F)

data_trim <- data[,c(2:3,5:6)]

# computing distance of the mutation to the center of the peak, in percentage

distance <- apply(data_trim, 1, function(x){
  # where is the center of the atac peak?
  
  center_peak <- round((x[3] + x[4])/2)
  
  # distance to the center of the peak of the mutation, in percentage
  
  distance2center <- ((center_peak - x[1]) / ((x[4] - x[3])/2))*100
  
})

ggdata <- data.frame(Distance = distance)

ggplot(ggdata, aes(x = ggdata$Distance ,y = (..count..)*100/sum(..count..))) +
  geom_histogram(bins = 200, fill = 'royalblue1') + theme_minimal() +
  ggtitle('Relative Position of Mutations \n in Differential ATAC peaks') + theme(plot.title = element_text(hjust = 0.5)) +
  xlab('Relative distance to center') + ylab('Frequency') +
  geom_vline(aes(xintercept=0), lwd = 1, lty = 2, col = 'grey30')+ 
  theme(plot.title = element_text(size = 22, face = 'bold'),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=16, face = 'bold'), 
        legend.text=element_text(size=14))

ggsave('results/Position_of_mutations_in_all_ATAC.pdf')


##### distribution of muts in upregulated and downregulated peaks ####

## downregulated genes

library(ggplot2)

setwd('~/cavefish/ATAC/new_genome/analysis_cave_vs_surface/')

## loading the data

data_down <- read.delim('results/Mut_in_DOWN_diff_peaks.bed', header = F, stringsAsFactors = F)

data_down_trim <- data_down[,c(2:3,5:6)]

# computing distance of the mutation to the center of the peak, in percentage

distance_down <- apply(data_down_trim, 1, function(x){
  # where is the center of the atac peak?
  
  center_peak <- round((x[3] + x[4])/2)
  
  # distance to the center of the peak of the mutation, in percentage
  
  distance2center <- ((center_peak - x[1]) / ((x[4] - x[3])/2))*100
  
})

### upregulated data

## loading the data

data_up<- read.delim('results/Mut_in_UP_diff_peaks.bed', header = F, stringsAsFactors = F)

data_up_trim <- data_up[,c(2:3,5:6)]

# computing distance of the mutation to the center of the peak, in percentage

distance_up <- apply(data_up_trim, 1, function(x){
  # where is the center of the atac peak?
  
  center_peak <- round((x[3] + x[4])/2)
  
  # distance to the center of the peak of the mutation, in percentage
  
  distance2center <- ((center_peak - x[1]) / ((x[4] - x[3])/2))*100
  
})

#### Non diff peaks

data <- read.delim('results/Mut_in_NO_diff_peaks.bed', header = F, stringsAsFactors = F)

data_trim <- data[,c(2:3,5:6)]

# computing distance of the mutation to the center of the peak, in percentage

distance <- apply(data_trim, 1, function(x){
  # where is the center of the atac peak?
  
  center_peak <- round((x[3] + x[4])/2)
  
  # distance to the center of the peak of the mutation, in percentage
  
  distance2center <- ((center_peak - x[1]) / ((x[4] - x[3])/2))*100
  
})


ggdata <- data.frame(Distance = c(distance_down, distance_up, distance), Group = factor(c(rep('Downregulated', length(distance_down)),
                                                                         rep('Upregulated', length(distance_up)),
                                                                         rep('Non Differential', length(distance))),
                                                                         levels = c('Downregulated',
                                                                                    'Upregulated',
                                                                                    'Non Differential')))

ggplot(ggdata, aes(x = ggdata$Distance ,y = (..count..*100) / sum(..count..), color = Group)) +
  geom_freqpoly(bins = 100, lwd = 2) + theme_minimal() +
  ggtitle('Relative Position of Mutations \n in Differential ATAC peaks') + theme(plot.title = element_text(hjust = 0.5)) +
  xlab('Relative distance to center') + ylab('Frequency') +
  geom_vline(aes(xintercept=0), lwd = 1, lty = 2, col = 'grey30')+ 
  theme(plot.title = element_text(size = 22, face = 'bold'),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=16, face = 'bold'), 
        legend.text=element_text(size=14)) 

ggplot(ggdata, aes(x = ggdata$Distance , y = ..density.., color = Group, fill = Group)) +
  geom_histogram(bins = 200) + theme_minimal() +
  ggtitle('Relative Position of Mutations \n in Differential ATAC peaks') + theme(plot.title = element_text(hjust = 0.5)) +
  xlab('Relative distance to center') +
  geom_vline(aes(xintercept=0), lwd = 1, lty = 2, col = 'grey30')+ 
  theme(plot.title = element_text(size = 22, face = 'bold'),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=16, face = 'bold'), 
        legend.text=element_text(size=14)) 

  # ggplot(ggdata, aes(x = ggdata$Distance, color = Group, fill = Group)) +
  # geom_density( alpha = 0.2) + theme_minimal() +
  # ggtitle('Relative Position of Mutations \n in Differential ATAC peaks') + theme(plot.title = element_text(hjust = 0.5)) +
  # xlab('Relative distance to center') +
  # geom_vline(aes(xintercept=0), lwd = 1, lty = 2, col = 'grey30')+ 
  # theme(plot.title = element_text(size = 22, face = 'bold'),
  #       axis.text=element_text(size=14),
  #       axis.title=element_text(size=18,face="bold"),
  #       legend.title=element_text(size=16, face = 'bold'), 
  #       legend.text=element_text(size=14)) 

##### Distribution of mutations in NO diff_peaks #########

data <- read.delim('results/Mut_in_NO_diff_peaks.bed', header = F, stringsAsFactors = F)

data_trim <- data[,c(2:3,5:6)]

# computing distance of the mutation to the center of the peak, in percentage

distance <- apply(data_trim, 1, function(x){
  # where is the center of the atac peak?
  
  center_peak <- round((x[3] + x[4])/2)
  
  # distance to the center of the peak of the mutation, in percentage
  
  distance2center <- ((center_peak - x[1]) / ((x[4] - x[3])/2))*100
  
})

ggdata <- data.frame(Distance = distance)

ggplot(ggdata, aes(x = ggdata$Distance ,y = (..count..)*100/sum(..count..))) +
  geom_histogram(bins = 200, fill = 'royalblue1') + theme_minimal() +
  ggtitle('Relative Position of Mutations \n in Non differential ATAC peaks') + theme(plot.title = element_text(hjust = 0.5)) +
  xlab('Relative distance to center') + ylab('Percentage') +
  geom_vline(aes(xintercept=0), lwd = 1, lty = 2, col = 'grey30')+ 
  theme(plot.title = element_text(size = 22, face = 'bold'),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=16, face = 'bold'), 
        legend.text=element_text(size=14))

#### KLF5 mut possition in mut diff peaks #####

library(ggplot2)

data <-  read.delim('results/motifs/Scan_of_mut_TFBS/klf5_100bp_diff_mut_peak.out', header = T, stringsAsFactors = F)

hist(data$Offset, breaks = 100)

