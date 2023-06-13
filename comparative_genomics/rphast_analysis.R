
# Rphast estimation of conserved elements in teleost species ####

#### library loading ####

require(rphast)
require(ape)

### Loading the alignment, the tree and the model ####

tree <- read.tree("data/tree_single.nw")

plot(tree)

tree_wo_sp_txt <- "((Oryzias_latipes:1,Esox_lucius:1)Euteleosteomorpha:1,(Clupea_harengus:1,((Danio_rerio:1,Sinocyclocheilus_grahami:1)Cyprinoidei:1,((Surface_Astyanax:1,(Colossoma_macropomum:1,Pygocentrus_nattereri:1)Serrasalmidae:1)Characoidei:1,(Electrophorus_electricus:1,(Pangasianodon_hypophthalmus:1,Ictalurus_punctatus:1)Siluroidei:1)Anc01:1)Characiphysae:1)Otophysi:1)Otomorpha:1)Clupeocephala;"
   
### Loading the CDS features:

chr_feats <- read.feat("data/Astyanax_mexicanus.Astyanax_mexicanus-2.0.103_mod.gtf")

chr_feats <-  add.introns.feat(chr_feats)

## For loop that generates the most cons regions of Surface astyanax

CDS_coverage <- NULL
chrom_coverage <- NULL

for (chr_maf in list.files(pattern = "no_cave.ss$",path = "data/mafs_chr_v2", full.names = F)) {
   
   print(chr_maf)
   
   chr_name <- gsub(".maf_no_cave.ss", "", chr_maf)
   
# loading the alignment
# chr19_surface.maf is a maf file which DONT have cavefish. We want this according to: http://compgen.cshl.edu/rphast/vignette2.pdf

align <- read.msa(list.files(path = "data/mafs_chr_v2", pattern = chr_maf, full.names = T),
                  ordered = T,format = "SS")


chr_feats_filtered <- chr_feats[chr_feats$seqname == chr_name,]

chr_feats_filtered$seqname <- "Surface_Astyanax"

wholeChrom <- feat(seqname =  "Surface_Astyanax", src=".", feature = "all", 
                   start = align$offset, end = align$offset+ncol.msa(align, "Surface_Astyanax"))

intergenicFeats <- inverse.feat(chr_feats_filtered, region.bounds=wholeChrom)

intergenicFeats$feature <- "intergenic"

chr_feats_filtered <- rbind.feat(chr_feats_filtered, intergenicFeats)

# loading model from file

if (file.exists(paste0("data/mafs_chr_v2/", chr_name, "_model_wo_cave.mod"))) {
   neutralMod <- read.tm(paste0("data/mafs_chr_v2/", chr_name, "_model_wo_cave.mod"))
   
} else {
   
   ### Loading only those features that are in our chromosome
   
   message("Getting 4d degenarated sites & model. This will take a while... ¯\\_(ツ)_/¯")
   
   table(chr_feats_filtered$feature)
   
   align4d <- get4d.msa(align, chr_feats_filtered)
   
   # generating the model
   
   neutralMod <- phyloFit(align4d, tree = tree_wo_sp_txt, subst.mod="REV")
   
   write.tm(neutralMod, paste0("data/mafs_chr_v2/", chr_name, "_model_wo_cave.mod"))
   
}


#using phastcons to estimate the most conserved elements of the alignment. The expected length and coverage have to be tuned.
# we have to look for a method to do it properly
# We can start with these values taken from https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1267528663_lo7zo4OFXa48Yody0aUhmi4hRiVT&db=fr3&c=chr21&g=cons8way
# --expected-length=45 --target-coverage=0.3 and --rho=0.3 

message("Running PhastCons...")

pc <- phastCons(align, neutralMod, viterbi = T, expected.length=45, target.coverage=0.3, rho = 0.3)

# pcEM <- phastCons(align, model, viterbi = T, estimate.transitions = T)

# most conserved elements
consElements <- pc$most.conserved

message("Done! Writting results.")

# this shows how many bases are predicted to be conserved
# coverage.feat(consElements)

# this shows the percentage of bases that are conserved. Since its a close tree, kind of makes sense

# coverage.feat(consElements)/coverage.feat(wholeChrom) # 38% of the chr is conserved. Makes kind of sense
chr_feats_cds <- chr_feats_filtered[chr_feats_filtered$feature == "CDS",]

CDS_coverage <- c(CDS_coverage,
                  coverage.feat(chr_feats_cds, consElements) * 100 /coverage.feat(chr_feats_cds)) # 60%% of covering of CDS is the target

chrom_coverage <- c(chrom_coverage, (coverage.feat(consElements)/coverage.feat(wholeChrom))*100)


cons_elements_bed <- consElements[,c("seqname","start","end","attribute","score","strand")]

cons_elements_bed$seqname <- chr_name

cons_elements_bed$attribute <- paste(cons_elements_bed$seqname, cons_elements_bed$start, cons_elements_bed$end, sep = "_")

write.table(cons_elements_bed, paste0("results/", format(Sys.Date(), format="%Y_%m_%d"),
                                      "_most_cons_regions_wo_Cave_raw.bed"), append = T,
            sep = "\t", quote = F, col.names = F, row.names = F)

}

stats <- data.frame(chr = list.files(pattern = "*maf$",path = "data/mafs_chr_v2", full.names = F),
                    CDS_coverage = CDS_coverage)

write.table(stats, "results/PhyloCons_cov_stats.txt", sep = "\t", row.names = F, quote = F)

 #### Now lets start with the accelerated element analysis ######

#see http://compgen.cshl.edu/rphast/vignette2.pdf

coverage_of_cons_by_inf_regions <- NULL

consElements <- read.delim("results/2022_07_19_most_cons_regions_wo_Cavefish_filtered_100bp_NonCoding.bed", header = F)


for (chr_maf in list.files(pattern = "*maf$",path = "data/mafs_chr_v2", full.names = F)) { # for loop
   
   print(chr_maf)

   gc() # clean the memory
   
   model <- rphast::read.tm(list.files(pattern = paste0(gsub("\\.maf", "", chr_maf),"_neutral_model.mod"),
                                       path = "data/mafs_chr_v2", full.names = T))
   
   
align_complete <- strip.gaps.msa(read.msa(list.files(path = "data/mafs_chr_v2", 
                                                     pattern = paste0(chr_maf, "_complete_aln.ss"), full.names = T),
                                          ordered = T,format = "SS"))
#unlink("data/AstMex_aln_chr19.maf")

# positions that dont have missing information in cavefish

hasCavefish <- informative.regions.msa(align_complete,1,spec = "Cave_Astyanax")

#### positions that have at least 4 species

hasAtLeastFour <- informative.regions.msa(align_complete, 5,
                                          spec = c("Surface_Astyanax", "Pygocentrus_nattereri", "Colossoma_macropomum",
                                                   "Electrophorus_electricus", "Ictalurus_punctatus", "Pangasianodon_hypophthalmus"))


### Overlap of these two set of regions and the conserved ones

consElements_filt <- consElements[consElements$V1 == gsub("\\.maf", "", chr_maf),]

consElements_filt <- feat(seqname = "Surface_Astyanax",
                          start = consElements_filt$V2,
                          end = consElements_filt$V3)

informativeElements <- coverage.feat(consElements_filt, hasCavefish, hasAtLeastFour, get.feats = T)

coverage_of_cons_by_inf_regions <- c(coverage_of_cons_by_inf_regions, 
                                      coverage.feat(informativeElements)*100/coverage.feat(consElements_filt))

# 90% of our most conserved sites have passed this filter, we can tighten it up, but we can do it later

# Regularize the size of cons elements to facilitate computation

splitLength <- 50

splitElements <- split.feat(informativeElements, f=splitLength, drop = T)

# now we run the phyloP program for computing the likelihood ratio for every feature and aproximates p-values

obsPhyloP <- phyloP(model, msa = align_complete, mode = "ACC", features = splitElements, subtree = "Cave_Astyanax")

# we need now to compute empirical pvalues. In order to do so we will generate random alignmets for features and compare
# the random dist vs the observed (simulating a lot of alignments)

elementAlign <- extract.feature.msa(copy.msa(align_complete), informativeElements)

nrep <- 100000 # this is the recomended in the vignete

simMsa <- sample.msa(elementAlign, nrep*splitLength, replace=TRUE)

# produce features allowing long alignment to be interpreted as 
# concatenation of shorter alignments

startIdx <- seq(from=1, by=splitLength, length.out=nrep)

features <- feat(seqname=names.msa(simMsa)[1], src="sim", feat=".",
                    start=startIdx,
                    end=startIdx+splitLength-1)
 
nonParaPhyloP <- phyloP(model, msa=simMsa, mode="ACC",
                           features=features, subtree="Cave_Astyanax")
### Okey now lets compare the two models
png(paste0("results/plots_stats/Cave_acc_nonCoding_", gsub("\\.maf", "", chr_maf),"_Stat_plots.png"),width = 500, height = 500)
par(mfrow=c(1,2), mar=c(5,2,4,2))

qqplot(nonParaPhyloP$lnlratio,obsPhyloP$lnlratio,
       xlim=c(0,15),ylim=c(0,15), xlab="Simulated likelihood ratio",
       ylab="Observed likelihood ratio")
abline(0, 1, lty=2)

# yeah well, We can see that they are different, at least :/
plot(density(obsPhyloP$lnlratio,adjust=3), lty=1,xlim=c(0,1),
      xlab="Likelihood Ratio",
      ylab="Density",main="", col="red")
lines(density(nonParaPhyloP$lnlratio,adjust=3), lty=1,
         col="black",xlim=c(0,1))

dev.off()
# computing of the empirical pval

empirical.pval <- function(x, dist) {
   sum(x <= dist)/length(dist)
}

nonParaPval <- sapply(obsPhyloP$lnlratio, empirical.pval,
                         nonParaPhyloP$lnlratio)
nonParaFDR <- p.adjust(nonParaPval, method="bonferroni")

# Extracting those that are significant from our set of accelerated elements

nonParaSigFeats <- splitElements[nonParaFDR < 0.05,]
nrow(nonParaSigFeats)

class(nonParaSigFeats)

nonParaSigFeats$feature <- "CaveAccRegions"

nonParaSigFeats$score <- obsPhyloP$lnlratio[nonParaFDR < 0.05]



# write.feat(nonParaSigFeats, "CaveAccRegions.gff")

## Sanity check that we have actually accelerated regions 

#### TURNED OFF: Too much time running, we have seen that the branch length is increased greatly in other runs

# consEleModel <- phyloFit(elementAlign, init.mod=model, no.opt=c("backgd", "ratematrix"))
# 
# CaveAccRegionsAln <- extract.feature.msa(copy.msa(align_complete), nonParaSigFeats)
# 
# AccModel <- phyloFit(CaveAccRegionsAln, init.mod=model, no.opt=c("backgd", "ratematrix"))
# 
# maxX <- depth.tree(AccModel$tree, "Cave_Astyanax")
# 
# pdf(paste0("results/", gsub("\\.maf", "", chr_maf),"_Stat_plots.pdf"))
# par(mfrow=c(1,2))
# 
# plot.phylo(read.tree(text=consEleModel$tree), x.lim=c(0, maxX),
#             main=paste0("Conserved Elements of ",  gsub("\\.maf", "", chr_maf)))
# 
# plot.phylo(read.tree(text=AccModel$tree), x.lim=c(0, maxX),
#            main=paste0("Accelerated Elements of ",  gsub("\\.maf", "", chr_maf)))
# dev.off()

nonParaSigFeats$seqname <- gsub("\\.maf", "", chr_maf) # this goes here in order to not disrupt the tree computation above



Acc_elements_bed <- data.frame(chr = nonParaSigFeats$seqname,
                               start = nonParaSigFeats$start,
                               end = nonParaSigFeats$end,
                               name = paste0("Acc_region_",gsub("\\.maf", "", chr_maf),"_", rep(1:length(nonParaSigFeats$seqname))),
                               score = nonParaSigFeats$score,
                               strand = "+")

# Acc_elements_bed$attribute <- paste(Acc_elements_bed$seqname, Acc_elements_bed$start, Acc_elements_bed$end, sep = "_")

write.table(Acc_elements_bed, paste0("results/", format(Sys.Date(), format="%Y_%m_%d"),
                                     "_Acc_non_coding_regions_of_cave.bed"), append = T,
            sep = "\t", quote = F, col.names = F, row.names = F)

}


########## Analysis of Astyanax accelerated regions ##########

## For loop that generates the most cons regions of Teleost without Astyanax mexicanus

## We follow this tutorial for doing this: http://compgen.cshl.edu/rphast/vignette1.pdf



feats <- read.feat("data/Danrer10_UCSC_ensgenes.gtf")

# we add the introns

feats <- add.introns.feat(feats)


# feats <- feats[!(feats$feature %in% c("start_codon, stop_codon")) ,]



tree_wo_ast <- collapse.singles(read.tree(file = "data/tree_wo_Astyanax.nw"))

# write.tree(tree_wo_ast, "data/tree_wo_Astyanax.nw")
tree_wo_ast_txt <- "((Clupea_harengus:1,((Sinocyclocheilus_grahami:1,Danio_rerio:1)Cyprinoidei:1,((Pygocentrus_nattereri:1,Colossoma_macropomum:1):NaN,(Electrophorus_electricus:1,(Ictalurus_punctatus:1,Pangasianodon_hypophthalmus:1)Siluroidei:1)Anc01:1)Characiphysae:1)Otophysi:1)Otomorpha:1,(Esox_lucius:1,Oryzias_latipes:1)Euteleosteomorpha:1)Clupeocephala;"

CDS_coverage <- NULL
chrom_coverage <- NULL

for (chr_maf in list.files(pattern = "*_wo_Ast.ss$",path = "data/mafs_chr_Danrer10_ref", full.names = F)) {
   
   print(chr_maf)
   
   chr_name <- gsub("_danrer10_reference_fixed_wo_Ast.ss", "", chr_maf)
   
   # loading the alignment
   # chr1.maf_wo_Ast.ss is an alignmentfile which DONT have Astyanax We want this according to: http://compgen.cshl.edu/rphast/vignette2.pdf
   
   align <- read.msa(list.files(path = "data/mafs_chr_Danrer10_ref", pattern = chr_maf, full.names = T),
                     ordered = T,format = "SS")
   
   ## If loop that will try to find the phylomodel, but if not, it will generate it
   
   if (file.exists(paste0("data/mafs_chr_Danrer10_ref/", chr_name, "_model_wo_Astyanax.mod"))) {
      neutralMod <- read.tm(paste0("data/mafs_chr_Danrer10_ref/", chr_name, "_model_wo_Astyanax.mod"))
      
   } else {
   
   ### Loading only those features that are in our chromosome
      
   message("Getting 4d degenarated sites & model. This will take a while... ¯\\_(ツ)_/¯")
   
   chr_feats <- feats[feats$seqname == chr_name,]
   
   chr_feats$seqname <- "Danio_rerio"
   
   wholeChrom <- feat(seqname =  "Danio_rerio", src=".", feature = "all", 
                      start = align$offset, end = align$offset+ncol.msa(align, "Danio_rerio"))
   
   intergenicFeats <- inverse.feat(chr_feats, region.bounds=wholeChrom)
   
   intergenicFeats$feature <- "intergenic"
   
   chr_feats <- rbind.feat(chr_feats, intergenicFeats)
   
   table(chr_feats$feature)
   
   align4d <- get4d.msa(align, chr_feats)
   
   # generating the model
   
   neutralMod <- phyloFit(align4d, tree = tree_wo_ast_txt, subst.mod="REV")

   write.tm(neutralMod, paste0("data/mafs_chr_Danrer10_ref/", chr_name, "_model_wo_Astyanax.mod"))
   
   }
   # we want to filter only those features that are from the same chr than the alignment
   chr_feats <- feats[feats$seqname == chr_name,]
   
   chr_feats$seqname <- "Danio_rerio"
   
   wholeChrom <- feat(seqname =  "Danio_rerio", src=".", feature = "all", 
                      start = align$offset, end = align$offset+ncol.msa(align, "Danio_rerio"))
   
   intergenicFeats <- inverse.feat(chr_feats, region.bounds=wholeChrom)
   
   intergenicFeats$feature <- "intergenic"
   
   chr_feats <- rbind.feat(chr_feats, intergenicFeats)
   
   #using phastcons to estimate the most conserved elements of the alignment. The expected length and coverage have to be tuned.
   # we have to look for a method to do it properly
   # We can start with these values taken from https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1267528663_lo7zo4OFXa48Yody0aUhmi4hRiVT&db=fr3&c=chr21&g=cons8way
   # --expected-length=45 --target-coverage=0.3 and --rho=0.3 
   
   pc <- phastCons(align, neutralMod, viterbi = T, expected.length=45, target.coverage=0.3, rho = 0.3)
   
   # pcEM <- phastCons(align, model, viterbi = T, estimate.transitions = T)
   
   # most conserved elements
   consElements <- pc$most.conserved
   
   # this shows how many bases are predicted to be conserved
   # coverage.feat(consElements)
   
   # this shows the percentage of bases that are conserved. Since its a close tree, kind of makes sense
   
   chrom_coverage <- c(chrom_coverage, (coverage.feat(consElements)/coverage.feat(wholeChrom))*100) # 38% of the chr is conserved. Makes kind of sense
   
   CDS_coverage <- c(CDS_coverage,
                     coverage.feat(chr_feats[chr_feats$feature == "CDS",], consElements) * 100 /coverage.feat(chr_feats[chr_feats$feature == "CDS",])) # 60%% of covering of CDS is the target
   
   cons_elements_bed <- consElements[,c("seqname","start","end","attribute","score","strand")]
   
   cons_elements_bed$seqname <- chr_name
   
   cons_elements_bed$attribute <- paste(cons_elements_bed$seqname, cons_elements_bed$start, cons_elements_bed$end, sep = "_")
   
   write.table(cons_elements_bed, "results/2022_07_07_most_cons_regions_wo_Astyanax_danrer10_coord_v2.bed", append = T,
               sep = "\t", quote = F, col.names = F, row.names = F)
   
}

stats <- data.frame(chr = list.files(pattern = "*maf$",path = "data/mafs_chr_v2", full.names = F),
                    CDS_coverage = CDS_coverage,
                    Chr_coverage = chrom_coverage)

write.table(stats, "results/PhyloCons_Cons_regions_WO_Ast_stats_v2.txt", sep = "\t", row.names = F, quote = F)

#### Now lets start with the accelerated element analysis ######

#see http://compgen.cshl.edu/rphast/vignette2.pdf

tree <- "((Clupea_harengus:1,((Sinocyclocheilus_grahami:1,Danio_rerio:1)Cyprinoidei:1,(((Pygocentrus_nattereri:1,Colossoma_macropomum:1)Serrasalmidae:1,(Surface_Astyanax:1,Cave_Astyanax:1)Anc02:1)Characoidei:1,(Electrophorus_electricus:1,(Ictalurus_punctatus:1,Pangasianodon_hypophthalmus:1)Siluroidei:1)Anc01:1)Characiphysae:1)Otophysi:1)Otomorpha:1,(Esox_lucius:1,Oryzias_latipes:1)Euteleosteomorpha:1)Clupeocephala;"

coverage_of_cons_by_inf_regions <- NULL

consElements <- read.delim("results/2022_07_07_most_cons_regions_wo_Astyanax_danrer10_coord_v2.bed", header = F)

consElements <- consElements[(consElements$V3 - consElements$V2)>100,]
# we include this again to compute the models with the complete alignment

feats <- read.feat("data/Danrer10_UCSC_ensgenes.gtf")

# we add the introns

feats <- add.introns.feat(feats)


# feats <- feats[!(feats$feature %in% c("start_codon, stop_codon")) ,]

for (chr_maf in list.files(pattern = "*complete_aln.ss",path = "data/mafs_chr_Danrer10_ref", full.names = F)) { # for loop
   
   print(chr_maf)
   
   chr_name <- gsub("_danrer10_reference_fixed_complete_aln.ss", "", chr_maf)
   
   gc() # clean the memory
   
   align_complete <- strip.gaps.msa(read.msa(list.files(path = "data/mafs_chr_Danrer10_ref", 
                                                        pattern = chr_maf,
                                                        full.names = T),
                                             ordered = T,format = "SS"))
   
   # if loop that will try to look for model file, if not, it will generate it
   
   if (file.exists(paste0("data/mafs_chr_Danrer10_ref/",chr_name, "_danrer10_ref_complete_aln.mod"))) {
      
      model <- read.tm(paste0("data/mafs_chr_Danrer10_ref/",chr_name, "_danrer10_ref_complete_aln.mod"))
                              
   } else {
   
      
      message("Getting 4d degenarated sites & model. This will take a while... ¯\\_(ツ)_/¯")
      
   ### we need to generate a model, since we dont have the model of the complete aln of the chromosome
   
   ### Loading only those features that are in our chromosome
   
   chr_feats <- feats[feats$seqname == chr_name,]
   
   chr_feats$seqname <- "Danio_rerio"
   
   wholeChrom <- feat(seqname =  "Danio_rerio", src=".", feature = "all", 
                      start = align_complete$offset, end = align_complete$offset+ncol.msa(align_complete, "Danio_rerio"))
   
   intergenicFeats <- inverse.feat(chr_feats, region.bounds=wholeChrom)
   
   intergenicFeats$feature <- "intergenic"
   
   chr_feats <- rbind.feat(chr_feats, intergenicFeats)
   
   table(chr_feats$feature)
   
  
   
   align4d <- get4d.msa(align_complete, chr_feats)
   
   # generating the model
   
   model <- phyloFit(align4d, tree = tree, subst.mod="REV")
   
   write.tm(model, paste0("data/mafs_chr_Danrer10_ref/",chr_name, "_danrer10_ref_complete_aln.mod"))
   }
   
   #unlink("data/AstMex_aln_chr19.maf")
   
   # positions that dont have missing information in cavefish
   
   hasCavefish <- informative.regions.msa(align_complete,1,spec = "Cave_Astyanax")
   
   #### positions that have at least 4 species
   
   hasAtLeastFour <- informative.regions.msa(align_complete, 5,
                                             spec = c("Surface_Astyanax", "Pygocentrus_nattereri", "Colossoma_macropomum",
                                                      "Electrophorus_electricus", "Ictalurus_punctatus", "Pangasianodon_hypophthalmus"))
   
   
   ### Overlap of these two set of regions and the conserved ones
   
   consElements_filt <- consElements[consElements$V1 == chr_name,]
   
   consElements_filt <- feat(seqname = "Danio_rerio",
                             start = consElements_filt$V2,
                             end = consElements_filt$V3)
   
   informativeElements <- coverage.feat(consElements_filt, hasCavefish, hasAtLeastFour, get.feats = T)
   
   coverage_of_cons_by_inf_regions <- c(coverage_of_cons_by_inf_regions, 
                                        coverage.feat(informativeElements)*100/coverage.feat(consElements_filt))
   
   # 90% of our most conserved sites have passed this filter, we can tighten it up, but we can do it later
   
   # Regularize the size of cons elements to facilitate computation
   
   splitLength <- 50
   
   splitElements <- split.feat(informativeElements, f=splitLength, drop = T)
   
   # now we run the phyloP program for computing the likelihood ratio for every feature and aproximates p-values
   
   obsPhyloP <- phyloP(model, msa = align_complete, mode = "ACC", features = splitElements, subtree = "Anc02") #Anc02 is the ancestor node of Astyanax
   
   # we need now to compute empirical pvalues. In order to do so we will generate random alignmets for features and compare
   # the random dist vs the observed (simulating a lot of alignments)
   
   elementAlign <- extract.feature.msa(copy.msa(align_complete), informativeElements)
   
   nrep <- 100000 # this is the recomended in the vignete
   
   simMsa <- sample.msa(elementAlign, nrep*splitLength, replace=TRUE)
   
   # produce features allowing long alignment to be interpreted as 
   # concatenation of shorter alignments
   
   startIdx <- seq(from=1, by=splitLength, length.out=nrep)
   
   features <- feat(seqname=names.msa(simMsa)[1], src="sim", feat=".",
                    start=startIdx,
                    end=startIdx+splitLength-1)
   
   nonParaPhyloP <- phyloP(model, msa=simMsa, mode="ACC",
                           features=features, subtree="Anc02")
   ### Okey now lets compare the two models
   png(paste0("results/plots_stats/", chr_name,"Acc_Astyanax_in_danrer10_v2_Stat_plots.png"),width = 500, height = 500)
   par(mfrow=c(1,2), mar=c(5,2,4,2))
   
   qqplot(nonParaPhyloP$lnlratio,obsPhyloP$lnlratio,
          xlim=c(0,15),ylim=c(0,15), xlab="Simulated likelihood ratio",
          ylab="Observed likelihood ratio")
   abline(0, 1, lty=2)
   
   # yeah well, We can see that they are different, at least :/
   plot(density(obsPhyloP$lnlratio,adjust=3), lty=1,xlim=c(0,1),
        xlab="Likelihood Ratio",
        ylab="Density",main="", col="red")
   lines(density(nonParaPhyloP$lnlratio,adjust=3), lty=1,
         col="black",xlim=c(0,1))
   
   dev.off()
   # computing of the empirical pval
   
   empirical.pval <- function(x, dist) {
      sum(x <= dist)/length(dist)
   }
   
   nonParaPval <- sapply(obsPhyloP$lnlratio, empirical.pval,
                         nonParaPhyloP$lnlratio)
   nonParaFDR <- p.adjust(nonParaPval, method="bonferroni")
   
   # Extracting those that are significant from our set of accelerated elements
   
   nonParaSigFeats <- splitElements[nonParaFDR < 0.05,]
   nrow(nonParaSigFeats)
   
   class(nonParaSigFeats)
   
   nonParaSigFeats$feature <- "AstyanaxAccRegions"
   
   nonParaSigFeats$score <- obsPhyloP$lnlratio[nonParaFDR < 0.05]
   
   
   
   # write.feat(nonParaSigFeats, "CaveAccRegions.gff")
   
   ## Sanity check that we have actually accelerated regions 
   
   #### TURNED OFF: Too much time running, we have seen that the branch length is increased greatly in other runs
   
   # consEleModel <- phyloFit(elementAlign, init.mod=model, no.opt=c("backgd", "ratematrix"))
   # 
   # CaveAccRegionsAln <- extract.feature.msa(copy.msa(align_complete), nonParaSigFeats)
   # 
   # AccModel <- phyloFit(CaveAccRegionsAln, init.mod=model, no.opt=c("backgd", "ratematrix"))
   # 
   # maxX <- depth.tree(AccModel$tree, "Cave_Astyanax")
   # 
   # pdf(paste0("results/", gsub("\\.maf", "", chr_maf),"_Stat_plots.pdf"))
   # par(mfrow=c(1,2))
   # 
   # plot.phylo(read.tree(text=consEleModel$tree), x.lim=c(0, maxX),
   #             main=paste0("Conserved Elements of ",  gsub("\\.maf", "", chr_maf)))
   # 
   # plot.phylo(read.tree(text=AccModel$tree), x.lim=c(0, maxX),
   #            main=paste0("Accelerated Elements of ",  gsub("\\.maf", "", chr_maf)))
   # dev.off()
   
   nonParaSigFeats$seqname <- chr_name# this goes here in order to not disrupt the tree computation above
   
   
   
   Acc_elements_bed <- data.frame(chr = nonParaSigFeats$seqname,
                                  start = nonParaSigFeats$start,
                                  end = nonParaSigFeats$end,
                                  name = paste0("Acc_region_",chr_name,"_", rep(1:length(nonParaSigFeats$seqname))),
                                  score = nonParaSigFeats$score,
                                  strand = "+")
   
   # Acc_elements_bed$attribute <- paste(Acc_elements_bed$seqname, Acc_elements_bed$start, Acc_elements_bed$end, sep = "_")
   
   write.table(Acc_elements_bed, "results/2022_07_08_Acc_regions_Astyanax_in_danrer10_coord_v2.bed", append = T,
               sep = "\t", quote = F, col.names = F, row.names = F)
   
}


########## Analysis of Ictalurus Accelerated regions ##########

## For loop that generates the most cons regions of Teleost without Astyanax mexicanus

## We follow this tutorial for doing this: http://compgen.cshl.edu/rphast/vignette1.pdf


### Conserved regions of Teleosts without Ictalurus ####
feats <- read.feat("data/Astyanax_mexicanus.Astyanax_mexicanus-2.0.103_mod.gtf")

# we add the introns

feats <- add.introns.feat(feats)

# feats <- feats[!(feats$feature %in% c("start_codon, stop_codon")) ,]

tree_wo_sp <- collapse.singles(read.tree(file = "data/tree_wo_Ipcoco.nw"))

# write.tree(tree_wo_sp, "data/tree_wo_spyanax.nw")
tree_wo_sp_txt <- "((Oryzias_latipes:1,Esox_lucius:1)Euteleosteomorpha:1,(Clupea_harengus:1,((Danio_rerio:1,Sinocyclocheilus_grahami:1)Cyprinoidei:1,((Pangasianodon_hypophthalmus:1,Electrophorus_electricus:1)Anc01:1,((Colossoma_macropomum:1,Pygocentrus_nattereri:1)Serrasalmidae:1,(Cave_Astyanax:1,Surface_Astyanax:1)Anc02:1)Characoidei:1)Characiphysae:1)Otophysi:1)Otomorpha:1)Clupeocephala;"

CDS_coverage <- NULL
chrom_coverage <- NULL

for (chr_maf in list.files(pattern = "*_wo_IpCoco.ss$",path = "data/mafs_chr_wo_IpCoco", full.names = F)) {
   
   print(chr_maf)
   
   chr_name <- gsub("_wo_IpCoco.ss", "", chr_maf)
   
   # loading the alignment
   # chr1.maf_wo_Ast.ss is an alignmentfile which DONT have Astyanax We want this according to: http://compgen.cshl.edu/rphast/vignette2.pdf
   
   align <- read.msa(list.files(path = "data/mafs_chr_wo_IpCoco", pattern = chr_maf, full.names = T),
                     ordered = T,format = "SS")
   
   ## If loop that will try to find the phylomodel, but if not, it will generate it
   
   if (file.exists(paste0("data/mafs_chr_wo_IpCoco/", chr_name, "_model_wo_IpCoco.mod"))) {
      neutralMod <- read.tm(paste0("data/mafs_chr_wo_IpCoco/", chr_name, "_model_wo_IpCoco.mod"))
      
   } else {
      
      ### Loading only those features that are in our chromosome
      
      message("Getting 4d degenarated sites & model. This will take a while... ¯\\_(ツ)_/¯")
      
      chr_feats <- feats[feats$seqname == chr_name,]
      
      chr_feats$seqname <- "Surface_Astyanax"
      
      wholeChrom <- feat(seqname =  "Surface_Astyanax", src=".", feature = "all", 
                         start = align$offset, end = align$offset+ncol.msa(align, "Surface_Astyanax"))
      
      intergenicFeats <- inverse.feat(chr_feats, region.bounds=wholeChrom)
      
      intergenicFeats$feature <- "intergenic"
      
      chr_feats <- rbind.feat(chr_feats, intergenicFeats)
      
      table(chr_feats$feature)
      
      align4d <- get4d.msa(align, chr_feats)
      
      # generating the model
      
      neutralMod <- phyloFit(align4d, tree = tree_wo_sp_txt, subst.mod="REV")
      
      write.tm(neutralMod, paste0("data/mafs_chr_wo_IpCoco/", chr_name, "_model_wo_IpCoco.mod"))
      
   }
   # we want to filter only those features that are from the same chr than the alignment
   chr_feats <- feats[feats$seqname == chr_name,]
   
   chr_feats$seqname <- "Surface_Astyanax"
   
   wholeChrom <- feat(seqname =  "Surface_Astyanax", src=".", feature = "all", 
                      start = align$offset, end = align$offset+ncol.msa(align, "Surface_Astyanax"))
   
   intergenicFeats <- inverse.feat(chr_feats, region.bounds=wholeChrom)
   
   intergenicFeats$feature <- "intergenic"
   
   chr_feats <- rbind.feat(chr_feats, intergenicFeats)
   
   #using phastcons to estimate the most conserved elements of the alignment. The expected length and coverage have to be tuned.
   # we have to look for a method to do it properly
   # We can start with these values taken from https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1267528663_lo7zo4OFXa48Yody0aUhmi4hRiVT&db=fr3&c=chr21&g=cons8way
   # --expected-length=45 --target-coverage=0.3 and --rho=0.3 
   
   pc <- phastCons(align, neutralMod, viterbi = T, expected.length=45, target.coverage=0.3, rho = 0.3)
   
   # pcEM <- phastCons(align, model, viterbi = T, estimate.transitions = T)
   
   # most conserved elements
   consElements <- pc$most.conserved
   
   # this shows how many bases are predicted to be conserved
   # coverage.feat(consElements)
   
   # this shows the percentage of bases that are conserved. Since its a close tree, kind of makes sense
   
   chrom_coverage <- c(chrom_coverage, (coverage.feat(consElements)/coverage.feat(wholeChrom))*100) # 38% of the chr is conserved. Makes kind of sense
   
   CDS_coverage <- c(CDS_coverage,
                     coverage.feat(chr_feats[chr_feats$feature == "CDS",], consElements) * 100 /coverage.feat(chr_feats[chr_feats$feature == "CDS",])) # 60%% of covering of CDS is the target
   
   cons_elements_bed <- consElements[,c("seqname","start","end","attribute","score","strand")]
   
   cons_elements_bed$seqname <- chr_name
   
   cons_elements_bed$attribute <- paste(cons_elements_bed$seqname, cons_elements_bed$start, cons_elements_bed$end, sep = "_")
   
   write.table(cons_elements_bed, paste0("results/", format(Sys.Date(), format="%Y_%m_%d"),
                                         "_most_cons_regions_wo_IpCoco.bed"), append = T,
               sep = "\t", quote = F, col.names = F, row.names = F)
   
}

stats <- data.frame(chr = gsub("_wo_IpCoco.ss", "", list.files(pattern = "*_wo_IpCoco.ss$",path = "data/mafs_chr_wo_IpCoco", full.names = F)),
                    CDS_coverage = CDS_coverage,
                    Chr_coverage = chrom_coverage)

write.table(stats, "results/PhyloCons_Cons_regions_WO_IpCoco_stats.txt", sep = "\t", row.names = F, quote = F)


#### Accelerated regions of IpCoco ######

#see http://compgen.cshl.edu/rphast/vignette2.pdf

coverage_of_cons_by_inf_regions <- NULL

consElements <- read.delim("results/2022_07_21_most_cons_regions_wo_IpCoco_filtered_100bp_NonCoding.bed", header = F)

for (chr_maf in list.files(pattern = "*maf$",path = "data/mafs_chr_v2", full.names = F)) { # for loop
   
   
   chr_name <- gsub("\\.maf", "", chr_maf)
   print(chr_name)
   
   gc() # clean the memory
   
   model <- rphast::read.tm(list.files(pattern = paste0(chr_name,"_neutral_model.mod"),
                                       path = "data/mafs_chr_v2", full.names = T))
   
   
   align_complete <- strip.gaps.msa(read.msa(list.files(path = "data/mafs_chr_v2", 
                                                        pattern = paste0(chr_maf, "_complete_aln.ss"), full.names = T),
                                             ordered = T,format = "SS"))
   #unlink("data/AstMex_aln_chr19.maf")
   
   # positions that dont have missing information in cavefish
   
   #### positions that have at least 4 species
   
   hasAtLeastFour <- informative.regions.msa(align_complete, 5,
                                             spec = c("Surface_Astyanax", "Pygocentrus_nattereri", "Colossoma_macropomum",
                                                      "Electrophorus_electricus", "Ictalurus_punctatus", "Pangasianodon_hypophthalmus"))
   
   
   ### Overlap of these two set of regions and the conserved ones
   
   consElements_filt <- consElements[consElements$V1 == chr_name,]
   
   consElements_filt <- feat(seqname = "Surface_Astyanax",
                             start = consElements_filt$V2,
                             end = consElements_filt$V3)
   
   informativeElements <- coverage.feat(consElements_filt, hasAtLeastFour, get.feats = T)
   
   coverage_of_cons_by_inf_regions <- c(coverage_of_cons_by_inf_regions, 
                                        coverage.feat(informativeElements)*100/coverage.feat(consElements_filt))
   
   # 90% of our most conserved sites have passed this filter, we can tighten it up, but we can do it later
   
   # Regularize the size of cons elements to facilitate computation
   
   splitLength <- 50
   
   splitElements <- split.feat(informativeElements, f=splitLength, drop = T)
   
   # now we run the phyloP program for computing the likelihood ratio for every feature and aproximates p-values
   
   message("Computing Observed PhyloP scores")
   
   obsPhyloP <- phyloP(model, msa = align_complete, mode = "ACC", features = splitElements, subtree = "Ictalurus_punctatus")
   
   # we need now to compute empirical pvalues. In order to do so we will generate random alignmets for features and compare
   # the random dist vs the observed (simulating a lot of alignments)
   
   elementAlign <- extract.feature.msa(align_complete, informativeElements)
   
   nrep <- 100000 # this is the recomended in the vignete
   
   simMsa <- sample.msa(elementAlign, nrep*splitLength, replace=TRUE)
   
   # produce features allowing long alignment to be interpreted as 
   # concatenation of shorter alignments
   
   startIdx <- seq(from=1, by=splitLength, length.out=nrep)
   
   features <- feat(seqname=names.msa(simMsa)[1], src="sim", feat=".",
                    start=startIdx,
                    end=startIdx+splitLength-1)
   message("Computing Simulated PhyloP")
   
   nonParaPhyloP <- phyloP(model, msa=simMsa, mode="ACC",
                           features=features, subtree="Ictalurus_punctatus")
   ### Okey now lets compare the two models
   png(paste0("results/plots_stats/", chr_name,"_Ictalurus_punctatus_NonCoding_Stat_plots.png"),width = 500, height = 500)
   par(mfrow=c(1,2), mar=c(5,2,4,2))
   
   qqplot(nonParaPhyloP$lnlratio,obsPhyloP$lnlratio,
          xlim=c(0,15),ylim=c(0,15), xlab="Simulated likelihood ratio",
          ylab="Observed likelihood ratio")
   abline(0, 1, lty=2)
   
   # yeah well, We can see that they are different, at least :/
   plot(density(obsPhyloP$lnlratio,adjust=3), lty=1,xlim=c(0,1),
        xlab="Likelihood Ratio",
        ylab="Density",main="", col="red")
   lines(density(nonParaPhyloP$lnlratio,adjust=3), lty=1,
         col="black",xlim=c(0,1))
   
   dev.off()
   # computing of the empirical pval
   
   empirical.pval <- function(x, dist) {
      sum(x <= dist)/length(dist)
   }
   
   message("Computing empirical p.values and writting results...")
   
   nonParaPval <- sapply(obsPhyloP$lnlratio, empirical.pval,
                         nonParaPhyloP$lnlratio)
   nonParaFDR <- p.adjust(nonParaPval, method="bonferroni")
   
   # Extracting those that are significant from our set of accelerated elements
   
   nonParaSigFeats <- splitElements[nonParaFDR < 0.05,]
   nrow(nonParaSigFeats)
   
   class(nonParaSigFeats)
   
   nonParaSigFeats$feature <- "IpCocoAccRegions"
   
   nonParaSigFeats$score <- obsPhyloP$lnlratio[nonParaFDR < 0.05]
   
   
   
   # write.feat(nonParaSigFeats, "CaveAccRegions.gff")
   
   ## Sanity check that we have actually accelerated regions 
   
   #### TURNED OFF: Too much time running, we have seen that the branch length is increased greatly in other runs
   
   # consEleModel <- phyloFit(elementAlign, init.mod=model, no.opt=c("backgd", "ratematrix"))
   # 
   # CaveAccRegionsAln <- extract.feature.msa(copy.msa(align_complete), nonParaSigFeats)
   # 
   # AccModel <- phyloFit(CaveAccRegionsAln, init.mod=model, no.opt=c("backgd", "ratematrix"))
   # 
   # maxX <- depth.tree(AccModel$tree, "Cave_Astyanax")
   # 
   # pdf(paste0("results/", chr_name,"_Stat_plots.pdf"))
   # par(mfrow=c(1,2))
   # 
   # plot.phylo(read.tree(text=consEleModel$tree), x.lim=c(0, maxX),
   #             main=paste0("Conserved Elements of ",  chr_name))
   # 
   # plot.phylo(read.tree(text=AccModel$tree), x.lim=c(0, maxX),
   #            main=paste0("Accelerated Elements of ",  chr_name))
   # dev.off()
   
   nonParaSigFeats$seqname <- chr_name # this goes here in order to not disrupt the tree computation above
   
   
   
   Acc_elements_bed <- data.frame(chr = nonParaSigFeats$seqname,
                                  start = nonParaSigFeats$start,
                                  end = nonParaSigFeats$end,
                                  name = paste0("Acc_region_",chr_name,"_", rep(1:length(nonParaSigFeats$seqname))),
                                  score = nonParaSigFeats$score,
                                  strand = "+")
   
   # Acc_elements_bed$attribute <- paste(Acc_elements_bed$seqname, Acc_elements_bed$start, Acc_elements_bed$end, sep = "_")
   
   write.table(Acc_elements_bed, "results/Acc_regions_Ictalurus_punctatus_NonCoding.bed", append = T,
               sep = "\t", quote = F, col.names = F, row.names = F)
   
}

########## Analysis of Piranha Accelerated regions ##########

## For loop that generates the most cons regions of Teleost without Astyanax mexicanus

## We follow this tutorial for doing this: http://compgen.cshl.edu/rphast/vignette1.pdf


### Conserved regions of Teleosts without Piranha ####
feats <- read.feat("data/Astyanax_mexicanus.Astyanax_mexicanus-2.0.103_mod.gtf")

# we add the introns

feats <- add.introns.feat(feats)

# feats <- feats[!(feats$feature %in% c("start_codon, stop_codon")) ,]



# write.tree(tree_wo_sp, "data/tree_wo_spyanax.nw")
tree_wo_sp_txt <- "((Oryzias_latipes:1,Esox_lucius:1)Euteleosteomorpha:1,(Clupea_harengus:1,((Danio_rerio:1,Sinocyclocheilus_grahami:1)Cyprinoidei:1,((Colossoma_macropomum:1,(Cave_Astyanax:1,Surface_Astyanax:1)Anc02:1)Characoidei:1,(Electrophorus_electricus:1,(Pangasianodon_hypophthalmus:1,Ictalurus_punctatus:1)Siluroidei:1)Anc01:1)Characiphysae:1)Otophysi:1)Otomorpha:1)Clupeocephala;"

CDS_coverage <- NULL
chrom_coverage <- NULL

for (chr_maf in list.files(pattern = "*_wo_PigNat.ss$",path = "data/mafs_chr_wo_PigNat", full.names = F)) {
   
   print(chr_maf)
   
   chr_name <- gsub("_wo_PigNat.ss", "", chr_maf)
   
   # loading the alignment
   # chr1.maf_wo_Ast.ss is an alignmentfile which DONT have Astyanax We want this according to: http://compgen.cshl.edu/rphast/vignette2.pdf
   
   align <- read.msa(list.files(path = "data/mafs_chr_wo_PigNat", pattern = chr_maf, full.names = T),
                     ordered = T,format = "SS")
   
   ## If loop that will try to find the phylomodel, but if not, it will generate it
   
   if (file.exists(paste0("data/mafs_chr_wo_PigNat/", chr_name, "_model_wo_PigNat.mod"))) {
      neutralMod <- read.tm(paste0("data/mafs_chr_wo_PigNat/", chr_name, "_model_wo_PigNat.mod"))
      
   } else {
      
      ### Loading only those features that are in our chromosome
      
      message("Getting 4d degenarated sites & model. This will take a while... ¯\\_(ツ)_/¯")
      
      chr_feats <- feats[feats$seqname == chr_name,]
      
      chr_feats$seqname <- "Surface_Astyanax"
      
      wholeChrom <- feat(seqname =  "Surface_Astyanax", src=".", feature = "all", 
                         start = align$offset, end = align$offset+ncol.msa(align, "Surface_Astyanax"))
      
      intergenicFeats <- inverse.feat(chr_feats, region.bounds=wholeChrom)
      
      intergenicFeats$feature <- "intergenic"
      
      chr_feats <- rbind.feat(chr_feats, intergenicFeats)
      
      table(chr_feats$feature)
      
      align4d <- get4d.msa(align, chr_feats)
      
      # generating the model
      
      neutralMod <- phyloFit(align4d, tree = tree_wo_sp_txt, subst.mod="REV")
      
      write.tm(neutralMod, paste0("data/mafs_chr_wo_PigNat/", chr_name, "_model_wo_PigNat.mod"))
      
   }
   # we want to filter only those features that are from the same chr than the alignment
   chr_feats <- feats[feats$seqname == chr_name,]
   
   chr_feats$seqname <- "Surface_Astyanax"
   
   wholeChrom <- feat(seqname =  "Surface_Astyanax", src=".", feature = "all", 
                      start = align$offset, end = align$offset+ncol.msa(align, "Surface_Astyanax"))
   
   intergenicFeats <- inverse.feat(chr_feats, region.bounds=wholeChrom)
   
   intergenicFeats$feature <- "intergenic"
   
   chr_feats <- rbind.feat(chr_feats, intergenicFeats)
   
   #using phastcons to estimate the most conserved elements of the alignment. The expected length and coverage have to be tuned.
   # we have to look for a method to do it properly
   # We can start with these values taken from https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1267528663_lo7zo4OFXa48Yody0aUhmi4hRiVT&db=fr3&c=chr21&g=cons8way
   # --expected-length=45 --target-coverage=0.3 and --rho=0.3 
   
   pc <- phastCons(align, neutralMod, viterbi = T, expected.length=45, target.coverage=0.3, rho = 0.3)
   
   # pcEM <- phastCons(align, model, viterbi = T, estimate.transitions = T)
   
   # most conserved elements
   consElements <- pc$most.conserved
   
   # this shows how many bases are predicted to be conserved
   # coverage.feat(consElements)
   
   # this shows the percentage of bases that are conserved. Since its a close tree, kind of makes sense
   
   chrom_coverage <- c(chrom_coverage, (coverage.feat(consElements)/coverage.feat(wholeChrom))*100) # 38% of the chr is conserved. Makes kind of sense
   
   CDS_coverage <- c(CDS_coverage,
                     coverage.feat(chr_feats[chr_feats$feature == "CDS",], consElements) * 100 /coverage.feat(chr_feats[chr_feats$feature == "CDS",])) # 60%% of covering of CDS is the target
   
   cons_elements_bed <- consElements[,c("seqname","start","end","attribute","score","strand")]
   
   cons_elements_bed$seqname <- chr_name
   
   cons_elements_bed$attribute <- paste(cons_elements_bed$seqname, cons_elements_bed$start, cons_elements_bed$end, sep = "_")
   
   write.table(cons_elements_bed, paste0("results/", format(Sys.Date(), format="%Y_%m_%d"),
                                         "_most_cons_regions_wo_PigNat.bed"), append = T,
               sep = "\t", quote = F, col.names = F, row.names = F)
   
}

stats <- data.frame(chr = gsub("_wo_PigNat.ss", "", list.files(pattern = "*_wo_PigNat.ss$",path = "data/mafs_chr_wo_PigNat", full.names = F)),
                    CDS_coverage = CDS_coverage,
                    Chr_coverage = chrom_coverage)

write.table(stats, "results/PhyloCons_Cons_regions_WO_PigNat_stats.txt", sep = "\t", row.names = F, quote = F)


#### Accelerated regions of IpCoco ######

#see http://compgen.cshl.edu/rphast/vignette2.pdf

coverage_of_cons_by_inf_regions <- NULL

consElements <- read.delim("results/2022_07_26_most_cons_regions_wo_PigNat_filt_100bp_merged_nonCoding.bed", header = F)

for (chr_maf in list.files(pattern = "*maf$",path = "data/mafs_chr_v2", full.names = F)) { # for loop
   
   
   chr_name <- gsub("\\.maf", "", chr_maf)
   print(chr_name)
   
   gc() # clean the memory
   
   model <- rphast::read.tm(list.files(pattern = paste0(chr_name,"_neutral_model.mod"),
                                       path = "data/mafs_chr_v2", full.names = T))
   
   
   align_complete <- strip.gaps.msa(read.msa(list.files(path = "data/mafs_chr_v2", 
                                                        pattern = paste0(chr_maf, "_complete_aln.ss"), full.names = T),
                                             ordered = T,format = "SS"))
   #unlink("data/AstMex_aln_chr19.maf")
   
   # positions that dont have missing information in cavefish
   
   #### positions that have at least 4 species
   
   hasAtLeastFour <- informative.regions.msa(align_complete, 5,
                                             spec = c("Surface_Astyanax", "Pygocentrus_nattereri", "Colossoma_macropomum",
                                                      "Electrophorus_electricus", "Ictalurus_punctatus", "Pangasianodon_hypophthalmus"))
   
   
   ### Overlap of these two set of regions and the conserved ones
   
   consElements_filt <- consElements[consElements$V1 == chr_name,]
   
   consElements_filt <- feat(seqname = "Surface_Astyanax",
                             start = consElements_filt$V2,
                             end = consElements_filt$V3)
   
   informativeElements <- coverage.feat(consElements_filt, hasAtLeastFour, get.feats = T)
   
   coverage_of_cons_by_inf_regions <- c(coverage_of_cons_by_inf_regions, 
                                        coverage.feat(informativeElements)*100/coverage.feat(consElements_filt))
   
   # 90% of our most conserved sites have passed this filter, we can tighten it up, but we can do it later
   
   # Regularize the size of cons elements to facilitate computation
   
   splitLength <- 50
   
   splitElements <- split.feat(informativeElements, f=splitLength, drop = T)
   
   # now we run the phyloP program for computing the likelihood ratio for every feature and aproximates p-values
   
   message("Computing Observed PhyloP scores")
   
   obsPhyloP <- phyloP(model, msa = align_complete, mode = "ACC", features = splitElements, subtree = "Pygocentrus_nattereri")
   
   # we need now to compute empirical pvalues. In order to do so we will generate random alignmets for features and compare
   # the random dist vs the observed (simulating a lot of alignments)
   
   elementAlign <- extract.feature.msa(align_complete, informativeElements)
   
   nrep <- 100000 # this is the recomended in the vignete
   
   simMsa <- sample.msa(elementAlign, nrep*splitLength, replace=TRUE)
   
   # produce features allowing long alignment to be interpreted as 
   # concatenation of shorter alignments
   
   startIdx <- seq(from=1, by=splitLength, length.out=nrep)
   
   features <- feat(seqname=names.msa(simMsa)[1], src="sim", feat=".",
                    start=startIdx,
                    end=startIdx+splitLength-1)
   message("Computing Simulated PhyloP")
   
   nonParaPhyloP <- phyloP(model, msa=simMsa, mode="ACC",
                           features=features, subtree="Pygocentrus_nattereri")
   ### Okey now lets compare the two models
   png(paste0("results/plots_stats/", chr_name,"_Pygocentrus_nattereri_NonCoding_Stat_plots.png"),width = 500, height = 500)
   par(mfrow=c(1,2), mar=c(5,2,4,2))
   
   qqplot(nonParaPhyloP$lnlratio,obsPhyloP$lnlratio,
          xlim=c(0,15),ylim=c(0,15), xlab="Simulated likelihood ratio",
          ylab="Observed likelihood ratio")
   abline(0, 1, lty=2)
   
   # yeah well, We can see that they are different, at least :/
   plot(density(obsPhyloP$lnlratio,adjust=3), lty=1,xlim=c(0,1),
        xlab="Likelihood Ratio",
        ylab="Density",main="", col="red")
   lines(density(nonParaPhyloP$lnlratio,adjust=3), lty=1,
         col="black",xlim=c(0,1))
   
   dev.off()
   # computing of the empirical pval
   
   empirical.pval <- function(x, dist) {
      sum(x <= dist)/length(dist)
   }
   
   message("Computing empirical p.values and writting results...")
   
   nonParaPval <- sapply(obsPhyloP$lnlratio, empirical.pval,
                         nonParaPhyloP$lnlratio)
   nonParaFDR <- p.adjust(nonParaPval, method="bonferroni")
   
   # Extracting those that are significant from our set of accelerated elements
   
   nonParaSigFeats <- splitElements[nonParaFDR < 0.05,]
   nrow(nonParaSigFeats)
   
   class(nonParaSigFeats)
   
   nonParaSigFeats$feature <- "PigNatAccRegions"
   
   nonParaSigFeats$score <- obsPhyloP$lnlratio[nonParaFDR < 0.05]
   
   
   
   # write.feat(nonParaSigFeats, "CaveAccRegions.gff")
   
   ## Sanity check that we have actually accelerated regions 
   
   #### TURNED OFF: Too much time running, we have seen that the branch length is increased greatly in other runs
   
   # consEleModel <- phyloFit(elementAlign, init.mod=model, no.opt=c("backgd", "ratematrix"))
   # 
   # CaveAccRegionsAln <- extract.feature.msa(copy.msa(align_complete), nonParaSigFeats)
   # 
   # AccModel <- phyloFit(CaveAccRegionsAln, init.mod=model, no.opt=c("backgd", "ratematrix"))
   # 
   # maxX <- depth.tree(AccModel$tree, "Cave_Astyanax")
   # 
   # pdf(paste0("results/", chr_name,"_Stat_plots.pdf"))
   # par(mfrow=c(1,2))
   # 
   # plot.phylo(read.tree(text=consEleModel$tree), x.lim=c(0, maxX),
   #             main=paste0("Conserved Elements of ",  chr_name))
   # 
   # plot.phylo(read.tree(text=AccModel$tree), x.lim=c(0, maxX),
   #            main=paste0("Accelerated Elements of ",  chr_name))
   # dev.off()
   
   nonParaSigFeats$seqname <- chr_name # this goes here in order to not disrupt the tree computation above
   
   
   
   Acc_elements_bed <- data.frame(chr = nonParaSigFeats$seqname,
                                  start = nonParaSigFeats$start,
                                  end = nonParaSigFeats$end,
                                  name = paste0("Acc_region_",chr_name,"_", rep(1:length(nonParaSigFeats$seqname))),
                                  score = nonParaSigFeats$score,
                                  strand = "+")
   
   # Acc_elements_bed$attribute <- paste(Acc_elements_bed$seqname, Acc_elements_bed$start, Acc_elements_bed$end, sep = "_")
   
   write.table(Acc_elements_bed, "results/Acc_regions_Pygocentrus_nattereri_NonCoding.bed", append = T,
               sep = "\t", quote = F, col.names = F, row.names = F)
   
}
