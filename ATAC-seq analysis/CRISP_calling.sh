#!/bin/bash

#SBATCH -J CRISP
#SBATCH -p day
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --constrain="SET2"

hostname

echo "Preparing files"

SCRATCH=/ssdscratch/agilgal_$SLURM_JOB_ID

genome=/home/ska/agilgal/genomes/AstMex2.0/AstMex_v2_fixed_names.fa
genome_index=/home/ska/agilgal/genomes/AstMex2.0/AstMex_v2_fixed_names.fa.fai

mkdir -p $SCRATCH || exit $?

cp *bam $SCRATCH/
#cp *bam.bai $SCRATCH/
cp $genome $SCRATCH/genome.fa
cp $genome_index $SCRATCH/genome.fa.fai
cp pools.txt $SCRATCH/
cp ../all_peaks_nonorm.bed $SCRATCH/

wd=$PWD

cd $SCRATCH
mkdir results

echo "Files ready. Starting to mark duplicates"

CRISP --bams pools.txt --ref genome.fa --bed all_peaks_nonorm.bed --VCF results/variant_calls.vcf 

dest=/home/ska/agilgal/cavefish/ATAC/new_genome/analysis_cave_vs_surface/results

echo $dest

mv results/* $dest

rm -rf $SCRATCH || exit $?