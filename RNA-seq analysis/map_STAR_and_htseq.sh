#!/bin/bash
#SBATCH -J STAR
##SBATCH -p day
#SBATCH -N 1
#SBATCH -n 6

hostname

## Change manually
indexDirectory=/home/ska/cpaliou/catshark/shark_genome/STAR_index_catshark
gtf_path=~/catshark/shark_genome/Annotation_new/20210914_sScyCan1.1_annot_FixChrNames_filtered_v2_1gene1transcript_FINAL.gtf


## Don't change (except cores)
SCRATCH=/home/ska/cpaliou/tmp/cpaliou_$SLURM_JOB_ID
mkdir -p $SCRATCH || exit $?
mkdir $SCRATCH/starIndex
mkdir $SCRATCH/results
wd=$PWD

folder=$1
namePrefix=${folder:0:(-1)}

echo $namePrefix

echo "Merging fq files"

# The & allow us to throw procesess in parallel, so we speed things up

zcat $folder/*_1.fq.gz | gzip > $folder/${namePrefix}_merged_1.fq.gz &

zcat $folder/*_2.fq.gz | gzip > $folder/${namePrefix}_merged_2.fq.gz &

# after using & is important to use the wait comand so everything is done before going on

wait %1 %2

echo "Copying files into scratch..."

cp $indexDirectory/* $SCRATCH/starIndex/

file1=$folder/*merged_1.fq.gz
file2=$folder/*merged_2.fq.gz

cp $file1 $SCRATCH/1.fq.gz &
cp $file2 $SCRATCH/2.fq.gz &
cp $gtf_path $SCRATCH/model.gtf &

wait %1 %2 %3

cd $SCRATCH

echo "Done"

echo "Starting alignment..."
#Add if this option if no strand specific and wanna use cufflinks
#--outSAMstrandField intronMotif
#gunzip *.fq.gz
#STAR 
# /home/ska/agilgal/programs/STAR-2.5.3a/bin/Linux_x86_64_static/STAR --readFilesCommand zcat --genomeDir ./starIndex --outSAMtype BAM SortedByCoordinate --sjdbGTFfile model.gtf --runThreadN 25 --readFilesIn 1.fq.gz --outFileNamePrefix results/$namePrefix 
/home/ska/agilgal/programs/STAR-2.5.3a/bin/Linux_x86_64_static/STAR --readFilesCommand zcat --genomeDir ./starIndex --outSAMtype BAM SortedByCoordinate --sjdbGTFfile model.gtf --runThreadN 6 --readFilesIn 1.fq.gz 2.fq.gz --outFileNamePrefix results/$namePrefix 

echo "Alignment Done!"
echo "Starting counting"

## Call to htseq with these options:

# -f bam because input is bam

# --nounique all to take all the reads that map to exons shared between trancripts

# -s no we dont take into count the strandess of the assay

# -r name BAM is sorted by name. If the BAM is sorted by position use -r pos

htseq-count -f bam --nonunique all -s no -r pos results/${namePrefix}*bam model.gtf > results/${namePrefix}.counts
                                                    

echo "Counting Done, cleaning up"

dest=${wd}/${folder}
cp -r results/* $dest

rm -rf $SCRATCH
