#!/bin/bash

# Usage Autoclosest_gene_ens94.sh FOLDER (folder contains the sub/sobexp.bed files)

folder=$1

## Change where you have the TSS file

tss=/home/ska/agilgal/genomes/AstMex2.0/Enseml/Ensembl_Ast_Mex2_hub_TSS_bedsort.bed


cd $folder

for i in *.bed 

do
	echo $i 

	sortBed -i $i > tmp_sorted.bed

	## Upstream genes

	closestBed -D ref -id -nonamecheck -k 2 -a tmp_sorted.bed -b $tss > ${i:0:(-4)}_2genes_upst.txt

	## Downstream genes

	closestBed -D ref -iu -nonamecheck -k 2 -a tmp_sorted.bed -b $tss > ${i:0:(-4)}_2genes_downst.txt

	## Merge of both lists

	cat ${i:0:(-4)}_2genes_* |cut -f 7 |uniq > ${i:0:(-4)}_2genes_1col.txt

	## generate updown files

	cat ${i:0:(-4)}_2genes_upst.txt ${i:0:(-4)}_2genes_downst.txt > ${i:0:(-4)}_2genes_updown.txt

done

rm tmp_sorted.bed

cd ..