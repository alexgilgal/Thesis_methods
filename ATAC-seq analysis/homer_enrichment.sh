#!/bin/bash

## Usage bash motifs.sh <peak/BED file> <output directory> <Background set of peaks in BED format>


echo $1
echo $2
findMotifsGenome.pl $1 ~/genomes/AstMex2.0/AstMex_v2_fixed_names.fa $2 -size given -mset vertebrates -p 12 -bg $3 -nomotif
