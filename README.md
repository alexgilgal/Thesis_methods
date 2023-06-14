# Thesis methods

The main goal of this repository is to show the persons that read my PhD dissertation how the computational analyses were done. The scripts have commented code inside them, so everyone can follow what the scripts are doing. Also, the Materials and Methods section explains what every script does.



***



## Installation
The use of the scripts requires several programs in order to run properly and reproduce the results. To save the user the tedious installation processes and the searching of the programs manually, we have generated conda recipes. These conda recipes can be used directly by users in order to generate ready-to-analyse conda environments. The only thing the user has to change from these recipes are the prefix variables.

If you don't have Conda installed, we recommend using Anaconda (or Miniconda) to start. Check [this link](https://docs.conda.io/en/latest/miniconda.html) to install conda from Miniconda. 

In each folder of the different analyses, there are YML files that contain these recipes.

Creating the environment for the ATAC-seq would be like this 
```
conda env create -f ATAC_env.yml
```
Once you have created the environment, do not forget to activate the environment with the following command:

```
conda activate ATAC_env
```

With these steps, you have loaded all the necessary programs to reproduce all our ATAC-seq analyses.

## Contents

#### **ATAC-seq analysis**
In this folder, you can find all the scripts necessary to reproduce ATAC-seq analyses. If you intend to use these scripts, please see above how to install the corresponding conda recipes.

These scripts should be used in the following order (this is a suggestion):

- **ATAC_pipe.pl:** This is our pipeline to map raw ATAC-seq and preprocess this data.
- **idr_mod_job.sh**: This script identifies reproducible peaks between samples.
- **ATAC_deseq2.R**: Using this script, the user can determine which peaks are differential between different conditions.
- **homer_enrichment.sh**: This script can be used with any BED file, but we used it mainly with the differential peaks from the previous script. This provides the user with the enrichment analysis of TFBS in the regions of interest.
- **combine_homer.R**: this script joints several Homer results in one dot plot for easy visualization of those results.

Independently, the user can use the script CRISP_calling.sh in order to call variants in the ATAC-seq signal. Then the user can run the scripts inside the folder Motifbreak_tobias in order to understand how these variants affect TFBS and Tf binding to the DNA in certain conditions.

#### **RNA-seq analysis**
This folder contains the scripts used for analysing RNA-seq and its results. Some of these analyses can be used with gene lists, independently if they have been generated with RNA-seq or not.

In general, we have followed the next steps:

- **map_STAR_and_htseq.sh**: This script maps the RNA-seq raw data against the indexed reference genome and transcriptome. Once the reads are mapped, the count tables are generated using htseq.
- **deseq_rnaseq.R**: Using this script, we detected the genes that changed their expression between conditions. The result of this is a gene list of upregulated and downregulated genes.
- **Amphi_zebra_fams.R**: In the case of amphioxus, differentially expressed genes were also translated into zebrafish genes via orthogroups. The input of this function is the  gene IDs of amphioxus and a table containing the orthology information between amphioxus and zebrafish genes. This orthology table contains orthogroups, or gene families, to which each amphioxus and zebrafish gene belongs. For each amphioxus gene, the function identifies to which orthogroup this amphioxus gene belongs and retrieves the entire group of zebrafish genes that belong to the same orthogroup. Since zebrafish have undergone three rounds of WGDs, one amphioxus gene can have several zebrafish orthologs. 
- **GO_enrichment.R**: Using this script, we determined the function of the gene lists used as input. For example, we determined the function of the differentially expressed genes. This script can be used with any gene list, not only differential genes.
- **GO_dotplots.R**: The result of the previous script can be visualized using this script. this script takes several TXT GO tables produced by the previous script and aggregates them in a unique plot.

#### **Integration_ATAC_RNA**
The contents of this file allow us to join the results coming from ATAC and RNA-seq.
We have used two main strategies:
- Associate gene to peaks. To do that, we have used the script **Gene_2_peaks.sh**. This script needs as input the folder where the BED files that will be analyzed (ATAC-seq peaks) are, and also a BED file of the TSS of the genes. This script will generate several TXT as output. the most important ones are the gene list of genes associated with peaks and the relationship of these peaks to the genes. Using this information, we could further analyze these integrated data.
- Using ananse. In this folder, there is a pdf guide of how we used ananse for cavefish and how to use it with other organisms. The developer team of Ananse has done great work documenting their tool [here](https://anansepy.readthedocs.io/en/master/), so check the page for more information on how it works and how to use it.

#### **Comparative genomics**

This folder contains the commands and scripts used in order to generate our comparative genomics analyses. If the user wants to use Cactus, we recommend reading its documentation carefully [here](https://github.com/ComparativeGenomicsToolkit/cactus). We have also deposited the genome sources used in our analyses and their phylogenetic relationship in the tree_single.nw.
The script rphast_analysis.R was used in order to compute conserved and accelerated regions between different species.


## Authors and acknowledgement
The main author of this repository is Alejandro Gil (agilgal AT upo DOT es agilgal@upo.es)
Some scripts and code come from diverse persons from my lab, like Juan J. Tena (ATAC_pipe.pl) or Alberto Perez Posadas.


## License
MIT License

Copyright (c) 2023 Alejandro Gil

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


