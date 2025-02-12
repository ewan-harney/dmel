# Data and scripts supporting InterChromaTE project
This repository contains the data and scripts that accompany the MS: _Transgenerational effects of heat shock on gene regulation and fitness-related traits are stronger in arid than temperate Drosophila populations_ which is available to view on bioRxiv at doi: . This work formed part of the InterChromaTE project (Interactions between chromatin and transposable elements in rapid adaptation to environmental stress) funded by the European Union’s Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 101030460.

Genomic data for this project are available from NCBI (https://www.ncbi.nlm.nih.gov/) Sequence Read Archive (SRA) under accession number PRJNA1002872.


## Project Outline
This project aimed to assess how natural variation in the epigenome and genome might interact to influence gene regulation within and across generations, using experimental data collected in the vinegar fly _Drosophila melanogaster_. 
We carried out laboratory experiments to measure gene expression and chromatin accessibility responses to heat shock in female _D. melanogaster_ from arid (Spanish) and temperate (Finnish) climates and their consequences for offspring phenotypes, and associated the expression and chromatin responses with population variation in TEs. We then measured the same molecular and phenotypic traits three generations later to explore transgenerational inheritance. 
The finding of these experiments are detailed in the manuscript on bioRxiv:

## Description of Data
The dataset consists of 11 directories. Within each directory is a seperate ##_README.txt with more information about the files and subdirectories within. Below is a brief description of the 10 directories:

* **01_kallisto** : 24 subdirectories, one for each RNA-seq sample. Within each are the kallisto results, including the abundance.tsv file, which is used for DESeq2 analysis.
* **02_nfcoreatac** : 1 subdirectory containing peaks files for postions of open chromatin. The .bam files used by interpareto were too large to be included in this data deposit.
* **03_intepareto** : 8 Rdata objects, containing matched and integrated RNA and ATAC results for all 8 pairwise comparisons considered in the study.
* **04_TEresults** : 3 subdirectories, one for Tlex3 results, one for Temp2 results, and one for PoPoolationTE2 results
* **05_TEconsensus** : 6 text files. Four .lst files, including reference (ref) and non-reference (den) insertions for P1 (population 1 = Akaa) and P2 (population 2 : Manzanares). Two .txt files, showing TE-gene associations for each population
* **06_downstream** : 3 Rdata objects and 7 summary .txt tables. Output from the _downstream_analysis_final_250203.R script 
* **07_GOenrich** : 10 summary. txt tables. Output from GO enrichment analyses in the _downstream_analysis_final_250203.R script
* **08_phenotypic** : 7 text files containing ctmax measures from the F2, counts of eggs, pupae and adults from the offspring of F3 and F6 flies, and age-to-pupation and age-to-eclosion data for these same F3 and F6 cohorts.
* **09_scripts** : 13 bash scripts, 3 R scripts and an R function used in the analysis of the aforementioned data
* **10_reference** : 10 text files (.txt/.bed/.tsv/.fasta) that provide information about the iso1 TE references or the identity, position or orientation of genes and transcripts within the reference genome.
* **11_meta** : 5 text files (.txt/.csv) that were used as meta data in the analysis of this dataset.