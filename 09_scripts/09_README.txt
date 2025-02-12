Within this 09_interchromate_scripts folder are 17 script files. The scripts are referred to in other README.txt files where they are required, They are listed here according to which part of the analysis they refer to.

RNA-Seq analysis with Kallisto
- InterChromaTE_rna_kallisto-index.sh
- InterChromaTE_rna_kallisto-quant.sh

ATAC-Seq analysis with nf-core/atac
- InterChromaTE_nf-core-atacseq.sh

Combining RNA-Seq and ATAC-Seq results
- InterChromaTE_intepareto_prepTranscriptTable_250204.sh
- InterChromaTE_intepareto_run_250204.R
- doMatch_4_function.R

TE Prep
- InterChromaTE_TEdisco_prep-fastp.sh
TE Tlex
- InterChromaTE_TEdisco_Tlex3-prep.sh
TE Temp2
- InterChromaTE_TEdisco_TEMP2.sh
TE PoPoolationTE2
- InterChromaTE_TEdisco_PopTE2-mpileup.sh
- InterChromaTE_TEdisco_PopTE2-prepReads.sh
- InterChromaTE_TEdisco_PopTE2-prepRepeatMask.sh
- InterChromaTE_TEdisco_PopTE2-stats.shCombining results
TE Combining Tlex3, Temp2 and PopTE2 results
- InterChromaTE_TEdisco_combResults_250203.sh

Combining results from Intepareto and TE discovery pipeline, including GO enrichment
- InterChromaTE_downstream_analysis_final_250203.R

Phenotypic analysis
- InterChromaTE_phenotypic_analysis_final_250203.R



