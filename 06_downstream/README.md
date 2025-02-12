# 06_downstream
Within this 06_downstream folder are 10 results files combining some of the main results from the InterChromaTE_downstream_analysis_final_250203.R script

- deseq2_rnaseq_lists_250203.Rdata
- deseq2_atacseq_lists_250203.Rdata
- interchromate_integrated_results_250203.Rdata
- summary_table_Fbgn.txt

- Akaa_TE_F3_tab.txt
- Akaa_TE_F6_tab.txt
- Akaa_TEs_F3F6transgen.txt
- Manz_TE_F3_tab.txt
- Manz_TE_F6_tab.txt
- Manz_TEs_F3F6transgen.txt


The first 2 Rdata objects provide the lists of differentially expressed genes (deseq2_rnaseq_lists_250203.Rdata) and differentially accessbile regions (deseq2_atacseq_lists_250203.Rdata) identified by the various DESeq2 analyses, as well as the results of set analysis comparing different DE and DA gene sets, identified with the setdiff() and interset() functions.

The thrid Rdata object, interchromate_integrated_results_250203.Rdata contains several summary table combining RNA-seq and ATAC-seq results (LFC and P-values) together with gene and TE-presence information. Most of these objects are split by population, but summary_table_Fbgn combines both populations and both generations.

This final summary file is therefore also included as a standalone table, summary_table_Fbgn.txt

The final 6 files are the precursors to tables 4, 5 and 6 from the manuscript "Transgenerational effects of heat shock on Drosophila melanogaster gene regulation and life history depend on adaptive environment", identifying genes proximate to TEs with significant patterns of expression and accessibility. Further filtering of the results in these text files was carried out to generate tables 4,5 and 6.