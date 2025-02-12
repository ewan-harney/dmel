Within this 10_interchromate_reference folder are 10 files that serve as reference files for the analyses: 

- ISO1_euchromatin.bed - Positions of euchromatin identified by Rech et al. 2022 (doi: 10.1038/s41467-022-29518-8)
- ISO1_famTE_FlyBase_Repet_correct_vpopte2_sortbyTE.bed - List of 2417 reference insertions in D. melanogaster ISO1 reference genome, together with chromosome, postion and TE family 
- TEcopies_2417TEs_ISO1_nochr.txt - List of 2417 reference insertions in D. melanogaster ISO1 reference genome, together with chromsome, position and strand information 
- te-hierarchy.txt - D. melanogaster TE hierarchy: IDs followed by family (e.g. gypsy, Doc) and order (e.g. LTR, LINE)
- TElist_2417TEs_ISO1_nochr.txt - List of 2417 reference insertions in D. melanogaster ISO1 reference genome
- consensuses_curated_v4.fasta - D. melanogaster consensus TE nucleotide sequences based on those available from Dfam, Berkeley Drosophila Genome Project (BDGP) and those identified in the population study of Rech et al, 2022  (doi: 10.1038/s41467-022-29518-8)
- FBgn_strand.txt - Flybase genes, they symbols and orientation
- fbtn_fbgn.tsv - Associations between flybase genes (FBgn IDs) and flybase transcripts (FBtr IDs)
- flybase_genes.bed - Chromosome and position information for flybase genes, as well as transcript and symbol information

Notes:
- Files containing the 2417 reference insertions contain two different ID numbers: FBti IDs for pre-existing insertions, ISO1 IDs for those identified in population study of Rech et al. 2022 (doi: 10.1038/s41467-022-29518-8)
- Throughout the analyses we used Dmel_FB_r6_v46.fa as a reference genome, which is available from flybase (https://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.46_FB2022_03/). We focused only on Chromsomes X, 2L, 2R, 3L, and 3R (Muller Elemetns, A, B, C, D and E). Some analyses required adding non-numeric characters to the start of the reference - these containing the prefix 'chr' to the chromosome name in the reference, i.e. "chr2L", "chr2R" etc (if the name of the reference doesn't feature the 'chr' string, or features the string 'nochr', it will be as in the original reference, i.e. "2L", "2R" etc). Some analyses also required the mitochondrial genome to be included as an addtional sequence in the reference fasta. This is indicated as Dmel_FB_r6_v46_chr_Mt