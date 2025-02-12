# 07_GOenrich
Within this 07_GOenrich folder are 10 results files for GO enrichment analysis produced by the InterChromaTE_downstream_analysis_final_250203.R script

- 1.  GO_AS_CUEU_out.txt (GO, Akaa Shared, Chromatin Up Expression Up)
- 2.  GO_AS_CDEU_out.txt (GO, Akaa Shared, Chromatin Down Expression Up)
- 3.  GO_AU_CDEU_out.txt (GO, Akaa Unique, Chromatin Down Expression Up)
- 4.  GO_AU_CDED_out.txt (GO, Akaa Unique, Chromatin Down Expression Down)
- 5.  GO_MU_CDEU_out.txt (GO, Manz Unique, Chromatin Down Expression Up)

- 6.  M_3up6up_out.txt (Manz, F3 expression up, F6 expression up)
- 7.  M_3dn6up_out.txt (Manz, F3 expression down, F6 expression up)

- 8.  Ak_TED_ED_out.txt (Akaa, TEs, Denovo (non reference), Expression Down)
- 9.  Ma_TED_CU_out.txt (Manz, TEs, Denovo (non reference), Chromatin Up)
- 10. Ma_TER_CU_out.txt (Akaa, TEs, Reference, Chromatin Up)

Tables 1-5 relate to sets of genes which showed both significant expression and significant accessibility differences in the F3. The 4 characters before "_out.txt" provide information about the direction of change between control and heat shock, so CUEU (Chromatin Up Expression Up) is functional enrichment for genes that had significantly more accessible chromatin and significantly higher transcript expression in the heat shocked treatment. AS (Akaa Shared) implies that these genes were the same in both populations (so can be read as 'Both populations'), while AU (Akaa Unique) and MU (Manz Unique). include genes that showed significant expression and accessibility only in one population.

Tables 6-7 relate to sets of genes that showed significant increases in expression following heat shock for Manz in both the F3 and F6 generations. Significant enrichment was found for genes that showed increased expression in both generations (3up6up), and also those that showed reduced expression in the F3 but increased expression in the F6 (3dn6up).

Tables 8-10 relate to sets of genes that were proximate to either non-reference (TED) or reference TEs (TER) in either Akaa (Ak) or Manz (Ma), and which showed particulart patterns of expression (e.g. ED = expression down) or chromatin accessibility (CU = Chromatin Up) following heat shock.

All these tables were generated with the script InterChromaTE_downstream_analysis_final_250203.R