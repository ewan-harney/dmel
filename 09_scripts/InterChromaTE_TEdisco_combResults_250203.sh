cd /path/to/SuppMat/interchromate_TEresults

#############################################################################
# STEP 1 - TEMP2 output replace annoying names (per sample)
for FILE in TEMP2/*bed
do
    TRIM=${FILE%_R*}
    sample=${TRIM##*/}
    awk '{ gsub("UnFUnClUnAlig_RXX-TRIM_B22_2","NewFam16",$0); print $0 }' TEMP2/${sample}_R.insertion.bed |\
    awk '{ gsub("UnFUnClUnAlig_RLX-incomp_TOM_3","Gypsy-1",$0); print $0 }' | awk '{ gsub("UnFUnClUnAlig_RLX-incomp_TOM_2","Gypsy-35",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAlig_RLX-incomp_TOM","Micropia",$0); print $0 }' | awk '{ gsub("UnFUnClUnAlig_RLX-incomp_SLO_4","BEL-4",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAlig_RLX-incomp_SLO","gypsy6A",$0); print $0 }' | awk '{ gsub("UnFUnClUnAlig_RLX-incomp_RAL177","Copia-1",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAlig_RLX-incomp_MUN_3","Gypsy-24",$0); print $0 }' | awk '{ gsub("UnFUnClUnAlig_RLX-incomp_MUN","Gypsy-26",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAlig_RLX-incomp_FIN","Pifo",$0); print $0 }' | awk '{ gsub("UnFUnClUnAlig_RLX-incomp_COR_3","Gypsy-4",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAlig_RLX-incomp_COR_2","Gypsy-1",$0); print $0 }' | awk '{ gsub("UnFUnClUnAlig_RLX-incomp_COR","Gypsy-1",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAlig_RIX-incomp_JUT_2","Doc3",$0); print $0 }' | awk '{ gsub("UnFUnClUnAlig_RIX-comp_TEN","I-element",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAlig_noCat_SLO","Rt1c",$0); print $0 }' | awk '{ gsub("UnFUnClUnAlig_DXX-MITE_COR","NewFam14",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAlig_DHX-incomp_COR","INE-1",$0); print $0 }' | awk '{ gsub("UnFUnClUnAligUnFmclCluster051_DTX-incomp_SLO","NewFam06",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAligUnFmclCluster039_RLX-incomp_TOM","Gypsy-6",$0); print $0 }' | awk '{ gsub("UnFUnClUnAligUnFmclCluster039_RLX-incomp_COR","Gypsy-1",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAligcon9_UnFUnCl001_DTX-incomp","Mariner",$0); print $0 }' | awk '{ gsub("UnFUnClUnAligcon8_UnFmcl025_RLX-comp","Copia1",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAligcon7_UnFmcl002_RXX-LARD","LARD",$0); print $0 }' | awk '{ gsub("UnFUnClUnAligcon6_UnFmcl028_RLX-incomp","Chimpo",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAligcon4_UnFUnCl002_DHX-incomp","INE-1",$0); print $0 }' | awk '{ gsub("UnFUnClUnAligcon4_UnFmcl034_RLX-incomp","Gypsy-2",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAligcon4_UnFmcl032_RLX-incomp","Copia2",$0); print $0 }' | awk '{ gsub("UnFUnClUnAligcon41_UnFmcl001_RLX-incomp","Gypsy-2",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAligcon2_UnFUnCl003_RLX-incomp","Stalker4",$0); print $0 }' | awk '{ gsub("UnFUnClUnAligcon24_UnFmcl003_noCat","Invader3",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAligcon13_UnFmcl013_RLX-incomp","Micropia",$0); print $0 }' | awk '{ gsub("UnFUnClUnAligcon12_UnFmcl011_RLX-incomp","Gypsy-24_Dya",$0); print $0 }' |\
    awk '{ gsub("UnFUnClUnAligcon10_UnFmcl019_RLX-incomp","Copia-3_Dbp",$0); print $0 }' >  interchromate_TEresults/${sample}_temp2_insertm5_prep.bed
# STEP 2 - TEMP2 output tidy up and cut down to necessary info (per sample)
    awk '{ gsub("[^,_]+_","",$4); print $1 "\t" $2 "\t" $3 "\t" $4 }' interchromate_TEresults/${sample}_temp2_insertm5_prep.bed | awk '{ gsub(":.+","",$4); print $1 "\t" $2 "\t" $3 "\t" $4 }' | tail -n +2 | sort -k1,1 -k2n,2 | awk 'BEGIN {OFS="\t"};{ gsub("-","",$4); print tolower($0)}' | awk '{ gsub("dbuz","",$0); print $0 }' | sed 's/^chr//g' | awk '{print $1 "\t" $2 "\t" $3 "\t" NR "\t" $4 }' | awk '{ $4=sprintf("%0.6i", $4); print $1 "\t" $2 "\t" $3 "\t" "F2_F_1_temp_" $4 "\t" $5 "\t" "denovo"}' > interchromate_TEresults/${sample}_temp2_insertm5_nochr_ids.bed
done

awk '{print tolower($0)}' ISO1_euchromatin.bed > iso1_euchromatin.bed

# STEP 3 - Combine tidied TEMP2 from different samples
cd interchromate_TEresults

POP_LIST=$(echo -e "F3_P1\nF3_P2")
for Pop in ${POP_LIST}
do
### Prep TEMP2 file
    cat ${Pop}1_temp2_insertm5_nochr_ids.bed  ${Pop}2_temp2_insertm5_nochr_ids.bed ${Pop}3_temp2_insertm5_nochr_ids.bed | sort -k1,1 -k2,2n | bedtools merge -i - -c 5, -o distinct | bedtools intersect -a - -b ${Pop}1_temp2_insertm5_nochr_ids.bed -wao | bedtools intersect -a - -b ${Pop}2_temp2_insertm5_nochr_ids.bed -wao | bedtools intersect -a - -b ${Pop}3_temp2_insertm5_nochr_ids.bed -wao | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $9 "\t" $16 "\t" $23}' > ${Pop}_any_three_samps_temp2_insertm5_nochr.bed
### Prep PopTE2 file
    awk '-F[\t,]' '{if ($9 > 0 && $10 > 0 && $11 > 0) print $0}' ../PopTE2/${Pop}_all_three_samps.flybaseref.strand.teinsertions | awk '{ gsub("_[^_]+","",$5); print $2 "\t" $3 "\t" $3+1 "\t" $5 "\t" $9 "\t" $10 "\t" $11}' | sort -k1,1 -k2n,2 | awk 'BEGIN {OFS="\t"};{ gsub("-","",$4); print tolower($0)}' | awk '{ gsub("dbuz","",$0); print $0 }' | sed 's/^chr//g' | grep -v 'reference' > ${Pop}_any_three_samps_pop2_insert_denovo_nochr.bed
### Combine TEMP2 and PopTE2 results
    cat interchromate_TEresults/${Pop}_any_three_samps_pop2_insert_denovo_nochr.bed interchromate_TEresults/${Pop}_any_three_samps_temp2_insertm5_nochr.bed | sort -k1,1 -k2,2n | bedtools merge -i - -c 4 -o distinct -d 25 -delim "|" | bedtools intersect -a - -b interchromate_TEresults/${Pop}_any_three_samps_pop2_insert_denovo_nochr.bed -wao | bedtools intersect -a - -b interchromate_TEresults/${Pop}_any_three_samps_temp2_insertm5_nochr.bed -wao | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $9 "\t" $10 "\t" $11 "\t" $17 "\t" $18 "\t" $19 }' | uniq  > interchromate_TEresults/Full_output_${Pop}_pop2_temp2_denovo_md25.bed
### filter - must be in 2 samples of each program, and family must match
    awk '{if (($4 !~ /\|/ && $4 !~ /,/ ) &&  (($5 > 0 && $6 > 0 ) || ( $5 > 0 &&  $7 > 0 ) || ( $6 > 0 && $7 > 0)) && (($8 != "." && $9 != "." ) || ( $8 != "." &&  $10 != "." ) || ( $9 != "." && $10!= "."))) print $0}' Full_output_${Pop}_pop2_temp2_denovo_md25.bed | bedtools intersect -b iso1_euchromatin.bed -a - > Both_prog_2samps_${Pop}_pop2_temp2_denovo_md25.bed
### Give the insertion a name
    awk '{print $1 "\t" $2 "\t" $3 "\t" NR "\t" $4}' Both_prog_2samps_${Pop}_pop2_temp2_denovo_md25.bed | awk -v Pop="${Pop}" '{ $4=sprintf("%0.3i", $4); print $1 "\t" $2 "\t" $3 "\t" Pop "_nonref_" $4 "\t" $5 "\t" "denovo"}' > ${Pop}_denovo_strict.bed
done


#############################################################################
# STEP 6 - Prepare and tidy tlex output (per sample)

#awk '{print tolower($0)}' TEcopies_2417TEs_ISO1_nochr_sorted.txt > TEcopies_2417TEs_iso1_nochr_sorted.txt

sort -k4,4 ISO1_famTE_FlyBase_Repet_correct_vpopte2_sortbyTE.bed > ISO1_famTE_FlyBase_Repet_correct_vpopte2_sortbyTE_sorted.bed

TEList=TEcopies_2417TEs_ISO1_nochr_sorted.txt
annotation=ISO1_famTE_FlyBase_Repet_correct_vpopte2_sortbyTE_sorted.bed

cd ../Tlex

for FILE in */Tresults
do
    TRIM=${FILE%/*}
    sample=${TRIM##*-}
    sort -k2,2 tlex_F3-${sample}/Tresults | awk '{if ($5 == "present" || $5 == "polymorphic") print $0 }' | join -1 1 -2 2 ../interchromate_TEresults/${TEList} - | awk '{print tolower ($2 "\t" $3 "\t" $4 "\t" $1)}' | join -1 4 -2 4 ../interchromate_TEresults/${annotation} - | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $5 "\t" "reference"}' > ../interchromate_TEresults/tlex_F3_${sample}_prespoly_fam.bed
done

#############################################################################
# STEP 7 - Combine tidied Tlex output

cd ../interchromate_TEresults

POP_LIST=$(echo -e "F3_P1\nF3_P2")
for Pop in ${POP_LIST}
do
### Tlex results all possible insertions (found in any one sample)
    sort -k4,4 ${annotation} | join -a 1 -e'.' -o '1.1,1.2,1.3,0,1.5,2.5' -1 4 -2 4 - tlex_${Pop}1_prespoly_fam.bed | join -a 1 -e'.' -o '1.1,1.2,1.3,0,1.5,1.6,2.5' -1 4 -2 4 -  tlex_${Pop}2_prespoly_fam.bed | join -a 1 -e'.' -o '1.1,1.2,1.3,0,1.5,1.6,1.7,2.5' -1 4 -2 4 -  tlex_${Pop}3_prespoly_fam.bed | awk '{if ($6 != "." || $7 != "." || $8 != ".") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "reference" "\t" $6 "\t" $7 "\t" $8}' | sort -k1,1 -k2,2n | bedtools intersect -b iso1_euchromatin.bed -a - > ${Pop}_any_three_samps_tlex_combined_nochr.bed
### Tlex results identified in all 3 samples
    sort -k4,4 ${annotation} | join -a 1 -e'.' -o '1.1,1.2,1.3,0,1.5,2.5' -1 4 -2 4 - tlex_${Pop}1_prespoly_fam.bed | join -a 1 -e'.' -o '1.1,1.2,1.3,0,1.5,1.6,2.5' -1 4 -2 4 -  tlex_${Pop}2_prespoly_fam.bed | join -a 1 -e'.' -o '1.1,1.2,1.3,0,1.5,1.6,1.7,2.5' -1 4 -2 4 -  tlex_${Pop}3_prespoly_fam.bed | awk '{if ($6 != "." && $7 != "." && $8 != ".") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "reference" "\t" $6 "\t" $7 "\t" $8}' | sort -k1,1 -k2,2n | bedtools intersect -b iso1_euchromatin.bed -a - > ${Pop}_all_samps_tlex_combined_nochr.bed
done

###############
##
## Summary files

## print tolower first column of flybase genes
awk '{print tolower ($1) "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 }' flybase_genes.bed > flybase_genes_lc.bed

# Need to use bedtools window to relate to the genes

promoterLength=1000
POP_LIST=$(echo -e "F3_P1\nF3_P2")
for Pop in ${POP_LIST}
do
### prep file for reference
    awk '{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ${Pop}_all_samps_tlex_combined_nochr.bed | sort -k1,1 -k2,2n > final/${Pop}_TEs_chr_ref.bed
### simple file to be included in downstream R analysis
    bedtools window -a flybase_genes_lc.bed -b ${Pop}_TEs_chr_ref.bed -w $promoterLength | cut -f 5 | sort -u  > final/genes_flybase_TE_${Pop}_strict_promoter_${promoterLength}_ref.lst
### full list for sup mat TEs with gene hits
    bedtools window -a ${Pop}_TEs_chr_ref.bed -b flybase_genes_lc.bed -w $promoterLength | awk '{print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $4 "\t" $5 "\t" $10 "\t" $11 "\t" $12 "\t" $8 "\t" $9}' > final/TE_ref_${Pop}_strict_promoter_${promoterLength}_gene_present.txt
### full list for sup mat TEs with no gene hits
    bedtools window -v -a ${Pop}_TEs_chr_ref.bed -b flybase_genes_lc.bed -w $promoterLength | awk '{print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $4 "\t" $5 "\t.\t.\t.\t.\t."}' > final/TE_ref_${Pop}_strict_promoter_${promoterLength}_gene_absent.txt
done

######
# Now for the denovo lists
promoterLength=1000
POP_LIST=$(echo -e "F3_P1\nF3_P2")
for Pop in ${POP_LIST}
do
### prep file for denovo
    awk '{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}'  ${Pop}_denovo_strict.bed | sort -k1,1 -k2,2n > final/${Pop}_TEs_chr_den.bed
### simple file to be included in downstream R analysis
    bedtools window -a flybase_genes_lc.bed -b ${Pop}_TEs_chr_den.bed -w $promoterLength | cut -f 5 | sort -u  > final/genes_flybase_TE_${Pop}_strict_promoter_${promoterLength}_den.lst
### full list for sup mat TEs with gene hits
    bedtools window -a ${Pop}_TEs_chr_den.bed -b flybase_genes_lc.bed -w $promoterLength | awk '{print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $4 "\t" $5 "\t" $10 "\t" $11 "\t" $12 "\t" $8 "\t" $9}' > final/TE_den_${Pop}_strict_promoter_${promoterLength}_gene_present.txt
### full list for sup mat TEs with no gene hits
    bedtools window -v -a ${Pop}_TEs_chr_den.bed -b flybase_genes_lc.bed -w $promoterLength | awk '{print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $4 "\t" $5 "\t.\t.\t.\t.\t."}' > final/TE_den_${Pop}_strict_promoter_${promoterLength}_gene_absent.txt
done

#####
cat final/TE_*P1*txt | sort -k1,1 -k2,2n > final/TE_P1_all_insertions_and_genes.txt
cat final/TE_*P2*txt | sort -k1,1 -k2,2n > final/TE_P2_all_insertions_and_genes.txt
