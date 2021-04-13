#!/bin/bash
#SBATCH -p general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbrotman@email.unc.edu
#SBATCH --mem 30000 
#SBATCH -t 2:00:00 
#SBATCH -J data_analysis
#SBATCH -o slurm.data_analysis_overallCounts.out 
#SBATCH -e slurm.data_analysis_overallCounts.err

date

# Filter sQTLs by FDR
#--------------------------------------
# awk 'BEGIN{FS="\t";OFS="\t"}{if ($12 < 7.73E-05) print $0}' ../plusBMI_k30_nominals_distToStart250kb.txt > plusBMI_k30_nominals_distToStart250kb_FDR01.txt
# awk 'BEGIN{FS="\t";OFS="\t"}{if ($12 < 5.46E-04) print $0}' ../plusBMI_k30_nominals_distToStart250kb.txt > plusBMI_k30_nominals_distToStart250kb_FDR05.txt
echo "finished filtering sQTLs by FDR"

# Subset molecules for lead variants within 250kb
#--------------------------------------
# awk 'BEGIN{FS="\t"}{if ($14 == 1) print $1}' ../METSIM_gene-level_cis-eQTL_nominal_adjPEERf40BMI_geneType_geneName_250kb_FDR01.tab > genes_wLeadVariant_within250kb_FDR01.txt
# awk 'BEGIN{FS="\t"}{if ($14 == 1) print $1}' ../METSIM_exon-level_cis-eQTL_nominal_adjPEERf50BMI_250kb_FDR0.01_geneType_geneName.tab > exons_wLeadVariant_within250kb_FDR01.txt
# awk 'BEGIN{FS="\t"}{if ($14 == 1) print $1}' plusBMI_k30_nominals_distToStart250kb_FDR01.txt > splices_wLeadVariant_within250kb_FDR01.txt
echo "finished extracting molecules with lead variant within 250kb FDR 1%"

# awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}{if ($1 in a) print $0}' genes_wLeadVariant_within250kb_FDR01.txt ../METSIM_gene-level_cis-eQTL_nominal_adjPEERf40BMI_geneType_geneName_250kb_FDR01.tab > METSIM_gene-level_250kb_FDR01_leadVariants250kb.tab
# awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}{if ($1 in a) print $0}' exons_wLeadVariant_within250kb_FDR01.txt ../METSIM_exon-level_cis-eQTL_nominal_adjPEERf50BMI_250kb_FDR0.01_geneType_geneName.tab > METSIM_exon-level_250kb_FDR01_leadVariants250kb.tab
# awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}{if ($1 in a) print $0}' splices_wLeadVariant_within250kb_FDR01.txt plusBMI_k30_nominals_distToStart250kb_FDR01.txt > METSIM_splice-level_250kb_FDR01_leadVariants250kb.tab
echo "finished subsetting files for only molecules with a lead variant within 250kb FDR 1%"

# awk 'BEGIN{FS="\t"}{if ($14 == 1) print $1}' ../METSIM_gene-level_cis-eQTL_nominal_adjPEERf40BMI_geneType_geneName_250kb_FDR05.tab > genes_wLeadVariant_within250kb_FDR05.txt
# awk 'BEGIN{FS="\t"}{if ($14 == 1) print $1}' ../METSIM_exon-level_cis-eQTL_nominal_adjPEERf50BMI_250kb_FDR0.05_geneType_geneName.tab > exons_wLeadVariant_within250kb_FDR05.txt
# awk 'BEGIN{FS="\t"}{if ($14 == 1) print $1}' plusBMI_k30_nominals_distToStart250kb_FDR05.txt > splices_wLeadVariant_within250kb_FDR05.txt
echo "finished extracting molecules with lead variant within 250kb FDR 5%"

# awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}{if ($1 in a) print $0}' genes_wLeadVariant_within250kb_FDR05.txt ../METSIM_gene-level_cis-eQTL_nominal_adjPEERf40BMI_geneType_geneName_250kb_FDR05.tab > METSIM_gene-level_250kb_FDR05_leadVariants250kb.tab
# awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}{if ($1 in a) print $0}' exons_wLeadVariant_within250kb_FDR05.txt ../METSIM_exon-level_cis-eQTL_nominal_adjPEERf50BMI_250kb_FDR0.05_geneType_geneName.tab > METSIM_exon-level_250kb_FDR05_leadVariants250kb.tab
# awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}{if ($1 in a) print $0}' splices_wLeadVariant_within250kb_FDR05.txt plusBMI_k30_nominals_distToStart250kb_FDR05.txt > METSIM_splice-level_250kb_FDR05_leadVariants250kb.tab
echo "finished subsetting files for only molecules with a lead variant within 250kb FDR 5%"

# Add corresponding genes to splice junctions
#--------------------------------------
# awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$0}' METSIM_splice-level_250kb_FDR01_leadVariants250kb.tab | sed 's/:[0-9]*:[0-9]*:/:/' | awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$2; next}{if ($1 in a) print a[$1], $0}' ../clustersToGenes_Filt5in25.txt2 - | cut -f 1,3- > METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed.tab
# awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$0}' METSIM_splice-level_250kb_FDR05_leadVariants250kb.tab | sed 's/:[0-9]*:[0-9]*:/:/' | awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$2; next}{if ($1 in a) print a[$1], $0}' ../clustersToGenes_Filt5in25.txt2 - | cut -f 1,3- > METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed.tab
echo "finished adding genes to the splice junctions"

# Choose one of the duplicate genes (we did this manually but here is the end list)
#--------------------------------------
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$2;next}{if ($1 in a) {print $0, a[$1]} else{print $0, $1}}' FDR05_dupGene_selection.txt METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed.tab > METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL.tab
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$2;next}{if ($1 in a) {print $0, a[$1]} else{print $0, $1}}' FDR05_dupGene_selection.txt METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed.tab > METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL.tab

awk 'BEGIN{FS="\t";OFS="\t"}{if ($16 != "NA") print $0}' METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL.tab > METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA.tab
awk 'BEGIN{FS="\t";OFS="\t"}{if ($16 != "NA") print $0}' METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL.tab > METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL_noNA.tab

# Number of genes from all analyses
#--------------------------------------
echo "# of FDR 1% genes from gene-level cis-eQTLs 1Mb"
sort -u -k1,1 ../METSIM_gene-level_cis-eQTL_nominal_adjPEERf40BMI_FDR01_geneType_geneName.tab | wc -l

echo "# of FDR 1% genes from the gene-level cis-eQTLs 250kb"
sort -u -k1,1 METSIM_gene-level_250kb_FDR01_leadVariants250kb.tab | wc -l

echo "# of FDR 1% genes from exon-level cis-eQTLs"
sed 's/_[0-9]*_[0-9]*//g' METSIM_exon-level_250kb_FDR01_leadVariants250kb.tab | sort -u -k1,1 | wc -l 

echo "# of FDR 1% genes from splice-level cis-eQTLs"
sort -u -k1,1 METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | wc -l
sort -u -k16,16 METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | wc -l

# Number of exons and splices
#--------------------------------------
echo "# of FDR 1% exons from exon-level cis-eQTLs"
sort -u -k1,1 METSIM_exon-level_250kb_FDR01_leadVariants250kb.tab | wc -l

echo "# of FDR 1% splices from splice-level cis-eQTLs"
sort -u -k2,2 METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | wc -l

# Number of associations
#--------------------------------------
echo "# of FDR 1% SNP-gene associations from gene-level cis-eQTLs 1Mb"
wc -l ../METSIM_gene-level_cis-eQTL_nominal_adjPEERf40BMI_FDR01_geneType_geneName.tab

echo "# of FDR 1% SNP-gene associations from gene-level cis-eQTLs 250kb"
wc -l METSIM_gene-level_250kb_FDR01_leadVariants250kb.tab

echo "# of FDR 1% SNP-exon associations from exon-level cis-eQTLs"
wc -l METSIM_exon-level_250kb_FDR01_leadVariants250kb.tab 

echo "# of FDR 1% SNP-gene associations from exon-level cis-eQTLs"
sed 's/_[0-9]*_[0-9]*//g' METSIM_exon-level_250kb_FDR01_leadVariants250kb.tab | awk -F "\t" '!seen[$1, $8]++' | wc -l 

echo "# of FDR 1% SNP-splice associations from splice-level cis-eQTLs"
wc -l METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA.tab

echo "# of FDR 1% SNP-gene associations from splice-level cis-eQTLs"
awk -F "\t" '!seen[$1, $9]++' METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | wc -l 
awk -F "\t" '!seen[$16, $9]++' METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | wc -l

# Number of eSNPs
#--------------------------------------
echo "# of FDR 1% eSNPs from gene-level cis-eQTLs"
sort -k8,8 -u METSIM_gene-level_250kb_FDR01_leadVariants250kb.tab | wc -l

echo "# of FDR 1% eSNPs from exon-level cis-eQTLs"
sort -k8,8 -u METSIM_exon-level_250kb_FDR01_leadVariants250kb.tab | wc -l

echo "# of FDR 1% eSNPs from splice-level cis-eQTLs"
sort -k9,9 -u METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | wc -l

date
