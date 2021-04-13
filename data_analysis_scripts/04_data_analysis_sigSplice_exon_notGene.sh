#!/bin/bash
#SBATCH -p general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbrotman@email.unc.edu
#SBATCH --mem 50000 
#SBATCH -t 2:00:00 
#SBATCH -J data_analysis
#SBATCH -o slurm.data_analysis_sigSplice_exon_notGene.out 
#SBATCH -e slurm.data_analysis_sigSplice_exon_notGene.err

date

##################################
# Significant in splice and exons but not genes
##################################

# significant in splices and exons only
#---------------------------------
echo "number of splice and exon eGenes only FDR 1%"
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' splice_sigALL_FDR01.txt splice_exon_sigBoth_FDR01.txt | sort -k1,1 -u > splice_exon_sigBothOnly_FDR01_unique.txt

# Match in the lead splice junction
#---------------------------------
awk 'BEGIN{FS="\t"}{if ($15 == 1) print $0}' METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA_ENSG.tab | sort -k13,13 -g | sort -k16,16 -u > splice-level_topSNPs_perGene_FDR01_allColumns.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1,$9]=$0; next} ($1,$2) in a {print $1, $2, a[$1,$2]}' splice-level_topSNPs_perGene_FDR01_allColumns.txt splice_exon_sigBothOnly_FDR01_unique.txt | awk 'BEGIN{FS="\t";OFS="\t"} {print $1, $18, $2, $4, $15}' > tmp.txt

# Match in the exon and p-value
#---------------------------------
awk 'BEGIN{FS="\t";OFS="\t"} {print $1, $0}' ../METSIM_exon-level_cis-eQTL_nominal_adjPEERf50BMI_250kb_FDR0.01_geneType_geneName.tab | sed 's/_[0-9]*_[0-9]*//1' > METSIM_exon-level_250kb_FDR01_leadVariants250kb_wExon.tab
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1,$9]=$0; next} ($1,$3) in a {print $0, a[$1,$3]}' METSIM_exon-level_250kb_FDR01_leadVariants250kb_wExon.tab tmp.txt | cut -f 1-5,7,18 > tmp2.txt

# Match in the gene p-value (NOTE: there are instances where the sQTL variant is not within 250kb of gene TSS therefore I used the 1Mb gene results here)
#---------------------------------
awk 'BEGIN{FS="\t"}{if ($14 == 1) print $1}' ../METSIM_gene-level_cis-eQTL_nominal_adjPEERf40BMI_geneType_geneName_250kb.tab > genes_wLeadVariant_within250kb.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}{if ($1 in a) print $0}' genes_wLeadVariant_within250kb.txt ../METSIM_gene-level_cis-eQTL_nominal_adjPEERf40BMI_geneType_geneName_250kb.tab > METSIM_gene-level_250kb_leadVariants250kb.tab

awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1,$3]=$0; next} ($1,$8) in a {print a[$1,$8], $0}' tmp2.txt ../METSIM_gene-level_cis-eQTL_nominal_adjPEERf40BMI_geneType_geneName.tab | cut -f 1-7,19 > tmp3.txt
#awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1,$3]=$0; next} ($1,$8) in a {print a[$1,$8], $0}' tmp2.txt METSIM_gene-level_250kb_leadVariants250kb.tab | cut -f 1-7,19 > tmp3.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}{if ($1 in a) print $0, a[$1]}' ../match_by_variant/gene-level_topSNPs_perGene.txt tmp3.txt | cut -f 1-8,10,11 > sig_splice_exon_noGene_250kb.txt

# rm tmp*.txt
