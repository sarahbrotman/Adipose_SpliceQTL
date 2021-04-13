#!/bin/bash
#SBATCH -p general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbrotman@email.unc.edu
#SBATCH --mem 30000 
#SBATCH -t 2:00:00 
#SBATCH -J data_analysis
#SBATCH -o slurm.data_analysis_scatterPlot.out 
#SBATCH -e slurm.data_analysis_scatterPlot.err

date

# ------------------------------------------------------------------
# Comparing gene-level and exon-level (top SNPs)
# ------------------------------------------------------------------

# Extracting best SNP per gene (all)
#---------------------------------------
#awk 'BEGIN{FS="\t"}{if ($14 == 1) print $1}' ../METSIM_gene-level_cis-eQTL_nominal_adjPEERf40BMI_geneType_geneName_250kb.tab > genes_wLeadVariant_within250kb.txt
#awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}{if ($1 in a) print $0}' genes_wLeadVariant_within250kb.txt ../METSIM_gene-level_cis-eQTL_nominal_adjPEERf40BMI_geneType_geneName_250kb.tab > METSIM_gene-level_250kb_leadVariants250kb.tab

#awk 'BEGIN{FS="\t"}{if ($14 == 1) print $1}' ../METSIM_exon-level_cis-eQTL_nominal_adjPEERf50BMI_250kb_geneType_geneName_exon.tab > exons_wLeadVariant_within250kb.txt
#awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}{if ($1 in a) print $0}' exons_wLeadVariant_within250kb.txt ../METSIM_exon-level_cis-eQTL_nominal_adjPEERf50BMI_250kb_geneType_geneName_exon.tab > METSIM_exon-level_250kb_leadVariants250kb.tab

#awk 'BEGIN{FS="\t"}{if ($14 == 1) print $0}' METSIM_gene-level_250kb_leadVariants250kb.tab | cut -f 1,8,12,13 > gene-level_topSNPs_perGene.txt
#sed 's/_[0-9]*_[0-9]*//1' METSIM_exon-level_250kb_leadVariants250kb.tab > METSIM_exon-level_250kb_leadVariants250kb_noExon.tab
#awk 'BEGIN{FS="\t"}{if ($14 == 1) print $0}' METSIM_exon-level_250kb_leadVariants250kb_noExon.tab | cut -f 1,8,12,13 | sort -k3,3 -g | sort -k1,1 -u > exon-level_topSNPs_perGene.txt



# Matching p-values and betas FDR 1% (gene-level top SNPs)
#---------------------------------------
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$8) in a {print $1, $8, $12, $13, a[$1,$8]}' OFS='\t' gene-level_topSNPs_perGene.txt METSIM_exon-level_250kb_leadVariants250kb_noExon.tab | cut -f 1-4,7,8 > exon_gene_topSNP_pval_beta_match.txt
sed -i '1i gene\tSNP\texon_pval\texon_beta\tgene_pval\tgene_beta' exon_gene_topSNP_pval_beta_match.txt
awk '{if ($3<=0.000401 || $5<=0.00115) print $0 }' exon_gene_topSNP_pval_beta_match.txt > exon_gene_topSNP_pval_beta_match.eitherOrFDR01.txt
sed -i '1i gene\tSNP\texon_pval\texon_beta\tgene_pval\tgene_beta' exon_gene_topSNP_pval_beta_match.eitherOrFDR01.txt

echo "finished comparing gene-level topSNP and exon-level"


# Matching p-values and betas FDR 1% (exon-level top SNPs)
#---------------------------------------
#awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$8) in a {print $1, $8, $12, $13, a[$1,$8]}' OFS='\t' exon-level_topSNPs_perGene.txt METSIM_gene-level_250kb_leadVariants250kb.tab | cut -f 1-4,7,8 > gene_exon_topSNP_pval_beta_match.txt
#sed -i '1i gene\tSNP\tgene_pval\tgene_beta\texon_pval\texon_beta' gene_exon_topSNP_pval_beta_match.txt
#awk '{if ($3<=0.00115 || $5<=0.000401) print $0 }' gene_exon_topSNP_pval_beta_match.txt > gene_exon_topSNP_pval_beta_match.eitherOrFDR01.txt
#sed -i '1i gene\tSNP\tgene_pval\tgene_beta\texon_pval\texon_beta' gene_exon_topSNP_pval_beta_match.eitherOrFDR01.txt

echo "finished comparing exon-level topSNP and gene-level"

# Extracting best exon per gene
#---------------------------------------
sort -k3,3 -g exon_gene_topSNP_pval_beta_match.eitherOrFDR01.txt | sort -r -k1,1 -u > bestExon_gene_topSNP_pval_beta_match.eitherOrFDR01.txt
sort -k5,5 -g gene_exon_topSNP_pval_beta_match.eitherOrFDR01.txt | sort -r -k1,1 -u > gene_bestExon_topSNP_pval_beta_match.eitherOrFDR01.txt

date
