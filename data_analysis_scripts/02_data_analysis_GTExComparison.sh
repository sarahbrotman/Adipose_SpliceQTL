#!/bin/bash
#SBATCH -p general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbrotman@email.unc.edu
#SBATCH --mem 30000 
#SBATCH -t 2:00:00 
#SBATCH -J data_analysis
#SBATCH -o slurm.data_analysis_GTExComparison.out 
#SBATCH -e slurm.data_analysis_GTExComparison.err

#---------------------------------
# GTEx analyses
#---------------------------------

date

# Count number of sQTLs
#---------------------------------
echo "# of sQTLs METSIM"
wc -l METSIM_splice-level_250kb_FDR05_leadVariants250kb.tab
echo "# of sQTLs GTEx"
zcat ../GTExV8_comparison/GTEx_Analysis_v8_sQTL/Adipose_Subcutaneous.v8.sqtl_signifpairs.txt.gz | cut -f 1,2 | sed '1d' | sort | uniq | wc -l

# Count number of sQTLs mapped back to annotation
#---------------------------------
echo "# of sQTLs mapped back to annotation METSIM"
wc -l METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL_noNA.tab
echo "# of sQTLs mapped back to annotation GTEx"
zcat ../GTExV8_comparison/GTEx_Analysis_v8_sQTL/Adipose_Subcutaneous.v8.sqtl_signifpairs.txt.gz | cut -f 1,2 | sed '1d' | sort | uniq | wc -l

# Count number of splice junctions with an sQTL
#---------------------------------
echo "# of splice junctions with an sQTL METSIM"
#sort -k13,13 -g METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | sort -k1,1 -u | sort -k2,2 -u | wc -l
sort -k13,13 -g METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | sort -k16,16 -u | sort -k2,2 -u | wc -l
echo "# of splice junctions with an sQTL GTEx"
zcat ../GTExV8_comparison/GTEx_Analysis_v8_sQTL/Adipose_Subcutaneous.v8.sqtl_signifpairs.txt.gz | cut -f 2 | sed '1d' | sed 's/:/\t/4' | sort -k1,1 -u | wc -l

# Count number of genes with an sQTL
#---------------------------------
echo "# of genes with an sQTL METSIM"
#cut -f 1 METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | sort | uniq | wc -l
cut -f 16 METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | sort | uniq | wc -l
echo "# of genes with an sQTL GTEx"
zcat ../GTExV8_comparison/GTEx_Analysis_v8_sQTL/Adipose_Subcutaneous.v8.sqtl_signifpairs.txt.gz | cut -f 2 | sed '1d'| sed 's/:/\t/4' | sort -k2,2 -u | wc -l

# Intersect METSIM and GTEx splice junctions
#---------------------------------
module load bedtools

#sort -k13,13 -g METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | sort -k1,1 -u | sort -k2,2 -u | awk 'BEGIN{FS="\t";OFS="\t"}{print $3, $4, $5, $2}' | sed 's/^/chr/1' > METSIM_sQTL_FDR05_oneSplice_perGene.bed
sort -k13,13 -g METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | sort -k16,16 -u | sort -k2,2 -u | awk 'BEGIN{FS="\t";OFS="\t"}{print $3, $4, $5, $2}' | sed 's/^/chr/1' > METSIM_sQTL_FDR05_oneSplice_perGene.bed
bedtools intersect -f 1E-9 -wo -a ../GTExV8_comparison/GTEx_Analysis_v8_sQTL/Adipose_Subcutaneous.v8.sqtl_signifpairs_afterLiftover.bed -b METSIM_sQTL_FDR05_oneSplice_perGene.bed > GTEx_METSIM_FDR05_intersect_oneSplice_perGene.txt

# Count number of splice junctions in other study
#--------------------------------- 
echo "# of GTEx splice junctions in METSIM"
sort -k4,4 -u GTEx_METSIM_FDR05_intersect_oneSplice_perGene.txt | wc -l
echo "# of METSIM splice junctions in GTEx"
sort -k8,8 -u GTEx_METSIM_FDR05_intersect_oneSplice_perGene.txt | wc -l

date
