#!/bin/bash
#SBATCH -p general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbrotman@email.unc.edu
#SBATCH --mem 30000 
#SBATCH -t 2:00:00 
#SBATCH -J data_analysis
#SBATCH -o slurm.data_analysis_eGeneOverlaps.out 
#SBATCH -e slurm.data_analysis_eGeneOverlaps.err

date

##################################
# e/sGene overlap analysis
##################################

#---------------------------------
# Extract lead SNPs per gene
#---------------------------------

awk 'BEGIN{FS="\t"}{if ($14 == 1) print $0}' METSIM_gene-level_250kb_FDR01_leadVariants250kb.tab | cut -f 1,8,12,13 > gene-level_topSNPs_perGene_FDR01.txt
sed 's/_[0-9]*_[0-9]*//1' METSIM_exon-level_250kb_FDR01_leadVariants250kb.tab > METSIM_exon-level_250kb_FDR01_leadVariants250kb_noExon.tab
awk 'BEGIN{FS="\t"}{if ($14 == 1) print $0}' METSIM_exon-level_250kb_FDR01_leadVariants250kb_noExon.tab | cut -f 1,8,12,13 | sort -k3,3 -g | sort -k1,1 -u > exon-level_topSNPs_perGene_FDR01.txt
awk 'BEGIN{FS="\t"}{if ($15 == 1) print $0}' METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | awk 'BEGIN{FS="\t";OFS="\t"}{print $16, $9, $13, $14}' | sort -k3,3 -g | sort -k1,1 -u > splice-level_topSNPs_perGene_FDR01.txt
echo "finished extracting top SNPs FDR 1%"

awk 'BEGIN{FS="\t"}{if ($14 == 1) print $0}' METSIM_gene-level_250kb_FDR05_leadVariants250kb.tab | cut -f 1,8,12,13 > gene-level_topSNPs_perGene_FDR05.txt
sed 's/_[0-9]*_[0-9]*//1' METSIM_exon-level_250kb_FDR05_leadVariants250kb.tab > METSIM_exon-level_250kb_FDR05_leadVariants250kb_noExon.tab
awk 'BEGIN{FS="\t"}{if ($14 == 1) print $0}' METSIM_exon-level_250kb_FDR05_leadVariants250kb_noExon.tab | cut -f 1,8,12,13 | sort -k3,3 -g | sort -k1,1 -u > exon-level_topSNPs_perGene_FDR05.txt
awk 'BEGIN{FS="\t"}{if ($15 == 1) print $0}' METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | awk 'BEGIN{FS="\t";OFS="\t"}{print $16, $9, $13, $14}' | sort -k3,3 -g | sort -k1,1 -u > splice-level_topSNPs_perGene_FDR05.txt
echo "finished extracting top SNPs FDR 5%"


#---------------------------------
# Match in ENSG IDs
#---------------------------------
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$10]=$7;next}{if ($1 in a) print a[$1], $2, $3, $4, $5}' ../gencode.v19.annotation.reextracted.20180724_addStrand_gene.txt splice-level_topSNPs_perGene_FDR01.txt > splice-level_topSNPs_perGene_FDR01_ENSG.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$10]=$7;next}{if ($1 in a) print a[$1], $2, $3, $4, $5}' ../gencode.v19.annotation.reextracted.20180724_addStrand_gene.txt splice-level_topSNPs_perGene_FDR05.txt > splice-level_topSNPs_perGene_FDR05_ENSG.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$10]=$7;next}{if ($16 in a) print a[$16], $0}' ../gencode.v19.annotation.reextracted.20180724_addStrand_gene.txt METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | cut -f 1,3- > METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA_ENSG.tab
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$10]=$7;next}{if ($16 in a) print a[$16], $0}' ../gencode.v19.annotation.reextracted.20180724_addStrand_gene.txt METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL_noNA.tab | cut -f 1,3- > METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL_noNA_ENSG.tab


#---------------------------------
# Splice lead variant overlaps 
#---------------------------------
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$8) in a {print $1, $8, a[$1,$8]}' OFS='\t' splice-level_topSNPs_perGene_FDR01_ENSG.txt METSIM_exon-level_250kb_FDR01_leadVariants250kb_noExon.tab | cut -f 1,2 > splice_exon_sigBoth_FDR01.txt
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$8) in a {print $1, $8, a[$1,$8]}' OFS='\t' splice-level_topSNPs_perGene_FDR01_ENSG.txt METSIM_gene-level_250kb_FDR01_leadVariants250kb.tab | cut -f 1,2  > splice_gene_sigBoth_FDR01.txt
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$8) in a {print $1, $8, a[$1,$8]}' OFS='\t' splice_exon_sigBoth_FDR01.txt METSIM_gene-level_250kb_FDR01_leadVariants250kb.tab | cut -f 1,2 > splice_sigALL_FDR01.txt
echo "finished matching splices, exons, and genes FDR 1%"

awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$8) in a {print $1, $8, a[$1,$8]}' OFS='\t' splice-level_topSNPs_perGene_FDR05_ENSG.txt METSIM_exon-level_250kb_FDR05_leadVariants250kb_noExon.tab | cut -f 1,2 > splice_exon_sigBoth_FDR05.txt
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$8) in a {print $1, $8, a[$1,$8]}' OFS='\t' splice-level_topSNPs_perGene_FDR05_ENSG.txt METSIM_gene-level_250kb_FDR05_leadVariants250kb.tab | cut -f 1,2  > splice_gene_sigBoth_FDR05.txt
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$8) in a {print $1, $8, a[$1,$8]}' OFS='\t' splice_exon_sigBoth_FDR05.txt METSIM_gene-level_250kb_FDR05_leadVariants250kb.tab | cut -f 1,2 > splice_sigALL_FDR05.txt
echo "finished matching splices, exons, and genes FDR 5%"

# significant in splices and exons only
#---------------------------------
echo "number of splice and exon eGenes only FDR 1%"
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' splice_sigALL_FDR01.txt splice_exon_sigBoth_FDR01.txt | sort -k1,1 -u | wc -l
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' splice_sigALL_FDR01.txt splice_exon_sigBoth_FDR01.txt > splice_exon_sigBothOnly_FDR01.txt

echo "number of splice and exon eGenes only FDR 5%"
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' splice_sigALL_FDR05.txt splice_exon_sigBoth_FDR05.txt | sort -k1,1 -u | wc -l
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' splice_sigALL_FDR05.txt splice_exon_sigBoth_FDR05.txt > splice_exon_sigBothOnly_FDR05.txt

# significant in splices and genes only
#---------------------------------
echo "number of splice and gene eGenes only FDR 1%"
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' splice_sigALL_FDR01.txt splice_gene_sigBoth_FDR01.txt | sort -k1,1 -u | wc -l
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' splice_sigALL_FDR01.txt splice_gene_sigBoth_FDR01.txt > splice_gene_sigBothOnly_FDR01.txt

echo "number of splice and gene eGenes only FDR 5%"
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' splice_sigALL_FDR05.txt splice_gene_sigBoth_FDR05.txt | sort -k1,1 -u | wc -l
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' splice_sigALL_FDR05.txt splice_gene_sigBoth_FDR05.txt > splice_gene_sigBothOnly_FDR05.txt

# significant in all three
#---------------------------------
echo "number of eGenes in all analyses FDR 1%"
sort -k1,1 -u splice_sigALL_FDR01.txt | wc -l

echo "number of eGenes in all analyses FDR 5%"
sort -k1,1 -u splice_sigALL_FDR05.txt | wc -l

# significant in only splices
#---------------------------------
echo "number of eGenes in only splices FDR 1%"
cat splice_exon_sigBothOnly_FDR01.txt splice_gene_sigBothOnly_FDR01.txt splice_sigALL_FDR01.txt | sort -k1,1 -u > tmp.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' tmp.txt splice-level_topSNPs_perGene_FDR01_ENSG.txt | cut -f 1,2 | sort -k1,1 -u | wc -l

echo "number of eGenes in only splices FDR 5%"
cat splice_exon_sigBothOnly_FDR05.txt splice_gene_sigBothOnly_FDR05.txt splice_sigALL_FDR05.txt | sort -k1,1 -u > tmp2.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' tmp2.txt splice-level_topSNPs_perGene_FDR05_ENSG.txt | cut -f 1,2 | sort -k1,1 -u | wc -l

rm tmp*.txt
echo "FINISHED SPLICE LEAD VARIANTS..."


#---------------------------------
# Exon lead variant overlaps 
#---------------------------------
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$8) in a {print $1, $8, a[$1,$8]}' OFS='\t' exon-level_topSNPs_perGene_FDR01.txt METSIM_gene-level_250kb_FDR01_leadVariants250kb.tab | cut -f 1,2 | sort -k1,1 -u > tmp.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$7]=$10;next}{if ($1 in a) print $0, a[$1]}' ../gencode.v19.annotation.reextracted.20180724_addStrand_gene.txt tmp.txt > exon_gene_sigBoth_FDR01.txt
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$9) in a {print $1, $9, a[$1,$9]}' OFS='\t' exon-level_topSNPs_perGene_FDR01.txt METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA_ENSG.tab | cut -f 1,2 > exon_splice_sigBoth_FDR01.txt
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$9) in a {print $1, $9, a[$1,$9]}' OFS='\t' exon_gene_sigBoth_FDR01.txt METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA_ENSG.tab | cut -f 1,2 > exon_sigALL_FDR01.txt
echo "finished matching splices, exons, and genes FDR 1%"

awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$8) in a {print $1, $8, a[$1,$8]}' OFS='\t' exon-level_topSNPs_perGene_FDR05.txt METSIM_gene-level_250kb_FDR05_leadVariants250kb.tab | cut -f 1,2 | sort -k1,1 -u > tmp2.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$7]=$10;next}{if ($1 in a) print $0, a[$1]}' ../gencode.v19.annotation.reextracted.20180724_addStrand_gene.txt tmp2.txt > exon_gene_sigBoth_FDR05.txt
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$9) in a {print $1, $9, a[$1,$9]}' OFS='\t' exon-level_topSNPs_perGene_FDR05.txt METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL_noNA_ENSG.tab | cut -f 1,2 > exon_splice_sigBoth_FDR05.txt
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$9) in a {print $1, $9, a[$1,$9]}' OFS='\t' exon_gene_sigBoth_FDR05.txt METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL_noNA_ENSG.tab | cut -f 1,2 > exon_sigALL_FDR05.txt
echo "finished matching exons, splices, and genes FDR 5%"

# significant in splices and exons only
#---------------------------------
echo "number of splice and exon eGenes only FDR 1%"
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' exon_sigALL_FDR01.txt exon_splice_sigBoth_FDR01.txt | sort -k1,1 -u | wc -l
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' exon_sigALL_FDR01.txt exon_splice_sigBoth_FDR01.txt > exon_splice_sigBothOnly_FDR01.txt

echo "number of splice and exon eGenes only FDR 5%"
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' exon_sigALL_FDR05.txt exon_splice_sigBoth_FDR05.txt | sort -k1,1 -u | wc -l
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' exon_sigALL_FDR05.txt exon_splice_sigBoth_FDR05.txt > exon_splice_sigBothOnly_FDR05.txt

# significant in exons and genes only
#---------------------------------
echo "number of splice and gene eGenes only FDR 1%"
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' exon_sigALL_FDR01.txt exon_gene_sigBoth_FDR01.txt | sort -k1,1 -u | wc -l
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' exon_sigALL_FDR01.txt exon_gene_sigBoth_FDR01.txt > exon_gene_sigBothOnly_FDR01.txt

echo "number of splice and gene eGenes only FDR 5%"
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' exon_sigALL_FDR05.txt exon_gene_sigBoth_FDR05.txt | sort -k1,1 -u | wc -l
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' exon_sigALL_FDR05.txt exon_gene_sigBoth_FDR05.txt > exon_gene_sigBothOnly_FDR05.txt

# significant in all three
#---------------------------------
echo "number of eGenes in all analyses FDR 1%"
sort -k1,1 -u exon_sigALL_FDR01.txt | wc -l

echo "number of eGenes in all analyses FDR 5%"
sort -k1,1 -u exon_sigALL_FDR05.txt | wc -l

# significant in only exons
#---------------------------------
echo "number of eGenes in only exons FDR 1%"
cat exon_splice_sigBothOnly_FDR01.txt exon_gene_sigBothOnly_FDR01.txt exon_sigALL_FDR01.txt | sort -k1,1 -u > tmp3.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' tmp3.txt exon-level_topSNPs_perGene_FDR01.txt | cut -f 1-2 | sort -k1,1 -u | wc -l

echo "number of eGenes in only exons FDR 5%"
cat exon_splice_sigBothOnly_FDR05.txt exon_gene_sigBothOnly_FDR05.txt exon_sigALL_FDR05.txt | sort -k1,1 -u > tmp4.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' tmp4.txt exon-level_topSNPs_perGene_FDR05.txt | cut -f 1-2 | sort -k1,1 -u | wc -l

rm tmp*.txt
echo "FINISHED EXON LEAD VARIANTS..."


#---------------------------------
# Gene lead variant overlaps 
#---------------------------------
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$8) in a {print $1, $8, a[$1,$8]}' OFS='\t' gene-level_topSNPs_perGene_FDR01.txt METSIM_exon-level_250kb_FDR01_leadVariants250kb_noExon.tab | cut -f 1,2 | sort -k1,1 -u > tmp.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$7]=$10;next}{if ($1 in a) print $0, a[$1]}' ../gencode.v19.annotation.reextracted.20180724_addStrand_gene.txt tmp.txt > gene_exon_sigBoth_FDR01.txt
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$9) in a {print $1, $9, a[$1,$9]}' OFS='\t' gene-level_topSNPs_perGene_FDR01.txt METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA_ENSG.tab | cut -f 1,2 > gene_splice_sigBoth_FDR01.txt
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$9) in a {print $1, $9, a[$1,$9]}' OFS='\t' gene_exon_sigBoth_FDR01.txt METSIM_splice-level_250kb_FDR01_leadVariants250kb_geneName_fixed_FINAL_noNA_ENSG.tab | cut -f 1,2 > gene_sigALL_FDR01.txt
echo "finished matching splices, exons, and genes FDR 1%"

awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$8) in a {print $1, $8, a[$1,$8]}' OFS='\t' gene-level_topSNPs_perGene_FDR05.txt METSIM_exon-level_250kb_FDR05_leadVariants250kb_noExon.tab | cut -f 1,2 | sort -k1,1 -u > tmp.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$7]=$10;next}{if ($1 in a) print $0, a[$1]}' ../gencode.v19.annotation.reextracted.20180724_addStrand_gene.txt tmp.txt > gene_exon_sigBoth_FDR05.txt
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$9) in a {print $1, $9, a[$1,$9]}' OFS='\t' gene-level_topSNPs_perGene_FDR05.txt METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL_noNA_ENSG.tab | cut -f 1,2 > gene_splice_sigBoth_FDR05.txt
awk 'NR==FNR {a[$1,$2]=$0; next} ($1,$9) in a {print $1, $9, a[$1,$9]}' OFS='\t' gene_exon_sigBoth_FDR05.txt METSIM_splice-level_250kb_FDR05_leadVariants250kb_geneName_fixed_FINAL_noNA_ENSG.tab | cut -f 1,2 > gene_sigALL_FDR05.txt
echo "finished matching exons, splices, and genes FDR 5%"

# significant in splices and genes only
#---------------------------------
echo "number of splice and gene eGenes only FDR 1%"
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' gene_sigALL_FDR01.txt gene_splice_sigBoth_FDR01.txt | sort -k1,1 -u | wc -l
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' gene_sigALL_FDR01.txt gene_splice_sigBoth_FDR01.txt > gene_splice_sigBothOnly_FDR01.txt

echo "number of splice and gene eGenes only FDR 5%"
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' gene_sigALL_FDR05.txt gene_splice_sigBoth_FDR05.txt | sort -k1,1 -u | wc -l
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' gene_sigALL_FDR05.txt gene_splice_sigBoth_FDR05.txt > gene_splice_sigBothOnly_FDR05.txt

# significant in exons and genes only
#---------------------------------
echo "number of exon and gene eGenes only FDR 1%"
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' gene_sigALL_FDR01.txt gene_exon_sigBoth_FDR01.txt | sort -k1,1 -u | wc -l
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' gene_sigALL_FDR01.txt gene_exon_sigBoth_FDR01.txt > gene_exon_sigBothOnly_FDR01.txt

echo "number of exon and gene eGenes only FDR 5%"
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' gene_sigALL_FDR05.txt gene_exon_sigBoth_FDR05.txt | sort -k1,1 -u | wc -l
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' gene_sigALL_FDR05.txt gene_exon_sigBoth_FDR05.txt > gene_exon_sigBothOnly_FDR05.txt

# significant in all three
#---------------------------------
echo "number of eGenes in all analyses FDR 1%"
sort -k1,1 -u gene_sigALL_FDR01.txt | wc -l

echo "number of eGenes in all analyses FDR 5%"
sort -k1,1 -u gene_sigALL_FDR05.txt | wc -l

# significant in only genes
#---------------------------------
echo "number of eGenes in only genes FDR 1%"
cat gene_splice_sigBothOnly_FDR01.txt gene_exon_sigBothOnly_FDR01.txt gene_sigALL_FDR01.txt | sort -k1,1 -u > tmp3.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' tmp3.txt gene-level_topSNPs_perGene_FDR01.txt | cut -f 1-2 | sort -k1,1 -u | wc -l

echo "number of eGenes in only genes FDR 5%"
cat gene_splice_sigBothOnly_FDR05.txt gene_exon_sigBothOnly_FDR05.txt gene_sigALL_FDR05.txt | sort -k1,1 -u > tmp3.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR==NR{a[$1]=$0;next}!a[$1]{print $0}' tmp3.txt gene-level_topSNPs_perGene_FDR05.txt | cut -f 1-2 | sort -k1,1 -u | wc -l

rm tmp*.txt
echo "FINISHED GENE LEAD VARIANTS..."

date
