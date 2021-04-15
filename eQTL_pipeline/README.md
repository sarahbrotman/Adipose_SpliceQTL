# eQTL Pipeline

## Table of Contents
1. Quantify genes and exons using QTLtools for each sample
2. Quantify splice junctions using LeafCutter and calculate PSI
3. Filter transcripts for 5 counts in 25% of samples
4. Inverse normal transform transcripts
5. Adjust for BMI and inverse normal transform BMI-adjusted residuals
6. Run PEER on inverse normalized BMI-adjusted residuals
7. Create input files for QTLtools
8. Run QTLtools for all sets of PEER factors
9. Choose optimal number of PEER factors from QTLtools output files
10. Filter results from 1Mb to 250kb distance from variant to start site

## 1.	Quantify genes and exons using QTLtools for each sample
```
QTLtools quan --bam sample482.bam --gtf gencode.v19.annotation.gtf.gz --samples 482 --out-prefix
482 --filter-mapping-quality 150 --filter-mismatch 5 --filter-mismatch-total 5 --rpkm
--no-merge

```
## 2. Quantify splice junctions using LeafCutter and calculate PSI

https://davidaknowles.github.io/leafcutter/index.html

## 3.	Filter transcripts (genes, exons, splice junctions) for 5 counts in 25% of samples, adjust for TMMs, and report in CPMs

```R
# Load edgeR
#---------------------------------------
library(edgeR)

# Load expression file (counts) and create DGEList 
#---------------------------------------
counts_exons <- read.delim("exon_quantification_matrix_counts.txt", sep = "\t", row.names = 1,
header = TRUE)
norm_exons <- DGEList(counts=counts_exons)

# Filter out exons that are not meeting the threshold (>5 counts in 25% individuals)
#---------------------------------------
sample_size <- 426
keep_exons_counts <- rowSums(norm_exons$counts>=5) >= round(0.25*sample_size)
exons_counts_filtered <- norm_exons[keep_exons_counts,]

# Calculate TMM values
#---------------------------------------
exons_counts_filtered_TMMadj <- calcNormFactors(exons_counts_filtered)
norm.factors<-exons_counts_filtered_TMMadj$samples$norm.factors
eff.lib.size=exons_counts_filtered_TMMadj$samples$lib.size * norm.factors

# Generate TMM normalized CPMs
#---------------------------------------
METSIM_5counts_adjTMM_25percent_exons=cpm(exons_counts_filtered_TMMadj)

# Write to table
#---------------------------------------
write.table(METSIM_5counts_adjTMM_25percent_exons, "METSIM_5counts_adjTMM_25percent_exons.txt",
quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
```

## 4.	Inverse normal transform transcripts

```R
input_data<-read.table("METSIM_5counts_adjTMM_25percent_exons.txt", sep="\t", header=FALSE)
attach(input_data)

#to get residuals
output <-matrix(NA, nrow=426, ncol=87217)  # nrow is the sample size, ncol is the number of
traits
for (i in 1:87217) {
	y<-input_data[,i]
	
	#inverse normal transform
	inv<-qnorm((rank(y,na.last="keep")-0.5)/sum(!is.na(y)))

  	output[,i]<-inv 
}

phenotypes.out <-data.frame(output)
write.table(phenotypes.out,file="METSIM_5counts_adjTMM_25percent_exons_invNorm.txt", col.names=F, row.names=F, quote=F,
sep="\t")
```

## 5.	Adjust for BMI and inverse normal transform BMI-adjusted residuals

```R
library(preprocessCore)

phenotypes<-read.table("METSIM_5counts_adjTMM_25percent_exons_invNorm_BMI.txt", sep="\t", header=TRUE, row.names=1)
attach(phenotypes)

#to get residuals
output<-matrix(NA, nrow=426, ncol=18563)  # nrow is the sample size, ncol is the number of traits
for (i in 11:18574) {
	y<-phenotypes[,i]


	#regression
	model=lm(y ~ BMI, na.action=na.exclude)
	resid=resid(model)
	
	#inverse normal transform
	inv<-qnorm((rank(resid,na.last="keep")-0.5)/sum(!is.na(resid)))

  	output[,i-11]<-inv 
  	#output[,i-10]<-resid
}

phenotypes.out<-data.frame(output)
write.table(phenotypes.out,file="METSIM_5counts_adjTMM_25percent_exons_invNorm_BMI_FINAL.txt", row.names=F,
col.names=F, quote=F, sep="\t")
```

## 6.	Run PEER on inverse normalized BMI-adjusted residuals
Run 0-100 PEER factors in sets of 10 (i.e. 0, 10, 20, 30…)

```R
# Necessary imports 
library(peer)

# Read in data, initialize model 
expr = read.delim(METSIM_5counts_adjTMM_25percent_exons_invNorm_BMI_FINAL.txt', header=FALSE, sep="\t")
model = PEER()
PEER_setPhenoMean(model,as.matrix(expr))
PEER_setNmax_iterations(model, 150)

# Set number of observed data and confounders 
PEER_setPhenoMean(model,as.matrix(expr))
PEER_setNk(model,80)

# Actually run the models 
PEER_update(model)

# Get output
factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)

write.table(factors, "METSIM.PEERfactors.k80.txt", sep="\t")
write.table(weights, "METSIM.PEERweights.k80.txt", sep="\t")
write.table(precision, "METSIM.PEERprecision.k80.txt", sep="\t")
write.table(residuals, "METSIM.PEERresiduals.k80.txt", sep="\t")
```

## 7.	Create input files for QTLtools:
  + An indexed phenotype data matrix in BED format. Use the inverse normal transformed gene expression (unadjusted for BMI or PEER factors)
  + An indexed genotype data matrix in VCF format
  + A covariate file with BMI and PEER factors

## 8.	Run QTLtools for all sets of PEER factors

```
QTLtools cis --vcf genotypes.vcf.gz --bed exons.bed.gz --cov known_covariates.txt.gz --nominal 1 --out nominals.txt –-window	1000000 --normal
```

## 9.	Choose optimal number of PEER factors from QTLtools output files
  + Calculate 1% FDR (use pvalues to calculate qvalues, then filter by 0.01) https://github.com/nfusi/qvalue
  + Count number of unique exons with an association at FDR1% threshold
  + Increase in PEER factors until we see less than a 1% increase in exon with an additional set of 10 PEER factors

## 10. Filter results from 1Mb to 250kb distance from variant to start site
```
sed 's/ /\t/g' METSIM_exon-level_cis-eQTL_nominal_adjPEERf80BMI.txt | awk 'BEGIN{FS="\t";OFS="\t"}{if ($7>-250000 && $7 < 250000) print $0}' > METSIM_exon-level_cis-eQTL_nominal_adjPEERf80BMI_250kb.txt
```

