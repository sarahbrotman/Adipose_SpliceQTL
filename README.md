# Adipose_SpliceQTL

## Description
Alternate splicing events can create isoforms that alter gene function, and genetic variants associated with alternate gene isoforms may reveal molecular mechanisms of disease. In this project, we use subcutaneous adipose tissue samples to identify splice junction quantitative trait loci (sQTL). We compare the results to both gene and exon eQTL in the same adipose tissue dataset. Finally, we identify colocalized sQTL signals with cardiometabolic traits to detect splicing events and candidate mechanisms that may contribute to gene function at GWAS loci. 

## Citation
TBD

## Table of Contents
1. METSIM data and summary statistics
2. Publicly available datasets and software
3. eQTL mapping pipeline
4. Scripts for follow-up data analysis

## 1. METSIM data and summary statistics
Please, see the METSIM_data_summaryStats directory to download the data.  
  
Included in the directory are the raw counts per gene or exon (quantified by QTLtools) for all 426 METSIM adipose tissue samples, PSI per splice junction (quantified by LeafCutter) for all METSIM samples, and the full summary statistics for the gene, exon, and splice junction QTL analyses.   


## 2. Publicly available datasets and software
Many publically available datasets and software were used throughout the completion of this project. Please see the following websites to download software and cite the corresponding papers if you use the data/software. 

### Datasets
Dataset | Website
--------|-------
GENCODE v19 (July 2013)| https://www.gencodegenes.org/human/release_19.html
GWAS Catalog (January 2020) | https://www.ebi.ac.uk/gwas/downloads
ENCODE | https://www.encodeproject.org/
GTEx v8 | https://gtexportal.org/home/datasets

### Software/Tools
Software/Tool | Website
--------------|--------
UCSC Genome Browser LiftOver | http://genome.ucsc.edu/cgi-bin/hgLiftOver
Sigora v3.0.5 | https://cran.r-project.org/web/packages/sigora/sigora.pdf
Bedtools | https://bedtools.readthedocs.io/en/latest/
ggplot2 | https://ggplot2.tidyverse.org
swiss | https://github.com/statgen/swiss
STAR v2.4.2a | https://github.com/alexdobin/STAR
LeafCutter v2.17.4 | https://davidaknowles.github.io/leafcutter/index.html
QTLtools v1.1 | https://qtltools.github.io/qtltools/
PEER | https://github.com/PMBio/peer/wiki
edgeR | https://rdrr.io/bioc/edgeR/
qvalue | https://github.com/nfusi/qvalue
locusZoom v1.4 | https://genome.sph.umich.edu/wiki/LocusZoom_Standalone

## 3. eQTL mapping pipeline
We used the same pipeline to identify gene, exon, and splice juntion QTLs throughout this project. Please, see the eQTL_pipeline directory for details.

## 4. Scripts for follow-up data analysis
