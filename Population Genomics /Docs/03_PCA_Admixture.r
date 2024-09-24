#Run a PCA

library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

options(bitmapType = "cairo")

setwd("~/Projects/eco_genomics_2/Population Genomics /")
vcfR <- read.vcfR("Outputs/vcf_final.filtered.vcf.gz")
#set wd, path picks up where working directory left off 
#we need to thin the SNPs for linkage disequilibrium before we run PCA and admixture to satisfy the assumptions of independence among loci

vcf.thin <- distance_thin(vcf,min.distance = 500)
#make a new vcf file based on a min distance of 500 to thin 

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
dim(meta)
meta2 <- meta[meta$id %in% colnames(vcf.thin@gt[,-1]),]
dim(meta2)
write.vcf(vcf.thin, "Outputs/vcf_final_filtered.thinned.vcf.gz")

#hide the uncompressed vcf too big for github outside of repo

system("gunzip -c ~/Projects/eco_genomics_2/Population\ Genomics\ /Outputs/vcf_final_filtered.thinned.vcf.gz > ~/vcf_final_filtered.thinned.vcf")

geno <- vcf2geno(input.file="/gpfs1/home/s/r/kgusek/vcf ")
centPCA <- LEA::pca("~/vcf_final_filtered.thinned.vcf", scale=TRUE)