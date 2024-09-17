library(vcfR)
install.packages("ape")
library(vcfR)
install.packages("vcfR")
library(vcfR)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

list.files()

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")

head(vcf)

dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")

gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t",quote="")
chr1 <- create.chromR(name="Chromosome 1", vcf=vcf, seq=dna, ann=gff)
plot(chr1)
pdf(file="~/Projects/eco_genomics_2/Population Genomics /Figures/ChromoPlot.pdf")
chromoqc(chr1, xlim=c(1e1, 1.1e8))
dev.off()
library(vcfR)
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
list.files("variants/")
vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")
vcf
head(vcf)
DP <- extract.gt(vcf, element="DP", as.numeric=T) #extract elements from vcf files 
DP[1:5,1:10] #number of reads total for that chromosome in that position 
dim(DP)
quantile(DP)
DP[DP==0] <- NA
quantile(DP)
quantile(DP, na.rm=T)
#visualize the matrix of depth and missingness in the vcf file
heatmap.bp(DP)
heatmap.bp(DP[1:1000,], rlabels=F, clabels=F)
library(SNPfiltR)
hard_filter(vcf) #histogram of DP values 
vcf.filt <- hard_filter(vcf, depth=3) 
#if you dont have at least three reads at a given SNP position for an individual set it at NA
#could explore values DP 5, 10...
max_depth(vcf.filt) # max depth filter 
vcf.file <- max_depth(vcf.filt, maxdepth=60)#filter out genotype with >60 reads/SNPs
meta <- read.csv("metadata/meta4vcf.csv", header=T)
meta2 <- meta[,c(1,4)]
names(meta2) <- c("id","pop")
meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta2$pop)
vcf.filt.indMiss <- missing_by_sample(vcf.filt,
                                      popmap = meta2,
                                      cutoff=0.75)
#filtering at 75%, drop individuals from 5/6 pops
vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)
vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, cutoff=0.5)
DP2 <- extract.gt(vcf.filt.indSNPMiss,
                  element="DP",
                  as.numeric=T)

heatmap.bp(DP2[1:5000,],
           rlabels=F, clabels=F)
write.vcf(vcf.filt.indSNPMiss,
          "~/Projects/eco_genomics_2/Population Genomics /Outputs/vcf_final.filtered.vcf.gz")
