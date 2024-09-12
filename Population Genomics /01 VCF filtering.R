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
