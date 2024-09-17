# Coding and data notes for the population genomics module

## Author: Kaeli Gusek

### 09/10/2024: Intro to Centaurea GBS data and working with VCF files

We'll be analyzing the GBS data from three regions (EU, NE, and PNW) starting today with variant call format files (VCFs)

Data was viewed in the command line in a fastq file. The letters represented the Q scores of the data (I was 1 in 10,000 Bp that there is an error), showed accurate data

cd: change directory

ls -l: long list (can be abbreviated as ll)

zcat: View a zip file

| : "pipe" send the output of one function to another

Head: view the first 10 lines of code

CCACA: bar code for that specific sequence

If you are in the correct directory with a specific file or directory you can use TAB to auto complete the file name or directory

### 9/12/2024 View VCF file and talk about visualizing

so moving forward all points will direct to /gpfs1/cl/pbio3990

sam: sequence alignment file

bam: binary version

Bai: index file

bcftools: takes each bam file and looks at where there is likely a SNP

setwd: set working directory

plot(chr1): plots for the data

chromoqc(chr1, xlim=c(1e1, 1.1e8)): plots from one telomere to the other, each dot is a SNP in the genome

pdf: pdf writer

### 9/17/24 Filtering Strategies 

Filtered vcf files based on missingness and depth, cut off outliers

0/0 reference homozygote

0/1 Heterozygote

./. N/A

Low depth: limited accuracy on the nucleotide in that position/what the genotype call should be

High depth: assembly error in the ref genome or gene paralogy (reads from duplicated genes mapping to the same position)

Missingness: individual level or SNP level (depending on site, may be poorly represented)

Low freq alleles: less than 1% in the population

The white on the heat map is missing data

row is individuals

columns are sequence data

df[ rows, column] leaving it blank will leave all of the rows and columns

? in console bring up an index to give information on a function

DP \<- extract.gt(vcf, element="DP", as.numeric=T): extract elements from vcf files DP[1:5,1:10]: number of reads total for that chromosome in that position

quantile(DP, na.rm=T): visualize the matrix of depth and missingness in the vcf file

Cutoff: set a cutoff value for individuals you do not want to include based on depth and missingness

if you dont have at least three reads at a given SNP position for an individual set it at NA

max_depth(vcf.filt): max depth filter

write.vcf(vcf.filt.indSNPMiss, "\~/Projects/eco_genomics_2/Population Genomics /Outputs/vcf_final.filtered.vcf.gz"): save files, make sure to include file types (vcf.gz)
