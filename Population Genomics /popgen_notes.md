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
