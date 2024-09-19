# Estimating diversity and genetig differentiation in the filtered centaurea data 

library(vcfR)
library(tidyverse)
library(qqman)

#helps solve plotting issues 
X11.options(type="cairo")

#read in our new vcf file from our repo outputs/directory
vcf <- read.vcfR("~/Projects/eco_genomics_2/Population Genomics /Outputs/vcf_final.filtered.vcf.gz")
#vcf <- read.vcfR("/gpfs1/cl/pbio3990/PopulationGenomics/variants/Centaurea_filtered.vcf.gz")
#read our metadata- info on pop origin and region 
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

head(meta) #vcf files has 595 samples
dim(meta) #meta has 629 inds

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),] #before comma is the rows and after the comma is the columns
#@ is a list 

dim(meta2)

vcf.div <- genetic_diff(vcf,
                        pops=as.factor(meta2$region),
                        method = "nei")
str(vcf.div)#tells you the structure of the data frame as apposed to head showing you the first few lines 

#CM is the real chromosome assembly 
#IDs that start with JARYM are scaffolds

chr.main <- unique(vcf.div$CHROM)[1:8] #assigning new variable
chr.main


chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1))) #assign numbers to each chromosome
chrnum

vcf.div.MHplot <- left_join(chrnum, vcf.div, join_by(chr.main==CHROM))

vcf.div.MHplot <- vcf.div.MHplot %>% 
                        filter(Gst>0) %>%
                        mutate(SNP=paste0(chr.main,"_",POS))
vcf.div.MHplot$V2 = as.numeric(vcf.div.MHplot$V2)
vcf.div.MHplot$POS = as.numeric(vcf.div.MHplot$POS)

manhattan(vcf.div.MHplot,
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue1","orange3"),
          logp=F,
          ylab="Fst among regions",
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.999))
#each dot is a SNP
#each color is an alternating chromosome
#amt of non random mating that is occurring 
# the V is most likely the centromere, they are difficult to sequence, dip in SNP diversity that is a tell tale sign that you are at a contromere
