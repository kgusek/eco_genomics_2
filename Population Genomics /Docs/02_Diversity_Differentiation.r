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
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.999)) #changes the quantline where the line rests (0.5 shows the median) the amount of Fst is extremely low 
#each dot is a SNP
#each color is an alternating chromosome
#amt of non random mating that is occurring 
# the V is most likely the centromere, they are difficult to sequence, dip in SNP diversity that is a tell tale sign that you are at a contromere

write.csv(vcf.div.MHplot, "~/Projects/eco_genomics_2/Population Genomics /Outputs/Genetic_Diff_byRegion.csv",
          quote=F,
          row.names=F) #saves the data table that made the manhattan plot to outputs 
#9/24/23
#look at the genetic diversity that is maintained between these groups 

#where the Hs values are contained
names(vcf.div.MHplot)
#Hs values are stored in rows 4-9
options(bitmapType = "cairo")

#use tidyverse to take the file that has each Hs value in a separate column and put it into one column for all of the Hs values 
vcf.div.MHplot %>%
  as_tibble()%>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) +
  geom_histogram(position="identity", alpha=0.5, bins=50)+
  labs(title="genmoe-wide expected heterozygosity (Hs)", fill="regions",
       x="Gene diversity within Regions", y="Counts of SNPs")

#make a ggplot with all of the Hs values and color them by the names of different columns 
#pipe the results from one step to the next %>%

ggsave("Histogram_Genome_Diversity_by_Region.pdf",
       path="~/Projects/eco_genomics_2/Population Genomics /Figures")
#ggsave saves the last plot made

vcf.div.MHplot %>%
  as_tibble()%>%
  pivot_longer(c(4:9)) %>%
  group_by(name) %>%
  filter(value!=0) %>%
  summarise(avgHs=mean(value), StdDev=sd(value), N_Hs=n())
#table with each of the different regions and its genome wide expected heterozygosity 
