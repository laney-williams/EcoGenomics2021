library(ggplot2)
library(gridExtra)

# Get the list of admixed individuals:
Admixed <- read.table("Admixed.Inds",header=F)

# Get the meta data:
meta <- read.table("Combined_Transect_Sampling_Data_2020.txt", sep="\t",header=T)

# merge them together:
meta_admx <- merge(meta, Admixed, by.x="ID", by.y="V1")
str(meta_admx)  

# Read in the Admixture coefficients for KBals that we made from the K=5 file:
KBals <- read.table("Admixed_KBals", sep="\t", header=F)
names(KBals) = c("ID","KBals")

# Second merge:
meta_admx_KBals <- merge(meta_admx,KBals,by="ID")


# Bring in phenotype data:
pheno <- read.table("climDat.txt",sep="\t",header=T)

# Merge pheno data with meta and KBals:
meta_admx_KBals_pheno <- merge(meta_admx_KBals,pheno,by="ID")

######  Bring in Association results from Plink   ######

#  Growing Degree Days
plotGrowingDegree <- ggplot(meta_admx_KBals_pheno,aes(x=KBals,y=mean_cGDDfreeze, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Mean No. of GrowingDegree Days") 

plotGrowingDegree

# Chilling Days
plotChillingDays <- ggplot(meta_admx_KBals_pheno,aes(x=KBals,y=med_DD0, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Mean No. of Chilling Degree Days") 

plotChillingDays

# Final Freeze
plotFinalFreeze <- ggplot(meta_admx_KBals_pheno,aes(x=KBals,y=mean_finalFreeze, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Average Date of Last Freeze") 

plotFinalFreeze

grid.arrange(plotGrowingDegree, plotChillingDays, plotFinalFreeze, nrow = 3)

# linear models testing trait ~ genome-wide admixture association
summary(lm(mean_cGDDfreeze~KBals + Transect.x, data=meta_admx_KBals_pheno))

summary(lm(med_DD0~KBals + Transect.x, data=meta_admx_KBals_pheno))

summary(lm(mean_finalFreeze~KBals + Transect.x, data=meta_admx_KBals_pheno))

### GLM using results from Plink

### Growing Degree ###
GrowingDegree <- read.table("plink2.mean_cGDDfreeze.glm.linear",skip=1,sep="\t",header=F)
names(GrowingDegree) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
GrowingDegree2 <- GrowingDegree[which(GrowingDegree$TEST=="ADD"),]

# Define association outliers as the upper 1% of p-values

#########  Growing Degree Days  #########
snps <- read.table("Chr05.kept.sites",sep="\t",header=T)

GrowingDegree2 <- cbind(snps, GrowingDegree2[,-c(1:2)])
GrowingDegree2$outlier = ifelse(GrowingDegree2$P<quantile(GrowingDegree2$P,0.01),2,1)

p1 <- ggplot(GrowingDegree2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=GrowingDegree2$outlier, color=GrowingDegree2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Growing Degree Days")

p1

####### Chilling Days  #########
ChillingDays <- read.table("plink2.med_DD0.glm.linear",skip=1,sep="\t",header=F)
names(ChillingDays) = c("CHROM",  "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
ChillingDays2 <- ChillingDays[which(ChillingDays$TEST=="ADD"),]
ChillingDays2 <- cbind(snps, ChillingDays2[,-c(1,2)])
ChillingDays2$outlier = ifelse(ChillingDays2$P<quantile(ChillingDays2$P,0.01),2,1)

p2 <- ggplot(ChillingDays2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=ChillingDays2$outlier, color=ChillingDays2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Chilling Degree Days")

p2

##### Date of Last Freeze ######

FinalFreeze <- read.table("plink2.mean_finalFreeze.glm.linear",skip=1,sep="\t",header=F)
names(FinalFreeze) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
FinalFreeze <- FinalFreeze[which(FinalFreeze$TEST=="ADD"),]
FinalFreeze2 <- cbind(snps, FinalFreeze[,-c(1,2)])
FinalFreeze2$outlier = ifelse(FinalFreeze2$P<quantile(FinalFreeze2$P,0.01),2,1)

p3 <- ggplot(FinalFreeze2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=FinalFreeze2$outlier, color=FinalFreeze2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Date of Last Freezing Event")

p3

grid.arrange(p1, p2, p3, nrow = 3)

#### Bud Flush #####

budflush <- read.table("plink2.FLUSH.glm.linear",skip=1,sep="\t",header=F)
names(budflush) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
budflush <- budflush[which(budflush$TEST=="ADD"),]
budflush2 <- cbind(snps, budflush[,-c(1,2)])
budflush2$outlier = ifelse(budflush2$P<quantile(budflush2$P,0.01),2,1)

p4 <- ggplot(budflush2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=budflush2$outlier, color=budflush2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Bud flush")

p4

grid.arrange(p3, p4, nrow = 2)

### Genomic Ranges ###

library(GenomicRanges)
library(GenomicFeatures)

# Define windows where peaks of association fall

CHR="Chr05"

### Bud Flush GR
budflushGR <- GRanges(CHR,IRanges(budflush2$POS-2.5e4,budflush2$POS+2.5e4),POS=budflush2$POS, P=budflush2$P, outlier=budflush2$outlier)

budflushGRout <- unlist(reduce(split(budflushGR, ~outlier)))
budflushGRout$outlier <- names(budflushGRout)
budflushGRCand <- subset(budflushGRout, outlier==2)

budflushGRCand # Print the candidate regions

### Growing Degree Days

GrowingDegreeGR <- GRanges(CHR,IRanges(GrowingDegree2$POS-2.5e4,GrowingDegree2$POS+2.5e4),POS=GrowingDegree2$POS, P=GrowingDegree2$P, outlier=GrowingDegree2$outlier)

GrowingDegreeGRout <- unlist(reduce(split(GrowingDegreeGR, ~outlier)))
GrowingDegreeGRout$outlier <- names(GrowingDegreeGRout)
GrowingDegreeGRCand <- subset(GrowingDegreeGRout, outlier==2)

GrowingDegreeGRCand # Print the candidate regions

### Chilling Degree Days

ChillingDaysGR <- GRanges(CHR,IRanges(ChillingDays2$POS-2.5e4,ChillingDays2$POS+2.5e4),POS=ChillingDays2$POS, P=ChillingDays2$P, outlier=ChillingDays2$outlier)

ChillingDaysGRout <- unlist(reduce(split(ChillingDaysGR, ~outlier)))
ChillingDaysGRout$outlier <- names(ChillingDaysGRout)
ChillingDaysGRCand <- subset(ChillingDaysGRout, outlier==2)

ChillingDaysGRCand # Print the candidate regions

### Date of Final Freeze

FinalFreezeGR <- GRanges(CHR,IRanges(FinalFreeze2$POS-2.5e4,FinalFreeze2$POS+2.5e4),POS=FinalFreeze2$POS, P=FinalFreeze2$P, outlier=FinalFreeze2$outlier)

FinalFreezeGRout <- unlist(reduce(split(FinalFreezeGR, ~outlier)))
FinalFreezeGRout$outlier <- names(FinalFreezeGRout)
FinalFreezeGRCand <- subset(FinalFreezeGRout, outlier==2)

FinalFreezeGRCand # Print the candidate regions

### Look at overlaps between climate variable and bud flush

# Bud Flush and Growing Degree Days
overlap_BF_GrowingDegree <- subsetByOverlaps(budflushGRCand, GrowingDegreeGRCand)
length(overlap_BF_GrowingDegree)

overlap_BF_GrowingDegree # Print the overlapping regions

### No overlap between Bud Flush and Growing Degree Days

# Bud Flush and Chilling Days
overlap_BF_ChillingDays <- subsetByOverlaps(budflushGRCand, ChillingDaysGRCand)
length(overlap_BF_ChillingDays)

overlap_BF_ChillingDays # Print the overlapping regions

## No overlap between bud flush and chilling degree days

# Bud Flush and date of Final Freeze
overlap_BF_FinalFreeze <- subsetByOverlaps(budflushGRCand, FinalFreezeGRCand)
length(overlap_BF_FinalFreeze)

overlap_BF_FinalFreeze

### 1 region of overlap between bud flush and date of final freeze!!

### Create a transcript database ###
# Import the GFF annotation file and make a transcript database
txdb <- makeTxDbFromGFF("Ptrichocarpa_533_v4.1.gene.gff3.gz", format="gff3")

txdb

# How many chromosomes are present?
head(seqlevels(txdb))

# Subset the database for just your chromosome of interest
seqlevels(txdb) <- CHR # subset for just your chromosome

# Reduce the transcript database to just the non-redundant gene names, instead of multiple entries for all the variant transcript types per gene
genes <- unlist(reduce(transcriptsBy(txdb, by="gene"))) 
genes$geneID <- names(genes)


candGenes <- subsetByOverlaps(genes, overlap_BF_FinalFreeze)

write.table(candGenes$geneID, paste0("candGenes",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")
