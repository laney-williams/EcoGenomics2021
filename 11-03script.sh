#!/bin/bash

cd /data/project_data/PopGenomics/shared/Chr05_sweeps/

CHR="Chr05"  # Be sure to customize to your chromosome number!!!

echo $CHR  # Does it look right?

# Then make some new folders to store future results in:

mkdir LAI

cd LAI/

mkdir Admixed

R

# Import K=5 Admixture run
Qscores <- read.table("/data/project_data/PopGenomics/shared/poplar_hybrids.LDpruned.5.Q", sep=" ",header=F)
names(Qscores) = c("K1","K2","K3","K4","K5")

# Import meta-data
meta <- read.table("/data/project_data/PopGenomics/Combined_Transect_Sampling_Data_2020.txt",sep="\t",header=T)

# Combine them 
merged <- cbind(meta,Qscores)
str(merged)

for(i in 1:nrow(merged)){
if(merged$K4[i]>=0.99){
    merged$Anc[i]="Bals"
    }else if (sum(c(merged$K1[i],merged$K2[i],merged$K5[i]))>=0.99){
    merged$Anc[i]="Tricho"
    } else if(merged$K3[i]>0.5){
    merged$Anc[i]="PopSpp"
    }else{
    merged$Anc[i]="Admx"
    }
}

table(merged$Anc)

# Subset individuals that just correspond to Balsam
Bals_Inds_Ref <- merged[merged$Anc=="Bals",1]  
length(Bals_Inds_Ref) # Should net you 46 individuals

write.table(Bals_Inds_Ref, "Balsam.Inds", quote=F, row.names=F, col.names=F)

Tricho_Inds_Ref <- merged[merged$Anc=="Tricho",1]
length(Tricho_Inds_Ref) # Should net you 80 individuals
write.table(Tricho_Inds_Ref, "Tricho.Inds", quote=F, row.names=F, col.names=F)

Admixed_Inds <- merged[merged$Anc=="Admx",1]
length(Admixed_Inds) # Should net you 442 individuals
write.table(Admixed_Inds, "Admixed.Inds", quote=F, row.names=F, col.names=F)

quit() # choose 'n' when prompted

# First we grab just the Balsam and Tricho reference individuals

vcftools --gzvcf /data/project_data/PopGenomics/shared/${CHR}_sweeps/${CHR}.recode.vcf --keep Balsam.Inds --recode --stdout | gzip -c >poplar_hybrids.maf05.${CHR}.BalsRef.vcf.gz

vcftools --gzvcf /data/project_data/PopGenomics/shared/${CHR}_sweeps/${CHR}.recode.vcf --keep Tricho.Inds --recode --stdout | gzip -c >poplar_hybrids.maf05.${CHR}.TrichoRef.vcf.gz

vcftools --gzvcf poplar_hybrids.maf05.${CHR}.BalsRef.vcf.gz --kept-sites --out ${CHR}

# Example while loop:

while read ID
do
  echo "$ID"
done < Admixed.Inds

# Once in a screen...

CHR="Chr05"
while read ID
do
  vcftools --gzvcf /data/project_data/PopGenomics/shared/${CHR}_sweeps/${CHR}.recode.vcf \
  --indv $ID \
  --recode \
  --stdout | gzip -c   >Admixed/poplar_hybrids.maf05.${CHR}.${ID}.vcf.gz
done < Admixed.Inds
