#!/bin/bash

cd /data/project_data/PopGenomics/shared/Chr05_sweeps/LAI

CHR="Chr05"

# First, make a new dir to store the results:
mkdir Loter_out

while read ID
do
for file in Admixed/poplar_hybrids.maf05.${CHR}.${ID}.vcf.gz
do
loter_cli -r poplar_hybrids.maf05.${CHR}.BalsRef.vcf.gz poplar_hybrids.maf05.${CHR}.TrichoRef.vcf.gz \
-a $file \
-f vcf \
-o Loter_out/${ID}_LAI.txt \
-n 1 \
-pc -v
done
done < Admixed.Inds