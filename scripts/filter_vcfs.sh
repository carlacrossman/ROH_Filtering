#!/bin/bash

#############################################
# Use VCFTOOLS to drop genotypes based on 
# depth and genotype quality.
#############################################

module load vcftools/0.1.16 

for SPECIES in cod rattlesnake rightwhale
do

vcftools --gzvcf merged_${SPECIES}_raw_allsites_high.vcf.gz --recode --recode-INFO-all --out merged_${SPECIES}_dp_gq_allsites_high --minDP 10 --minGQ 30

vcftools --gzvcf merged_${SPECIES}_raw_allsites_low.vcf.gz --recode --recode-INFO-all --out merged_${SPECIES}_dp_gq_allsites_low --minDP 5 --minGQ 20

gzip merged_${SPECIES}_dp_gq_allsites_high.recode.vcf
gzip merged_${SPECIES}_dp_gq_allsites_low.recode.vcf

done


####################################################
# Loop through each species 
# 1) Use BCFTOOLS to drop sites based on missingness
#    and mapping quality.
# 2) Restrict the output to bi-allelic SNPs.
# 3) Update the INFO fields with +fill-tags.
####################################################

module load bcftools/1.22

for SPECIES in cod rattlesnake rightwhale
do

bcftools filter -i 'COUNT(FORMAT/GT="mis")<2 && INFO/MQ>30' -Oz -o merged_${SPECIES}_filtered_allsites_high.recode.vcf.gz merged_${SPECIES}_dp_gq_allsites_high.recode.vcf.gz

bcftools filter -i 'COUNT(FORMAT/GT="mis")<5 && INFO/MQ>20' -Oz -o merged_${SPECIES}_filtered_allsites_low.recode.vcf.gz merged_${SPECIES}_dp_gq_allsites_low.recode.vcf.gz

bcftools view -m2 -Oz -o merged_${SPECIES}_filtered_variants_low.recode.vcf.gz merged_${SPECIES}_filtered_allsites_low.recode.vcf.gz  

bcftools view -m2 -M2 -v snps -Oz -o merged_${SPECIES}_filtered_variants_high.recode.vcf.gz merged_${SPECIES}_filtered_allsites_high.recode.vcf.gz 

bcftools +fill-tags merged_${SPECIES}_filtered_variants_high.recode.vcf.gz -Oz -o merged_${SPECIES}_filtered_variants_updatetags_high.vcf.gz -- -t all

bcftools +fill-tags merged_${SPECIES}_filtered_variants_low.recode.vcf.gz -Oz -o merged_${SPECIES}_filtered_variants_updatetags_low.vcf.gz -- -t all

done
