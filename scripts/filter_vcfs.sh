#!/bin/bash

#############################################
# Use VCFTOOLS to drop genotypes based on 
# depth and genotype quality.
#############################################

module load vcftools/0.1.16 

for SPECIES in cod rattlesnake rightwhale
do

vcftools --gzvcf merged_${SPECIES}_raw_allsites_high.vcf.gz --recode --recode-INFO-all --out merged_${SPECIES}_dp_gq_allsites_high --minDP 10 --minGQ 30

vcftools --gzvcf merged_${SPECIES}_raw_allsites_mid.vcf.gz --recode --recode-INFO-all --out merged_${SPECIES}_dp_gq_allsites_mid --minDP 8 --minGQ 25

vcftools --gzvcf merged_${SPECIES}_raw_allsites_low.vcf.gz --recode --recode-INFO-all --out merged_${SPECIES}_dp_gq_allsites_low --minDP 5 --minGQ 20

gzip merged_${SPECIES}_dp_gq_allsites_high.recode.vcf
gzip merged_${SPECIES}_dp_gq_allsites_mid.recode.vcf
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

bcftools filter -i 'COUNT(FORMAT/GT="mis")<3 && INFO/MQ>25' -Oz -o merged_${SPECIES}_filtered_allsites_mid.recode.vcf.gz merged_${SPECIES}_dp_gq_allsites_mid.recode.vcf.gz

bcftools filter -i 'COUNT(FORMAT/GT="mis")<5 && INFO/MQ>20' -Oz -o merged_${SPECIES}_filtered_allsites_low.recode.vcf.gz merged_${SPECIES}_dp_gq_allsites_low.recode.vcf.gz

bcftools view -m2 -M2 -v snps -Oz -o merged_${SPECIES}_filtered_variants_low.recode.vcf.gz merged_${SPECIES}_filtered_allsites_low.recode.vcf.gz  

bcftools view -m2 -M2 -v snps -Oz -o merged_${SPECIES}_filtered_variants_mid.recode.vcf.gz merged_${SPECIES}_filtered_allsites_mid.recode.vcf.gz 

bcftools view -m2 -M2 -v snps -Oz -o merged_${SPECIES}_filtered_variants_high.recode.vcf.gz merged_${SPECIES}_filtered_allsites_high.recode.vcf.gz 

bcftools +fill-tags merged_${SPECIES}_filtered_variants_high.recode.vcf.gz -Oz -o merged_${SPECIES}_filtered_variants_updatetags_high.vcf.gz -- -t all

bcftools +fill-tags merged_${SPECIES}_filtered_variants_mid.recode.vcf.gz -Oz -o merged_${SPECIES}_filtered_variants_updatetags_mid.vcf.gz -- -t all

bcftools +fill-tags merged_${SPECIES}_filtered_variants_low.recode.vcf.gz -Oz -o merged_${SPECIES}_filtered_variants_updatetags_low.vcf.gz -- -t all

done

######################################################
# Performing MAF filtering on MID samples with MAC >=2
######################################################

# Note this was changed to MAC>=2 to account for missing data

module load vcftools/0.1.16 

for SPECIES in cod rattlesnake rightwhale
do

vcftools --gzvcf merged_${SPECIES}_filtered_variants_updatetags_mid.vcf.gz --recode --recode-INFO-all --out merged_${SPECIES}_filtered_variants_mid_maf --maf 0.1

gzip merged_${SPECIES}_filtered_variants_mid_maf.recode.vcf

done

####################################################
# Repeat Masking
####################################################

module load seqkit/2.5.1 bedtools/2.31.0

seqkit locate --bed -r -P -p "[a-z]+" ../cod/reference/GCF_902167405.1_gadMor3.0_genomic.fna | bedtools merge -i - > cod_repeats.bed

seqkit locate --bed -r -P -p "[a-z]+" GCA_039880765.1_rSisCat1_p1.0_genomic.fna | bedtools merge -i - > rattlesnake_repeats.bed

sed -e '/^>/! s/[[:lower:]]/Z/g' GCF_028564815.1_mEubGla1.1.hap2._XY_genomic.fna > hardmask_Z_rightwhale_reference.fasta

seqkit locate --bed -r -P -p "Z" ../rightwhale/reference/hardmask_Z_rightwhale_reference.fasta | bedtools merge -i - > rightwhale_repeats.bed


bedtools intersect -v -a ../raw_vcf/merged_cod_filtered_variants_updatetags_mid.vcf.gz -b cod_repeats.bed -wa -header > merged_cod_filtered_variants_mid_repeatmasked.vcf

bedtools intersect -v -a ../raw_vcf/merged_rattlesnake_filtered_variants_updatetags_mid.vcf.gz -b rattlesnake_repeats.bed -wa -header > merged_rattlesnake_filtered_variants_mid_repeatmasked.vcf

bedtools intersect -v -a ../raw_vcf/merged_rightwhale_filtered_variants_updatetags_mid.vcf.gz -b rightwhale_repeats.bed -wa -header > merged_rightwhale_filtered_variants_mid_repeatmasked.vcf

gzip merged_cod_filtered_variants_mid_repeatmasked.vcf
gzip merged_rattlesnake_filtered_variants_mid_repeatmasked.vcf
gzip merged_rightwhale_filtered_variants_mid_repeatmasked.vcf
