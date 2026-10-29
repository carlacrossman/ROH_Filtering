#!/bin/bash

module load bcftools/1.22
module load plink/2.00-20231024-avx2

####################
# Define variables #
####################

PREFIX=$(basename $1)

bcftools tabix 00_data/$PREFIX.vcf.gz

##########################################################################
# Rename chromosome and remove scaffolds, INDELS and triallelic variants #
##########################################################################

zcat "00_data/$PREFIX.vcf.gz" | head -n100 |  grep '^##contig' | grep -Eo '(NC_|CM)[0-9]+\.1' | \
  awk  'BEGIN { OFS="\t" } !seen[$1]++ {print $1, "chr" ++n}' > 00_data/$PREFIX.chrmap

bcftools annotate --rename-chrs 00_data/$PREFIX.chrmap \
    -r "$(cut -f1 00_data/$PREFIX.chrmap | paste -sd,)" \
    -Ou 00_data/"$PREFIX".vcf.gz \
    | bcftools view -v snps -m2 -M2 \
    -Oz -o 00_data/"$PREFIX"_chr.vcf.gz 
    
bcftools index -t 00_data/"$PREFIX"_chr.vcf.gz

NB_CHR=$(wc -l 00_data/$PREFIX.chrmap | cut -f1 -d ' ')

######################
# Convert VCF to PED #
######################

plink2 --vcf 00_data/"$PREFIX"_chr.vcf.gz \
    --export ped \
    --allow-extra-chr 0 \
    --out 00_data/"$PREFIX"_chr \
    --chr-set $NB_CHR
