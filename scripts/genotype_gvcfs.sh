#!/bin/bash

###################################################
# Joint genotype calling to be done on each species
###################################################

######################################
# Load chromosome list for the species
######################################

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" chrom_list)

########################
# Load necessary modules
########################

module load gatk/4.6.1.0

##########################################
# Specify memory allocation for JAVA tools
##########################################

export JAVA_TOOL_OPTIONS="-Xmx24G"

#######################################
# Run GenotypeGVCFs for each chromosome
#######################################

gatk --java-options "-Xmx24g" GenotypeGVCFs \
   -R ../reference/GCF_902167405.1_gadMor3.0_genomic.fna \
   -V gendb://genomicsDB_low_${CHR}/ \
   -O cod_${CHR}_low.vcf.gz
   
gatk --java-options "-Xmx24g" GenotypeGVCFs \
   -R ../reference/GCF_902167405.1_gadMor3.0_genomic.fna \
   -V gendb://genomicsDB_mid_${CHR}/ \
   -O cod_${CHR}_mid.vcf.gz

gatk --java-options "-Xmx24g" GenotypeGVCFs \
   -R ../reference/GCF_902167405.1_gadMor3.0_genomic.fna \
   -V gendb://genomicsDB_high_${CHR}/ \
   -O cod_${CHR}_high.vcf.gz
