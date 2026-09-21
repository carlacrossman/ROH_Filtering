#!/bin/bash
#SBATCH --job-name=genotype_gvcf
#SBATCH --account=def-pbentzen
#SBATCH --array=1-23
#SBATCH --output=%x-%a-%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --time=12:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" chrom_list)

export JAVA_TOOL_OPTIONS="-Xmx24G"

module load gatk/4.6.1.0

gatk --java-options "-Xmx24g" GenotypeGVCFs \
   -R ../reference/GCF_902167405.1_gadMor3.0_genomic.fna \
   -V gendb://genomicsDB_low_${CHR}/ \
   -O cod_${CHR}_low.vcf.gz

gatk --java-options "-Xmx24g" GenotypeGVCFs \
   -R ../reference/GCF_902167405.1_gadMor3.0_genomic.fna \
   -V gendb://genomicsDB_high_${CHR}/ \
   -O cod_${CHR}_high.vcf.gz
