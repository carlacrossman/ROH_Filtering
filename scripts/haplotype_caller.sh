#!/bin/bash
#SBATCH --job-name=haplotype-gvcf
#SBATCH --account=def-pbentzen
#SBATCH --array=1-230%50
#SBATCH --output=%x-%a-%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=25G
#SBATCH --time=1:00:00

################################################################
# Load data file and relevant information needed for each sample
################################################################

LIST=$(sed -n "${SLURM_ARRAY_TASK_ID}p" data_list)
SAMPLE=$(echo $LIST | awk '{print $1}')
SCAFFOLD=$(echo $LIST | awk '{print $2}')

########################
# Load necessary modules
########################

module load gatk/4.6.1.0

##########################################
# Specify memory allocation for JAVA tools
##########################################

export JAVA_TOOL_OPTIONS="-Xmx24G"

######################################################
# Call haplotypes with GATK for each scaffold
# Repeat for each species with proper reference genome
######################################################


gatk --java-options "-Xmx24g" HaplotypeCaller \
  -R ../reference/GCF_902167405.1_gadMor3.0_genomic.fna \
  -I cod_${SAMPLE}_paired-nodups_high.bam \
  -O ../vcf/${SAMPLE}_${SCAFFOLD}_high_whole.g.vcf.gz \
  -ERC GVCF \
  -L ${SCAFFOLD}

gatk --java-options "-Xmx24g" HaplotypeCaller \
  -R ../reference/GCF_902167405.1_gadMor3.0_genomic.fna \
  -I cod_${SAMPLE}_paired-nodups_low.bam \
  -O ../vcf/${SAMPLE}_${SCAFFOLD}_low_whole.g.vcf.gz \
  -ERC GVCF \
  -L ${SCAFFOLD}
