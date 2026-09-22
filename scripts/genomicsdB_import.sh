#!/bin/bash
#SBATCH --job-name=genomicsdb
#SBATCH --account=def-pbentzen
#SBATCH --array=2-23
#SBATCH --output=%x-%a-%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --time=8:00:00

####################################
# This is to be run for each species
####################################

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

########################################
# Run GenomicsDBImport for each scaffold
# Repeat for each species
########################################

gatk --java-options "-Xmx24g" GenomicsDBImport \
    --genomicsdb-workspace-path genomicsDB_low_${CHR} \
    --sample-name-map ${CHR}_low_sample_map \
    -L ${CHR} \
    --tmp-dir ${SLURM_TMPDIR}

gatk --java-options "-Xmx24g" GenomicsDBImport \
    --genomicsdb-workspace-path genomicsDB_high_${CHR} \
    --sample-name-map ${CHR}_high_sample_map \
    -L ${CHR} \
    --tmp-dir ${SLURM_TMPDIR}
