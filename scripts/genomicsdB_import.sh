#!/bin/bash

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
    --genomicsdb-workspace-path genomicsDB_mid_${CHR} \
    --sample-name-map ${CHR}_mid_sample_map \
    -L ${CHR} \
    --tmp-dir ${SLURM_TMPDIR}

gatk --java-options "-Xmx24g" GenomicsDBImport \
    --genomicsdb-workspace-path genomicsDB_high_${CHR} \
    --sample-name-map ${CHR}_high_sample_map \
    -L ${CHR} \
    --tmp-dir ${SLURM_TMPDIR}
