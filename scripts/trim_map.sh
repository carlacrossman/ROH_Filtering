#!/bin/bash

################################################################
# Load data file and relevant information needed for each sample
################################################################

FILES=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_info)
SAMPLE=$(echo $FILES | awk '{print $1}')
SPECIES=$(echo $FILES | awk '{print $2}')
REFERENCE=$(echo $FILES | awk '{print $3}')

#######################################
# Copy relevant files to temp directory
#######################################

cp ${SPECIES}/data/${SAMPLE}*.fastq.gz ${SLURM_TMPDIR}
cp ${SPECIES}/reference/GC* ${SLURM_TMPDIR}
cp adapters.fa ${SLURM_TMPDIR}

cd ${SLURM_TMPDIR}

########################
# Load necessary modules
########################

module load bwa/0.7.18 samtools/1.22.1 trimmomatic/0.41

##########################################
# Specify memory allocation for JAVA tools
##########################################

export JAVA_TOOL_OPTIONS="-Xmx9g"

#########################################
# Run Trimmomatic for loose filters 'low'
#########################################

java -Xmx9G -jar $EBROOTTRIMMOMATIC/trimmomatic-0.41.jar PE -threads 14 -phred33 \
    ${SAMPLE}_1.fastq.gz \
    ${SAMPLE}_2.fastq.gz \
    ${SPECIES}_${SAMPLE}_FP_low.fq.gz \
    ${SPECIES}_${SAMPLE}_FUP_low.fq.gz \
    ${SPECIES}_${SAMPLE}_RP_low.fq.gz \
    ${SPECIES}_${SAMPLE}_RUP_low.fq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 SLIDINGWINDOW:4:20 \
    MINLEN:36
    
#########################################
# Run Trimmomatic for mid filters 'mid'
#########################################

java -Xmx9G -jar $EBROOTTRIMMOMATIC/trimmomatic-0.41.jar PE -threads 14 -phred33 \
    ${SAMPLE}_1.fastq.gz \
    ${SAMPLE}_2.fastq.gz \
    ${SPECIES}_${SAMPLE}_FP_mid.fq.gz \
    ${SPECIES}_${SAMPLE}_FUP_mid.fq.gz \
    ${SPECIES}_${SAMPLE}_RP_mid.fq.gz \
    ${SPECIES}_${SAMPLE}_RUP_mid.fq.gz \
    ILLUMINACLIP:adapters.fa:2:30:15 LEADING:20 SLIDINGWINDOW:5:25 \
    AVGQUAL:25 MINLEN:36

##############################################
# Run Trimmomatic for stringent filters 'high'
##############################################

java -Xmx9G -jar $EBROOTTRIMMOMATIC/trimmomatic-0.41.jar PE -threads 14 -phred33 \
    ${SAMPLE}_1.fastq.gz \
    ${SAMPLE}_2.fastq.gz \
    ${SPECIES}_${SAMPLE}_FP_high.fq.gz \
    ${SPECIES}_${SAMPLE}_FUP_high.fq.gz \
    ${SPECIES}_${SAMPLE}_RP_high.fq.gz \
    ${SPECIES}_${SAMPLE}_RUP_high.fq.gz \
    ILLUMINACLIP:adapters.fa:2:40:15 LEADING:25 SLIDINGWINDOW:5:30 \
    AVGQUAL:30 MINLEN:36


########################################
# Map trimmed reads to reference genomes
########################################

bwa mem -M -t 15 \
 -R '@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SPECIES}\tCN:NCBI\tPL:ILLUMINA' \
  ${REFERENCE} \
  ${SPECIES}_${SAMPLE}_FP_low.fq.gz \
  ${SPECIES}_${SAMPLE}_RP_low.fq.gz > \
  ${SPECIES}_${SAMPLE}_paired_low.sam

samtools sort -o ${SPECIES}_${SAMPLE}_paired_low.bam ${SPECIES}_${SAMPLE}_paired_low.sam
rm ${SPECIES}_${SAMPLE}_paired_low.sam

bwa mem -M -t 15 \
 -R '@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SPECIES}\tCN:NCBI\tPL:ILLUMINA' \
  ${REFERENCE} \
  ${SPECIES}_${SAMPLE}_FP_mid.fq.gz \
  ${SPECIES}_${SAMPLE}_RP_mid.fq.gz > \
  ${SPECIES}_${SAMPLE}_paired_mid.sam

samtools sort -o ${SPECIES}_${SAMPLE}_paired_mid.bam ${SPECIES}_${SAMPLE}_paired_mid.sam
rm ${SPECIES}_${SAMPLE}_paired_mid.sam

bwa mem -M -t 15 \
 -R '@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SPECIES}\tCN:NCBI\tPL:ILLUMINA' \
  ${REFERENCE} \
  ${SPECIES}_${SAMPLE}_FP_high.fq.gz \
  ${SPECIES}_${SAMPLE}_RP_high.fq.gz > \
  ${SPECIES}_${SAMPLE}_paired_high.sam

samtools sort -o ${SPECIES}_${SAMPLE}_paired_high.bam ${SPECIES}_${SAMPLE}_paired_high.sam
rm ${SPECIES}_${SAMPLE}_paired_high.sam

####################
# Load Picard module
####################

# Picard requires an older java version 
# that forces an older version of trimmomatic
# Load here to avoid problems

module load picard/3.1.0

######################
# Mark duplicate reads
######################

java -Xmx9G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  --REMOVE_DUPLICATES true --CREATE_INDEX true \
  --INPUT ${SPECIES}_${SAMPLE}_paired_low.bam \
  --OUTPUT ${SPECIES}_${SAMPLE}_paired-nodups_low.bam \
  --METRICS_FILE ${SPECIES}_${SAMPLE}_paired-nodups_low.metrics
  
java -Xmx9G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  --REMOVE_DUPLICATES true --CREATE_INDEX true \
  --INPUT ${SPECIES}_${SAMPLE}_paired_mid.bam \
  --OUTPUT ${SPECIES}_${SAMPLE}_paired-nodups_mid.bam \
  --METRICS_FILE ${SPECIES}_${SAMPLE}_paired-nodups_mid.metrics

java -Xmx9G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  --REMOVE_DUPLICATES true --CREATE_INDEX true \
  --INPUT ${SPECIES}_${SAMPLE}_paired_high.bam \
  --OUTPUT ${SPECIES}_${SAMPLE}_paired-nodups_high.bam \
  --METRICS_FILE ${SPECIES}_${SAMPLE}_paired-nodups_high.metrics


###########################################
# Copy finished files to relevant directory 
###########################################

cp ${SPECIES}_${SAMPLE}_paired-nodups* ${SLURM_SUBMIT_DIR}/${SPECIES}/bam_files/
