#!/bin/bash

#############################################
# Calculate the mean coverage for each sample 
# across each chromosome.
# Calculate the mean of the autosomes scaled 
# by chromosome length in excel.
#############################################


for ITERATOR in 2..20
LIST=$(sed -n "${ITERATOR}p" cod_bamlist)
SAMPLE=$(echo $LIST | awk '{print $1}')
FILE=$(echo $LIST | awk '{print $2}')

samtools coverage ${FILE} > ${SAMPLE}_coverage
