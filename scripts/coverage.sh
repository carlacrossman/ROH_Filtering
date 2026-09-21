#!/bin/bash

################################################################
# Load data file and relevant information needed for each sample
################################################################


for ITERATOR in 2..20
LIST=$(sed -n "${ITERATOR}p" cod_bamlist)
SAMPLE=$(echo $LIST | awk '{print $1}')
FILE=$(echo $LIST | awk '{print $2}')

samtools coverage ${FILE} > ${SAMPLE}_coverage
