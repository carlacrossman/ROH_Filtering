#!/bin/bash

module load bcftools/1.22

PREFIX=$(basename $1)

bcftools roh 00_data/"$PREFIX"_chr.vcf.gz \
    -Orz -o 03_bcftools/"$PREFIX".bcftoolsROH.gz
