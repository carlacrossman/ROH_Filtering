############################
### --- Analysis Log --- ###
############################

This is a running log of what I am doing. It began as a bit of chicken scratch as a proof of concept, but evolved into more strategic thinking around line 500.



**Trimming and Mapping**

Set two thresholds high and low.

```
#!/bin/bash
#SBATCH --account=def-pbentzen
#SBATCH --job-name=trim_map
#SBATCH --array=1-10
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=15
#SBATCH --time=4:00:00

PREFIX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" fastq_list)
REFERENCE=~/links/scratch/ROH/cod/reference/GCF


cp *.fastq.gz ${SLURM_TMPDIR}
cp ${REFERENCE}* ${SLURM_TMPDIR}
cp ~/links/projects/def-pbentzen/awf/reference/adapters.fa ${SLURM_TMPDIR}

cd ${SLURM_TMPDIR}

module load bwa/0.7.18 samtools/1.20 trimmomatic/0.39 picard

export JAVA_TOOL_OPTIONS="-Xmx9g"

java -Xmx9G -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 14 -phred33 \
    ${PREFIX}_1.fastq.gz \
    ${PREFIX}_2.fastq.gz \
    ${PREFIX}_FP_low.fq.gz \
    ${PREFIX}_FUP_low.fq.gz \
    ${PREFIX}_RP_low.fq.gz \
    ${PREFIX}_RUP_low.fq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 SLIDINGWINDOW:5:15 \
    MINLEN:36


java -Xmx9G -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 14 -phred33 \
    ${PREFIX}_1.fastq.gz \
    ${PREFIX}_2.fastq.gz \
    ${PREFIX}_FP_high.fq.gz \
    ${PREFIX}_FUP_high.fq.gz \
    ${PREFIX}_RP_high.fq.gz \
    ${PREFIX}_RUP_high.fq.gz \
    ILLUMINACLIP:adapters.fa:2:30:15 LEADING:20 SLIDINGWINDOW:5:20 \
    AVGQUAL:30 MINLEN:36

REFERENCE_GENOME=GCF_902167405.1_gadMor3.0_genomic.fna


bwa mem -M -t 15 \
  ${REFERENCE_GENOME} \
  ${PREFIX}_FP_low.fq.gz \
  ${PREFIX}_RP_low.fq.gz > \
  ${PREFIX}_paired_low.sam

bwa mem -M -t 15 \
  ${REFERENCE_GENOME} \
  ${PREFIX}_FP_high.fq.gz \
  ${PREFIX}_RP_high.fq.gz > \
  ${PREFIX}_paired_high.sam

samtools sort -o ${PREFIX}_paired_high.bam ${PREFIX}_paired_high.sam

samtools sort -o ${PREFIX}_paired_low.bam ${PREFIX}_paired_low.sam

export JAVA_TOOL_OPTIONS="-Xmx9g"

java -Xmx9G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  --REMOVE_DUPLICATES true --CREATE_INDEX true \
  --INPUT ${PREFIX}_paired_low.bam \
  --OUTPUT ${PREFIX}_paired-nodups_low.bam \
  --METRICS_FILE ${PREFIX}_paired-nodups_low.metrics

java -Xmx9G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  --REMOVE_DUPLICATES true --CREATE_INDEX true \
  --INPUT ${PREFIX}_paired_high.bam \
  --OUTPUT ${PREFIX}_paired-nodups_high.bam \
  --METRICS_FILE ${PREFIX}_paired-nodups_high.metrics


cp ${PREFIX}_paired-nodups* ${SLURM_SUBMIT_DIR}

```
Should have added readgroup here in BWA. 

### May 28, 2026 ###

**Fixing RG and indexing bam**

```
SAMPLE=SRR13374377
samtools addreplacerg -r 'ID:L001\tSM:SRR13374377\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_high_header.bam ${SAMPLE}_paired-nodups_high.bam 
SAMPLE=SRR13374381
samtools addreplacerg -r 'ID:L001\tSM:SRR13374381\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_high_header.bam ${SAMPLE}_paired-nodups_high.bam 
SAMPLE=SRR13374380
samtools addreplacerg -r 'ID:L001\tSM:SRR13374380\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_high_header.bam ${SAMPLE}_paired-nodups_high.bam
SAMPLE=SRR13374384
samtools addreplacerg -r 'ID:L001\tSM:SRR13374384\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_high_header.bam ${SAMPLE}_paired-nodups_high.bam
SAMPLE=SRR13374385
samtools addreplacerg -r 'ID:L001\tSM:SRR13374385\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_high_header.bam ${SAMPLE}_paired-nodups_high.bam
SAMPLE=SRR13374386
samtools addreplacerg -r 'ID:L001\tSM:SRR13374386\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_high_header.bam ${SAMPLE}_paired-nodups_high.bam 

samtools index SRR13374374_paired-nodups_high_header.bam
samtools index SRR13374377_paired-nodups_high_header.bam
samtools index SRR13374380_paired-nodups_high_header.bam
samtools index SRR13374381_paired-nodups_high_header.bam
samtools index SRR13374384_paired-nodups_high_header.bam
samtools index SRR13374385_paired-nodups_high_header.bam
samtools index SRR13374386_paired-nodups_high_header.bam


SAMPLE=SRR13374346
samtools addreplacerg -r 'ID:L001\tSM:SRR13374346\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_low_header.bam ${SAMPLE}_paired-nodups_low.bam 
SAMPLE=SRR13374347
samtools addreplacerg -r 'ID:L001\tSM:SRR13374347\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_low_header.bam ${SAMPLE}_paired-nodups_low.bam 
SAMPLE=SRR13374348
samtools addreplacerg -r 'ID:L001\tSM:SRR13374348\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_low_header.bam ${SAMPLE}_paired-nodups_low.bam 
SAMPLE=SRR13374374
samtools addreplacerg -r 'ID:L001\tSM:SRR13374374\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_low_header.bam ${SAMPLE}_paired-nodups_low.bam 

SAMPLE=SRR13374377
samtools addreplacerg -r 'ID:L001\tSM:SRR13374377\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_low_header.bam ${SAMPLE}_paired-nodups_low.bam 
SAMPLE=SRR13374381
samtools addreplacerg -r 'ID:L001\tSM:SRR13374381\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_low_header.bam ${SAMPLE}_paired-nodups_low.bam 
SAMPLE=SRR13374380
samtools addreplacerg -r 'ID:L001\tSM:SRR13374380\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_low_header.bam ${SAMPLE}_paired-nodups_low.bam
SAMPLE=SRR13374384
samtools addreplacerg -r 'ID:L001\tSM:SRR13374384\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_low_header.bam ${SAMPLE}_paired-nodups_low.bam
SAMPLE=SRR13374385
samtools addreplacerg -r 'ID:L001\tSM:SRR13374385\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_low_header.bam ${SAMPLE}_paired-nodups_low.bam
SAMPLE=SRR13374386
samtools addreplacerg -r 'ID:L001\tSM:SRR13374386\tLB:cod\tCN:NCBI\tPL:ILLUMINA' -o ${SAMPLE}_paired-nodups_low_header.bam ${SAMPLE}_paired-nodups_low.bam 

samtools index SRR13374346_paired-nodups_low_header.bam
samtools index SRR13374347_paired-nodups_low_header.bam
samtools index SRR13374348_paired-nodups_low_header.bam
samtools index SRR13374374_paired-nodups_low_header.bam
samtools index SRR13374377_paired-nodups_low_header.bam
samtools index SRR13374380_paired-nodups_low_header.bam
samtools index SRR13374381_paired-nodups_low_header.bam
samtools index SRR13374384_paired-nodups_low_header.bam
samtools index SRR13374385_paired-nodups_low_header.bam
samtools index SRR13374386_paired-nodups_low_header.bam
```

**Variant Calling**

```
#!/bin/bash
#SBATCH --job-name=high-haplotype-gvcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-pbentzen
#SBATCH --array=2-10
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --time=8:00:00

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sra_list)

export JAVA_TOOL_OPTIONS="-Xmx24G"

gatk --java-options "-Xmx24g" HaplotypeCaller \
  -R reference/GCF_902167405.1_gadMor3.0_genomic.fna \
  -I ${SAMPLE}_paired-nodups_high_header.bam \
  -O ${SAMPLE}_high_whole.g.vcf.gz \
  -ERC GVCF


#!/bin/bash
#SBATCH --job-name=low-haplotype-gvcf
#SBATCH --output=%x-%j.out
#SBATCH --account=def-pbentzen
#SBATCH --array=1-10
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --time=8:00:00

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sra_list)

export JAVA_TOOL_OPTIONS="-Xmx24G"

gatk --java-options "-Xmx24g" HaplotypeCaller \
  -R reference/GCF_902167405.1_gadMor3.0_genomic.fna \
  -I ${SAMPLE}_paired-nodups_low_header.bam \
  -O ${SAMPLE}_low_whole.g.vcf.gz \
  -ERC GVCF
```

Here I will have gvcf files for each individual. This is used for PLINK I believe. Some of the others will use them combined with genomicsDB and genotypeGVCFs

```
#!/bin/bash
#SBATCH --job-name=genomicsdb_low
#SBATCH --output=%x-%j.out
#SBATCH --account=def-pbentzen
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --time=4:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" chr_list)

export JAVA_TOOL_OPTIONS="-Xmx24G"

gatk --java-options "-Xmx24g" GenomicsDBImport \
        --genomicsdb-update-workspace-path genomicsDB_low_${CHR} \
        -V ~/links/scratch/ROH/cod/gvcf/SRR13374346_low_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374347_low_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374348_low_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374374_low_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374377_low_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374380_low_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374381_low_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374384_low_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374385_low_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374386_low_whole.g.vcf.gz \
    -L ${CHR} \
    --tmp-dir ${SLURM_TMPDIR}


#!/bin/bash
#SBATCH --job-name=genomicsdb_high
#SBATCH --output=%x-%j.out
#SBATCH --account=def-pbentzen
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --time=4:00:00

export JAVA_TOOL_OPTIONS="-Xmx24G"

gatk --java-options "-Xmx24g" GenomicsDBImport \
        --genomicsdb-update-workspace-path ~/links/scratch/ROH/cod/gvcf/genomicsDB_high/ \
        -V ~/links/scratch/ROH/cod/gvcf/SRR13374346_high_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374347_high_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374348_high_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374374_high_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374377_high_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374380_high_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374381_high_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374384_high_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374385_high_whole.g.vcf.gz \
		-V ~/links/scratch/ROH/cod/gvcf/SRR13374386_high_whole.g.vcf.gz \
        --tmp-dir ${SLURM_TMPDIR}
```

### June 1, 2026 ###

Well I forgot that the servers were down today and tomorrow for upgrading. When they come back on I will first merge the high and low vcf files. These include all sites. 

### June 3, 2026 ###
Running to merge and filter the vcf files.
```
module load bcftool/1.22 vcftools/0.1.16


bcftools concat cod_NC_044048.1_high.vcf.gz  cod_NC_044060.1_high.vcf.gz cod_NC_044049.1_high.vcf.gz  cod_NC_044061.1_high.vcf.gz cod_NC_044050.1_high.vcf.gz  cod_NC_044062.1_high.vcf.gz cod_NC_044051.1_high.vcf.gz  cod_NC_044063.1_high.vcf.gz cod_NC_044052.1_high.vcf.gz  cod_NC_044064.1_high.vcf.gz cod_NC_044053.1_high.vcf.gz  cod_NC_044065.1_high.vcf.gz cod_NC_044054.1_high.vcf.gz  cod_NC_044066.1_high.vcf.gz cod_NC_044055.1_high.vcf.gz  cod_NC_044067.1_high.vcf.gz cod_NC_044056.1_high.vcf.gz  cod_NC_044068.1_high.vcf.gz cod_NC_044057.1_high.vcf.gz  cod_NC_044069.1_high.vcf.gz cod_NC_044058.1_high.vcf.gz  cod_NC_044070.1_high.vcf.gz cod_NC_044059.1_high.vcf.gz -Oz -o merged_cod_raw_allsites_high.vcf.gz

bcftools concat cod_NC_044048.1.vcf.gz  cod_NC_044060.1.vcf.gz cod_NC_044049.1.vcf.gz  cod_NC_044061.1.vcf.gz cod_NC_044050.1.vcf.gz  cod_NC_044062.1.vcf.gz cod_NC_044051.1.vcf.gz  cod_NC_044063.1.vcf.gz cod_NC_044052.1.vcf.gz  cod_NC_044064.1.vcf.gz cod_NC_044053.1.vcf.gz  cod_NC_044065.1.vcf.gz cod_NC_044054.1.vcf.gz  cod_NC_044066.1.vcf.gz cod_NC_044055.1.vcf.gz  cod_NC_044067.1.vcf.gz cod_NC_044056.1.vcf.gz  cod_NC_044068.1.vcf.gz cod_NC_044057.1.vcf.gz  cod_NC_044069.1.vcf.gz cod_NC_044058.1.vcf.gz  cod_NC_044070.1.vcf.gz cod_NC_044059.1.vcf.gz -Oz -o ../vcf_files/merged_cod_raw_allsites_low.vcf.gz

mv merged_cod_raw_allsites_high.vcf.gz ../vcf_files/merged_cod_raw_allsites_high.vcf.gz
```

```
module load vcftools/0.1.16

##Filtering quality and missingness
vcftools --gzvcf merged_cod_raw_allsites_high.vcf.gz --recode --recode-INFO-all --out merged_cod_dp_gq_allsites_high --minDP 10 --minGQ 30

vcftools --gzvcf merged_cod_raw_allsites_low.vcf.gz --recode --recode-INFO-all --out merged_cod_dp_gq_allsites_low --minDP 5 --minGQ 20

gzip merged_cod_dp_gq_allsites_high.recode.vcf
gzip merged_cod_dp_gq_allsites_low.recode.vcf

module load bcftools/1.22

bcftools filter -i 'COUNT(FORMAT/GT="mis")<2 && INFO/MQ>30' -Oz -o merged_cod_filtered_allsites_high.recode.vcf.gz merged_cod_dp_gq_allsites_high.recode.vcf.gz

bcftools filter -i 'COUNT(FORMAT/GT="mis")<5 && INFO/MQ>20' -Oz -o merged_cod_filtered_allsites_low.recode.vcf.gz merged_cod_dp_gq_allsites_low.recode.vcf.gz


##Filtering variants only
bcftools view -m2 -Oz -o merged_cod_filtered_variants_low.recode.vcf.gz merged_cod_filtered_allsites_low.recode.vcf.gz  

bcftools view -m2 -M2 -v snps -Oz -o merged_cod_filtered_variants_high.recode.vcf.gz merged_cod_filtered_allsites_high.recode.vcf.gz 

bcftools roh -G 30 -Orz --threads 8 -o cod_filtered_allsites_low_ROH.txt.gz merged_cod_filtered_allsites_low.recode.vcf.gz

bcftools roh -G 30 -Orz --threads 8 -o cod_filtered_variants_high_ROH.txt.gz merged_cod_filtered_variants_high.recode.vcf.gz

bcftools roh -G 30 -Orz --threads 8 -o cod_filtered_allsites_high_ROH.txt.gz merged_cod_filtered_allsites_high.recode.vcf.gz

bcftools roh -G 30 -Orz --threads 8 -o cod_filtered_variants_low_ROH.txt.gz merged_cod_filtered_variants_low.recode.vcf.gz
```
I am going to start with just a simple bcftools -ROH -G30


### June 4, 2026 ###

These worked great. 

```
zcat cod_filtered_variants_low_ROH.txt.gz | tail -n+4 \
| sed 's/# //g' | sed -E 's/\[[0-9]\]//g' | sed 's/ (bp)//g' \
| sed 's/ (average fwd-bwd phred score)//g' \
| tr ' ' '_'> LOW_cod_filtered_variants_ROH.txt

DATA=LOW_cod_filtered_variants_ROH.txt
while read -r SAMPLE ; do 
cat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $6>=1e5){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/653355439, num, sum/num}'
done < sra_list

###

zcat cod_filtered_variants_high_ROH.txt.gz | tail -n+4 \
| sed 's/# //g' | sed -E 's/\[[0-9]\]//g' | sed 's/ (bp)//g' \
| sed 's/ (average fwd-bwd phred score)//g' \
| tr ' ' '_'> HIGH_cod_filtered_variants_ROH.txt

DATA=HIGH_cod_filtered_variants_ROH.txt
while read -r SAMPLE ; do 
cat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $6>=1e5){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/653355439, num, sum/num}'
done < sra_list

###

zcat cod_filtered_allsites_high_ROH.txt.gz | tail -n+4 \
| sed 's/# //g' | sed -E 's/\[[0-9]\]//g' | sed 's/ (bp)//g' \
| sed 's/ (average fwd-bwd phred score)//g' \
| tr ' ' '_'> HIGH_cod_filtered_allsites_ROH.txt

DATA=HIGH_cod_filtered_allsites_ROH.txt
while read -r SAMPLE ; do 
cat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $6>=1e5){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/653355439, num, sum/num}'
done < sra_list > high_allsites_cod_ROH_ind

###

zcat cod_filtered_allsites_low_ROH.txt.gz | tail -n+4 \
| sed 's/# //g' | sed -E 's/\[[0-9]\]//g' | sed 's/ (bp)//g' \
| sed 's/ (average fwd-bwd phred score)//g' \
| tr ' ' '_'> LOW_cod_filtered_allsites_ROH.txt

DATA=LOW_cod_filtered_allsites_ROH.txt
while read -r SAMPLE ; do 
cat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $6>=1e5){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/653355439, num, sum/num}'
done < sra_list > low_allsites_cod_ROH_ind

```

Now I am going to try and make plink file from my four datasets.
```{r}
plink --vcf merged_cod_filtered_allsites_low.recode.vcf.gz --recode --allow-extra-chr --out cod_filtered_allsites_low
plink --vcf merged_cod_filtered_variants_low.recode.vcf.gz --recode --allow-extra-chr --out cod_filtered_variants_low
plink --vcf merged_cod_filtered_allsites_high.recode.vcf.gz --recode --allow-extra-chr --out cod_filtered_allsites_high
plink --vcf merged_cod_filtered_variants_high.recode.vcf.gz --recode --allow-extra-chr --out cod_filtered_variants_high

```
Now I had .ped & .map files for the four datasets.

Ideas from other papers:
```
#generate basic plink map file
module load plink/1.9b_6.21-x86_64

plink --vcf ~/scratch/IBDNe/IBDSeq/narw_all_filters_no_missing_aug22.vcf.gz --recode --allow-extra-chr --out narw_ibdne

plink --vcf ~/scratch/IBDNe/IBDSeq/srw_all_filters_no_missing_aug22.vcf.gz --recode --allow-extra-chr --out srw_ibdne

sed -i 's,\(HiC_scaffold_.*\)\t\.\(\t0\t\)\(.*\),\1\t\1_\3\2\3,' narw_ibdne.map
```
Past thoughts:
"I ran scikit-allel with phet_roh=0 which would mean that a run of homozygosity should not contain any heterozygous sites. The FROH values were very low and there were very few ROHs over 1000000 bp in length.

I don't necessarily agree with the assumptions that are made. I can see making assumptions about potential variant sites inbetween heteozygous base calls, but with robust variant calling, I don't think that the assumption or a probabilty that a heterozygote is found in a ROH is a good one. That said, if that is the assumption all tests are making, it doesn't make a lot of sense for us to not do that, because you would be comparing apples to oranges at that point. 
I am going to implement the PLINK, VCFTOOLS and BCFTOOLS per Robinson et al 2022 and Foote et al 2021."

Foote:
window- based approach implemented in plink v1.07
Sliding window size was set to 300 kb, with a minimum of 50 SNPs at a minimum density of 1 SNP per 50 kb required to call a ROH. To we allowed up to 3 heterozygote sites per 300 kb window within called ROHs, as per Ceballos et al. (2018).

```
plink --tfile chr1 --homozyg-snp 50 --homozyg-kb 300 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 3 --homozyg-window-missing 10 --homozyg-window-threshold 0.05 --out chr1
```

Robinson: 
We identified runs of homozygosity (ROH) ≥ 1 Mb with three methods: the --homozyg function in Plink v1.90b6.12 (40), the --LROH function in Vcftools v0.1.16, (41), and the roh function in Bcftools v1.10.2

Plink and Vcftools were run with default parameters. VCF files containing all genotyped sites were provided for Plink, whereas files containing just variant sites were provided for Vcftools. For Bcftools, we provided variants only, used option -G 30, and specified that the function should use allele frequencies from the 15 unrelated vaquitas in identifying ROH
```
#BCFTOOLS:
KEEP=samples_keep.list
bcftools roh -e ${KEEP} -G 30 -Orz -o ${VCF}_bcftoolsROH.txt.gz ${VCF}

zcat ${VCF}_bcftoolsROH.txt.gz | tail -n+4 \
| sed 's/# //g' | sed -E 's/\[[0-9]\]//g' | sed 's/ (bp)//g' \
| sed 's/ (average fwd-bwd phred score)//g' \
| tr ' ' '_'> ${VCF}_bcftoolsROH.txt

DATA=vaquita_20_simple_PASS_autos_variants.vcf.gz_bcftoolsROH.txt.gz
while read -r SAMPLE ; do 
zcat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $6>=1e6){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/2237397028, num, sum/num}'
done < samples.list

#vcftools

VCF=vaquita_20_simple_PASS_autos_variants.vcf.gz
for CHR in chr{1..21} ; do
vcftools --gzvcf ${VCF} --LROH --chr ${CHR} --out ${VCF}_vcftools_${CHR}
done

OUT=${VCF}_vcftools.LROH
cp ${VCF}_vcftools_chr1.LROH  ${OUT}
cat vaquita_20_simple_PASS_autos_variants.vcf.gz_vcftools_chr{2..21}.LROH | grep -v "^CHROM" >> ${OUT}

#plink - including invariamt sites

# Convert simplified, concatenated VCF file to plink format
export PLINK=~/project/programs/plink1.9/plink
export VCF=vaquita_20_simple_PASS_autos_variants.vcf.gz
${PLINK} --make-bed --keep-allele-order --double-id --chr-set 21 --vcf ${VCF} --out ${VCF%.vcf.gz}

# Identify ROH and reformat output
FILE=vaquita_20_simple_PASS_autos_variants
${PLINK} --bfile ${VCF%.vcf.gz} --out ${VCF%.vcf.gz}_ROH --homozyg
sed -i -E 's/^\ +//g' ${FILE}_ROH.hom
sed -i -E 's/\ +/\t/g' ${FILE}_ROH.hom 
```

```
# ATTEMPTED CODE FOR Bcftools

VCF=~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz
KEEP=../narw
bcftools roh -e ${KEEP} -G 30 -Orz --threads 8 -o narw_bcftoolsROH.txt.gz ${VCF}

zcat narw_bcftoolsROH.txt.gz | tail -n+4 \
| sed 's/# //g' | sed -E 's/\[[0-9]\]//g' | sed 's/ (bp)//g' \
| sed 's/ (average fwd-bwd phred score)//g' \
| tr ' ' '_'> narw_bcftoolsROH.txt

DATA=narw_bcftoolsROH.txt
while read -r SAMPLE ; do 
cat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $6>=1e6){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/2170759585, num, sum/num}'
done < ../narw

VCF=~/projects/def-frasiert/RW_WGS/vcf/locked/srw_unrelated_filtered_aug5.vcf.gz
KEEP=../srw
bcftools roh -e ${KEEP} -G 30 -Orz --threads 8 -o srw_bcftoolsROH.txt.gz ${VCF}

zcat srw_bcftoolsROH.txt.gz | tail -n+4 \
| sed 's/# //g' | sed -E 's/\[[0-9]\]//g' | sed 's/ (bp)//g' \
| sed 's/ (average fwd-bwd phred score)//g' \
| tr ' ' '_'> srw_bcftoolsROH.txt

DATA=srw_bcftoolsROH.txt
while read -r SAMPLE ; do 
cat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $6>=1e6){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/2296311778, num, sum/num}'
done < ../srw

# I don't have gvcf file for PLINK

# VCFTOOLS
VCF=~/projects/def-frasiert/RW_WGS/vcf/locked/narw_unrelated_filtered_aug5.vcf.gz
for CHR in HiC_Scaffold_{1..22} ; do
vcftools --gzvcf ${VCF} --LROH --chr ${CHR} --out ${VCF}_vcftools_${CHR}
done

OUT=narw_vcftools.LROH
cp narw_vcftools_roh_HiC_scaffold_1.LROH  ${OUT}
cat narw_vcftools_roh_HiC_scaffold_{2..22}.LROH | grep -v "^CHROM" >> ${OUT}

VCF=~/projects/def-frasiert/RW_WGS/vcf/locked/srw_unrelated_filtered_aug5.vcf.gz
for CHR in HiC_scaffold_{1..21} ; do
vcftools --gzvcf ${VCF} --LROH --chr ${CHR} --out srw_vcftools_roh_${CHR}
done

OUT=srw_vcftools.LROH
cp srw_vcftools_roh_HiC_scaffold_1.LROH  ${OUT}
cat srw_vcftools_roh_HiC_scaffold_{2..21}.LROH | grep -v "^CHROM" >> ${OUT}

```


### August 25, 2026 ###

Downloading Eastern Massasauga Rattlesnake genomes
```
module load sra-toolkit

wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15274910/SRR15274910
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15274909/SRR15274909
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15274908/SRR15274908
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15274907/SRR15274907
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15274906/SRR15274906
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15274905/SRR15274905
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15274904/SRR15274904
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15274903/SRR15274903
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15274902/SRR15274902
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15274901/SRR15274901


fastq-dump --split-files --gzip -v SRR15274910
fastq-dump --split-files --gzip -v SRR15274909
fastq-dump --split-files --gzip -v SRR15274908
fastq-dump --split-files --gzip -v SRR15274907
fastq-dump --split-files --gzip -v SRR15274906
fastq-dump --split-files --gzip -v SRR15274905
fastq-dump --split-files --gzip -v SRR15274904
fastq-dump --split-files --gzip -v SRR15274903
fastq-dump --split-files --gzip -v SRR15274902
fastq-dump --split-files --gzip -v SRR15274901
```

### September 5, 2026 ###

Downloading the NARW files:
```
module load sra-toolkit

wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR29386324/SRR29386324
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR29386336/SRR29386336
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR29386335/SRR29386335
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR29386321/SRR29386321
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR29386320/SRR29386320
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR29386319/SRR29386319
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR29386318/SRR29386318
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR29386317/SRR29386317
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR29386316/SRR29386316
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR29386315/SRR29386315

#!/bin/bash
#SBATCH --job-name=unpack_fastq
#SBATCH --output=%x-%j.out
#SBATCH --account=def-pbentzen
#SBATCH --array=1-10
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G
#SBATCH --time=3:00:00

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples)

fastq-dump --split-files --gzip -v ${SAMPLE}


#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-pbentzen
#SBATCH --job-name=index_reference
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=8

module load bwa/0.7.18 samtools/1.22.1

bwa index GCA_039880765.1_rSisCat1_p1.0_genomic.fna

samtools faidx GCA_039880765.1_rSisCat1_p1.0_genomic.fna

samtools dict GCA_039880765.1_rSisCat1_p1.0_genomic.fna -o GCA_039880765.1_rSisCat1_p1.0_genomic.dict
```

```
#!/bin/bash
#SBATCH --account=def-pbentzen
#SBATCH --job-name=trim_map
#SBATCH --array=1-30
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=15
#SBATCH --time=12:00:00

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

bwa mem -M -t 15 \
 -R '@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SPECIES}\tCN:NCBI\tPL:ILLUMINA' \
  ${REFERENCE} \
  ${SPECIES}_${SAMPLE}_FP_high.fq.gz \
  ${SPECIES}_${SAMPLE}_RP_high.fq.gz > \
  ${SPECIES}_${SAMPLE}_paired_high.sam

samtools sort -o ${SPECIES}_${SAMPLE}_paired_high.bam ${SPECIES}_${SAMPLE}_paired_high.sam

samtools sort -o ${SPECIES}_${SAMPLE}_paired_low.bam ${SPECIES}_${SAMPLE}_paired_low.sam

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
  --INPUT ${SPECIES}_${SAMPLE}_paired_high.bam \
  --OUTPUT ${SPECIES}_${SAMPLE}_paired-nodups_high.bam \
  --METRICS_FILE ${SPECIES}_${SAMPLE}_paired-nodups_high.metrics


###########################################
# Copy finished files to relevant directory 
###########################################

cp ${SPECIES}_${SAMPLE}_paired-nodups* ${SLURM_SUBMIT_DIR}/${SPECIES}/bam_files/
```

Testing this on array job 30. Looks good ran on all.

### September 8, 2026 ###

This ran great on the cod and rattlesnake genomes. Kept having issues with some of the right whale based on the size of the data. Changed the code about for handling larger files and deleting some intermediates along the way to free up space. 

Array job 6 and 10 ran in full.
Array 2 and 9 are running with fewer resources but still going with a few tweaks. They will need to be checked.

Array jobs 1,3-5,7-8 are running this:

```
#!/bin/bash
#SBATCH --account=def-pbentzen
#SBATCH --job-name=trim_map
#SBATCH --array=1,3,4,5,7,8
#SBATCH --output=%x-%a-%j.out
#SBATCH --mem-per-cpu=45G
#SBATCH --cpus-per-task=30
#SBATCH --time=12:00:00

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

rm *RUP*.fq.gz
rm *FUP*.fq.gz

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

rm *RUP*.fq.gz
rm *FUP*.fq.gz


########################################
# Map trimmed reads to reference genomes
########################################

rm ${SAMPLE}*.fastq.gz

bwa mem -M -t 30 \
 -R '@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SPECIES}\tCN:NCBI\tPL:ILLUMINA' \
  ${REFERENCE} \
  ${SPECIES}_${SAMPLE}_FP_low.fq.gz \
  ${SPECIES}_${SAMPLE}_RP_low.fq.gz > \
  ${SPECIES}_${SAMPLE}_paired_low.sam

samtools sort -o ${SPECIES}_${SAMPLE}_paired_low.bam ${SPECIES}_${SAMPLE}_paired_low.sam
rm ${SPECIES}_${SAMPLE}_paired_low.sam

bwa mem -M -t 30 \
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
  --INPUT ${SPECIES}_${SAMPLE}_paired_high.bam \
  --OUTPUT ${SPECIES}_${SAMPLE}_paired-nodups_high.bam \
  --METRICS_FILE ${SPECIES}_${SAMPLE}_paired-nodups_high.metrics


###########################################
# Copy finished files to relevant directory 
###########################################

cp ${SPECIES}_${SAMPLE}_paired-nodups* ${SLURM_SUBMIT_DIR}/${SPECIES}/bam_files/

```

### Septemeber 11, 2026 ###

I have been working through the code above in chunks for the right whale files because they are so big.

Looking at coverage.

```
ls *.bam > rattlesnake_bams
sed -i 's,\(rattlesnake_\)\(.*\)\(_paired-nodups\)\(_[a-z]\+\)\(.bam\),\2\4\t\1\2\3\4\5,' rattlesnake_bams 
```

```
samtools coverage ${BAM_FILE}
```

**Depth Problem**

Looks like the rattlesnake data were top ups to existing datasets to increase sequencing depth. I am going to redo with 10 new samples, that were run once with sifficient data:

SRR36656545
SRR36656544
SRR36656543
SRR36656542
SRR36656541
SRR36656540
SRR36656539
SRR36656538
SRR36656537
SRR36656536

```
module load sra-toolkit

wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR36656545/SRR36656545
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR36656544/SRR36656544
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR36656543/SRR36656543
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR36656542/SRR36656542
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR36656541/SRR36656541
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR36656540/SRR36656540
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR36656539/SRR36656539
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR36656538/SRR36656538
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR36656537/SRR36656537
wget  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR36656536/SRR36656536


fastq-dump --split-files --gzip -v SRR36656545
fastq-dump --split-files --gzip -v SRR36656544
fastq-dump --split-files --gzip -v SRR36656542
fastq-dump --split-files --gzip -v SRR36656541
fastq-dump --split-files --gzip -v SRR36656540
fastq-dump --split-files --gzip -v SRR36656539
fastq-dump --split-files --gzip -v SRR36656538
fastq-dump --split-files --gzip -v SRR36656537
fastq-dump --split-files --gzip -v SRR36656536
fastq-dump --split-files --gzip -v SRR36656543
```

I will be rerunning the trim-map pipeline when they are done unpacking.

### September 14, 2026 ###

I have marked bam files for all samples now.
```
for SAMPLE in SRR13374346 SRR13374347 SRR13374348 SRR13374374 SRR13374377 SRR13374380 SRR13374381 SRR13374384 SRR13374385 SRR13374386
do
 sed "s,\(.*\),$SAMPLE\t\1," chrom_list >> data_list
done


for SAMPLE in SRR29386315 SRR29386316 SRR29386317 SRR29386318 SRR29386319 SRR29386320 SRR29386321 SRR29386324 SRR29386335 SRR29386336
do
 sed "s,\(.*\),$SAMPLE\t\1," chrom_list >> data_list
done

for SAMPLE in SRR36656536 SRR36656537 SRR36656538 SRR36656539 SRR36656540 SRR36656541 SRR36656542 SRR36656543 SRR36656544 SRR36656545
do
 sed "s,\(.*\),$SAMPLE\t\1," chrom_list >> data_list
done
```


**Variant Calling**

```
#!/bin/bash
#SBATCH --job-name=haplotype-gvcf
#SBATCH --account=def-pbentzen
#SBATCH --array=1-10
#SBATCH --output=%x-%a-%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=25G
#SBATCH --time=2:00:00

LIST=$(sed -n "${SLURM_ARRAY_TASK_ID}p" data_list)
SAMPLE=$(echo $LIST | awk '{print $1}')
SCAFFOLD=$(echo $LIST | awk '{print $2}')

export JAVA_TOOL_OPTIONS="-Xmx24G"

gatk --java-options "-Xmx24g" HaplotypeCaller \
  -R ../reference/GCF_902167405.1_gadMor3.0_genomic.fna \
  -I ${SAMPLE}_paired-nodups_high.bam \
  -O ../vcf/${SAMPLE}_${SCAFFOLD}_high_whole.g.vcf.gz \
  -ERC GVCF \
  -L ${SCAFFOLD}

gatk --java-options "-Xmx24g" HaplotypeCaller \
  -R ../reference/GCF_902167405.1_gadMor3.0_genomic.fna \
  -I ${SAMPLE}_paired-nodups_low.bam \
  -O ../vcf/${SAMPLE}_${SCAFFOLD}_low_whole.g.vcf.gz \
  -ERC GVCF \
  -L ${SCAFFOLD}
```

### Septemeber 15, 2026 ###

The sample name didn't population properly in the bam file. I think I can correct it with a sample map.

Creating sample map files for GenomicsDB Imports. First for cod:

SRR13374346  SRR13374346_CHR_low_whole.g.vcf.gz
SRR13374347  SRR13374347_CHR_low_whole.g.vcf.gz
SRR13374348  SRR13374348_CHR_low_whole.g.vcf.gz
SRR13374374  SRR13374374_CHR_low_whole.g.vcf.gz
SRR13374377  SRR13374377_CHR_low_whole.g.vcf.gz
SRR13374380  SRR13374380_CHR_low_whole.g.vcf.gz
SRR13374381  SRR13374381_CHR_low_whole.g.vcf.gz
SRR13374384  SRR13374384_CHR_low_whole.g.vcf.gz
SRR13374385  SRR13374385_CHR_low_whole.g.vcf.gz
SRR13374386  SRR13374386_CHR_low_whole.g.vcf.gz


SRR36656536 SRR36656536_CHR_low_whole.g.vcf.gz
SRR36656541 SRR36656541_CHR_low_whole.g.vcf.gz
SRR36656537 SRR36656537_CHR_low_whole.g.vcf.gz
SRR36656542 SRR36656542_CHR_low_whole.g.vcf.gz
SRR36656538 SRR36656538_CHR_low_whole.g.vcf.gz
SRR36656543 SRR36656543_CHR_low_whole.g.vcf.gz
SRR36656539 SRR36656539_CHR_low_whole.g.vcf.gz
SRR36656544 SRR36656544_CHR_low_whole.g.vcf.gz
SRR36656540 SRR36656540_CHR_low_whole.g.vcf.gz
SRR36656545 SRR36656545_CHR_low_whole.g.vcf.gz

SRR29386315 SRR29386315_CHR_low_whole.g.vcf.gz
SRR29386316 SRR29386316_CHR_low_whole.g.vcf.gz
SRR29386317 SRR29386317_CHR_low_whole.g.vcf.gz
SRR29386318 SRR29386318_CHR_low_whole.g.vcf.gz
SRR29386319 SRR29386319_CHR_low_whole.g.vcf.gz
SRR29386320 SRR29386320_CHR_low_whole.g.vcf.gz
SRR29386321 SRR29386321_CHR_low_whole.g.vcf.gz
SRR29386324 SRR29386324_CHR_low_whole.g.vcf.gz
SRR29386335 SRR29386335_CHR_low_whole.g.vcf.gz
SRR29386336 SRR29386336_CHR_low_whole.g.vcf.gz


```
for CHR in NC_083716.1 NC_083717.1 NC_083718.1 NC_083719.1 NC_083720.1 NC_083721.1 NC_083722.1 NC_083723.1 NC_083724.1 NC_083725.1 NC_083726.1 NC_083727.1 NC_083728.1 NC_083729.1 NC_083730.1 NC_083731.1 NC_083732.1 NC_083733.1 NC_083734.1 NC_083735.1 
do
  sed "s,CHR,${CHR},g" list > ${CHR}_low_sample_map
  sed "s,low,high,g" ${CHR}_low_sample_map > ${CHR}_high_sample_map
done

for CHR in NC_083716.1 NC_083717.1 NC_083718.1 NC_083719.1 NC_083720.1 NC_083721.1 NC_083722.1 NC_083723.1 NC_083724.1 NC_083725.1 NC_083726.1 NC_083727.1 NC_083728.1 NC_083729.1 NC_083730.1 NC_083731.1 NC_083732.1 NC_083733.1 NC_083734.1 NC_083735.1 
do
  sed -i "s,\s,\t,g" ${CHR}_low_sample_map
  sed -i "s,\s,\t,g" ${CHR}_high_sample_map
done
```

**GenomicsDB Import**

```
#!/bin/bash
#SBATCH --job-name=genomicsdb
#SBATCH --account=def-pbentzen
#SBATCH --array=1-23
#SBATCH --output=%x-%a-%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --time=8:00:00

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" chrom_list)

export JAVA_TOOL_OPTIONS="-Xmx24G"

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
```

**Joint genotype calling**

```
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
```


### September 17, 2026 ###

Merging VCF files

```
module load bcftool/1.22


bcftools concat cod_NC_044048.1_high.vcf.gz cod_NC_044049.1_high.vcf.gz cod_NC_044050.1_high.vcf.gz cod_NC_044051.1_high.vcf.gz cod_NC_044052.1_high.vcf.gz cod_NC_044053.1_high.vcf.gz cod_NC_044054.1_high.vcf.gz cod_NC_044055.1_high.vcf.gz cod_NC_044056.1_high.vcf.gz cod_NC_044057.1_high.vcf.gz cod_NC_044058.1_high.vcf.gz cod_NC_044059.1_high.vcf.gz cod_NC_044060.1_high.vcf.gz cod_NC_044061.1_high.vcf.gz cod_NC_044062.1_high.vcf.gz cod_NC_044063.1_high.vcf.gz cod_NC_044064.1_high.vcf.gz cod_NC_044065.1_high.vcf.gz cod_NC_044066.1_high.vcf.gz cod_NC_044067.1_high.vcf.gz cod_NC_044068.1_high.vcf.gz cod_NC_044069.1_high.vcf.gz cod_NC_044070.1_high.vcf.gz -Oz -o merged_cod_raw_allsites_high.vcf.gz

bcftools concat cod_NC_044048.1_low.vcf.gz cod_NC_044049.1_low.vcf.gz cod_NC_044050.1_low.vcf.gz cod_NC_044051.1_low.vcf.gz cod_NC_044052.1_low.vcf.gz cod_NC_044053.1_low.vcf.gz cod_NC_044054.1_low.vcf.gz cod_NC_044055.1_low.vcf.gz cod_NC_044056.1_low.vcf.gz cod_NC_044057.1_low.vcf.gz cod_NC_044058.1_low.vcf.gz cod_NC_044059.1_low.vcf.gz cod_NC_044060.1_low.vcf.gz cod_NC_044061.1_low.vcf.gz cod_NC_044062.1_low.vcf.gz cod_NC_044063.1_low.vcf.gz cod_NC_044064.1_low.vcf.gz cod_NC_044065.1_low.vcf.gz cod_NC_044066.1_low.vcf.gz cod_NC_044067.1_low.vcf.gz cod_NC_044068.1_low.vcf.gz cod_NC_044069.1_low.vcf.gz cod_NC_044070.1_low.vcf.gz -Oz -o merged_cod_raw_allsites_low.vcf.gz

bcftools concat rattlesnake_CM078115.1_high.vcf.gz rattlesnake_CM078116.1_high.vcf.gz rattlesnake_CM078117.1_high.vcf.gz rattlesnake_CM078118.1_high.vcf.gz rattlesnake_CM078119.1_high.vcf.gz rattlesnake_CM078120.1_high.vcf.gz rattlesnake_CM078121.1_high.vcf.gz rattlesnake_CM078122.1_high.vcf.gz rattlesnake_CM078123.1_high.vcf.gz rattlesnake_CM078124.1_high.vcf.gz rattlesnake_CM078125.1_high.vcf.gz rattlesnake_CM078126.1_high.vcf.gz rattlesnake_CM078127.1_high.vcf.gz rattlesnake_CM078128.1_high.vcf.gz rattlesnake_CM078129.1_high.vcf.gz rattlesnake_CM078130.1_high.vcf.gz rattlesnake_CM078131.1_high.vcf.gz -Oz -o merged_rattlesnake_raw_allsites_high.vcf.gz


bcftools concat rattlesnake_CM078115.1_low.vcf.gz rattlesnake_CM078116.1_low.vcf.gz rattlesnake_CM078117.1_low.vcf.gz rattlesnake_CM078118.1_low.vcf.gz rattlesnake_CM078119.1_low.vcf.gz rattlesnake_CM078120.1_low.vcf.gz rattlesnake_CM078121.1_low.vcf.gz rattlesnake_CM078122.1_low.vcf.gz rattlesnake_CM078123.1_low.vcf.gz rattlesnake_CM078124.1_low.vcf.gz rattlesnake_CM078125.1_low.vcf.gz rattlesnake_CM078126.1_low.vcf.gz rattlesnake_CM078127.1_low.vcf.gz rattlesnake_CM078128.1_low.vcf.gz rattlesnake_CM078129.1_low.vcf.gz rattlesnake_CM078130.1_low.vcf.gz rattlesnake_CM078131.1_low.vcf.gz -Oz -o merged_rattlesnake_raw_allsites_low.vcf.gz


bcftools concat rightwhale_NC_083716.1_low.vcf.gz rightwhale_NC_083717.1_low.vcf.gz rightwhale_NC_083718.1_low.vcf.gz rightwhale_NC_083719.1_low.vcf.gz rightwhale_NC_083720.1_low.vcf.gz rightwhale_NC_083721.1_low.vcf.gz rightwhale_NC_083722.1_low.vcf.gz rightwhale_NC_083723.1_low.vcf.gz rightwhale_NC_083724.1_low.vcf.gz rightwhale_NC_083725.1_low.vcf.gz rightwhale_NC_083726.1_low.vcf.gz rightwhale_NC_083727.1_low.vcf.gz rightwhale_NC_083728.1_low.vcf.gz rightwhale_NC_083729.1_low.vcf.gz rightwhale_NC_083730.1_low.vcf.gz rightwhale_NC_083731.1_low.vcf.gz rightwhale_NC_083732.1_low.vcf.gz rightwhale_NC_083733.1_low.vcf.gz rightwhale_NC_083734.1_low.vcf.gz rightwhale_NC_083735.1_low.vcf.gz  -Oz -o merged_rightwhale_raw_allsites_low.vcf.gz


bcftools concat rightwhale_NC_083716.1_high.vcf.gz rightwhale_NC_083717.1_high.vcf.gz rightwhale_NC_083718.1_high.vcf.gz rightwhale_NC_083719.1_high.vcf.gz rightwhale_NC_083720.1_high.vcf.gz rightwhale_NC_083721.1_high.vcf.gz rightwhale_NC_083722.1_high.vcf.gz rightwhale_NC_083723.1_high.vcf.gz rightwhale_NC_083724.1_high.vcf.gz rightwhale_NC_083725.1_high.vcf.gz rightwhale_NC_083726.1_high.vcf.gz rightwhale_NC_083727.1_high.vcf.gz rightwhale_NC_083728.1_high.vcf.gz rightwhale_NC_083729.1_high.vcf.gz rightwhale_NC_083730.1_high.vcf.gz rightwhale_NC_083731.1_high.vcf.gz rightwhale_NC_083732.1_high.vcf.gz rightwhale_NC_083733.1_high.vcf.gz rightwhale_NC_083734.1_high.vcf.gz rightwhale_NC_083735.1_high.vcf.gz -Oz -o merged_rightwhale_raw_allsites_high.vcf.gz
```

**Filtering**
```
module load vcftools/0.1.16 

for SPECIES in cod rattlesnake rightwhale
do

vcftools --gzvcf merged_${SPECIES}_raw_allsites_high.vcf.gz --recode --recode-INFO-all --out merged_${SPECIES}_dp_gq_allsites_high --minDP 10 --minGQ 30

vcftools --gzvcf merged_${SPECIES}_raw_allsites_low.vcf.gz --recode --recode-INFO-all --out merged_${SPECIES}_dp_gq_allsites_low --minDP 5 --minGQ 20

gzip merged_${SPECIES}_dp_gq_allsites_high.recode.vcf
gzip merged_${SPECIES}_dp_gq_allsites_low.recode.vcf

done


module load bcftools/1.22

for SPECIES in cod rattlesnake rightwhale
do

bcftools filter -i 'COUNT(FORMAT/GT="mis")<2 && INFO/MQ>30' -Oz -o merged_${SPECIES}_filtered_allsites_high.recode.vcf.gz merged_${SPECIES}_dp_gq_allsites_high.recode.vcf.gz

bcftools filter -i 'COUNT(FORMAT/GT="mis")<5 && INFO/MQ>20' -Oz -o merged_${SPECIES}_filtered_allsites_low.recode.vcf.gz merged_${SPECIES}_dp_gq_allsites_low.recode.vcf.gz

bcftools view -m2 -Oz -o merged_${SPECIES}_filtered_variants_low.recode.vcf.gz merged_${SPECIES}_filtered_allsites_low.recode.vcf.gz  

bcftools view -m2 -M2 -v snps -Oz -o merged_${SPECIES}_filtered_variants_high.recode.vcf.gz merged_${SPECIES}_filtered_allsites_high.recode.vcf.gz 

bcftools +fill-tags merged_${SPECIES}_filtered_variants_high.recode.vcf.gz -Oz -o merged_${SPECIES}_filtered_variants_updatetags_high.vcf.gz -- -t all

bcftools +fill-tags merged_${SPECIES}_filtered_variants_low.recode.vcf.gz -Oz -o merged_${SPECIES}_filtered_variants_updatetags_low.vcf.gz -- -t all

done
```


### September 21, 2026 ###


Updating tags to be sure they contain the correct AF and depth etc. 

```
for SPECIES in cod rattlesnake rightwhale
do

bcftools +fill-tags merged_${SPECIES}_filtered_variants_high.recode.vcf.gz -Oz -o merged_${SPECIES}_filtered_variants_updatetags_high.vcf.gz -- -t all

bcftools +fill-tags merged_${SPECIES}_filtered_variants_low.recode.vcf.gz -Oz -o merged_${SPECIES}_filtered_variants_updatetags_low.vcf.gz -- -t all

done
```


I wanted to peak at what the data look like. Simple transformation of the bcftools on output of bcftools roh -G 30 -Or output input.gz

```
tail -n+4 rightwhale_low_bcftools_roh \
| sed 's/# //g' | sed -E 's/\[[0-9]\]//g' | sed 's/ (bp)//g' \
| sed 's/ (average fwd-bwd phred score)//g' \
| tr ' ' '_'> low_rightwhale_filtered_variants_ROH.txt

DATA=low_rightwhale_filtered_variants_ROH.txt
while read -r SAMPLE ; do 
cat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $6>=1e5){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/2512368044, num, sum/num}'
done < rightwhale_list > low_rightwhale_ROH_ind

###

tail -n+4 rightwhale_high_bcftools_roh \
| sed 's/# //g' | sed -E 's/\[[0-9]\]//g' | sed 's/ (bp)//g' \
| sed 's/ (average fwd-bwd phred score)//g' \
| tr ' ' '_'> high_rightwhale_filtered_variants_ROH.txt

DATA=high_rightwhale_filtered_variants_ROH.txt
while read -r SAMPLE ; do 
cat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $6>=1e5){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/2512368044, num, sum/num}'
done < rightwhale_list > high_rightwhale_ROH_ind
```

### Sepetmber 28, 2026 ###

**Intermdeiate thresholds**
```
#!/bin/bash
#SBATCH --account=def-pbentzen
#SBATCH --job-name=trim_map_mid
#SBATCH --array=1,3,4,5,7,8
#SBATCH --output=%x-%a-%j.out
#SBATCH --mem-per-cpu=45G
#SBATCH --cpus-per-task=30
#SBATCH --time=8:00:00

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

rm *RUP*.fq.gz
rm *FUP*.fq.gz

########################################
# Map trimmed reads to reference genomes
########################################

rm ${SAMPLE}*.fastq.gz


bwa mem -M -t 30 \
 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SPECIES}\tCN:NCBI\tPL:ILLUMINA" \
  ${REFERENCE} \
  ${SPECIES}_${SAMPLE}_FP_mid.fq.gz \
  ${SPECIES}_${SAMPLE}_RP_mid.fq.gz > \
  ${SPECIES}_${SAMPLE}_paired_mid.sam

samtools sort -o ${SPECIES}_${SAMPLE}_paired_mid.bam ${SPECIES}_${SAMPLE}_paired_mid.sam
rm ${SPECIES}_${SAMPLE}_paired_mid.sam

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
  --INPUT ${SPECIES}_${SAMPLE}_paired_mid.bam \
  --OUTPUT ${SPECIES}_${SAMPLE}_paired-nodups_mid.bam \
  --METRICS_FILE ${SPECIES}_${SAMPLE}_paired-nodups_mid.metrics


###########################################
# Copy finished files to relevant directory 
###########################################

cp ${SPECIES}_${SAMPLE}_paired-nodups* ${SLURM_SUBMIT_DIR}/${SPECIES}/bam_files/

```


### September 30, 2026 ###
Prepping files for genomics DB Import
```
# RightWhale

for CHR in NC_083716.1 NC_083717.1 NC_083718.1 NC_083719.1 NC_083720.1 NC_083721.1 NC_083722.1 NC_083723.1 NC_083724.1 NC_083725.1 NC_083726.1 NC_083727.1 NC_083728.1 NC_083729.1 NC_083730.1 NC_083731.1 NC_083732.1 NC_083733.1 NC_083734.1 NC_083735.1 
do
  sed "s,low,mid,g" ${CHR}_low_sample_map > ${CHR}_mid_sample_map
  sed -i "s,\s,\t,g" ${CHR}_mid_sample_map
done

# Cod

for CHR in NC_044048.1 NC_044049.1 NC_044050.1 NC_044051.1 NC_044052.1 NC_044053.1 NC_044054.1 NC_044055.1 NC_044056.1 NC_044057.1 NC_044058.1 NC_044059.1 NC_044060.1 NC_044061.1 NC_044062.1 NC_044063.1 NC_044064.1 NC_044065.1 NC_044066.1 NC_044067.1 NC_044068.1 NC_044069.1 NC_044070.1
do   
  sed "s,low,mid,g" ${CHR}_low_sample_map > ${CHR}_mid_sample_map
  sed -i "s,\s,\t,g" ${CHR}_mid_sample_map
done

# Rattlesnake

for CHR in CM078115.1 CM078116.1 CM078117.1 CM078118.1 CM078119.1 CM078120.1 CM078121.1 CM078122.1 CM078123.1 CM078124.1 CM078125.1 CM078126.1 CM078127.1 CM078128.1 CM078129.1 CM078130.1 CM078131.1
do   
  sed "s,low,mid,g" ${CHR}_low_sample_map > ${CHR}_mid_sample_map
  sed -i "s,\s,\t,g" ${CHR}_mid_sample_map
done
```
```
for CHR in CM078115.1 CM078116.1 CM078117.1 CM078118.1 CM078119.1 CM078120.1 CM078121.1 CM078122.1 CM078123.1 CM078124.1 CM078125.1 CM078126.1 CM078127.1 CM078128.1 CM078129.1 CM078130.1 CM078131.1
do   
  sed -i "s,mid,nid,g" ${CHR}_mid_sample_map
done
```

### October 1, 2026 ###

```
bcftools concat cod_NC_044048.1_mid.vcf.gz cod_NC_044049.1_mid.vcf.gz cod_NC_044050.1_mid.vcf.gz cod_NC_044051.1_mid.vcf.gz cod_NC_044052.1_mid.vcf.gz cod_NC_044053.1_mid.vcf.gz cod_NC_044054.1_mid.vcf.gz cod_NC_044055.1_mid.vcf.gz cod_NC_044056.1_mid.vcf.gz cod_NC_044057.1_mid.vcf.gz cod_NC_044058.1_mid.vcf.gz cod_NC_044059.1_mid.vcf.gz cod_NC_044060.1_mid.vcf.gz cod_NC_044061.1_mid.vcf.gz cod_NC_044062.1_mid.vcf.gz cod_NC_044063.1_mid.vcf.gz cod_NC_044064.1_mid.vcf.gz cod_NC_044065.1_mid.vcf.gz cod_NC_044066.1_mid.vcf.gz cod_NC_044067.1_mid.vcf.gz cod_NC_044068.1_mid.vcf.gz cod_NC_044069.1_mid.vcf.gz cod_NC_044070.1_mid.vcf.gz -Oz -o merged_cod_raw_allsites_mid.vcf.gz

bcftools concat rattlesnake_CM078115.1_mid.vcf.gz rattlesnake_CM078116.1_mid.vcf.gz rattlesnake_CM078117.1_mid.vcf.gz rattlesnake_CM078118.1_mid.vcf.gz rattlesnake_CM078119.1_mid.vcf.gz rattlesnake_CM078120.1_mid.vcf.gz rattlesnake_CM078121.1_mid.vcf.gz rattlesnake_CM078122.1_mid.vcf.gz rattlesnake_CM078123.1_mid.vcf.gz rattlesnake_CM078124.1_mid.vcf.gz rattlesnake_CM078125.1_mid.vcf.gz rattlesnake_CM078126.1_mid.vcf.gz rattlesnake_CM078127.1_mid.vcf.gz rattlesnake_CM078128.1_mid.vcf.gz rattlesnake_CM078129.1_mid.vcf.gz rattlesnake_CM078130.1_mid.vcf.gz rattlesnake_CM078131.1_mid.vcf.gz -Oz -o merged_rattlesnake_raw_allsites_mid.vcf.gz

bcftools concat rightwhale_NC_083716.1_mid.vcf.gz rightwhale_NC_083717.1_mid.vcf.gz rightwhale_NC_083718.1_mid.vcf.gz rightwhale_NC_083719.1_mid.vcf.gz rightwhale_NC_083720.1_mid.vcf.gz rightwhale_NC_083721.1_mid.vcf.gz rightwhale_NC_083722.1_mid.vcf.gz rightwhale_NC_083723.1_mid.vcf.gz rightwhale_NC_083724.1_mid.vcf.gz rightwhale_NC_083725.1_mid.vcf.gz rightwhale_NC_083726.1_mid.vcf.gz rightwhale_NC_083727.1_mid.vcf.gz rightwhale_NC_083728.1_mid.vcf.gz rightwhale_NC_083729.1_mid.vcf.gz rightwhale_NC_083730.1_mid.vcf.gz rightwhale_NC_083731.1_mid.vcf.gz rightwhale_NC_083732.1_mid.vcf.gz rightwhale_NC_083733.1_mid.vcf.gz rightwhale_NC_083734.1_mid.vcf.gz rightwhale_NC_083735.1_mid.vcf.gz  -Oz -o merged_rightwhale_raw_allsites_mid.vcf.gz
```

### October 2, 2026 ###

```
module load vcftools/0.1.16 

for SPECIES in cod rattlesnake rightwhale
do

vcftools --gzvcf merged_${SPECIES}_raw_allsites_mid.vcf.gz --recode --recode-INFO-all --out merged_${SPECIES}_dp_gq_allsites_mid --minDP 8 --minGQ 25

gzip merged_${SPECIES}_dp_gq_allsites_mid.recode.vcf

done


module load bcftools/1.22

for SPECIES in cod rattlesnake rightwhale
do

bcftools filter -i 'COUNT(FORMAT/GT="mis")<3 && INFO/MQ>25' -Oz -o merged_${SPECIES}_filtered_allsites_mid.recode.vcf.gz merged_${SPECIES}_dp_gq_allsites_mid.recode.vcf.gz

bcftools view -m2 -M2 -v snps -Oz -o merged_${SPECIES}_filtered_variants_mid.recode.vcf.gz merged_${SPECIES}_filtered_allsites_mid.recode.vcf.gz 

bcftools +fill-tags merged_${SPECIES}_filtered_variants_mid.recode.vcf.gz -Oz -o merged_${SPECIES}_filtered_variants_updatetags_mid.vcf.gz -- -t all

done

module load vcftools/0.1.16 

for SPECIES in cod rattlesnake rightwhale
do

vcftools --gzvcf merged_${SPECIES}_filtered_variants_updatetags_mid.vcf.gz --recode --recode-INFO-all --out merged_${SPECIES}_filtered_variants_mid_maf --maf 0.1

gzip merged_${SPECIES}_filtered_variants_mid_maf.recode.vcf

done
```

**Repeat masking**

```
#!/bin/bash
#SBATCH --account=def-pbentzen
#SBATCH --job-name=repeat
#SBATCH --output=%x-%j.out
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00

module load seqkit/2.5.1 bedtools/2.31.0

seqkit locate --bed -r -P -p "[a-z]+" ../cod/reference/GCF_902167405.1_gadMor3.0_genomic.fna | bedtools merge -i - > cod_repeats.bed
seqkit locate --bed -r -P -p "[a-z]+" GCA_039880765.1_rSisCat1_p1.0_genomic.fna | bedtools merge -i - > rattlesnake_repeats.bed
seqkit locate --bed -r -P -p "[a-z]+" GCF_028564815.1_mEubGla1.1.hap2._XY_genomic.fna | bedtools merge -i - > rightwhale_repeats.bed


bedtools intersect -v -a ../raw_vcf/merged_cod_filtered_variants_updatetags_mid.vcf.gz -b cod_repeats.bed -wa -header > merged_cod_filtered_variants_mid_repeatmasked.vcf

gzip merged_cod_filtered_variants_mid_repeatmasked.vcf

bedtools intersect -v -a ../raw_vcf/merged_rattlesnake_filtered_variants_updatetags_mid.vcf.gz -b rattlesnake_repeats.bed -wa -header > merged_rattlesnake_filtered_variants_mid_repeatmasked.vcf

gzip merged_rattlesnake_filtered_variants_mid_repeatmasked.vcf

bedtools intersect -v -a ../raw_vcf/merged_rightwhale_filtered_variants_updatetags_mid.vcf.gz -b rightwhale_repeats.bed -wa -header > merged_rightwhale_filtered_variants_mid_repeatmasked.vcf

gzip merged_rightwhale_filtered_variants_mid_repeatmasked.vcf
```

```
sed -e '/^>/! s/[[:lower:]]/Z/g' GCF_028564815.1_mEubGla1.1.hap2._XY_genomic.fna > hardmask_Z_rightwhale_reference.fasta
```
