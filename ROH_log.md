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

Array jobs 1,3-5,7-8 are runnin this:

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
