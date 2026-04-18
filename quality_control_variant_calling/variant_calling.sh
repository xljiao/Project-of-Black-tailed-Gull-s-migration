#!/bin/bash

REF="ref.genome.fasta"
SAMPLE="ind"                  
THREADS=10                    
PICARD_JAR="picard.jar"        

FQ1="${SAMPLE}.R1.fp.fastq.gz"
FQ2="${SAMPLE}.R2.fp.fastq.gz"

bwa index -a bwtsw ${REF}

bwa mem -t ${THREADS} -R "@RG\tID:${SAMPLE}\tLB:${SAMPLE}\tPL:DNBSEQ\tSM:${SAMPLE}" \
    ${REF} ${FQ1} ${FQ2} | samtools view -@ ${THREADS} -b -S > ${SAMPLE}.bam

# Sort BAM File

java -jar ${PICARD_JAR} SortSam \
    -I ${SAMPLE}.bam \
    -O ${SAMPLE}.sort.bam \
    --REFERENCE_SEQUENCE ${REF} \
    --TMP_DIR ./bam_tmp \
    --SORT_ORDER coordinate

# Mark PCR Duplicates

java -jar ${PICARD_JAR} MarkDuplicates \
    -I ${SAMPLE}.sort.bam \
    -O ${SAMPLE}.mkdup.bam \
    -M ${SAMPLE}_metrics.txt \
    -REMOVE_DUPLICATES false

# Alignment Statistics & Coverage Analysis

java -jar ${PICARD_JAR} CollectAlignmentSummaryMetrics \
    -R ${REF} \
    -I ${SAMPLE}.mkdup.bam \
    -O ${SAMPLE}.map_rate.txt

bedtools genomecov -ibam ${SAMPLE}.mkdup.bam > ${SAMPLE}.total.txt

grep "genome" ${SAMPLE}.total.txt | awk '
    $2==0 {print "Breadth of Coverage =", 1-$5}; 
    {sum+=$2*$3} END {print "Depth of Coverage =",sum/$4}; 
    $2>=10 {total+=$5} END {print "Percents of reads more than 10X =", total}; 
    $2>=15 {all+=$5} END {print "Percents of reads more than 15X =", all}
' > ${SAMPLE}_stat.txt

# Variant Calling (GATK HaplotypeCaller)

gatk HaplotypeCaller \
    -R ${REF} \
    -I ${SAMPLE}.mkdup.bam \
    -O ${SAMPLE}.g.vcf.gz \
    -ERC GVCF \
    --genotyping-mode DISCOVERY \
    --pcr-indel-model CONSERVATIVE \
    --sample-ploidy 2 \
    --min-base-quality-score 5 \
    --kmer-size 10 --kmer-size 25 \
    --native-pair-hmm-threads ${THREADS}

# Cohort Joint Genotyping (Requires multiple gVCFs)
gatk CombineGVCFs -R ${REF} -V gvcf.list -O com.g.vcf.gz
gatk GenotypeGVCFs -R ${REF} -V com.g.vcf.gz -O com.vcf.gz

# GATK Hard Filtering for SNPs

gatk SelectVariants \
    -V com.vcf.gz \
    -select-type SNP \
    -O snps.vcf.gz

gatk VariantFiltration \
    -V snps.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "SNP_FILTER" \
    -O snps_filter_tmp.vcf.gz

# Keep only the variants that passed the hard filter
zcat snps_filter_tmp.vcf.gz | grep -v "SNP_FILTER" > snps_filter_gatk_fin.vcf

# Strict Filtering using VCFtools

vcftools --vcf com.vcf.gz \
    --min-alleles 2 --max-alleles 2 \
    --maf 0.05 --max-maf 0.95 \
    --minQ 30 --minDP 5 --maxDP 100 \
    --max-missing 0.95 \
    --remove-indels \
    --recode \
    --out filtered_final_snps
