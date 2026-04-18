#!/bin/bash

VCF_IN="filter.vcf"                    
REMOVE_LIST="remove_relate_outg.txt"   
PREFIX="snp_structure"                 
THREADS=10                             
CHR_SET=34                             

STRUCTURE_DIR="./structure"
PLINK_DIR="${STRUCTURE_DIR}/plink_format"
LD_DIR="${STRUCTURE_DIR}/ld_prune"
PCA_DIR="${STRUCTURE_DIR}/pca"
ADMIX_DIR="${STRUCTURE_DIR}/admixture"

# Filter VCF: Remove Related Individuals and Outgroups

vcftools --vcf ${VCF_IN} \
         --remove ${REMOVE_LIST} \
         --recode \
         --out ${STRUCTURE_DIR}/${PREFIX}_unrelated

# Extract Autosomes

# Build an array of chromosome arguments dynamically (e.g., --chr chr1 --chr chr2 ...)
CHR_ARGS=()
for chr in $(seq 1 ${CHR_SET}); do
    CHR_ARGS+=("--chr" "chr${chr}")
done

vcftools --vcf ${STRUCTURE_DIR}/${PREFIX}_unrelated.recode.vcf \
         "${CHR_ARGS[@]}" \
         --recode \
         --out ${STRUCTURE_DIR}/${PREFIX}_auto

# Convert VCF to PLINK Format
vcftools --vcf ${STRUCTURE_DIR}/${PREFIX}_auto.recode.vcf \
         --plink \
         --out ${PLINK_DIR}/${PREFIX}_auto

plink --file ${PLINK_DIR}/${PREFIX}_auto \
      --make-bed \
      --out ${PLINK_DIR}/${PREFIX}_auto \
      --threads ${THREADS} \
      --allow-extra-chr \
      --chr-set ${CHR_SET}

# Linkage Disequilibrium (LD) Pruning
plink --bfile ${PLINK_DIR}/${PREFIX}_auto \
      --indep-pairwise 50kb 1 0.2 \
      --out ${LD_DIR}/${PREFIX}_ld \
      --threads ${THREADS} \
      --allow-extra-chr \
      --chr-set ${CHR_SET}

# Extract the pruned SNPs to create a new pruned dataset
plink --bfile ${PLINK_DIR}/${PREFIX}_auto \
      --extract ${LD_DIR}/${PREFIX}_ld.prune.in \
      --make-bed \
      --out ${LD_DIR}/${PREFIX}_ld \
      --allow-extra-chr \
      --chr-set ${CHR_SET}

# Principal Component Analysis (PCA)
plink --bfile ${LD_DIR}/${PREFIX}_ld \
      --threads 20 \
      --pca 10 \
      --out ${PCA_DIR}/pca_results \
      --allow-extra-chr \
      --chr-set ${CHR_SET}

# ADMIXTURE Analysis

# ADMIXTURE requires the older .ped/.map format (recode 12)
plink --bfile ${LD_DIR}/${PREFIX}_ld \
      --recode 12 \
      --out ${ADMIX_DIR}/${PREFIX}_admix \
      --allow-extra-chr \
      --chr-set ${CHR_SET}

# Loop through K=1 to 7
for K in 1 2 3 4 5 6 7; do 
    echo "Running ADMIXTURE for K=${K}..."
    admixture --cv -j${THREADS} ${PREFIX}_admix.ped ${K} | tee log${K}.out 
done
