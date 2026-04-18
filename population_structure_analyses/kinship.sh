#!/bin/bash
VCF_IN="snps_filter.recode.vcf" 
PREFIX="snps_filter_plink"         
THREADS=20                         
CHR_SET=34                          

PLINK_DIR="./plink_format"
KING_DIR="./king"


# Convert VCF to PLINK Format

vcftools --vcf ${PLINK_DIR}/${PREFIX}.vcf \
         --plink \
         --out ${PLINK_DIR}/${PREFIX}

plink --file ${PLINK_DIR}/${PREFIX} \
      --make-bed \
      --out ${PLINK_DIR}/${PREFIX} \
      --threads ${THREADS} \
      --allow-extra-chr \
      --chr-set ${CHR_SET}

# Kinship Calculation (KING)

# Run KING to calculate kinship coefficients
king -b ${PREFIX}.bed \
     --fam ${PREFIX}.fam \
     --bim ${PREFIX}.bim \
     --kinship \
     --prefix snps_king

# Filter for pairs with a kinship coefficient > 0.0442 (typically 3rd-degree relatives or closer)
# Note: Column 8 ($8) in KING's .kin0 output represents the estimated Kinship coefficient
awk '{if($8 > 0.0442) print $0}' snps_king.kin0 > related_individuals_filtered.txt
