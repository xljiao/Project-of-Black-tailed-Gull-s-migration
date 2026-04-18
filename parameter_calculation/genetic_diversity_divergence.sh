#!/bin/bash
# Script Purpose: Calculate Tajima's D, Pi, FST, and Dxy.

VCF="filtered_autosomal.vcf"
WINDOW=50000
STEP=25000

# Calculate Tajima's D and Pi using VCFtools (50kb non-overlapping windows)
vcftools --vcf ${VCF} --keep pop.list --TajimaD ${WINDOW} --out pop_tajima
vcftools --vcf ${VCF} --keep pop.list --window-pi ${WINDOW} --out pop_pi

# Calculate Sliding-Window FST (50kb window, 25kb step)
# This provides windowed FST to match Dxy and Pi scale.
vcftools --vcf ${VCF} \
         --weir-fst-pop pop1.list \
         --weir-fst-pop pop2.list \
         --fst-window-size ${WINDOW} \
         --fst-window-step ${STEP} \
         --out pop1_pop2_window_fst

# Calculate Dxy using genomics_general (Simon Martin's scripts)
# Convert VCF to .geno format
python ${parsevcf} -i ${VCF} \
       --minQual 30 \
       -o geno.gz \
       --ploidyMismatchToMissing \
       --ploidy 2 

# Calculate windowed statistics
python ${popgenWindows} \
       -T 10 \
       -w ${WINDOW} \
       -s ${STEP} \
       -m 10 \
       -g geno.gz \
       -o raw_dxy_results.csv.gz \
       -f phased \
       --windType coordinate \
       -p JA -p RU -p KDQ -p ZS -p FJ \
       --popsFile pops.txt 

# Correct Dxy values for SNP-only bias
python correct.py