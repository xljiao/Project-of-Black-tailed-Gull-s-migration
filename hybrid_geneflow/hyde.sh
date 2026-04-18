#!/bin/bash

VCF_AUTO="snp_tree_auto.recode.vcf" 
POP_MAP="group.txt"                   
PREFIX="hyde_run"                       

# IMPORTANT HyDe Parameters:
NUM_TAXA=130                            # Number of individuals/taxa in the alignment
NUM_SITES=4726711                       # EXACT number of sites (SNPs) in the alignment
TRIPLETS=6                              # Number of taxa to test as outgroups (adjust if needed)
THREADS=12                              # Number of CPU threads

# External Scripts
VCF2PHYLIP="vcf2phylip.py"             

# Output directory
HYDE_DIR="./hyde_output"

# Convert VCF to PHYLIP Format
python ${VCF2PHYLIP} \
       -i ${VCF_AUTO} \
       -o ${HYDE_DIR}/${PREFIX}

# Run HyDe Analysis

run_hyde_mp.py \
    -i ${PHY_OUT} \
    -m ${POP_MAP} \
    -o ${HYDE_DIR}/${PREFIX}_OUT \
    -n ${NUM_TAXA} \
    -s ${NUM_SITES} \
    -t ${TRIPLETS} \
    -j ${THREADS} \
    --prefix ${PREFIX}_group