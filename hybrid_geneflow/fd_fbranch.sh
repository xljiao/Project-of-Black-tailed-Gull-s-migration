#!/bin/bash


VCF_AUTO="snp_tree_auto.recode.vcf"    
POP_MAP="all.map"                      
TREE="tree.nwk"                        
OUTGROUP="OUT"                         

# Executable paths
DSUITE="Dsuite"                        
DTOOLS="dtools.py"                     

# Output directories
OUT_DIR="output"

# Calculate Genome-wide D-statistics (Dtrios)
${DSUITE} Dtrios ${VCF_AUTO} ${POP_MAP} \
          -t ${TREE} \
          -c \
          -o ${OUT_DIR}/all_auto

# Heuristic f-branch Calculation (Explaining f4-ratio results)
${DSUITE} Fbranch ${TREE} ${OUT_DIR}/all_auto_tree.txt > ${OUT_DIR}/fbranch_auto.out

# Plotting Fbranch Results
python ${DTOOLS} ${OUT_DIR}/fbranch_auto.out ${TREE} \
       --outgroup ${OUTGROUP} \
       --dpi 300 \
       --tree-label-size 10 \
       -n fbranch_auto_plot

