struct_d=/data/jiaoxiaolu/larus/vcf_145/structure/
filer_d=/data/jiaoxiaolu/larus/vcf_145/
vcftools=/home/software/vcftools_0.1.13/bin/vcftools
plipro_d=/data/jiaoxiaolu/larus/vcf_145/plink_format/
plink=/home/software/plink1.9/plink

${vcftools} --vcf ${filer_d}snps_filter_gatk_fin.vcf \
--recode --remove remove_relate --out ${filer_d}snp_tree


## vcf文件提取常染色体
echo "${vcftools} --vcf ${filer_d}snp_tree.recode.vcf \\" > ${filer_d}auto.sh
for chr in `seq 1 34`
do
echo "--chr chr${chr} \\" >> ${filer_d}auto.sh
done
echo "--recode --out ${filer_d}snp_tree_auto" >> ${filer_d}auto.sh

vcf2phylip=/home/jiaoxiaolu_org/software/vcf2phylip/vcf2phylip.py

python ${vcf2phylip} -i ${plipro_d}snp_tree_auto.vcf -o snp_tree_auto -f

### fasttree
fasttree=/home/jiaoxiaolu_org/software/FastTreeMP
fasttree_d=/data/jiaoxiaolu/larus/vcf_145/structure/fasttree/
${fasttree} -nt -gtr < ${njtree_d}snp_tree_auto.min4.fasta > ${fasttree_d}snp_tree_auto.tree

