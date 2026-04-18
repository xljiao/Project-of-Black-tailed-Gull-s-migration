# Fastp v0.21 to trim poor-quality reads and adapters

fastp -i ind.R1.fastq.gz -o ind.R1.fp.fastq.gz -I ind.R2.fastq.gz -O ind.fp.fastq.gz
-q 15 -u 50 -n 5 -l 150 -t 0 \
--adapter_sequence AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
--adapter_sequence_r2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
-h ind.fastp.html