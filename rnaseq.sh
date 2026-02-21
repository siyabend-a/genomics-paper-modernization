hisat2-build -p 4 data/ref/mm39.fa output/index/rnaseq/mm39
hisat2 -p 4 -x output/index/rnaseq/mm39 -U data/fastq/rnaseq/SRR006489.fastq.gz -S output/SRR006489.sam
hisat2 -p 4 -x output/index/rnaseq/mm39 -U data/fastq/rnaseq/SRR189589.fastq.gz -S output/SRR189589.sam
samtools view -@ 4 -bS output/SRR006489.sam | samtools sort -@ 4 -o SRR006489_sorted.bam
samtools view -@ 4 -bS output/SRR189589.sam | samtools sort -@ 4 -o SRR189589_sorted.bam
samtools merge -@ 4 merged_rnaseq.bam SRR189589_sorted.bam SRR006489_sorted.bam
samtools index merged_rnaseq.bam
