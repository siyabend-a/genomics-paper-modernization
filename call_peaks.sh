# TEST
# Call peaks
macs3 callpeak -t output/process/E2f1/E2f1_merged.bam -c output/process/GFP/GFP_merged.bam \
    -f BAM -g mm -n E2f1_vs_GFP --outdir output/peaks/E2f1 -q 0.01 --keep-dup all
wc -l output/peaks/E2f1/E2f1_vs_GFP_peaks.narrowPeak

# Check overlaps
bedtools window -a data/ref/mm39_TSS.bed -b output/peaks/E2f1/E2f1_vs_GFP_summits.bed -w 1000 > output/pre_TFAS.txt # For calculating TFAS(Transcription Association Factor) as mentioned in the study