tfs=("Oct4" "Sox2" "Nanog" "Esrrb" "Kfl4" "Tcfcp2l1" "Stat3" "Smad1" "Myc" "Mycn" "Zfx" "E2f1")
# Call peaks
for p in ${tfs[@]}; do
	macs3 callpeak -t output/process/$p/${p}_merged.bam -c output/process/GFP/GFP_merged.bam \
    	-f BAM -g mm -n ${p}_vs_GFP --outdir output/peaks/$p -q 0.05
	wc -l output/peaks/$p/${p}_vs_GFP_peaks.narrowPeak

	# Check overlaps
	bedtools window -a data/ref/mm39_ref_flat.bed -b output/peaks/$p/${p}_vs_GFP_peaks.narrowPeak -w 50000 > output/${p}_pretfas.txt # For calculating TFAS(Transcription Association Factor) as mentioned in the study
done
