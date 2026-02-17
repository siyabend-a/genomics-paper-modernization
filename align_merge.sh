GFP_runs=("SRR001998" "SRR001999" "SRR001996" "SRR001997") # Negative control
E2f1_runs=("SRR001990" "SRR001991" "SRR001988" "SRR001989")
Oct4_runs=("SRR002013" "SRR002012" "SRR002014" "SRR002015")
Sox2_runs=("SRR002023" "SRR002024" "SRR002025" "SRR002026")
Nanog_runs=("SRR002005" "SRR002009" "SRR002004" "SRR002007" "SRR002011" "SRR002006" "SRR002008" "SRR002010")
Esrrb_runs=("SRR001995" "SRR001993" "SRR001994" "SRR001992")
Kfl4_runs=("SRR002001" "SRR002003" "SRR002000" "SRR002002")
Tcfcp2l1_runs=("SRR002033" "SRR002032" "SRR002031" "SRR002034")
Stat3_runs=("SRR002019" "SRR002016" "SRR002017" "SRR002018")
Smad1_runs=("SRR002020" "SRR002021" "SRR002022")
Myc_runs=("SRR002042" "SRR002039" "SRR002040" "SRR002041")
Mycn_runs=("SRR002046" "SRR002043" "SRR002044" "SRR002045")
Zfx_runs=("SRR002035" "SRR002037" "SRR002036" "SRR002038")

tfs=("GFP" "Oct4" "Sox2" "Nanog" "Esrrb" "Kfl4" "Tcfcp2l1" "Stat3" "Smad1" "Myc" "Mycn" "Zfx" "E2f1")

mkdir -p output/process/{Oct4,Sox2,Nanog,Esrrb,Kfl4,Tcfcp2l1,Stat3,Smad1,Myc,Mycn,Zfx,E2f1,GFP}
# Align to generate SAM and convert to BAM
for p in ${tfs[@]}; do
    target_p="${p}_runs[@]" 
    for align in "${!target_p}"; do
        
        bowtie2 -x output/index/mm39_index -U data/fastq/$align.fastq.gz -S output/process/$p/$align.sam --threads 6
        sam="output/process/$p/${align}.sam"
        samtools view -bS $sam > ${sam/sam/bam}
        rm $sam
    done

    # Sort BAM
    for bam in output/process/$p/*bam; do
        samtools sort -@ 4 -o ${bam/bam/_sorted.bam} $bam
    done

    merged="${p}_merged.bam"
    # Merge BAM
    samtools merge -@ 6 output/process/$p/$merged output/process/$p/*sorted.bam

    # Index the merged BAM
    samtools index output/process/$p/$merged
done