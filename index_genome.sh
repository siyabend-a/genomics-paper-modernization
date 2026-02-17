mkdir -p output/index
cp data/ref/mm39.fa.gz output/index
cd output/index
gunzip mm39.fa.gz
bowtie2-build --threads 8 -f mm39.fa mm39_index