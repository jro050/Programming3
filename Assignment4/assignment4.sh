#!/bin/bash
chmod u+x assignment4.sh
chmod a+x assignment4.py

# EXPORT NGS=/data/dataprocessing/MinIONData/MG5267/
export NGS1=/data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R1_001_BC24EVACXX.filt.fastq
export NGS2=/data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R2_001_BC24EVACXX.filt.fastq
export OUTPUT=/students/2021-2022/master/Jan_DSLS/output
mkdir -p /students/2021-2022/master/Jan_DSLS/output
mkdir -p output

seq 25 2 31 | parallel -j16 'velveth $OUTPUT/{} {} -longPaired -fastq $NGS1 $NGS1 && velvetg $OUTPUT/{} && cat $OUTPUT/{}/contigs.fa | python3 assignment4.py -kmers {} >> output/output.csv'
python3 best_kmer.py

rm -rf $OUTPUT
