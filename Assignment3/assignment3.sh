#!/bin/bash

#SBATCH --time 7:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=BLASTPassignmentJanRombouts
#SBATCH --partition=assemblix

chmod u+x assignment3.sh
chmod a+x assignment3.py
export BLASTDB=~/local-fs/datasets/refseq_protein/refseq_protein
mkdir output
for x in {1..16} ; do /usr/bin/time -ao output/timings.txt -f "${x} %e" blastp -query MCRA.faa -db /local-fs/datasets/refseq_protein/refseq_protein -num_threads 16 -outfmt 6 >> output/blastoutput_num_threads_${x}.txt ; done

assignment3.py
