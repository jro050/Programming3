"""
Assignment 4 Programming 3
Data Sciences for Life Sciences
by Jan Rombouts
"""
#!/usr/bin/python
import numpy as np
import argparse as ap
import sys


argparser = ap.ArgumentParser(description="Kmer size")
argparser.add_argument("-kmers", action="store",
                           dest="kmers", required=True, type=int,
                           help="number of kmers used to build this assembly.")
                           
args = argparser.parse_args()
kmer = str(args.kmers)

results = []
input = sys.stdin.readlines()
for line in input:
    if line.startswith('>'):
        start = line.split('N')[0]
        results.append(start)
    else:
        results.append(line.rstrip('\n'))

records = ''.join(results)
sequence = records.split('>')[1:]

#get lengths of contigs
contigs_len = []
for record in sequence:
    seq = len(record)
    contigs_len.append(seq)
sorted_contigs_len =  sorted(contigs_len)[::-1]

def calculate_N50(list_of_lengths):
    '''calculate N50 from the list of lengths
    input: list of lengths
    output: n50'''
    a = 0 
    n50 = np.inf
    for i in range(len(list_of_lengths)):

        if list_of_lengths[i] < n50:
            n50 = list_of_lengths[i] 


        a +=  list_of_lengths[i]

        if a >= sum(list_of_lengths)/2:
            return n50


n50  = calculate_N50(sorted_contigs_len)
sys.stdout.write(f"{kmer},{n50}\n") 
