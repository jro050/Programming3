import pandas as pd
import shutil

df = pd.read_csv('/homes/jarombouts/Programming3/Assignment4/output/output.csv', header = None)

kmers = list(df[0].values)
n50 = list(df[1].values)

index = 0 
best_index = 0 
best_n50 = 0
for i in n50:

    if i > best_n50:
        best_index = index
        best_n50 = i
    index += 1
    
best_kmer = kmers[best_index]

original = f'/students/2021-2022/master/Jan_DSLS/output/{best_kmer}/contigs.fa'
target = '/homes/jarombouts/Programming3/Assignment4/output/contigs.fa'

shutil.move(original, target)