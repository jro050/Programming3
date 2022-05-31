#!/usr/bin/python3
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('output/timings.txt', sep=' ', header=None).rename({0:'num_threads',1:'time in seconds'}, axis=1)
df = df.drop_duplicates(subset=['num_threads'], keep='last')
print(df)
fig = plt.figure()
plt.plot(df['num_threads'], df['time in seconds'], 'k.-')
plt.title("Duration of Blastp with increasing number of threads")
plt.xlabel("Number of Threads")
plt.ylabel("Time taken (s)")
plt.savefig('output/timings.png')
