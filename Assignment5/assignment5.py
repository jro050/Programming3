'''
Assignment 5 for Programming 3
Data Science for Life Sciences
Author: Jan Rombouts
'''

from pyspark import SparkFiles
sc = SparkContext('local[8]')
sc.addFile(url)
sqlContext = SQLContext(sc)
SparkFiles.get("us_pop.csv")

df = sqlContext.read.csv(SparkFiles.get('us_pop.csv'), header=True, inferSchema= True)
df.printSchema()
# Input: InterPROscan output file
# test data:
# /data/dataprocessing/interproscan/all_bacilli.tsv

# Qs
# 1. How many distinct protein annotations are found in the dataset? I.e. how many distinc InterPRO numbers are there?
# 2. How many annotations does a protein have on average?
# 3. What is the most common GO Term found?
# 4. What is the average size of an InterPRO feature found in the dataset?
# 5. What is the top 10 most common InterPRO features?
# 6. If you select InterPRO features that are almost the same size (within 90-100%) as the protein itself, what is the top10 then?
# 7. If you look at those features which also have textual annotation, what is the top 10 most common word found in that annotation?
# 8. And the top 10 least common?
# 9. Combining your answers for Q6 and Q7, what are the 10 most commons words found for the largest InterPRO features?
# 10. What is the coefficient of correlation ($R^2$) between the size of the protein and the number of features found?


# Your output should be a CSV file with 3 columns;
# 1. in the first column the question number
# 2. in the second column the answer(s) to the question
# 3. in the third column the output of the scheduler's physical plan (using the `.explain()` PySpark method) as a string

# NB1: Make sure you use the /commons/conda environment
# NB2: Use only 16 threads maximum; `sc = SparkContext('local[16]')`
# NB3: Use the `csv` Python module to make the CSV file in "excel" format; this makes it easier to deal with the different answer types (number, string, list etc.)