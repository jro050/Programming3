'''
Assignment 5 for Programming 3
Data Science for Life Sciences
Author: Jan Rombouts
'''

from pyspark import SparkFiles
from pyspark import SparkContext
from pyspark.sql import Row
from pyspark.sql import SQLContext
# from interface_assignment5 import interface
from CreateSparkDF import create_spark_df
    
# df = sqlContext.read.csv(SparkFiles.get('us_pop.csv'), header=True, inferSchema= True)
# df.printSchema()
# Input: InterPROscan output file
# test data:

class InspectInterPRO:

    def __init__(self, spark_df):
        self.df = spark_df


    def unique_pro_annots(self):
        return self.df.select('interPRO_accession').distinct().count()
# 1. How many distinct protein annotations are found in the dataset? I.e. how many distinc InterPRO numbers are there?
# df[column].unique()


    def avg_pro_annots(self):

# 2. How many annotations does a protein have on average?
# df[column].mean()
        pass


    def max_go_term(self):
# 3. What is the most common GO Term found?
# argmax(df[column]) / # df[column].sum().max()
        pass


    def avg_feature_size(self):
# 4. What is the average size of an InterPRO feature found in the dataset?
# df[column].mean() > COMBINE WITH Q2 INTO 1 FUNCITON,S THEN CREATE AN EXTRA FUNCTION
# TO ANSWER THE ASSIGNMENT QUESTIONS
        pass


    def top_10_features(self):
# 5. What is the top 10 most common InterPRO features?
# df.value_counts(ascending=True)[0:10]
        pass


    def top_10_percent(self):
# 6. If you select InterPRO features that are almost the same size (within 90-100%) as the protein itself, what is the top10 then?
# top10percent = len(df)*0.1 > top10percent_df = df.iloc[length > top10percent] > df.value_counts(ascending=True)[0:10]
        pass


    # def COMBINEQ5Q6
# 7. If you look at those features which also have textual annotation, what is the top 10 most common word found in that annotation?
# df.value_counts(ascending=True)[0:10]


    def top_10_features(self):
# 8. And the top 10 least common?
# df.value_counts(ascending=False)[0:10] / df.value_counts(ascending=True)[-10:]
        pass


    # def COMBINEANSWERSPREVIOUSQUESTIONS
# 9. Combining your answers for Q6 and Q7, what are the 10 most commons words found for the largest InterPRO features?
# use top10percent_df > df.value_counts(ascending=True)[0:10]


    def cor_coef(self):
# 10. What is the coefficient of correlation ($R^2$) between the size of the protein and the number of features found?
# coeffcor = (stuff)^2
        pass

# Your output should be a CSV file with 3 columns;
# 1. in the first column the question number
# 2. in the second column the answer(s) to the question
# 3. in the third column the output of the scheduler's physical plan (using the `.explain()` PySpark method) as a string

# NB3: Use the `csv` Python module to make the CSV file in "excel" format; this makes it easier to deal with the different answer types (number, string, list etc.)


if __name__ == "__main__":
    path = '/data/dataprocessing/interproscan/all_bacilli.tsv'
    df = create_spark_df(path)
    inspect = InspectInterPRO(df)
    inspect.unique_pro_annots()
