'''
Assignment 5 for Programming 3
Data Science for Life Sciences
Author: Jan Rombouts
'''
import pandas as pd
from pyspark.sql import functions as F
from pyspark.sql.types import IntegerType
from interface_assignment5 import interface
from pyspark.ml.stat import Correlation
from CreateSparkDF import create_spark_df

class InspectInterPRO:
    def __init__(self, spark_df):
        self.df = spark_df

# 1. How many distinct protein annotations are found in the dataset? I.e. how many distinc InterPRO numbers are there?
    def unique_pro_annots(self):
        unique_pro = self.df.select('interPRO_accession').distinct().count()
        explain = spark_explain(self.df.select('interPRO_accession').distinct())
        return  unique_pro, explain

# 2. How many annotations does a protein have on average?
    def avg_pro_annots(self):
        avg_pro =  round(self.df.select('interPRO_accession').count() / self.unique_pro_annots())
        explain = spark_explain(self.df.select('interPRO_accession'))
        return avg_pro, explain

# 3. What is the most common GO Term found?
    def max_go_term(self):
        max_go = self.df.agg({"GO_annots": "max"}).collect()[0][0]
        explain = spark_explain(self.df.agg({"GO_annots": "max"}))
        return max_go, explain

# 4. What is the average size of an InterPRO feature found in the dataset?
    def avg_feature_size(self):
        avg_feat = self.df.agg({"seq_length": "mean"}).collect()[0][0]
        explain = spark_explain(self.df.agg({"seq_length": "mean"}))
        return avg_feat, explain

# 5. What is the top 10 most common InterPRO features?
    def top_10_features(self):
        sub_df = self.df[self.df.interPRO_accession.contains('IPR0')]
        top10list = sub_df.groupby('interPRO_accession').count().orderBy("count", ascending=False).take(10)
        top10 = []
        for nr in range(0,len(top10list)):
            top10.append(top10list[nr][0])
        explain = spark_explain(sub_df.groupby('interPRO_accession').count().orderBy("count", ascending=False))
        return top10, explain

# 6. If you select InterPRO features that are almost the same size (within 90-100%) as the protein itself, what is the top10 then?
    def top_10_percent(self):
        sub_df = self.df[self.df.interPRO_accession.contains('IPR0')]\
                    .withColumn('relative_seq_length',( 
                    self.df['stop_location'] - self.df['start_location'] ) / df['seq_length'])
        top10percent = sub_df[sub_df['relative_seq_length'] > 0.9]
        top10percentlist =  top10percent.groupby('interPRO_accession').count().orderBy("count", ascending=False).take(10)
        top10 = []
        for nr in range(0,len(top10percentlist)):
            top10.append(top10percentlist[nr][0])
        explain = spark_explain(sub_df.groupby('interPRO_accession').count().orderBy("count", ascending=False))
        return top10, explain

# 7. If you look at those features which also have textual annotation, what is the top 10 most common word found in that annotation?
    def top_10_words(self):
        sub_df = self.df[self.df.interPRO_accession.isin('-') == False]
        top10words = sub_df.withColumn('word', F.explode(F.split(F.col('interPRO_description'), ' ')))\
                    .groupBy('word')\
                    .count()\
                    .sort('count', ascending=False)\
                    .take(10)
        top10 = []
        for word in range(0,len(top10words)):
            top10.append(top10words[word][0])
        explain = spark_explain(sub_df.withColumn('word', F.explode(F.split(F.col('interPRO_description'), ' ')))\
                                .groupBy('word')\
                                .count()\
                                .sort('count', ascending=False))
        return top10, explain

# 8. And the top 10 least common?
    def bottom_10_words(self):
        sub_df = self.df[self.df.interPRO_accession.isin('-') == False]
        bottom10words = sub_df.withColumn('word', F.explode(F.split(F.col('interPRO_description'), ' ')))\
                    .groupBy('word')\
                    .count()\
                    .sort('count', ascending=True)\
                    .take(10)
        bottom10 = []
        for word in range(0,len(bottom10words)):
            bottom10.append(bottom10words[word][0])
        explain = spark_explain(sub_df.withColumn('word', F.explode(F.split(F.col('interPRO_description'), ' ')))\
                                .groupBy('word')\
                                .count()\
                                .sort('count', ascending=True))
        return bottom10, explain

# 9. Combining your answers for Q6 and Q7, what are the 10 most commons words found for the largest InterPRO features?
    def com_words_large_features(self):
        sub_df = self.df[self.df.interPRO_accession.contains('IPR0')]\
                    .withColumn('relative_seq_length',( 
                    self.df['stop_location'] - self.df['start_location'] ) / df['seq_length'])
        top10percent = sub_df[sub_df['relative_seq_length'] > 0.9]
        top10words = top10percent.withColumn('word', F.explode(F.split(F.col('interPRO_description'), ' ')))\
                    .groupBy('word')\
                    .count()\
                    .sort('count', ascending=False)\
                    .take(10)
        top10 = []
        for word in range(0,len(top10words)):
            top10.append(top10words[word][0])
        explain = spark_explain(top10percent.withColumn('word', F.explode(F.split(F.col('interPRO_description'), ' ')))\
                                .groupBy('word')\
                                .count()\
                                .sort('count', ascending=False))
        return top10, explain

# 10. What is the coefficient of correlation ($R^2$) between the size of the protein and the number of features found?
    def cor_coef(self):
        self.df = self.df.withColumn("seq_length", self.df["seq_length"].cast(IntegerType()))
        sub_df=self.df.select(self.df.protein_accession,self.df.interPRO_accession,self.df.seq_length)\
                        .filter(self.df.interPRO_accession != "-")\
                        .groupby(self.df.protein_accession,"seq_length")\
                        .count()
        corcoef = sub_df.corr('seq_length', 'count')
        explain = spark_explain(self.df.select(self.df.protein_accession,self.df.interPRO_accession,self.df.seq_length)\
                        .filter(self.df.interPRO_accession != "-")\
                        .groupby(self.df.protein_accession,"seq_length"))
        return corcoef, explain

def spark_explain(spark_command):
    return spark_command._sc._jvm.PythonSQLUtils.explainString(spark_command._jdf.queryExecution(), 'simple')

if __name__ == "__main__":
    # path = '/data/dataprocessing/interproscan/all_bacilli.tsv'
    path = interface()
    df = create_spark_df(path)
    inspect = InspectInterPRO(df)
    q1, e1 = inspect.unique_pro_annots()
    q2, e2 = inspect.avg_pro_annots()
    q3, e3 = inspect.max_go_term()
    q4, e4 = inspect.avg_feature_size()
    q5, e5 = inspect.top_10_features()
    q6, e6 = inspect.top_10_percent()
    q7, e7 = inspect.top_10_words()
    q8, e8 = inspect.bottom_10_words()
    q9, e9 = inspect.com_words_large_features()
    q10, e10 = inspect.cor_coef()

    output = {'question' : range(1,11),
                'answer' : [q1, q2, q3, q4, q5, q6, q7, q8, q9, q10],
                'explain' : [e1, e2, e3, e4, e5, e6, e7, e8, e9, e10]}
output = pd.DataFrame(output)
output.to_csv('assignment5.csv', sep=',', index=False)


    

