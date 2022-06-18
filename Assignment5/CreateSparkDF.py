'''
Module to create a SparkDF
as explained during Lecture 5
Assignment 5 for Programming 3
Data Science for Life Sciences
Author: Jan Rombouts
'''
from pyspark import SparkFiles
from pyspark import SparkContext
from pyspark.sql import SQLContext


def create_spark_df(path):
    sc = SparkContext('local[16]')
    sqlContext = SQLContext(sc)
    spark_df = sqlContext.read.options(delimiter="\t", header=False).csv(SparkFiles.get(path))
    spark_df = spark_df.toDF("protein_accession",
                            "seq_MD5_digest",
                            "seq_length",
                            "analysis_method",
                            "sig_accession",
                            "sig_description",
                            "start_location",
                            "stop_location",
                            "score",
                            "status_match",
                            "date",
                            "interPRO_accession",
                            "interPRO_description",
                            "GO_annots",
                            "pathway_annots")

    return spark_df
