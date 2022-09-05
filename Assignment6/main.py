"""
Main script to run Assignment 6 Programming 3
Data Sciences for Life Sciences
by Jan Rombouts
Date: 13-07-2022
"""
import time
import glob
import os
import multiprocessing as mp
import pandas as pd
from pubmed_xml_parser import PubMedXMLParser
from interface import Interface
from network_comp_assignment_6 import Connect2Network
from answer_questions import DaskQuestions

def create_df_csv(pubmed_file):
    '''
    Function to extract information from a PubMed XML file
    Uses PubMedXMLParser to extract the information
    The information is put into a Pandas DF and saved as CSV in csv_path
    Parameters:
    Input:
        - pubmed_file: a file path to PubMed XML file
    Returns:
        - record_df: created Pandas DF with extracted info of PubMed file
    '''

    csv_path = '/students/2021-2022/master/Jan_DSLS/prog3_final_assignment'
    fname = (pubmed_file.split('.')[0])
    if os.path.exists(f'{csv_path}/{fname}.csv'):
        return ''
    else:
        records = PubMedXMLParser(pubmed_file)
        pmids = list(records.get_pmid())
        main_authors = list(records.get_main_author())
        co_authors = list(records.get_co_authors())
        keyword_list = list(records.get_keywords())
        publish_date = list(records.get_publish_date())
        article_titles = list(records.get_title())
        references = list(records.get_references())
        journals = list(records.get_journal())
        languages = list(records.get_language())
        record_df = pd.DataFrame(
                        {'pmid': pmids,
                        'main_author': main_authors,
                        'co_author': co_authors,
                        'publish_date': publish_date,
                        'keyword_list': keyword_list,
                        'article_title': article_titles,
                        'journal': journals,
                        'language': languages,
                        'references': references
                        })
        record_df.to_csv(f'{csv_path}/{fname}.csv')
        print(f'pcreated a CSV for {pubmed_file}')
        return record_df


if __name__ == "__main__":
    AUTHKEY = b'whathasitgotinitspocketsesss?'
    args = Interface()
    fpath = args.folder_path()
    IP = args.host_name()
    PORTNUM = args.port_num()
    number_peons = args.number_peons()
    run_mode = args.run_mode()

    os.chdir(fpath)
    pubmed_files = list(glob.glob("*.xml"))

    host = Connect2Network(PORTNUM,AUTHKEY,IP)
    if run_mode == "c":
        client = mp.Process(target=host.runclient, args=(number_peons,))
        client.start()
        client.join()

    if run_mode == "s":
        server = mp.Process(target=host.runserver, args=(create_df_csv, pubmed_files))
        server.start()
        time.sleep(1)
        server.join()

    if run_mode == 'local':
        server = mp.Process(target=host.runserver, args=(create_df_csv, pubmed_files))
        server.start()
        time.sleep(1)
        client = mp.Process(target=host.runclient, args=(number_peons,))
        client.start()

    print("Succesfully created CSVs for all PubMed files :)")

    dask_qs = DaskQuestions('/students/2021-2022/master/Jan_DSLS/prog3_final_assignment')
    avg_co = dask_qs.avg_co_authors()
    year_avg_co = dask_qs.yearly_avg_co_authors()
    year_papers = dask_qs.yearly_papers()
    year_avg_refs = dask_qs.yearly_avg_refs()
    perc_languages = dask_qs.language_of_articles()
    art_per_auth = dask_qs.articles_per_main_author()

    answer_dict = {"Question":[
            "What is the average amount of co authors?",
            "What is the average amount of co authors per year?",
            "How many papers are published per year?",
            "What is the average amount of references per year per article?",
            "What is the distribution of article languages?",
            "What is the average amount of publications a main author has?"],
            "Answer":[avg_co, year_avg_co, year_papers,
                    year_avg_refs, perc_languages, art_per_auth]
        }

    output = pd.DataFrame(answer_dict)
    output.to_csv('~/homes/jarombouts/Programming3/Assignment6/assignment6.csv', sep=',', index=False)