'''
Script to answer questions for Assignment 6 Programming 3
Data Sciences for Life Sciences
by Jan Rombouts
Date: 22-07-2022
'''
import numpy as np
import pandas as pd
import dask.dataframe as dd


class DaskQuestions:
    '''
    Class to answer questions about PudMed data for course PROG3
    Author: Jan Rombouts
    '''

    def __init__(self, csv_path):
        '''
        Initiates DaskQuestions class
        Parameters:
        Input:
        - File path for all PubMed csv folders as created by create_df_csv
          in main.py
        Returns:
        - ddf: dask dataframe with all csv files
        - years: list of all publishing years found in ddf
        - len_df: length of the full ddf
        '''
        self.ddf = dd.read_csv(f'{csv_path}/*.csv')
        self.ddf['publish_date'] = self.ddf['publish_date'].astype(int)
        self.years = np.unique(self.ddf['publish_date'])
        self.len_df = len(self.ddf)


    def avg_co_authors(self):
        '''
        Calculates the average amount of co authors per article
        Returns:
            - average amount of co authors per article (int)
            also save outcome in avg_co_authors.npy as backup
        '''
        print("started avg_co_authors")
        co_list = list(self.ddf["co_author"])
        print(f'list is {len(co_list)} long')
        co_list = [len(co_list[i][1:-1].replace("'",'').split(',')) for i in range(len(co_list))]
        avg_co_auths = sum(co_list) / self.len_df
        np.save("avg_co_authors", avg_co_auths)
        return avg_co_auths


    def yearly_avg_co_authors(self):
        '''
        Calculates the average amount of co authors per article per year
        Returns:
            - year_co_auths (dict): key =  year, value  = average amount of co authors per article that year
            also save outcome in year_co_auths.npy as backup
        '''
        year_co_auths = {}
        for year in self.years:
            sub_df = self.ddf[self.ddf["publish_date"] == year]
            len_sub_df = len(sub_df)
            print(f"length of {year} is {len_sub_df}")
            co_list = list(sub_df["co_author"])
            co_list = [len(co_list[i][1:-1].replace("'",'').split(',')) for i in range(len_sub_df)]
            avg_co_author = sum(co_list) / len_sub_df
            year_co_auths[year] = avg_co_author
            print(f"finished year {year}")
        np.save("year_co_auths", year_co_auths)
        return year_co_auths


    def yearly_papers(self):
        '''
        Calculates the total amount of papers per year
        Returns:
            - paper_per_year (dict): key =  year, value  = total amount of articles that year
            also save outcome in paper_per_year.npy as backup
        '''
        paper_per_year = {}
        
        for year in self.years:
            sub_df = self.ddf[self.ddf["publish_date"] == year]
            paper_per_year[year] = len(sub_df)
            print(f"finished year {year}")
        np.save("paper_per_year", paper_per_year)
        return paper_per_year


    def yearly_avg_refs(self):
        '''
        Calculates the average amount of references per year
        Returns:
            - year_refs (dict): key =  year, value  = average amount of references per year
            also save outcome in yearly_avg_refs.npy as backup
        '''
        year_refs = {}
        for year in self.years:
            sub_df = self.ddf[self.ddf["publish_date"] == year]
            refs = list(sub_df["references"])
            len_sub_df = len(sub_df)
            print(f"length of {year} is {len_sub_df}")
            refs = [len(refs[i][1:-1].replace("'",'').split(',')) for i in range(len_sub_df)]
            avg_refs = sum(refs) / len_sub_df
            year_refs[year] = avg_refs
            print(f"finished year {year}")
        np.save("yearly_avg_refs", year_refs)
        return year_refs


    def language_of_articles(self):
        '''
        Calculates the percentage of articles in a languages
        Takes the first languages found in language if there are multiple
        Returns:
            - year_refs (dict): key =  language, value  = percentage of articles in that language
            also save outcome in language_dict.npy as backup
        '''
        language_dict = {}
        languages = np.unique([lang[2:5] for lang in np.unique(self.ddf['language'])])
        for lang in languages:
            print(f"started {lang}")
            lang = "['" + lang + "']"
            sub_df = self.ddf[self.ddf["language"] == lang]
            perc_lang = (len(sub_df) / self.len_df) *100
            language_dict[lang] = perc_lang
            print(f"finished {lang}")
        print(language_dict)
        np.save("language_dict", language_dict)
        return language_dict


    def articles_per_main_author(self):
        '''
        Calculates the average amount of articles each main author has
        Returns:
            - average amount of articles per main author (int)
            also save outcome in art_per_main_auth.npy as backup
        '''
        print("started calculating articles per main author")
        main_auths = len(np.unique(self.ddf['main_author']))
        print("finished calculating articles per main author")
        art_per_main_auth = self.len_df / main_auths
        np.save("art_per_main_auth", art_per_main_auth)
        return art_per_main_auth


if __name__ == "__main__":
    csv_path = '/students/2021-2022/master/Jan_DSLS/prog3_final_assignment'
    dask_qs = DaskQuestions(csv_path)
    avg_co = dask_qs.avg_co_authors()
    print("answered q1")
    year_papers = dask_qs.yearly_papers()
    print("answered q2")
    perc_languages = dask_qs.language_of_articles()
    print("answered q3")
    art_per_auth = dask_qs.articles_per_main_author()
    print("answered q4")
    year_avg_refs = dask_qs.yearly_avg_refs()
    print("answered q5")
    year_avg_co = dask_qs.yearly_avg_co_authors()
    print("answered q6")


    answer_dict = {"Question":[
        "What is the average amount of co authors?",
        "What is the average amount of co authors per year?",
        "How many papers are published per year?",
        "What is the average amount of references per year per article?",
        "What is the distribution of article languages?",
        "What is the average amount of publications a main author has?"],
        "Answer":[dask_qs, avg_co, year_avg_co, year_papers,
                year_avg_refs, perc_languages, art_per_auth]
    }

    output = pd.DataFrame(answer_dict)
    output.to_csv('~/homes/jarombouts/Programming3/Assignment6/assignment6.csv', sep=',', index=False)