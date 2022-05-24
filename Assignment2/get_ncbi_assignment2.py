"""
NCBI Handler class for Assignment 2 Programming 3
Data Sciences for Life Sciences
by Jan Rombouts
"""

import os
import time
import pickle
from Bio import Entrez, Medline

class NCBIHandler:
    '''
    Class to fetch data from ncbi using Biopython
    '''
    def __init__(self, pmid, n=10):
        self.pmid = pmid
        self.n_refs = n

    def ncbi_query(self):
    # try implementing retmax for max references instead of ref
        Entrez.email = "j.a.rombouts@hotmail.com"
        results = Entrez.read(Entrez.elink(dbfrom="pubmed",
                                    db="pmc",
                                    LinkName="pubmed_pmc_refs",
                                    id=self.pmid,
                                    api_key='747992ecd1ef352c0c878188cd4804449a08'))
        references = [f'{link["Id"]}' for link in results[0]["LinkSetDb"][0]["Link"]]
        if len(references) < self.n_refs:
            return references
        else:
            return references[0:self.n_refs] #only return first n references


    def make_dir(self):
        current_directory = os.getcwd()
        final_directory = os.path.join(current_directory, "output")
        if not os.path.exists(final_directory):
            os.makedirs(final_directory)
        print('Created output folder in current working directory')


    def download_ncbi_refs(self, pmid):
        Entrez.email = "j.a.rombouts@hotmail.com"
        # for pmid in ref_pmid:
        handle = Entrez.efetch(db="pmc", id=pmid, rettype="XML", retmode="xml",
                            api_key='747992ecd1ef352c0c878188cd4804449a08')
        with open(f'output/{pmid}.xml', 'wb') as file:
            file.write(handle.read())
        auth_list = self.download_article_authors(pmid)
        time.sleep(1/10) # limit to max 10 queries per second
        print(f"Succesfully fetched references from {pmid}")
        return auth_list


    def download_article_authors(self, pmid):
        Entrez.email = "j.a.rombouts@hotmail.com"
        handle = Entrez.efetch(db="pmc", LinkName="pubmed_pmc_refs",
                                id=pmid, rettype='medline',
                                retmode="text", api_key='747992ecd1ef352c0c878188cd4804449a08')
        record = list(Medline.parse(handle))
        auth_list = record[0]["AU"]
        pickle.dump(auth_list, open(f'output/{pmid}.authors.pkl', 'wb'))
        time.sleep(1/10) # limit to max 10 queries per second
        print(f"Succesfully fetched authors from {pmid}")
        return auth_list
