'''
DOCSTRING
'''

import time
import pickle
import os
import networkx as nx
from Bio import Entrez, Medline

# CREATE CLASS TO FEED ALL DICT ELEMENTS BELOW
# Create network with 1st authors as nodes (Authorlist[0]) CHECK
# Attributes:
# MAKE TRY-EXCEPT FOR EACH ITEM IN CASE ONE IS MISSING
# - From MedlineCitation
# - Rest of autorlist (co-authors) CHECK (INITIAL + LASTNAME)
# - Journal (check if electronic is more cited?) CHECK
# - PMID CHECK
# - DateCompleted OR PubDate DONE
# - Article[Language] TODO
# - ArticleTitle CHECK
# - KeywordList CHECK

# Out nodes:
# Parse ref 
# > if 1 item > parse the string
# > if 2 items > get [ArticleIdList][0] == PMID

class PubMedXMLParser:
    '''
    DOCSTRING
    '''
    def __init__(self, fpath):
        self.fpath = fpath
        self.records = self.read_xml()


    def make_dir(self):
        current_directory = os.getcwd()
        final_directory = os.path.join(current_directory, "output")
        if not os.path.exists(final_directory):
            os.makedirs(final_directory)
        print('Created output folder in current working directory')


    def read_xml(self):
        Entrez.email = "j.a.rombouts@hotmail.com"
        with open(self.fpath, "br") as handle:
            records = Entrez.read(handle)
        return records['PubmedArticle']


    def get_pmid(self):
        for record in self.records:
            try:
                yield record["MedlineCitation"]['PMID']
            except:
                yield ' '

    def get_main_author(self):
            for record in self.records:
                try:
                    last_name = record["MedlineCitation"]['Article']['AuthorList'][0]['LastName']
                    initials = record["MedlineCitation"]['Article']['AuthorList'][0]['Initials']
                    yield str(last_name + ' ' + initials)
                except:
                    yield ' '


    def get_co_authors(self):
        for record in self.records:
                try:
                    author_list = []
                    for author in record["MedlineCitation"]['Article']['AuthorList'][1:]:
                        last_name = author['LastName']
                        initials = author['Initials']
                        author_list.append(last_name + ' ' + initials)
                    yield author_list
                except:
                    yield []


    def get_journal(self):
        for record in self.records:
            try:
                journal = record["MedlineCitation"]['Article']['Journal']['Title']
                yield journal
            except:
                yield ' '


    def get_publish_date(self):
        for record in self.records:
            try:
                yield record["PubmedData"]['History'][0]['Year']
            except:
                yield ' '

    def get_title(self):
        for record in self.records:
            try:
                yield record["MedlineCitation"]['Article']['ArticleTitle']
            except:
                yield ' '

    def get_keywords(self):
        for record in self.records:
            try:
                yield record["MedlineCitation"]['KeywordList']
            except:
                yield []


    def get_language(self):
        for record in self.records:
            try:
                yield record["MedlineCitation"]['Article']['Language']
            except:
                yield ' '


    def get_references(self):
        for record in self.records:
            # print(record["PubmedData"]["ReferenceList"])
            ref_list = []
            for ref in  record["PubmedData"]["ReferenceList"][0]:
                print(ref)
                ref_list.append(ref)
                # if reference["ArticleIdList"]:
                #     references.append(reference["ArticleIdList"])
                # else:
                #     references.append(reference["Citation"])
            yield ref_list


#     def ncbi_query(self):
#     # try implementing retmax for max references instead of ref
#         Entrez.email = "j.a.rombouts@hotmail.com"
#         results = Entrez.read(Entrez.elink(dbfrom="pubmed",
#                                     db="pmc",
#                                     LinkName="pubmed_pmc_refs",
#                                     id=self.pmid,
#                                     api_key='747992ecd1ef352c0c878188cd4804449a08'))
#         references = [f'{link["Id"]}' for link in results[0]["LinkSetDb"][0]["Link"]]
#         if len(references) < self.n_refs:
#             return references
#         else:
#             return references[0:self.n_refs] #only return first n references



#     def download_ncbi_refs(self, pmid):
#         Entrez.email = "j.a.rombouts@hotmail.com"
#         # for pmid in ref_pmid:
#         handle = Entrez.efetch(db="pmc", id=pmid, rettype="XML", retmode="xml",
#                             api_key='747992ecd1ef352c0c878188cd4804449a08')
#         with open(f'output/{pmid}.xml', 'wb') as file:
#             file.write(handle.read())
#         auth_list = self.download_article_authors(pmid)
#         time.sleep(1/10) # limit to max 10 queries per second
#         print(f"Succesfully fetched references from {pmid}")
#         return auth_list


#     def download_article_authors(self, pmid):
#         Entrez.email = "j.a.rombouts@hotmail.com"
#         handle = Entrez.efetch(db="pmc", LinkName="pubmed_pmc_refs",
#                                 id=pmid, rettype='medline',
#                                 retmode="text", api_key='747992ecd1ef352c0c878188cd4804449a08')
#         record = list(Medline.parse(handle))
#         auth_list = record[0]["AU"]
#         pickle.dump(auth_list, open(f'output/{pmid}.authors.pkl', 'wb'))
#         time.sleep(1/10) # limit to max 10 queries per second
#         print(f"Succesfully fetched authors from {pmid}")
#         return auth_list
