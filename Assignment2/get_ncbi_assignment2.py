"""
NCBI Handler class for Assignment 2 Programming 3
Data Sciences for Life Sciences
by Jan Rombouts
"""

import os
import time
from Bio import Entrez, Medline

# Your script needs to analyze the XML of each of the references further to extract all the authors of the article. CHECK
# It should save the authors in a Python tuple and use the Pickle module to save it to the disk as CHECK
# output/PUBMED_ID.authors.pickle where PUBMEDID is of course the pubmed ID of the article in question CHECK
# NB3: you need to both save the downloaded xml files and communicate the reference PUBMEDs back to the server !!!

class NCBIHandler:
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
        self.download_article_authors(pmid)
        time.sleep(1/10) # limit to max 10 queries per second
        print(f"Succesfully fetched references from {pmid}")


    def download_article_authors(self, pmid):
        Entrez.email = "j.a.rombouts@hotmail.com"
        handle = Entrez.efetch(db="pmc", LinkName="pubmed_pmc_refs",
                                id=pmid, rettype='medline',
                                retmode="text", api_key='747992ecd1ef352c0c878188cd4804449a08')
        record = list(Medline.parse(handle))
        auth_list = record[0]["AU"]
        with open(f'output/{pmid}.authors.pickle', 'w') as file:
                file.write(str(auth_list))
        time.sleep(1/10) # limit to max 10 queries per second
        print(f"Succesfully fetched authors from {pmid}")

# if __name__ == "__main__":
#     args = interface()
#     make_dir("output")
#     ref_ids = ncbi_query(args.pubmed_id[0])

#     cpus = mp.cpu_count()
#     with mp.Pool(cpus) as pool:
#         results = pool.map(download_ncbi_refs, ref_ids)
#     print(f"Succesfully fetched references from {args.pubmed_id[0]}")

