"""
Assignment 1 Programming 3
Data Sciences for Life Sciences
by Jan Rombouts
"""

import multiprocessing as mp
import argparse as ap
from Bio import Entrez
import os
import time

def ncbi_query(pmid, n=10):
    Entrez.email = "j.a.rombouts@hotmail.com"
    results = Entrez.read(Entrez.elink(dbfrom="pubmed",
                                   db="pmc",
                                   LinkName="pubmed_pmc_refs",
                                   id=pmid,
                                   api_key='747992ecd1ef352c0c878188cd4804449a08'))
    references = [f'{link["Id"]}' for link in results[0]["LinkSetDb"][0]["Link"]]
    return references[0:n] #only return first n references


def make_dir(folder_name):
    current_directory = os.getcwd()
    final_directory = os.path.join(current_directory, folder_name)
    if not os.path.exists(final_directory):
        os.makedirs(final_directory)


def download_ncbi_refs(pmid):
    Entrez.email = "j.a.rombouts@hotmail.com"
    # for pmid in ref_pmid:
    handle = Entrez.efetch(db="pmc", id=pmid, rettype="XML", retmode="text",
                           api_key='747992ecd1ef352c0c878188cd4804449a08')
    with open(f'output/{pmid}.xml', 'wb') as file:
        file.write(handle.read())
    time.sleep(1/10) # limit to max 10 queries per second


def interface():
    argparser = ap.ArgumentParser(description="Script that downloads (default) 10 articles referenced by the given PubMed ID concurrently.")
    argparser.add_argument("-n", action="store",
                           dest="n", required=False, type=int,
                           help="Number of references to download concurrently.")
    argparser.add_argument("pubmed_id", action="store", type=str, nargs=1, help="Pubmed ID of the article to harvest for references to download.")
    args = argparser.parse_args()
    print("Getting: ", args.pubmed_id)
    return args


if __name__ == "__main__":
    args = interface()
    make_dir("output")
    ref_ids = ncbi_query(args.pubmed_id[0])

    cpus = mp.cpu_count()
    with mp.Pool(cpus) as pool:
        results = pool.map(download_ncbi_refs, ref_ids)
    print(f"Succesfully fetched references from {args.pubmed_id[0]}")

