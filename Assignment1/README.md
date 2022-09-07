# 1. Assignment 1
## 1.1. Goal
The goal of this assignment is to get used to programmatically querying NCBI and use the multiprocessing.Pool construct for concurrent programming in Python.
You will need to use the Biopython Querying facilities (see: Chapter 9 Biopython Tutorial to download the XML data for 10 articles. In particular, I recommend using the esearch, elink, efetch functions.

## 1.2. Hints
To do this succesfully, be sure to make an NCBI account and set up API tokens! See: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ Since the beginning of this year, you can only make a new NCBI account by linking an existing Microsoft/Google/Facebook etc. account to it. On the NCBI sign up page you can choose an existing login to connect to NCBI (e.g. your Hanze school account, which is a Microsoft account).

NB: Rate-limit your script to below the required requests per second using time.sleep() if necessary. (May require some experimentation.)

## 2. Deliverables
A script called "assignment1.py" that given 1 starting article's Pubmed ID, downloads 10 articles referenced by that first article. It should do this concurrently from PubMed using the Biopython API.

## 3. Command line
The only command-line argument you need is a query PubMed ID to ask Entrez about an article. So your script should be called as follows:

'python3 assignment1.py <pubmed_id>'

Your script should then download 10 articles cited by the article referenced by the Pubmed ID concurrently using the multiprocessing.Pool() and map() constructs! The 10 articles should be saved as PUBMEDID.xml in the directory "output" in the directory your script is run from.
