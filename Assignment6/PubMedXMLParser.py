'''
DOCSTRING
'''

import time
import pickle
import os
import regex as re
from weakref import ref
import networkx as nx
from Bio import Entrez

# Out nodes:
# Parse ref 
# > if 1 item > parse the string
# > if 2 items > get [ArticleIdList][0] == PMID

# Beide op een andere manier verwerken om te kunnen linken
# aan daadwerkelijke artikelen.
# Niet via Efetch (duurt lang!)
# Je hebt aan het eind van het parsen van alle xml files
# namelijk dezelfde informatie lokaal beschikbaar als NCBI heeft.
# Dus als je een mapping maakt tijdens de eerste ronde parsen
# van auteurs/artikel -> PMID en andersom, dan kun je die
# in een tweede ronde voor het opzetten van de networkx graph eruit halen.

# Hoe de graph eruit ziet hangt af van hoe je de vragen wilt beantwoorden.
# Maar hetzij een artikel als node en referenties tussen artikels als edge
# of auteurs als node en co-auteurschap als edges zijn een goede aanpak
# om de vragen te beantwoorden. De rest van de info (die relevant is voor de vragen!)
# zet je idd in attributes
# (Maar: attributes opvragen is relatief langzaam vergeleken met node/edge operaties
# dus encode de informatie die je echt vaak aanspreekt om de vragen te beantwoorden
# als nodes en edges, en bijkomstigheden in de attributen. Maak desnoods twee graphs!)


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
            # ref_list = []
            for temp_ref in record["PubmedData"]["ReferenceList"]:
                temp_ref = temp_ref["Reference"]
                # try:
                if len(temp_ref[0]) > 1:
                    yield self.get_ref_pmid(temp_ref)
                else:
                    yield self.get_ref_authors(temp_ref)
                # except:
                    # yield []

            # yield ref_list    

    def get_ref_pmid(self, references):
    ### another exception if the attribute of ArticleIdList = pmcid > try to parse with
    ### get_ref_authors, otherwise skip
        ref_list = []
        for ref in references:
            try:
                if (ref["ArticleIdList"][0].attributes)["IdType"] == "pubmed":
                    ref_list.append(ref["ArticleIdList"])
                else:
                    ref_list.append(self.get_ref_authors(ref))
            except:
                    pass
        return ref_list


    def get_ref_authors(self, references):
        regex_patterns = [
            r'[A-Z]{1}\.\ [A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑ]{1}[A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑa-zàáâäãåąčćęèéêëėįìíîïłńòóôöõøùúûüųūÿýżźñçčšž\-]+',
            r'[A-Z]{1}\.\ [A-Z]{1}\.\ [A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑ]{1}[A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑa-zàáâäãåąčćęèéêëėįìíîïłńòóôöõøùúûüųūÿýżźñçčšž\-]+',
            r'[A-Z]{1}\.\ [A-Z]{1}\.\ [A-Z]{1}\.\ [A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑ]{1}[A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑa-zàáâäãåąčćęèéêëėįìíîïłńòóôöõøùúûüųūÿýżźñçčšž\-]+',
            r'[A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑ]{1}[A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑa-zàáâäãåąčćęèéêëėįìíîïłńòóôöõøùúûüųūÿýżźñçčšž\-]+\, [A-Z]{1}\.',
            r'[A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑ]{1}[A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑa-zàáâäãåąčćęèéêëėįìíîïłńòóôöõøùúûüųūÿýżźñçčšž\-]+\, [A-Z]{1}\. [A-Z]{1}\.',
            r'[A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑ]{1}[A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑa-zàáâäãåąčćęèéêëėįìíîïłńòóôöõøùúûüųūÿýżźñçčšž\-]+\, [A-Z]{1}\. [A-Z]{1}\. [A-Z]{1}\.',
            r'[A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑ]{1}[A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑa-zàáâäãåąčćęèéêëėįìíîïłńòóôöõøùúûüųūÿýżźñçčšž\-]+\ [A-Z]{1,3}(?![a-z])'
        ]
        # ref_auth_list = []
        for ref in references:
            print(ref)
            ref_auth_list = [re.findall(pattern, ref["Citation"]) for pattern in regex_patterns]
            if len(ref_auth_list[-1]) > 0:
                ref_auth_list = ref_auth_list[-1]
            print(ref_auth_list)
        # except:
            pass
        pass

# Loop through all items starting at the shortest, if exact match > replace the shorter one 
# Exception if it starts with Name, A. B., C. Name << 

# Use re.findall
# START WITH LONGEST REGEX AND THEN GO DOWN IN SPECIFICITY
# For A. Name, B. Name: [A-Z]{1}\.\ [A-Z]{1}[a-z]+
# For A. B. Name, C. D. Name: [A-Z]{1}\.\ [A-Z]{1}\.\ [A-Z]{1}[a-z]+
# For A. B. C. Name: [A-Z]{1}\.\ [A-Z]{1}\.\ [A-Z]{1}\.\ [A-Z]{1}\.\ [A-Z]{1}[a-z]+[A-Z]{1}[a-z]+
# For Name, A., Name, B.,: [A-Z]{1}[a-z]+\, [A-Z]{1}\.
# For Name, A. B., Name C. D.,: [A-Z]{1}[a-z]+\, [A-Z]{1}\. [A-Z]{1}\.
# For Name AB, Name C: [A-Z]{1}[a-z]+\ [A-Z]{1,3}
# IF Citation = "\xa0" >> drop!
# >>>> Convert all to Surname Initials, put in 1 list and then drop dupes!
# Citation, 2 types
# with A. Name, B. Name, journal >> check if last parsed name is et al; drop it OR Regex for all letter + . + [letters]
# OR A. B. Name, C. D. Name
# OR Name, A., Name, B.,
# with Name A, Name BC. >> regex >> Senstive to CAPS to filter out Journal abbrevs!
# >> Again 


# for ref in ref_dict:
#                     temp_ref_list
#                     try:
#                         if len(ref) > 1:
#                             temp_ref_list.append(ref["ArticleIdList"])
#                             # print(ref["ArticleIdList"])
#                         else:
#                             temp_ref_list.append(ref["Citation"])
#                             # print(ref["Citation"])
#                     except:
#                         temp_ref_list.append(' ')


# Author Format: Surname initials >> Jones B
# Types of CItation:
# ArticleIdList >> PMID

# After that > link PMID to PMID and fetch the appropriate authors
# Also > link author list to author list
