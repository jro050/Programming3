'''
DOCSTRING
'''

from dataclasses import replace
import time
import pickle
import os
import regex as re
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
        self.capital = r'[A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑ]{1}'
        self.lowercase = r'[A-ZÀÁÂÄÃÅĄĆČĖĘÈÉÊËÌÍÎÏĮŁŃÒÓÔÖÕØÙÚÛÜŲŪŸÝŻŹÑa-zàáâäãåąčćęèéêëėįìíîïłńòóôöõøùúûüųūÿýżźñçčšž\-\']+'
        self.regex_patterns = [
                str(self.capital + r'\.\ ' + self.capital + self.lowercase),
                str(self.capital + r'\.\ ' + self.capital + r'\.\ ' + self.capital + self.lowercase),
                str(self.capital + r'\.\ ' + self.capital + r'\.\ ' + self.capital + r'\.\ ' + self.capital + self.lowercase),
                str(self.capital + self.lowercase + r'\,\ ' + self.capital + r'\.'),
                str(self.capital + self.lowercase + r'\,\ ' + self.capital + r'\.\ ' + self.capital + r'\.'),
                str(self.capital + self.lowercase + r'\,\ ' + self.capital + r'\.\ ' + self.capital + r'\.\ ' + self.capital + r'\.'),
                str(self.capital + self.lowercase + r'\ ' + self.capital[0:-1] + r',3}(?![a-z])')
        ]

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
            try:
                for temp_ref in record["PubmedData"]["ReferenceList"]:
                    temp_ref = temp_ref["Reference"]
                    if len(temp_ref[0]) > 1:
                        ref_list = self.get_ref_pmid(temp_ref)
                    else:
                        ref_list = self.get_ref_authors(temp_ref)
                yield ref_list
            except:
                yield []

            # yield ref_list    

    def get_ref_pmid(self, references):
        ref_list = []
        for ref in references:
            try:
                if (ref["ArticleIdList"][0].attributes)["IdType"] == "pubmed":
                    ref_list.append(ref["ArticleIdList"])
                else:
                    ref_list.append(self.get_ref_authors(ref))
            except:
                    return ref_list
        return ref_list


    def get_ref_authors(self, references):
        replace_char_dict = {",":"",".":""}
        ref_list = []
        for ref in references:
            ref_auth_lists = [re.findall(pattern, ref["Citation"]) for pattern in self.regex_patterns]
            if len(ref_auth_lists[-1]) > 0:
                ref_list.append(ref_auth_lists[-1])
            else:
                reg1_dict = {}
                for auth_list in ref_auth_lists[0:3]:
                    auth = [auth.translate(str.maketrans(replace_char_dict)) for auth in auth_list]
                    reg1_lname = self.get_auth_lastname(auth)
                    reg1_init = self.get_auth_initials(auth)
                    for i in range(len(reg1_lname)):
                        if reg1_lname[i] in reg1_dict.keys():
                            if len(reg1_lname[i]) > len(reg1_dict[reg1_lname[i]]):
                                reg1_dict[reg1_lname[i]] = reg1_init[i]
                        else:
                            reg1_dict[reg1_lname[i]] = reg1_init[i]
                reg2_dict = {}
                for auth_list in ref_auth_lists[3:6]:
                    auth = [auth.translate(str.maketrans(replace_char_dict)) for auth in auth_list]
                    reg2_lname = self.get_auth_lastname(auth)
                    reg2_init = self.get_auth_initials(auth)
                    for i in range(len(reg2_lname)):
                        if reg2_lname[i] in reg2_dict.keys():
                            if len(reg2_lname[i]) > len(reg2_dict[reg2_lname[i]]):
                                reg2_dict[reg2_lname[i]] = reg2_init[i]
                        else:
                            reg2_dict[reg2_lname[i]] = reg2_init[i]
                ref_list.append(self.get_correct_authors(reg1_dict, reg2_dict))
        return ref_list

    def get_correct_authors(self, reg1_dict, reg2_dict):
        if len(reg1_dict) == 0:
            auth_dict = reg2_dict
            auth_list = [''.join(map(str, i)) for i in zip(auth_dict.keys(), auth_dict.values())]
            return auth_list
        if len(reg2_dict) == 0:
            auth_dict = reg1_dict
            auth_list = [''.join(map(str, i)) for i in zip(auth_dict.keys(), auth_dict.values())]
            return auth_list
        if len(reg1_dict) == len(reg2_dict):
            auth_dict = {**reg1_dict, **reg2_dict}
            auth_list = [''.join(map(str, i)) for i in zip(auth_dict.keys(), auth_dict.values())]
            return auth_list
        if len(reg1_dict) > len(reg2_dict):
            co_auth_dict = reg1_dict
            main_auth_dict = reg2_dict
        else:
            co_auth_dict = reg2_dict
            main_auth_dict = reg1_dict

        main_auth = list(set(main_auth_dict) - set(co_auth_dict))[0]
        co_auth_dict[main_auth] = main_auth_dict[main_auth]
        auth_list = [''.join(map(str, i)) for i in zip(co_auth_dict.keys(), co_auth_dict.values())]
        return auth_list

    def get_auth_initials(self, author_list):
        initials = [re.findall(str(self.capital+ r'(?!' + self.lowercase[:-6] + r'])'), name) for name in author_list]
        initials = [''.join(letter) for letter in initials]
        return initials
        
    def get_auth_lastname(self, author_list):
        lastname = [re.findall(self.capital + self.lowercase, name) for name in author_list]
        lastname = [f'{it[0]} ' for it in lastname]
        return lastname

# Author Format: Surname initials >> Jones B
# Types of CItation:
# ArticleIdList >> PMID

# After that > link PMID to PMID and fetch the appropriate authors
# Also > link author list to author list
