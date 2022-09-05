'''
Script to Parse a PubMed XML file downloaded with Biopython
Author: Jan Rombouts
Date: 13-07-2022
'''

import regex as re
from Bio import Entrez

class PubMedXMLParser:
    '''
    Class with modules to extract information from PubMed XML file
    Author: Jan Rombouts
    '''
    def __init__(self, fpath):
        '''
        Initiates PubMedXMLParser class
        Parameters:
        Input:
            - File path for a single PubMed XML file
        Returns:
            - fpath: given file path
            - records: all unparsed articles in fpath
            - regex_patterns: regex used in multiple class methods
        '''
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


    def read_xml(self):
        '''
        Opens and reads PubMedXML file
        Returns:
            - list: all PubMed articles in fpath
        '''
        Entrez.email = "j.a.rombouts@hotmail.com"
        with open(self.fpath, "br") as handle:
            records = Entrez.read(handle)
        return records['PubmedArticle']


    def get_pmid(self):
        '''
        Extracts PMID from PubMed record
        Returns:
            - 1 record PMID
            - If no PMID found, returns ''
        '''
        for record in self.records:
            try:
                yield record["MedlineCitation"]['PMID']
            except:
                yield ' '


    def get_main_author(self):
        '''
        Extracts main author from PubMed record
        Returns:
            - str: record main author in format Lastname AB
            - If no main author found, returns ''
        '''
        for record in self.records:
            try:
                last_name = record["MedlineCitation"]['Article']['AuthorList'][0]['LastName']
                initials = record["MedlineCitation"]['Article']['AuthorList'][0]['Initials']
                yield str(last_name + ' ' + initials)
            except:
                yield ''


    def get_co_authors(self):
        '''
        Extracts co authors from PubMed record
        Returns:
            - list: record co authors in format [Lastname AB, Lastname CD]
            - If no co authors found, returns []
        '''
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
        '''
        Extracts journal from PubMed record
        Returns:
            - str: record journal
            - If no journal found, returns ''
        '''
        for record in self.records:
            try:
                journal = record["MedlineCitation"]['Article']['Journal']['Title']
                yield journal
            except:
                yield ' '


    def get_publish_date(self):
        '''
        Extracts publishing year from PubMed record
        Returns:
            - str: record publishing year
            - If no publshing year found, returns ''
        '''
        for record in self.records:
            try:
                yield record["PubmedData"]['History'][0]['Year']
            except:
                yield ' '


    def get_title(self):
        '''
        Extracts title from PubMed record
        Returns:
            - str: record title
            - If no article title found, returns ''
        '''
        for record in self.records:
            try:
                yield record["MedlineCitation"]['Article']['ArticleTitle']
            except:
                yield ' '


    def get_keywords(self):
        '''
        Extracts keyword list from PubMed record
        Returns:
            - list: record keywords
            - If no keyword list found, returns []
        '''
        for record in self.records:
            try:
                yield record["MedlineCitation"]['KeywordList']
            except:
                yield []


    def get_language(self):
        '''
        Extracts language from PubMed record
        Returns:
            - str: record language
            - If no language found, returns ''
        '''
        for record in self.records:
            try:
                yield record["MedlineCitation"]['Article']['Language']
            except:
                yield ' '


    def get_references(self):
        '''
        Extracts references from PubMed record
        Returns:
            - list: record references
            - If no PMID found, returns []
        '''
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


    def get_ref_pmid(self, references):
        '''
        Extracts PMIDs from reference list
        Parameters:
        Input:
            - reference list
        Returns:
            - list: record PMID references
            - If no PMIDs found, returns []
        '''
        ref_list = []
        for ref in references:
            try:
                if (ref["ArticleIdList"][0].attributes)["IdType"] == "pubmed":
                    pmid = int(ref["ArticleIdList"][0])
                    ref_list.append(pmid)
                else:
                    ref_list.append(self.get_ref_authors(ref))
            except:
                ref_list.append('')
        return ref_list


    def get_ref_authors(self, references):
        '''
        Extracts author names from reference list
        Parameters:
        Input:
            - reference list
        Returns:
            - list: author reference list
        '''
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
        '''
        Extracts unique author names from regex dictionaries
        Parameters:
        Input:
            - regex dictionaries
        Returns:
            - list: author reference list in format [Lastname AB, Lastname CD]
        '''
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
        auth_list = [''.join(map(str, i)) for i in zip(co_auth_dict.keys(), co_auth_dict.values())][::-1]
        return auth_list


    def get_auth_initials(self, author_list):
        '''
        Extracts author initials from reference author list
        Parameters:
        Input:
            - reference author list
        Returns:
            - list: author initials
        '''
        initials = [re.findall(str(self.capital+ r'(?!' + self.lowercase[:-6] + r'])'), name) for name in author_list]
        initials = [''.join(letter) for letter in initials]
        return initials

  
    def get_auth_lastname(self, author_list):
        '''
        Extracts author lastnames from PubMed record
        Parameters:
        Input:
            - author reference list
        Returns:
            - list: author lastnames
        '''
        lastname = [re.findall(self.capital + self.lowercase, name) for name in author_list]
        lastname = [f'{it[0]} ' for it in lastname]
        return lastname
