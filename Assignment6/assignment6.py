"""
Assignment 6 Programming 3
Data Sciences for Life Sciences
by Jan Rombouts
"""
import time
import socket
import networkx as nx
import multiprocessing as mp
from ParseInterface import ParseInterface
from ParsePubmedXML import PubMedXMLParser
from interface_assignment_6 import Interface
from network_comp_assignment_6 import Connect2Network

if __name__ == "__main__":
#  to test :python3 assignment6.py -s --port 4 --host '' /data/dataprocessing/NCBI/PubMed/pubmed21n0001.xml 
    AUTHKEY = b'whathasitgotinitspocketsesss?'
    args = Interface()
    parse_interface = ParseInterface(args.args)
    fpath = parse_interface.folder_path()
    IP = parse_interface.host_name()
    PORTNUM = parse_interface.port_num()
    number_peons = parse_interface.number_peons()
    run_mode = parse_interface.run_mode()
    
    records = PubMedXMLParser(fpath)
    pmids = [pmid for pmid in records.get_pmid()]
    keyword_list = [keywords for keywords in records.get_keywords()]
    publish_date = [dates for dates in records.get_publish_date()]
    article_titles = [title for title in records.get_title()]
    main_authors = [author for author in records.get_main_author()]
    co_authors = [co_author for co_author in records.get_co_authors()]
    journals = [journal for journal in records.get_journal()]
    languages = [language for language in records.get_language()]
    references = [reference for reference in records.get_references()]
    print(len(references))
    print(len(languages))
    print(len(journals))
    print(len(co_authors))
    print(len(main_authors))
    print(len(article_titles))
    print(len(publish_date))
    print(len(keyword_list))
    print(len(pmids))
    # for record in pmid:
    #     print(record)

    # auth_list = xml_file.read_xml()


    # host = Connect2Network(PORTNUM,AUTHKEY,IP)

    # if run_mode == "c":
    #     client = mp.Process(target=host.runclient, args=(number_peons,))
    #     client.start()
    #     client.join()

    # if run_mode == "s":
    #     server = mp.Process(target=host.runserver, args=(main_pmid.download_ncbi_refs, ref_ids))
    #     server.start()
    #     time.sleep(1)
    #     server.join()

    # if run_mode == 'local':
    #     server = mp.Process(target=host.runserver, args=(main_pmid.download_ncbi_refs, ref_ids))
    #     server.start()
    #     time.sleep(1)
    #     client = mp.Process(target=host.runclient, args=(number_peons,))
    #     client.start()


