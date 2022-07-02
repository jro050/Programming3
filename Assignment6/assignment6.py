"""
Assignment 6 Programming 3
Data Sciences for Life Sciences
by Jan Rombouts
"""
import time
import socket
import networkx as nx
import multiprocessing as mp
from interface_assignment_2 import Interface
from network_comp_assignment_2 import Connect2Network

class Graphstuff:
    '''
    Function to test for valid args for assignment 2
    Parameters:
        Input: arguments
        Output: validated or standard arguments
    '''
    def __init__(self, path):
        self.path = path


    def parse_xml(self):
        parse_graph = nx.read_graphml(self.path)



if __name__ == "__main__":
#  to test :python3 assignment6.py -p testfile.xml -s --port 4 --host ''
    AUTHKEY = b'whathasitgotinitspocketsesss?'
    interface = Interface()

#   input_args = TestArgs(interface.args)
    fpath = interface.args.p
    IP = input_args.test_host()
    PORTNUM = input_args.test_port()
    number_peons = input_args.test_peons()
    pmid = input_args.test_pmid()
    run_mode = input_args.test_mode()

    main_pmid = NCBIHandler(pmid, number_refs)
    ref_ids = main_pmid.ncbi_query()
    main_pmid.make_dir()

    host = Connect2Network(PORTNUM,AUTHKEY,IP)
    if run_mode == "c":
        client = mp.Process(target=host.runclient, args=(number_peons,))
        client.start()
        client.join()

    if run_mode == "s":
        server = mp.Process(target=host.runserver, args=(main_pmid.download_ncbi_refs, ref_ids))
        server.start()
        time.sleep(1)
        server.join()

    if run_mode == 'local':
        server = mp.Process(target=host.runserver, args=(main_pmid.download_ncbi_refs, ref_ids))
        server.start()
        time.sleep(1)
        client = mp.Process(target=host.runclient, args=(number_peons,))
        client.start()


#   def __init__(self, arguments):
#         self.args = arguments


#     def test_peons(self):
#         if not self.args.n:
#             return 4
#         else:
#             return self.args.n


#     def test_mode(self):
#         if self.args.c:
#             return 'c'
#         elif self.args.s:
#             return 's'
#         else:
#             return 'local'


#     def test_port(self):
#         if not self.args.port:
#             sock = socket.socket()
#             sock.bind(('', 0))
#             return sock.getsockname()[1]
#         else:
#             return self.args.port


#     def test_host(self):
#         if not self.args.host:
#             return ''
#         else:
#             return self.args.host


#     def test_pmid(self):
#         if 1 < len(self.args.pubmed_id[0]) > 8:
#             raise ValueError('Not a valid PMID')
#         else:
#             return self.args.pubmed_id[0]