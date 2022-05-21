"""
Assignment 2 Programming 3
Data Sciences for Life Sciences
by Jan Rombouts
"""
import multiprocessing as mp
# from multiprocessing.managers import BaseManager, SyncManager
import time
# import argparse as ap
# from Bio import Entrez

from interface_assignment_2 import Interface
from network_comp_assignment_2 import Connect2Network
from get_ncbi_assignment2 import NCBIHandler

# Download and analyze the NCBI articles from different "remote" clients
# You can test all this using '' as your IP (remember this means : localhost)
# For the time being, start clients by hand using the "ssh" command to connect terminals to different computers on the BIN network.

# # assignment2.py -n <number_of_peons_per_client> [-c | -s] --port <portnumber> --host <serverhost> -a <number_of_articles_to_download> STARTING_PUBMED_ID

# NB2: if you want to run truly networked, specify at least two hosts: the first one is assumed to run the server process, all hosts after that are clients
# NB3: you need to both save the downloaded xml files and communicate the reference PUBMEDs back to the server

if __name__ == "__main__":
    IP = ''
    PORTNUM = 5381
    AUTHKEY = b'whathasitgotinitspocketsesss?'
    # input_args = Interface()
    pmid = "30049270"
    main_pmid = NCBIHandler(pmid)
    ref_ids = main_pmid.ncbi_query()

    server_start = Connect2Network(PORTNUM,AUTHKEY,IP)
    server = mp.Process(target=server_start.runserver, args=(NCBIHandler.download_ncbi_refs, ref_ids))
    server.start()
    time.sleep(1)
    client = mp.Process(target=server_start.runclient, args=(4,))
    client.start()
    server.join()
    client.join()


# 1.2.1. In this first test iteration, you need to start clients manually:
# this means using ssh to log in to the hosts where the clients should run and starting them by hand
# 1.2.2. In this first iteration, only specify the first host (the server) with the "–host": CHECK
# this is the host where the server runs and which the clients should try to contact CHECK
# 1.2.3. You need to specify an option "-c" or "-s" which specifies whether your script starts as a server, or as a client. CHECK
# You can use the argparse.addmutuallyexclusivegroup option for this! CHECK
# 1.2.4. Perhaps superflously: if you start as a server the script needs to run on the host
# specified in "–host" argument (otherwise you can't bind the address in QueueManager)
# 1.2.5. Don't forget to run clients in "tmux" sessions!