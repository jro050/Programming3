"""
Interface class for Assignment 2 Programming 3
Data Sciences for Life Sciences
by Jan Rombouts
"""

import argparse as ap
import sys

# assignment2.py -n <number_of_peons_per_client> [-c | -s] --port <portnumber> --host <serverhost> -a <number_of_articles_to_download> STARTING_PUBMED_ID

class Interface:
    def __init__(self):
        self.args = self.interface()

    def interface(self):
        argparser = ap.ArgumentParser(description="Script that downloads (default) 10 articles referenced by the given PubMed ID concurrently.")
        argparser.add_argument("-n", action="store",
                            dest="n", required=False, type=int,
                            help="Number of peons per client.")
        modus = argparser.add_mutually_exclusive_group(required=False)
        modus.add_argument("-c", action="store", help="Starts the Client modus")
        modus.add_argument("-s", action="store", help="Starts the Server modus")
        argparser.add_argument("--port", action="store", required=False, type=int,
                            help="Portnumber")
        argparser.add_argument("--host", action="store", required=True, type=int,
                            help="Serverhost")
        argparser.add_argument("-a", action="store", required=False, type=int,
                            help="Number of articles to download.")
        argparser.add_argument("pubmed_id", action="store", type=str, nargs=1, help="Pubmed ID of the article to harvest for references to download.")
        # args = argparser.parse_args()
        
        if len(sys.argv) == 1:
            argparser.print_help()
        args = argparser.parse_args()
        return args

        # "-n", action="store",
        #                     dest="n", required=True, type=str,
        #                     help="Server or Client modus."
if __name__ == "__main__":
    Interface()