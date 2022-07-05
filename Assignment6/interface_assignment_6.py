"""
Interface class for Assignment 2 Programming 3
Data Sciences for Life Sciences
by Jan Rombouts
"""

import argparse as ap

class Interface:
    '''
    Interface class for assignment6.py for course PROG3
    '''
    def __init__(self):
        self.args = self.interface()

    def interface(self):
        argparser = ap.ArgumentParser(description="Script that parses all PubMed XML files in the given folder path.")
        argparser.add_argument("-n", action="store",
                            dest="n", required=False, type=int,
                            help="Number of peons per client.")
        modus = argparser.add_mutually_exclusive_group(required=False)
        modus.add_argument("-c", action="store_true", help="Starts the Client modus")
        modus.add_argument("-s", action="store_true", help="Starts the Server modus")
        argparser.add_argument("--port", action="store", required=False, type=int,
                            help="Portnumber")
        argparser.add_argument("--host", action="store", required=True,
                            help="Serverhost")
        argparser.add_argument(action="store",
                            dest="path", type=str,
                            help="path to folder containing all PubMed XML files")
        args = argparser.parse_args()
        return args


if __name__ == "__main__":
    Interface()
