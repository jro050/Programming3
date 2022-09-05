"""
Interface class for Assignment 6 Programming 3
Data Sciences for Life Sciences
by Jan Rombouts
Date: 13-07-2022
"""
import socket
import argparse as ap

class Interface:
    '''
    Interface class for main.py for course PROG3
    Author: Jan Rombouts
    '''
    def __init__(self):
        self.args = self.run_interface()

    def run_interface(self):
        '''
        Takes input arguments from commandline
        '''
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

    def number_peons(self):
        '''
        Return number of peons to use given by user
        Default = 4
        '''
        if not self.args.n:
            return 4
        else:
            return self.args.n


    def run_mode(self):
        '''
        Return mode to run server in given by user
        Default = local
        '''
        if self.args.c:
            return 'c'
        elif self.args.s:
            return 's'
        else:
            return 'local'


    def port_num(self):
        '''
        Return port to use to use given by user
        Default = socket
        '''
        if not self.args.port:
            sock = socket.socket()
            sock.bind(('', 0))
            return sock.getsockname()[1]
        else:
            return self.args.port


    def host_name(self):
        '''
        Return host to use given by user
        Default = local
        '''
        if not self.args.host:
            return ''
        else:
            return self.args.host


    def folder_path(self):
        '''
        Return file path given by user
        '''
        return self.args.path


if __name__ == "__main__":
    Interface()
