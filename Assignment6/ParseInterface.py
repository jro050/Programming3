'''
DOCSTRING
'''
import socket

class ParseInterface:
    '''
    DOCSTRING
    '''

    def __init__(self, arguments):
        self.args = arguments


    def number_peons(self):
        if not self.args.n:
            return 4
        else:
            return self.args.n


    def run_mode(self):
        if self.args.c:
            return 'c'
        elif self.args.s:
            return 's'
        else:
            return 'local'


    def port_num(self):
        if not self.args.port:
            sock = socket.socket()
            sock.bind(('', 0))
            return sock.getsockname()[1]
        else:
            return self.args.port


    def host_name(self):
        if not self.args.host:
            return ''
        else:
            return self.args.host


    def folder_path(self):
        return self.args.path
        # if 1 < len(self.args.pubmed_id[0]) > 8:
        #     raise ValueError('Not a valid PMID')
        # else:
        #     return self.args.pubmed_id[0]

