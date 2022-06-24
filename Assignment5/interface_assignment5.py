'''
Interface for Assignment 5 for Programming 3
Data Science for Life Sciences
Author: Jan Rombouts
'''

import argparse as ap

def interface():
    argparser = ap.ArgumentParser(description="Script that creates and analyzes InterPRO data loaded in a Spark DF.")
    argparser.add_argument("InterPRO_file_path", action="store", type=str,
                           nargs=1,
                           help="path of InterPRO file to analyze")
    args = argparser.parse_args()
    print("Analyzing: ", args.InterPRO_file_path[0])
    return args.InterPRO_file_path[0]


if __name__ == "__main__":
    interface()

