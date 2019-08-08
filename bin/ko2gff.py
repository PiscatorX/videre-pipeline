#!/usr/bin/env python
from Bio.KEGG import REST
import collections
import pandas as pd
import argparse
import sqlite3
import pprint 
import sys
import csv
import os



class ParseKEGG(object):

    def __init__(self):
        parser = argparse.ArgumentParser(description="""Trimm KEGG annotation file and remain entries with KO numbers""")
        parser.add_argument("kegg_file_taxonomy", metavar="[Kegg file taxonomy]", type=argparse.FileType('r'))
        #parser.add_argument("kegg_file_function", metavar="[Kegg file function]", type=argparse.FileType('r'))
        parser.add_argument("-m","--mapping", default = "ko_map.tsv", help="outfile to save kegg data")
        parser.add_argument('-d','--db_name', help ='SQLite database filename', required = True)
        args, unknown = parser.parse_known_args()
        db_name  = args.db_name
        self.conx  = sqlite3.connect(db_name)
        self.cur = self.conx.cursor()
        #self.kegg_function_fobj = csv.reader(args.kegg_file_function, delimiter="\t")
        self.kegg_taxonomy_fobj = csv.reader(args.kegg_file_taxonomy, delimiter="\t")
        self.mapping_file = args.mapping


        
    def processs_taxonomy(self, ):
        
        for data_row in self.kegg_taxonomy_fobj:    
            if (len(data_row) > 1):
                seq_id, KO = data_row
                try:
                     kegg_list = REST.kegg_list(KO).read()
                     identifier, definition = kegg_list.strip().split("\t", 1)[1:][0].split(";", 1)
                     definition = definition.rstrip("]").split("[EC:")
                     name = definition[0]
                     ec_number = ''
                     if len(definition) == 2:
                         ec_number = definition[1]
                     print(name, ec_number)
                except Exception as e:
                     sys.stderr.write("\t".join([KO, str(e)]))
                     

 
if __name__ == "__main__":
    ko = ParseKEGG()
    ko.processs_taxonomy()
    #ko. processs_KO()
    #ko.kegg_server
                                               
    
