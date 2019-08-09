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
        parser.add_argument("kegg_file", metavar="[Kegg file taxonomy]", type=argparse.FileType('r'))
        parser.add_argument('-d','--db_name', help ='SQLite database filename', required = True)
        args, unknown = parser.parse_known_args()
        db_name  = args.db_name
        self.conx  = sqlite3.connect(db_name)
        self.cur = self.conx.cursor()
        #self.kegg_function_fobj = csv.reader(args.kegg_file_function, delimiter="\t")
        self.kegg_file_fobj = csv.reader(args.kegg_file, delimiter="\t")
    
        
    def processs_function(self):
        
        for data_row in self.kegg_file_fobj:
            
            if (len(data_row) > 1):
                seq_id, K_number = data_row
                self.cur.execute("""SELECT K_number from kegg_function where K_number = ?""",(K_number,))
                name = definition = identifier = ec_number = None
                if not self.cur.fetchone():
                    try:
                        kegg_list = REST.kegg_list(K_number).read()
                        identifier, definition = kegg_list.strip().split("\t", 1)[1:][0].split(";", 1)
                        definition = definition.rstrip("]").split("[EC:")
                        name = definition[0]
                        ec_number = None
                        if len(definition) == 2:
                            ec_number = definition[1]
                    except Exception as e:
                        sys.stderr.write("\t".join([K_number, str(e)]))
                    finally:
                        self.cur.execute("""INSERT OR IGNORE INTO kegg_function(K_number, name, definition, ec)
                                             VALUES(?,?,?,?)""",(K_number, name, identifier, ec_number))
                else:
                    self.cur.execute("""INSERT OR IGNORE INTO kegg_map(seq_id, K_number)
                                             VALUES(?,?)""",(seq_id, K_number))
                print(seq_id, K_number, identifier, name, ec_number)
                self.conx.commit()


                      
if __name__ == "__main__":
    ko = ParseKEGG()
    ko.processs_function()
    
    #ko.kegg_server
                                               
    
