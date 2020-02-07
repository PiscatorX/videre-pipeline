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
        parser.add_argument("-f","--kegg_function", metavar="[Kegg file taxonomy]", type=argparse.FileType('r'))
        parser.add_argument("-t","--kegg_taxonomy", metavar="[Kegg file taxonomy]", type=argparse.FileType('r'))
        parser.add_argument('-d','--db_name', help ='SQLite database filename', required = True)
        self.args, unknown = parser.parse_known_args()
        db_name  = self.args.db_name
        self.conx  = sqlite3.connect(db_name)
        self.cur = self.conx.cursor()
        #self.kegg_function_fobj = csv.reader(seargs.kegg_file_function, delimiter="\t")
        if self.args.kegg_function:
            self.kegg_function_fobj = csv.reader(self.args.kegg_function, delimiter="\t")
        if self.args.kegg_taxonomy:
            self.kegg_taxonomy_fobj = csv.reader(self.args.kegg_taxonomy, delimiter="\t")

            
    def processs_function(self):
        for data_row in self.kegg_function_fobj:
            print(data_row)
            if (len(data_row) > 1):
                contig_id, K_number = data_row
                #print(data_row)
                self.get_kegg(K_number)
                self.cur.execute("""INSERT INTO contig2ko(contig_id, K_number)
                                    VALUES(?,?)""",(contig_id, K_number))
                
        self.conx.commit()
                

                
    def get_kegg(self, K_number):
        
        print(K_number)
        self.cur.execute("""SELECT K_number from kegg_reference  where K_number = ?""",(K_number,))
        if self.cur.fetchone():
            return 
        name = definition = identifier = ec_number = None
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
            self.cur.execute("""INSERT OR IGNORE INTO kegg_reference(K_number, name, identifier, ec_number)
                                 VALUES(?,?,?,?)""",(K_number, name, identifier, ec_number))
            sys.stderr.write(str(self.cur.lastrowid)+"\n")
            
        self.conx.commit()
        
            
    
    def processs_taxonomy(self):
       
        #['user:k141_151420', 'K01903', 'Bacteria', 'Actinobacteria', 'Ilumatobacter', 'aym:YM304_36100', '48.9062']
        ['conti_id', 'K_number', 'Level_1', 'Level_2', 'Level_3', 'tax_code', 'score']
        check_value  =  lambda x: x if x else None 
        for data_row in self.kegg_taxonomy_fobj:
            data_row[0] = data_row[0].split(":")[1]
            data_row = [ check_value(i) for i in data_row ]
            #print(data_row)
            if data_row[1]:
                self.get_kegg(data_row[1])
            self.cur.execute("""INSERT INTO contig2kegg_taxonomy(contig_id, K_number, Level_1, Level_2, Level_3, organism_code, score)
                                VALUES(?,?,?,?,?,?,?)""",(data_row))
        
        self.conx.commit()
        

        
                      
if __name__ == "__main__":
    ko = ParseKEGG()
    if ko.args.kegg_function:
         ko.processs_function()
    if ko.args.kegg_taxonomy:
        ko.processs_taxonomy()
                                               
    
