#!/usr/bin/env python
from Bio.KEGG import REST
import collections
import pandas as pd
import argparse
import pprint 
import sys
import csv
import os



class ParseKO(object):

    def __init__(self):
        parser = argparse.ArgumentParser(description="""Trimm KEGG annotation file and remain entries with KO numbers""")
        parser.add_argument("kegg_file", metavar="Kegg file", type=argparse.FileType('r'))
        parser.add_argument("-m","--mapping", default = "ko_map.tsv", help="outfile to save kegg data")
        args, unknown = parser.parse_known_args()
        self.kegg_csv_fobj = csv.reader(args.kegg_file, delimiter="\t")
        self.mapping_file = args.mapping

        
    def processs_KO(self):

        self.seq_ko_map  = collections.defaultdict(list)
        KO_counter  = collections.defaultdict(int)
        for count,  seq_id, KO in self.kegg_annotation():
            KO_counter[KO] += 1
            self.seq_ko_map[seq_id].append(KO)
        #df_ko = pd.DataFrame.from_dict(KO_counter, orient='index')
        #df_ko.to_csv("ko_counts.tsv", sep="\t")
        

    def kegg_annotation(self):
        count = 0
        for row in self.kegg_csv_fobj:
            if (len(row) > 1):
                count += 1
                yield [count] + row
                
        
    def kegg_server(self):
        
        self.ko_names_def = {}
        if os.path.exists(self.mapping_file):
            with open(self.mapping_file, newline='') as map_fobj:
                tsv_reader = csv.reader(map_fobj, delimiter="\t")
                print("*Entries found in existing map file")
                for row in tsv_reader:
                    seq_id,ko,name,definition = row
                    self.ko_names_def[ko] = [name, definition]
                    print("*{}".format("\t".join([seq_id,ko,name,definition])))
                    
        map_fobj =  open(self.mapping_file,"a", newline='')
        map_writer = csv.writer(map_fobj, delimiter="\t", )
        for seq_id, KO_list in self.seq_ko_map.items():
            for ko in KO_list:
                if ko in self.ko_names_def:
                    #print("{} data available skipping".format(ko))
                    continue
                try:
                    name, definition = REST.kegg_list(ko).read().strip().split("\t", 1)[1:][0].split(";", 1)
                    self.ko_names_def[ko] = [name, definition]
                    print("{}".format("\t".join([seq_id,ko,name,definition])))
                    map_writer.writerow([seq_id,ko,name,definition])
                    map_fobj.flush()
                except Exception as e:
                    sys.stderr.write("\t".join([ko, str(e)]))
                    self.ko_names_def[ko] = [ko, str(e)]
                    
        
        

        
        
                

                

if __name__ == "__main__":
    ko = ParseKO()
    ko. processs_KO()
    ko.kegg_server()

        

               
    
