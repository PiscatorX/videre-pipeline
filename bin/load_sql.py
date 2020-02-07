#!/usr/bin/env python
import argparse
import sqlite3



parser = argparse.ArgumentParser(description="""Trimm KEGG annotation file and remain entries with KO numbers""")
parser.add_argument("ko_map", type=argparse.FileType('r'))
parser.add_argument('-d','--db_name', help ='SQLite database filename', required = True)
args, unknown = parser.parse_known_args()
db_name  = args.db_name
conx  = sqlite3.connect(db_name)
cur = conx.cursor()

for line in args.ko_map:
    #k141_1042899', 'K01084', 'G6PC\tglucose-6-phosphatase \tEC:3.1.3.9\n']
    contig, K_number, annotation =  line.split("\t",2)
    annotation = annotation.strip().split("\t")
    name, definition = annotation[:2]
    ec_number  = None
    if ("[EC:" in definition):
        definition, ec_number = definition.rstrip("]").split("[EC:")
       
    elif len(annotation) == 3:
         ec_number = annotation[2].replace("EC:",'')
    print(K_number,  name, definition,  ec_number)    
    sql = """INSERT OR IGNORE INTO kegg_reference(K_number, name, identifier, ec_number)
               VALUES(?,?,?,?)"""
    cur.execute(sql, (K_number, name, definition,  ec_number))
conx.commit()
