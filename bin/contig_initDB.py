#! /usr/bin/env python
from Bio import SeqIO
import collections
import  sqlite3
import argparse
import pprint
import time
import sys
import os



class DB_Connect(object):
    
    def __init__(self):

        """
          load blast results in SQLite database

        """
        
        parser = argparse.ArgumentParser(description="""Connect to SQlite database""")
        parser.add_argument('-d','--db_name', help ='SQLite database filename', required = True)
        parser.add_argument('-c','--contigs', help ='fasta file with contigs')
        parser.add_argument('-f','--fix-id', dest = "fix_id", action="store_true", help ='remove description from fasta def line')
        parser.add_argument('-o','--outfile', help = "Output filename for fixed contig and tsv filename", required = True)
        self.args, unknown = parser.parse_known_args()
        self.db_name  = self.args.db_name
        self.db_path, self.db_fname  = os.path.split(self.db_name)
        self.conx  = sqlite3.connect(self.db_name)  
        self.contigs  = self.args.contigs
        
        if self.args.fix_id:
            if not self.args.contigs:
                raise argparse.ArgumentTypeError('fix-id requires fasta contig file. See -h/--help')
            fname, ext = os.path.splitext(self.args.outfile)
            if not ext:
                raise argparse.ArgumentTypeError('Contig file should have an *.fasta or *.fa extension')
            self.contigs_newfname_fp  = open(self.args.outfile, 'w')
            self.tsv_fp =  open('.'.join([fname, "tsv"]), 'w')
            
            
    def get_contig(self): 

        self.init_contigTable()
        contig_records = SeqIO.parse(self.contigs, "fasta")
        id_counts =  collections.defaultdict(int)
        contig_data = {}
        for contig in contig_records:
            id_counts[contig.id]=+1
            contig_id = contig.id+str(id_counts[contig.id]) if id_counts[contig.id] != 1 else contig.id 
            if self.args.fix_id:
                contig.id = contig_id
                descr = contig.description.split(" ",1)[1:][0]
                contig.description = ''
                contig_data = dict([field.split("=") for field in descr.split(' ')])
                SeqIO.write(contig, self.contigs_newfname_fp,  "fasta")
                self.contigs_newfname_fp.flush()
                print("{}\t{}".format(contig.id, descr),file=self.tsv_fp, flush = True) 
            contig_data.update({'id': contig_id , 'description': descr , 'length': str(len(contig))})
            cols = ', '.join(contig_data.keys())
            values = ', '.join([ "'"+val+"'" for val in contig_data.values() ])
            sql = "INSERT INTO contigs({}) VALUES ({})".format(cols, values )
            self.conx.execute(sql)
        self.conx.commit()
        if self.args.fix_id:
            self.contigs_newfname_fp.close()
            self.tsv_fp.close()
           

        
    def init_contigTable(self):
        
        contigs_table = """
        CREATE TABLE IF NOT EXISTS contigs(
        id VARCHAR(20) PRIMARY KEY,
        description VARCHAR,
        flag INT,   
        len  INT, 
        multi REAL,
        length INT)"""
        self.conx.execute(contigs_table)

        
    def init_blastTables(self):
        taxa_assignment  = """CREATE TABLE taxa_assignment(
                              contig_id  VARCHAR(20),
                              outcome VARCHAR(20),
                              hit_id  VARCHAR(20),
                              FOREIGN KEY(contig_id) REFERENCES contigs(id))""" 


        hit_ids = """CREATE TABLE hit_def(
                     accession VARCHAR(20)  NOT NULL,
                     DNA_END INTNOT NULL,
                     DNA_ID VARCHAR(20)  NOT NULL,
                     DNA_ORIENTATION VARCHAR(1) NOT NULL,
                     DNA_START INT   NOT NULL,
                     LENGTH INT NOT NULL,
                     ORGANISM VARCHAR,
                     TAXON_ID INT NOT NULL)"""        
        
        hsps ="""CREATE TABLE hsps(
                 hit_id VARCHAR(20),
                 align_length INT NOT NULL,
                 bits	REAL NOT NULL,
                 contig_id VARCHAR(20),
                 expect	REAL NOT NULL,
                 frame	VARCHAR NOT NULL,
                 gaps	INT NOT NULL,
                 identities	INT NOT NULL,
                 mismatches	INT NOT NULL,
                 perc_id	INT NOT NULL, 
                 positives	INT NOT NULL,
                 query_end	INT NOT NULL,
                 query_start	INT NOT NULL,
                 sbjct_end	INT NOT NULL, 
                 sbjct_start	INT NOT NULL,
                 score	REAL NOT NULL,
                 FOREIGN KEY(hit_id) REFERENCES hit_ids(accession))"""

        
        for sql_cmd in [hit_ids, taxa_assignment, hsps]:
            self.conx.execute(sql_cmd)



        
if __name__ == '__main__':
    db_connect = DB_Connect()
    if db_connect.contigs:
        db_connect.get_contig()
    db_connect.init_blastTables()
