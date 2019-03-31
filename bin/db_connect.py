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
        parser.add_argument('db_name', help ='SQLite database filename')
        parser.add_argument('-c','--contigs', help ='fasta file with contigs')
        parser.add_argument('-f','--fix-id', dest = "fix_id", action="store_true", help ='remove description from fasta def line')
        self.args, unknown = parser.parse_known_args()
        self.db_name  = self.args.db_name
        self.db_path, self.db_fname  = os.path.split(self.db_name)
        self.conx  = sqlite3.connect(self.db_name)  
        self.contigs  = self.args.contigs
        
        if self.args.fix_id:
            if not self.args.contigs:
                raise argparse.ArgumentTypeError('fix-id requires fasta contig file. See -h/--help')
            path, fname = os.path.split(self.args.contigs)
            fname, ext = os.path.splitext(fname)
            new_fname  = ''.join([fname+'_fx', ext])
            contigs_newfname = os.path.join(path, new_fname)
            self.contigs_newfname_fp  = open(contigs_newfname, 'w')
            tsv = '.'.join([fname, 'tsv'])
            tsv_fname = os.path.join(self.db_path, tsv)
            self.tsv_fp =  open(tsv_fname, 'w')
            
            
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
                contig_data = dict([field.split("=") for field in descr.split(' ')])
                SeqIO.write(contig, self.contigs_newfname_fp,  "fasta")
                self.contigs_newfname_fp.flush(contig_data)
                print("{}\t{}".format(contig.id, descr),file=self.tsv_fp, flush = True) 
            contig_data.update({'id': contig_id , 'description': contig.description, 'length': str(len(contig))})
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


        
if __name__ == '__main__':
    db_connect = DB_Connect()
    db_connect.get_contig()






#         TABLES = {        
#         TABLES['primers'] = (
#              " CREATE TABLE `primers` ("
#              " `gene` CHAR(255) NOT NULL,"
#              " `Fwd_id`  CHAR(30) NOT NULL,"
#              " `Fwd_Primer` CHAR(255) NOT NULL,"
#              " `Rev_id` CHAR(30) NOT NULL,"
#              " `Rev_Primer` CHAR(255) NOT NULL,"
#              " `Amplicon_length` CHAR(255),"
#              " `technology` CHAR(255),"
#              " `Reference`  CHAR(255),"
#              " `Cross_ref`  CHAR(255),"
#              " `notes` LONGTEXT,"
#              " primary key (`Fwd_id`, `Rev_id`)"
#              " ) ENGINE=InnoDB")
        
        
#         TABLES['primerprop'] = (
#              " CREATE TABLE `primerprop` ("
#              " `primer` CHAR(30) PRIMARY KEY,"
#              " `TmProd` LONGTEXT NOT NULL," 
#              " `DeltaS` LONGTEXT NOT NULL," 
#              " `Length` LONGTEXT NOT NULL," 
#              " `GC`     LONGTEXT NOT NULL," 
#              " `DeltaG` LONGTEXT NOT NULL," 
#              " `DeltaH` LONGTEXT NOT NULL," 
#              " `Tm`     LONGTEXT NOT NULL"
#              " ) ENGINE=InnoDB")

#         TABLES['amplicons'] = (
#              " CREATE TABLE `amplicons` ("
#              " `primer_ID` CHAR(255) NOT NULL,"	
#              " `amplimer`  CHAR(255) NOT NULL,"	
#              " `fwd`        INT NOT NULL,"	
#              " `fwd_mis`    INT NOT NULL,"	
#              " `rev`       INT NOT NULL,"	
#              " `rev_mis`   INT NOT NULL,"	
#              " `len`	   INT NOT NULL,"
#              " `seq_id`    VARCHAR(1000) NOT NULL"
#              " ) ENGINE=InnoDB")

        
#         TABLES['taxa_data'] = (
#              " CREATE TABLE `taxa_data` ("
#              " `seq_id` CHAR(30)  PRIMARY KEY NOT NULL,"
#              " `taxonomy` LONGTEXT NOT NULL"
#              " ) ENGINE=InnoDB")

#         self.TABLES = TABLES
#         self.DB_tables = list(TABLES.keys())


#     def create_tables(self):
        
#         try:
#             self.cnx.database = self.DB_NAME  
#         except mysql.connector.Error as err:
#             if err.errno == errorcode.ER_BAD_DB_ERROR:
#                 print('Database Error')
#                 raise Exception(err)
            
#         for table, ddl in list(self.TABLES.items()):
#             try:
#                 self.cursor.execute(ddl)
#             except mysql.connector.Error as err:
                
#                 if err.errno == errorcode.ER_TABLE_EXISTS_ERROR:
#                     print('Table "{}" already exists.'.format(table))
                    
#                 else:
#                     raise Exception(err)

# #!/usr/bin/env python

# from init_primerDB import PrimerDB
# import mysql.connector
# import argparse
# import csv


# class LoadDB(PrimerDB):

#     def __init__(self):
        
#         super(LoadDB, self).__init__()
#         parser = argparse.ArgumentParser(description="""Compile High-Throughput Sequencing (HTS) Primer database using a csv primer file.""",
#         epilog='NOTE: This utility expects database to be created and cvs files to have headers')
#         parser.add_argument('-p','--primers-file', dest='primers_file', action='store', 
#                             required=True, type=str)
#         parser.add_argument('-c','--custom-header', dest='headers', default=False, action='store_true', 
#                             required=False, help="""use cvs headers (default=False). Used when order of column header is different from expected. Column/header are expected in the following order: "Fwd_id, Fwd_Primer, Rev_id, Rev_Primer, Target, gene, Amplicon_length, technology, Reference, Cross_ref, notes". Important! header reading is case sensitive, if provided headers must match the case those provided here""")

#         self.args, unknown = parser.parse_known_args()
#         self.cnx.database = self.DB_NAME  
#         self.headers = self.args.headers     
#         self.csv_reader = csv.reader(open(self.args.primers_file))
#         self.table_cols = ["Fwd_id","Fwd_Primer","Rev_id","Rev_Primer",
#                           "gene","Amplicon_length","technology",
#                           "Reference","Cross_ref","notes"]
        
#     def LoadCSV(self):
        
#         #must skip the first line of csv file
#         top_line = next(self.csv_reader)
#         csv_headers = top_line if self.headers else self.table_cols
        
#         insert_dict = {}
#         for row in self.csv_reader:
#             row_dict = dict(zip(csv_headers, map(lambda w: w.strip(), row)))
#             if self.headers:
#                 for k in row_dict:
#                     if k in self.table_cols:
#                         insert_dict[k] = row_dict[k]
#                     else:
#                         print("column {} removed from data not inserted into DB (use -h/--help for help)".format(k))
#                 if  insert_dict:
#                     row_dict = insert_dict
#                     insert_dict = {}
#             self.DB_Insert(row_dict)
#         self.cnx.commit()
        
#     def DB_Insert(self, row_dict):
        
#         cols = ','.join(row_dict.keys()) 
#         values = ','.join( '"'+val.translate(None,""""'""")+'"'  for val in  row_dict.values())
#         sql = """INSERT INTO primers ({}) VALUES ({})""".format(cols, values)
        
#         try:
#             self.cursor.execute(sql)
#         except mysql.connector.errors.IntegrityError as err:
#             print(err)
#             pass









# LoadDB = LoadDB()
# LoadDB.LoadCSV()
                
    #     primer_db = PrimerDB()
#     primer_db.init_tables()
#     primer_db.create_DB()
#     primer_db.create_tables()

    
