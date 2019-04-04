#!/usr/bin/env python

import gffutils
import argparse
from collections import  OrderedDict
import pprint
import sys
import os

class ParseGff(object):

    def __init__(self, args):


        self.outfile = open(args.outfile,"w") if args.outfile else sys.stdout
        self.dfn = args.dfn if args.dfn else '.'.join([os.path.splitext(os.path.basename(args.gff))[0],'db'])
        self.gff_db = gffutils.create_db(args.gff, self.dfn, force=True, keep_order=True)  
        self.source_version = self.gff_db.directives[1].split()[2]
        self.genetable_partial = 1
        self.table2gff_map = OrderedDict([
            ("gene_callers_id", "file_order"),
            ("contig", "seqid"),
            ("start", "start"),
            ("stop", "end"),
            ("direction", "strand"),
            ("partial", None),
            ("source", "source"),
            ("version", None)])
        

        
    def gff2table(self):
        
        format2table = "{"+"}\t{".join(self.table2gff_map.keys())+"}"
        header =  "\t".join(self.table2gff_map.keys())
        print(header,file=self.outfile, flush=True)
        for feature in self.gff_db.all_features():
            feature_dict = vars(feature)
            table_dict = OrderedDict()
            for table,gff in self.table2gff_map.items():
                if gff is None:
                    if table == 'partial':
                        table_dict[table] = self.genetable_partial
                        continue
                    elif table == 'version':
                        table_dict[table] = self.source_version 
                        continue
                table_dict[table] = feature_dict[gff]
                
            print(format2table.format(**table_dict), file=self.outfile, flush=True)


            
if  __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="""Parse  GFF file""")
    parser.add_argument('gff')
    parser.add_argument('-d',"--database-fname", dest="dfn", help="filename of sqlite database to be used. Default: gff baseneme with '.db' extension" )
    parser.add_argument('--no-table', dest ='no_table', action='store_true', help="do not generate table. Default:False")
    parser.add_argument('-o','--outfile', help="gene table for Anvi'o")
    args = parser.parse_args()
    gff =  ParseGff(args)
    if not args.no_table:
        gff.gff2table()
    
