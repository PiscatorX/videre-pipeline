#! /usr/bin/env python

from Bio.Blast import  NCBIXML
import argparse


class  Bunker(self):

   """
      Parse and iterater through Blast results 
     
     
   """
   def __init__(self,blast_xml):
       
       self.blast_xml =  blast_xml
       print()












if  __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Parse  BLAST xml output""")
    group.add_argument('blast_xml')
    parser.add_argument('-m','--retmax', dest='retmax', action='store', required=False, default=15,type=int)
    parser.add_argument('-b','--batchsize', dest='batchsize', action='store', required=False, default=5,type=int)
    args = parser.parse_args()
    bunker = Bunker(args.blast_xml)
    








