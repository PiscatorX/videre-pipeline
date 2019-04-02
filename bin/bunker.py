#! /usr/bin/env python

from Bio.Blast import  NCBIXML
from collections import Counter, OrderedDict
import itertools
import argparse
import pprint
import sys



class  Bunker(object):

    """

      Parse and iterater through Blast results 
     
     
    """
    def __init__(self, blast_xml, cut_off, assign_on ='TAXON_ID'):
       self.Blast_XML_Output = NCBIXML.parse(open(blast_xml))
       self.assign_on = assign_on
       self.cut_off = cut_off
       assert (1 >= self.cut_off and self.cut_off >=  0.5), "0.5 =< cut_off =< 1"  

       
    def blast_xml(self):
               
       for  i  in itertools.count():

           fields = ['query',
                     'hit_id',
                     'hsps'
                        'identities',
                        'align_length',
                     
           ]
           try:
               blast_record =  next(self.Blast_XML_Output)
               query = blast_record.query
               align_ref = {}
               assign_data  = {}
               pprint.pprint(vars(blast_record))
               for alignment in blast_record.alignments:
                   
                   align_ref[alignment.hit_id] = vars(alignment)
                   hit_def = self.get_hitdef(alignment)
                   assign_data[alignment.hit_id] = hit_def[self.assign_on]
                   #only one HSP!
                   align_ref[alignment.hit_id]['hit_def'] = hit_def
                   align_ref[alignment.hit_id]['hsps'] = vars(alignment.hsps[0])
                   pprint.pprint(align_ref)
                   exit(1)
               pprint.pprint(assign_data)
               outcome, field  = self.get_consesus(assign_data)
               print(outcome, field)
               exit(1)
#Query ID, Subject ID, Percentage of identical matches, Alignment length, Number of mismatches, Number of gap openings, Start of alignment in query, End of alignment in query, Start of alignment in subject, End of alignment in subject, Expected value, Bit score
               if (outcome == "consesus"):
                   pass
                  #print(outcome)
                   
               #print("{}\t{}\t{}".format(i,blast_record.query, outcome)
           except StopIteration:
               break    



           
    def get_hitdef(self, alignment):
        
        descr_data = {}
        try:
            for field  in  alignment.hit_def.split(" /")[1:]:
                k,v = field.strip().split("=")[:2]
                descr_data[k] = self.clean(v)
        except Exception as e:
            print(e)
            
        return descr_data
                        
            
    def clean(self, value):
        #extensible 
        #remove stray strings to list
        for artifact in ["ORGANISM"]:
            value = value.replace(artifact,'')
        return value.strip()


    
    def get_consesus(self, assign_data):
        
         n = len(assign_data)
         field_counts =  OrderedDict(sorted(Counter(assign_data.values()).items(),
                                key=lambda items: (items[1], items[0]), reverse=True))
         
         for field, count in field_counts.items():
              if (count/float(n) >  self.cut_off):
                  return ("consesus", field)
              
         return ("mutiple_hits", None)
         

     
if  __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Parse  BLAST xml output""")
    parser.add_argument('blast_xml')
    parser.add_argument('-c','--consesus-cut-off', dest='cut_off', required=False, default=0.5,type=float)
    
    parser.add_argument('-b','--batchsize', dest='batchsize', action='store', required=False, default=5,type=int)
    args = parser.parse_args()
    bunker = Bunker(args.blast_xml, args.cut_off)
    bunker.blast_xml()

#assign_taxid is scripted for  MMETSP defline
#/NCGR_PEP_ID=MMETSP1392-20130828|78981_1 /ASSEMBLY_ACC=CAM_ASM_000868 /TAXON_ID=225041 /ORGANISM="Chlamydomonas chlamydogama, Strain SAG 11-48b" /LENGTH=446 /DNA_ID=CAMNT_0049618179 /DNA_START=1 /DNA_END=1339 / DNA_ORIENTATION=+

# {'DNA_END': '662',
#   'DNA_ID': 'CAMNT_0041872075',
#   'DNA_ORIENTATION': '+',
#   'DNA_START': '60',
#   'LENGTH': '202',
#   'ORGANISM': '"Emiliania huxleyi, Strain 374"',
#   'TAXON_ID': '2903',
#   'align_length': 48,
#   'bits': 105.5,
#   'expect': 9.9e-25,
#   'frame': (-2, 0),
#   'gaps': 0,
#   'identities': 47,
#   'match': 'PIVSIANIFEG CFVMNTCSHMTMGCVSTFWQSCGFKLNVCSTSVSPP',
#   'num_alignments': None,
#   'positives': 47,
#   'query': 'PIVSIANIFEGRCFVMNTCSHMTMGCVSTFWQSCGFKLNVCSTSVSPP',
#   'query_end': 322,
#   'query_start': 178,
#   'sbjct': 'PIVSIANIFEGSCFVMNTCSHMTMGCVSTFWQSCGFKLNVCSTSVSPP',
#   'sbjct_end': 48,
#   'sbjct_start': 1,
#   'score': 262.0,
#   'strand': (None, None)}

# 'score',
# 'bits',
# 'expect',
# 'num_alignments',
# 'identities',
# 'positives',
# 'gaps',
# 'align_length',
# 'strand',
# 'frame',
# 'query',
# 'query_start',
# 'query_end',
# 'match',
# 'sbjct',
# 'sbjct_start',
# 'sbjct_end'
