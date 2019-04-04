#! /usr/bin/env python

from Bio.Blast import  NCBIXML
from collections import Counter, OrderedDict
from db_connect import DB_Connect 
import itertools
import argparse
import pprint
import sys



class  Bunker(DB_Connect):

    """

      Parse and iterater through Blast results 
     
   """
   
    
    def __init__(self, blast_xml, cut_off, verbose, assign_on ='TAXON_ID'):
        
       super().__init__()
       self.Blast_XML_Output = NCBIXML.parse(open(blast_xml))
       self.assign_on = assign_on
       self.cut_off = cut_off
       assert (1 >= self.cut_off and self.cut_off >=  0.5), "0.5 =< cut_off =< 1"  
       self.verbose = verbose
       self.hsps_remove = [ 'match', 'query', 'sbjct', 'strand', 'num_alignments']
       self.hit_def_remove = ['ASSEMBLY_ACC', 'assembly_acc','NCGR_SAMPLE_ID']
       
    def blast_xml(self):
        
       print_fields = OrderedDict.fromkeys(['query',
                                            'hit_id',
                                            'perc_id',
                                            'align_length',
                                            'mismatches',
                                            'gaps',
                                            'query_end',
                                            'query_start',
                                            'sbjct_start',
                                            'sbjct_end',
                                            'expect',
                                            'bits',
                                            'outcome',
                                            'organism'])
       header  = "\t".join(print_fields.keys())
       print_format =  "{"+"}\t{".join(print_fields.keys())+"}"
       if self.verbose:
           print(header)
            
       for  i  in itertools.count():
           try:
               blast_record =  next(self.Blast_XML_Output)
               query = blast_record.query
               align_ref = {'contig_id': query }
               assign_data  = {}
               for alignment in blast_record.alignments:
                   #get alignment data
                   align_ref[alignment.hit_id] = vars(alignment)
                   #get_hitdef, implemented MMETSP data
                   hit_def = self.get_hitdef(alignment)
                   #only one HSP!
                   align_ref[alignment.hit_id]['hit_def'] = hit_def
                   align_ref[alignment.hit_id]['hsps'] = vars(alignment.hsps[0])
                   identities = align_ref[alignment.hit_id]['hsps']['identities']
                   align_length = align_ref[alignment.hit_id]['hsps']['align_length']
                   align_ref[alignment.hit_id]['hsps'].update({'perc_id' : round(100*identities/float(align_length), 2),
                                                            'mismatches' : align_length - identities})
                   assign_data[alignment.hit_id] = hit_def[self.assign_on]
                   if self.verbose:
                       #updating printing fields
                        print_fields.update(align_ref[alignment.hit_id]['hsps'])
                        print_fields.update(align_ref[alignment.hit_id])
                        #replace query sequence with query id
                        print_fields['query'] = query.split()[0]
               outcome, hit_ref, taxa_id = self.get_consesus(assign_data)
               print_fields['outcome'] = outcome
               if (outcome == "consesus"):
                   align_ref['outcome']  = outcome
                   align_ref["consesus"] = hit_ref
                   self.loadDB(hit_ref, align_ref)
                   organism, taxon_id = self.get_organism(align_ref, hit_ref, taxa_id)
                   if self.verbose:
                       print_fields['organism'] = organism
                       print(print_format.format(**print_fields))
               else:
                   
                   for hit_ref, taxa_id  in assign_data.items():
                       align_ref['outcome']  = outcome
                       align_ref["consesus"] = hit_ref
                       organism, taxon_id = self.get_organism(align_ref, hit_ref, taxa_id)
                       self.loadDB(hit_ref, align_ref)
                       if self.verbose:
                           print_fields['organism'] = organism
                           print(print_format.format(**print_fields))
           except StopIteration:
               break    
           self.conx.commit()


           
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
                  hit_ref, tax_id = [ (hit, tax_id)  for hit,tax_id in assign_data.items()\
                             if tax_id == field ][0]
                  return ("consesus", hit_ref, tax_id)   
              
         return ("mutiple hits", None, None)

     
    def get_organism(self, align_ref, hit_id, taxon_ref):

        organism = align_ref[hit_id]['hit_def'].get('ORGANISM', taxon_ref)            
        taxon_id = align_ref[hit_id]['hit_def']['TAXON_ID']
        assert taxon_id == taxon_ref,"hit_def error taxon ids not matching!"

        return organism, taxon_id

    def loadDB(self, hit_ref, align_ref):
         contig_id = align_ref['contig_id']
         accession = align_ref[hit_ref]['accession']
         taxa_sql = """INSERT INTO taxa_assignment(
                              contig_id, 
                              outcome, 
                              hit_id) 
                              VALUES 
                              ('{contig_id}',
                                '{outcome}',
                                '{consesus}')""".format(**align_ref)
         
         self.conx.execute(taxa_sql)
         hit_def = align_ref[hit_ref]['hit_def']
         
         #update hsp keys and values for database insert
         hsps = align_ref[hit_ref]['hsps']
         
         _ = [ hsps.pop(k, None) for k in self.hsps_remove ]
         hsps.update({'contig_id' : align_ref['contig_id'], 'hit_id' : accession })
         #update hit_def key and values for database insert
         _ = [ hit_def.pop(k, None) for k in self.hit_def_remove ]
         hit_def.update({'accession': accession }) 
         table_entries = {'hit_def': hit_def, 'hsps': hsps}
         for table,entry in table_entries.items():
             cols = ', '.join(entry.keys())
             values = ', '.join([ "'"+str(val)+"'" for val in entry.values() ])
             sql = "INSERT INTO {0} ({1}) VALUES ({2})".format(table, cols, values )
             pprint.pprint(sql)
             self.conx.execute(sql)
         
         
     
if  __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Parse  BLAST xml output""")
    parser.add_argument('blast_xml')
    parser.add_argument('-c','--consesus-cut-off', dest='cut_off', required=False, default=0.5,type=float)
    parser.add_argument('-b','--batchsize', dest='batchsize', action='store', required=False, default=5,type=int)
    parser.add_argument('-v','--verbose', action="store_true")
    args, unknown = parser.parse_known_args()
    bunker = Bunker(args.blast_xml, args.cut_off, args.verbose)
    bunker.blast_xml()
