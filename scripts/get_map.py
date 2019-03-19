#! /usr/bin/env python
from collections import defaultdict
from Bio import SeqIO
import argparse
import pprint
import sys




def parse_fasta(fasta_seqs, outfile):

    out = None
    seq_records = SeqIO.parse(fasta_seqs,"fasta")    
    out = open(outfile, "w") if outfile else sys.stdout
    
    print("accession\taccession.version\ttaxid\tgi",file=out, flush=True)
    for seq in seq_records:
        try:
            seq_data = defaultdict(lambda : "")
            for field  in  seq.description.split(" /")[1:]:
                k,v = field.strip().split("=")[:2]
                seq_data[k] = clean(v)
        except Exception as e:
            print(seq.id,seq.description)
            print(e)
            continue 
        print("{0}\t{1}\t{2}\t0".format(seq.id,seq_data['DNA_ID'],seq_data['TAXON_ID']),file=out, flush=True)

        
def clean(value):
    #extensible 
    #add stray strings to list
    
    for artifact in ["ORGANISM"]:
      value = value.replace(artifact,'')
    
    return value.strip()
        
if __name__ ==  '__main__':
    parser = argparse.ArgumentParser("generate taxa map file")
    parser.add_argument('fasta_seqs_file', help = "MMETSP sequence file")
    parser.add_argument("-o","--outfile", help = "MMETSP sequence file")
    args = parser.parse_args()
    fasta_seqs = args.fasta_seqs_file
    outfile = args.outfile
    parse_fasta(fasta_seqs, outfile)

#some data fail the parsing     
#sed -i -e  's|Strain 621\/1 \/|Strain621\/1\/|g'  mmetsp_pep/MMETSP0151.pep.fa 





























