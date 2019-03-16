#! /usr/bin/env python
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
        seq_data = dict([ field.split("=") for field  in  seq.description.split("/")[1:]])
        print("{0}\t{NCGR_PEP_ID}\t{TAXON_ID}\t0".format(seq.id,**seq_data))
        
        
if __name__ ==  '__main__':
    parser = argparse.ArgumentParser("generate taxa map file")
    parser.add_argument('fasta_seqs', help = "MMETSP sequence file")
    parser.add_argument("-o","--outfile", help = "MMETSP sequence file")
    args = parser.parse_args()
    fasta_seqs = args.fasta_seqs
    outfile = args.outfile
    parse_fasta(fasta_seqs, outfile)
    
