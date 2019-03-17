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
        seq_data = dict([ field.strip().split("=") for field  in  seq.description.split("/")[1:]])
        print("{0}\t{DNA_ID}\t{TAXON_ID}\t0".format(seq.id,**seq_data),file=out, flush=True)

        
if __name__ ==  '__main__':
    parser = argparse.ArgumentParser("generate taxa map file")
    parser.add_argument('fasta_seqs_file', help = "MMETSP sequence file")
    parser.add_argument("-o","--outfile", help = "MMETSP sequence file")
    args = parser.parse_args()
    fasta_seqs = args.fasta_seqs_file
    outfile = args.outfile
    parse_fasta(fasta_seqs, outfile)
    
