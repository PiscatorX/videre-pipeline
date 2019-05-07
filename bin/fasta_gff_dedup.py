#!/usr/bin/env python
from Bio import SeqIO
import collections
import argparse
import pprint



class FastaGFF(object):
    
    def __init__(self):
        
        parser = argparse.ArgumentParser(description="""Dedupe gff and fasta""")
        parser.add_argument('fasta')
        parser.add_argument('gff')
        parser.add_argument('-f','--fasta_out', required = True)
        parser.add_argument('-g','--gff_out',  required = True)
        self.args, unknown = parser.parse_known_args()
        self.fasta_fp  = open(self.args.fasta)
        self.gff_fp    = open(self.args.gff)
        self.fasta_out = open(self.args.fasta_out,'w')
        self.gff_out   = open(self.args.gff_out,'w')
        self.duplicates = collections.defaultdict(list)

        
    def dedup(self):
         fasta_records = SeqIO.parse(self.fasta_fp, "fasta")
         id_counts =  collections.defaultdict(int)
         fid_counts =  collections.defaultdict(int)
         dedup_i = 0
         for seq  in fasta_records:
             id_counts[seq.id]+=1
             seq_id = seq.id+"_"+str(id_counts[seq.id]) if id_counts[seq.id] != 1 else seq.id
             if id_counts[seq.id] != 1:
                 self.duplicates[seq.id].append(seq_id)
                 print("{}\tfasta duplicate: {}\t{}".format(dedup_i, seq.id, seq_id))
                 dedup_i += 1
                 seq.id = seq_id
             SeqIO.write(seq, self.fasta_out, "fasta")
             self.fasta_out.flush()
             
         print()
         dedup_j = 0   
         for track in self.gff_fp:
            track_data  = track.strip()
            if not track_data or track.startswith("#"):
                print(track, file=self.gff_out, end="")
                continue
            feature_id = track_data.split()[0]
            id_counts[seq.id]+=1 
            if self.duplicates.get(feature_id, False):
                fid_counts[feature_id] += 1
                if fid_counts[feature_id] != 1:
                    deduped_id = self.duplicates[feature_id].pop(0)
                    track = track.replace(feature_id,deduped_id)
                    print(track, file=self.gff_out, end="")
                    print("{}\tgff duplicate: {}\t{}".format(dedup_j,feature_id, deduped_id))
                    dedup_j += 1
                    self.gff_out.flush()
                    continue
                print(track, file=self.gff_out, end="")
        
         print("\nDedup:\n{} fasta\n{} gff\n".format(dedup_i, dedup_j))   
            
if __name__ == '__main__':
    fasta2dedup = FastaGFF()
    fasta2dedup.dedup()
