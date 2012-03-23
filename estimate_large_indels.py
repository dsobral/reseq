#!/usr/bin/env python
"""
This is a simple script to estimate large indels 
based on resequencing data comparing to a reference
"""
from Bio import SeqIO
import sys

def estimate_large_indels(reseq_contigs,mappings):
  sizes = {}
  f_reseq=open(reseq_contigs,"rU")
  for record in SeqIO.parse(f_reseq, "fasta") :
    sizes[record.id] = len(record.seq)
  f_reseq.close()
  f_maps=open(mappings,"r")
  for line in f_maps:
    print line,
  f_maps.close()
  return

instructions = "e.g. python estimate_large_deletions.py contigs.fasta mappings.map"

if __name__=="__main__":
    if len(sys.argv) == 3:
        estimate_large_indels(sys.argv[1],sys.argv[2])
    else:
        print instructions