#!/usr/bin/env python
"""
This is a simple script to mask a genome based on self-mapping info
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

def mask_genome(genome, self_maps_file):
  
  f_gseq=open(genome,"rU")
  genome_record=list(SeqIO.parse(f_gseq, "fasta"))[0]
  f_gseq.close()
  
  #We need to be able to change the sequence
  genome_seq = genome_record.seq.tomutable()
  
  f_maps=open(self_maps_file,"r")
  #first line is header
  f_maps.readline()
  for line in f_maps:
    line=line.strip()   
    flist=line.split("\t")
    st1=int(flist[1])
    en1=int(flist[2])
    st2=int(flist[5])
    en2=int(flist[6])
    #print flist[1] + "\t"+ flist[2] +"\t"+ flist[5] + "\t"+ flist[6]
    genome_seq[st1:en1] = "N" * (en1-st1)
    genome_seq[st2:en2] = "N" * (en2-st2)

  #Now write the masked genome
  masked_record = SeqRecord(genome_seq,id=genome_record.id)
  output_handle = open(genome+".masked", "w")
  SeqIO.write(masked_record, output_handle, "fasta")
  output_handle.close()
  
  return


instructions = "e.g. python mask_genome.py genome.fasta self_mappings.map"

if __name__=="__main__":
    if len(sys.argv) == 3:
        mask_genome(sys.argv[1],sys.argv[2])
    else:
        print instructions