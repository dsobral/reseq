#!/usr/bin/env python
"""
This is a simple script to extract a sequence from a fasta database
"""
from Bio import SeqIO
import sys

def add_coordinates(fasta_file, bedfile, outfile):
  
  starts = {}
  f_gseq=open(fasta_file,"rU")
  for record in SeqIO.parse(f_gseq, "fasta"):
    #for backwards compatibility (apparently), fasta descriptors contain the id!!
    start=int(record.description.split("-")[0].replace(record.id,""))
    starts[record.id]=start-1
  f_gseq.close()

  output_handle = open(outfile, "w")
  bed_handle=open(bedfile,"r")
  for bedline in bed_handle:
      list=bedline.split()
      output_handle.write(list[0]+"\t"+str(int(list[1])+starts[list[0]])+"\t"+str(int(list[2])+starts[list[0]])+"\n")
  bed_handle.close()
  output_handle.close()

  return


instructions = "e.g. python add_coordinates_in_bed.py fastafile bedfile outfile"

if __name__=="__main__":
    if len(sys.argv) == 4:
        add_coordinates(sys.argv[1],sys.argv[2],sys.argv[3])
    else:
        print instructions