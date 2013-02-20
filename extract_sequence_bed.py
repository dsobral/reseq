#!/usr/bin/env python
"""
This is a simple script to extract a sequence from a fasta database
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

class Interval:
  """A simple interval class"""
  def __init__(self, begin, end):
    if(begin > end):
      swap=begin
      begin=end
      end=swap
    self.begin=begin
    self.end=end

def extract_sequence_bed(fasta_file, outfile, bedfile):
  
  loci = {}
  fbed = open(bedfile,"r")
  for line in fbed:
    line=line.strip()
    fields = line.split("\t")
    region=fields[0]
    start=int(fields[1])
    end=int(fields[2])
    if(region in loci):
      loci[region].append(Interval(start,end))
    else :
      loci[region] = [Interval(start,end)]
  fbed.close()

  output_handle = open(outfile, "w")
  f_gseq=open(fasta_file,"rU")
  for record in SeqIO.parse(f_gseq, "fasta"):
    if(record.id in loci):    
      #Now write the masked genome
      for interval in loci[record.id]:
	seq_r = SeqRecord(record[interval.begin+1:interval.end].seq,id=record.id,description=str(interval.begin)+"-"+str(interval.end))
	SeqIO.write(seq_r, output_handle, "fasta")
	#SeqIO.write(record[interval.begin+1:interval.end], output_handle, "fasta")

  output_handle.close()
  f_gseq.close()
   
  return


instructions = "e.g. python extract_sequence_bed.py fastafile outputfile bedfile"

if __name__=="__main__":
    if len(sys.argv) == 4:
        extract_sequence_bed(sys.argv[1],sys.argv[2],sys.argv[3])
    else:
        print instructions