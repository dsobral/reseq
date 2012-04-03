#!/usr/bin/env python
"""
This is a simple script to estimate large indels 
based on resequencing data comparing to a reference
"""
from Bio import SeqIO
import sys, copy

class Interval:
  """A simple interval class"""
  def __init__(self, begin, end):
    if(begin > end):
      swap=begin
      begin=end
      end=swap
    self.begin=begin
    self.end=end
   
  #Merges two intervals: If interval is not mergeable, then return None
  def merge(self, interval):
    if((self.end < interval.begin) or (self.begin> interval.end)):
      return None
    else: 
      return Interval(min(self.begin,interval.begin), max(self.end, interval.end))

  #Intersect two intervals: If interval is not intersectable, then return None
  def intersect(self, interval):
    if((self.end < interval.begin) or (self.begin> interval.end)):
      return None
    else: 
      return Interval(max(self.begin,interval.begin), min(self.end, interval.end))

  #Returns true if interval is fully covered by parameter
  def is_fully_covered(self, interval):
    if((self.end <= interval.end) and (self.begin >= interval.begin)):
      return True
    else: 
      return False

  #Returns true if interval fully covers parameter
  def fully_covers(self, interval):
    return interval.is_fully_covered(self)


class GenomicMap:
  """A Class to represent genomic mappings"""
  def __init__(self, interval_a, interval_b):
    self.mappings = {}
    self.mappings[interval_a] = interval_b

  def addMapping(self, interval_a, interval_b):
    recon = copy.copy(self.mappings)
    for inter_a, inter_b in recon.iteritems():
      #If new mapping is covered by existing intervals, then ignore new mapping      
      if(interval_a.is_fully_covered(inter_a)):
	return self
      #If new mapping completly covers existing intervals then remove them
      if(interval_a.fully_covers(inter_a)):
	del self.mappings[inter_a]
    self.mappings[interval_a]=interval_b
    return self

def infer_genomic_event(src1,dst1,src2,dst2):
  #Determines genomic event based on "contiguous" sets of 
  # mappings between denovo assembly and reference genome

def infer_genomic_events(reseq_contigs_file,mappings_file):
  #Gets the resequencing contigs and reads mappings to reference
  sizes = {}
  f_reseq=open(reseq_contigs_file,"rU")
  for record in SeqIO.parse(f_reseq, "fasta") :
    sizes[record.id] = len(record.seq)
  f_reseq.close()
  
  contig_maps = {}
  f_maps=open(mappings_file,"r")
  #first line is header
  f_maps.readline()
  for line in f_maps:
    line=line.strip()   
    flist=line.split("\t")
    contig_id=flist[0]
    contig_interval = Interval(int(flist[1]),int(flist[2]))
    reference_interval = Interval(int(flist[5]),int(flist[6]))
    try:
      contig_maps[contig_id].addMapping(contig_interval,reference_interval)
    except KeyError:
      contig_maps[contig_id] = GenomicMap(contig_interval,reference_interval)	
  f_maps.close()
  
  #Infer events based on assembly...
  for contig, cont_maps in sorted(contig_maps.iteritems()):
    print contig+": "+str(sizes[contig])
    sorted_mappings = sorted(cont_maps.mappings.iterkeys(), key=lambda x: x.begin)
    for in range()
    #for inter_a, inter_b in sorted(cont_maps.mappings.iteritems(), key=lambda x: x[0].begin):
    #  print "\t"+str(inter_a.begin)+"-"+str(inter_a.end)+"->"+str(inter_b.begin)+"-"+str(inter_b.end)
  return

#Tests... the script will not continue if these fail
assert(Interval(10,20).merge(Interval(30,15)).begin==10)
assert(Interval(10,20).intersect(Interval(30,15)).begin==15)
assert(Interval(10,20).is_fully_covered(Interval(5,30)))

instructions = "e.g. python estimate_large_deletions.py contigs.fasta mappings.map"

if __name__=="__main__":
    if len(sys.argv) == 3:
        estimate_genomic_events(sys.argv[1],sys.argv[2])
    else:
        print instructions