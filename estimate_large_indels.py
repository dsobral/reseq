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

  def size(self):
    return self.end - self.begin

  #Subtract two intervals: If interval is not subtractable, then return both
  #Not done yet...
  #def intersect(self, interval):
  #  if((self.end < interval.begin) or (self.begin> interval.end)):
  #    return [self, interval]
  #  else: 
  #  if(self.is_fully_covered(interval):
  #  elif(self.fully_covers(interval):
  #  else:
  #    return [ Interval(max(self.begin,interval.begin), min(self.end, interval.end)) ]

  #Returns true if interval is fully covered by parameter
  def is_fully_covered(self, interval):
    if((self.end <= interval.end) and (self.begin >= interval.begin)):
      return True
    else: 
      return False

  #Returns true if interval fully covers parameter
  def fully_covers(self, interval):
    return interval.is_fully_covered(self)

class GenomicAlignment:
  """A Class to represent gapped genomic mappings"""
  def __init__(self, interval_a, interval_b, genome_repeats):
    self.mappings = {}
    self.mappings[interval_a] = interval_b
    self.genome_repeats = genome_repeats
  
  def addMapping(self, interval_a, interval_b):
    recon = copy.copy(self.mappings)
    
    #does not seem to work...
    #Ignore if genome portion overlaps a repeat area...
    #for repeat in self.genome_repeats:
    #  if(interval_b.intersect(repeat) != None):
    #	return self
      
    for inter_a, inter_b in recon.iteritems():
      #If new mapping is covered by existing intervals, then ignore new mapping (probably short duplication)     
      if(interval_a.is_fully_covered(inter_a)):
	return self
      #If new mapping completly covers existing intervals then remove them (probably short duplication(s))
      elif(interval_a.fully_covers(inter_a)):
      	del self.mappings[inter_a]
      #Now we need to check for repeats inside the regions...
      #elif(interval_a.intersect(inter_a) != None):
      #	#ignore the smaller portion, as it is probably a repeat..?
      #	if(interval_a.size < inter_a.size):
      #	  return self
      #	else:
      #	  del self.mappings[inter_a]
	
    self.mappings[interval_a]=interval_b
    return self

def infer_genomic_event(src1,dst1,src2,dst2, potential_repeats):
  #Determines genomic event based on "contiguous" sets of 
  # mappings between denovo assembly and reference genome
  
  #ignore mappings to places that overlap repeating areas...?
  #missing_block=?
  
  src_diff = src2.begin - src1.end
  dst_diff = dst2.begin - dst1.end
  if(src_diff == dst_diff):
    print str(dst1.end+1)+"\tvariation\t"+str(src_diff-1)+"bp"
  elif (src_diff > dst_diff):
    print str(dst1.end+1)+"\tinsertion\t"+str(src_diff-dst_diff)+"bp"
  else: 
    print str(dst1.end+1)+"\tdeletion\t"+str(dst_diff-src_diff)+"bp"
  #print str(src1.begin)+":"+str(src1.end)+" > "+str(dst1.begin)+":"+str(dst1.end)+" to "+str(src2.begin)+":"+str(src2.end)+" > "+str(dst2.begin)+":"+str(dst2.end)

def infer_genomic_events(reseq_contigs_file,mappings_file,self_genome_mappings):
  #Gets the resequencing contigs and reads mappings to reference
  sizes = {}
  f_reseq=open(reseq_contigs_file,"rU")
  for record in SeqIO.parse(f_reseq, "fasta") :
    sizes[record.id] = len(record.seq)
  f_reseq.close()
  
  #Make a more generic function...
  potential_repeats = []
  f_maps=open(self_genome_mappings,"r")
  f_maps.readline()
  for line in f_maps:
    line=line.strip()   
    flist=line.split("\t")
    potential_repeats.append(Interval(int(flist[1]),int(flist[2])))
    potential_repeats.append(Interval(int(flist[5]),int(flist[6])))
  f_maps.close()
    
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
      #I'm not sure if what I'm doing here is correct...
      contig_maps[contig_id].addMapping(contig_interval,reference_interval)
    except KeyError:
      contig_maps[contig_id] = GenomicAlignment(contig_interval,reference_interval, potential_repeats)	
  f_maps.close()
  
  #Infer events based on assembly...
  for contig, cont_maps in sorted(contig_maps.iteritems()):
    sorted_mappings = sorted(cont_maps.mappings.iterkeys(), key=lambda x: x.begin)
    for i in range(len(sorted_mappings)-1):
      infer_genomic_event(sorted_mappings[i],cont_maps.mappings[sorted_mappings[i]],sorted_mappings[i+1],cont_maps.mappings[sorted_mappings[i+1]], potential_repeats)
    #for inter_a, inter_b in sorted(cont_maps.mappings.iteritems(), key=lambda x: x[0].begin):
    #  print "\t"+str(inter_a.begin)+"-"+str(inter_a.end)+"->"+str(inter_b.begin)+"-"+str(inter_b.end)
  return

#Tests... the script will not continue if these fail
assert(Interval(10,20).merge(Interval(30,15)).begin==10)
assert(Interval(10,20).intersect(Interval(30,15)).begin==15)
assert(Interval(10,20).is_fully_covered(Interval(5,30)))

instructions = "e.g. python estimate_large_deletions.py contigs.fasta mappings.map genome_self_mappings.map"

if __name__=="__main__":
    if len(sys.argv) == 4:
        infer_genomic_events(sys.argv[1],sys.argv[2],sys.argv[3])
    else:
        print instructions