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

class GenomicMap:
  """A Class to represent gapped genomic mappings"""
  def __init__(self, interval_a, interval_b, strand):
    self.source = interval_a
    self.destination = interval_b
    self.strand = strand

class GenomicEvent:
  """A Class to represent an event"""
  def __init__(self, type, begin, length):
    self.type = type
    self.begin = begin
    self.length = length

def infer_genomic_event(gen_map1,gen_map2):
  
  #Ignore if source is not in the same strand...
  if(gen_map1.strand != gen_map2.strand): return None
  
  #Determines genomic event based on "contiguous" sets of 
  # mappings between denovo assembly and reference genome
  src_diff = gen_map2.source.begin - gen_map1.source.end
  dst_diff = gen_map2.destination.begin - gen_map1.destination.end
    
  if(src_diff == dst_diff):
    return GenomicEvent("variation", gen_map1.destination.end+1, src_diff-1)
  elif (src_diff > dst_diff):
    return GenomicEvent("insertion", gen_map1.destination.end+1, src_diff-dst_diff)
  else: 
    return GenomicEvent("deletion", gen_map1.destination.end+1, dst_diff-src_diff)

def infer_genomic_events(reseq_contigs_file,mappings_file,self_genome_mappings, outfile):
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
  
  covered_areas = []
  contig_maps = {}
  f_maps=open(mappings_file,"r")
  #first line is header
  f_maps.readline()
  for line in f_maps:
    line=line.strip()   
    flist=line.split("\t")
    contig_id=flist[0]
    contig_strand=flist[3]
    contig_interval = Interval(int(flist[1]),int(flist[2]))
    reference_interval = Interval(int(flist[5]),int(flist[6]))        
    
    #Ignore mappings that completely fall within "known" repeats
    #or ignore those that are mostly repeats... ? more than % repeats?
    repeated = False
    for ref_repeat in potential_repeats:
      if(ref_repeat.fully_covers(reference_interval)):
      	repeated=True
      	break
    #  intersect = ref_repeat.intersect(reference_interval)
    #  if(intersect is not None):
    #	if((float(intersect.size()) / float(reference_interval.size())) > 0.75):
    #	  repeated=True
    #	  break
      
    if(not repeated): 
    #	print contig_id+"\t"+str(contig_interval.begin)+"\t"+str(contig_interval.end)+"\t"+str(reference_interval.begin)+"\t"+str(reference_interval.end)
      try:
	recon=copy.copy(contig_maps[contig_id])
	overlapped = False
	for gmap in recon:
	  if(contig_interval.is_fully_covered(gmap.source)):
	    overlapped=True
	if(not overlapped):
	  for gmap in recon:
	    if(contig_interval.fully_covers(gmap.source)):
	      contig_maps[contig_id].remove(gmap)
	  contig_maps[contig_id].append(GenomicMap(contig_interval,reference_interval,contig_strand))	
      except KeyError:
	contig_maps[contig_id] = []
	contig_maps[contig_id].append(GenomicMap(contig_interval,reference_interval,contig_strand))	
  
  f_maps.close()
    
  events = []
  #Infer events based on assembly...
  for contig, cont_maps in sorted(contig_maps.iteritems()):
    sorted_mappings = sorted(cont_maps, key=lambda x: x.source.begin)
    for i in range(len(sorted_mappings)-1):
      new_event = infer_genomic_event(sorted_mappings[i],sorted_mappings[i+1])
      if(new_event != None): events.append(new_event)
  
  output=open(outfile,"w")
  #Use mappings to eliminate events that cannot happen because there is evidence against it...
  for event in sorted(events, key=lambda x: x.begin):
    output.write(str(event.begin)+"-"+str(event.begin+event.length)+"\t"+event.type+"\t("+str(event.length)+"bp)\n")
    event_interval=Interval(event.begin,event.begin+event.length-1)
    credible=True
    for contig, cont_maps in contig_maps.iteritems():
      for map in cont_maps:
    	#if(map.destination.fully_covers(event_interval)):
    	if(map.destination.intersect(event_interval) is not None):	  
    	  output.write("\t region "+str(map.destination.begin)+"-"+str(map.destination.end)+" is covered by "+contig+"\n")
    	  credible=False
    if(not credible): 
      output.write("\t Caution: event may be false\n")
  output.close()
  return

#Tests... the script will not continue if these fail
assert(Interval(10,20).merge(Interval(30,15)).begin==10)
assert(Interval(10,20).intersect(Interval(30,15)).begin==15)
assert(Interval(10,20).is_fully_covered(Interval(5,30)))

instructions = "e.g. python estimate_genome_events.py contigs.fasta mappings.map genome_self_mappings.map outfile.tab"

if __name__=="__main__":
    if len(sys.argv) == 5:
        infer_genomic_events(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    else:
        print instructions