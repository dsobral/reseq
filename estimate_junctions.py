#!/usr/bin/env python
"""
This script estimates junctions from a map file
"""
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

  #Returns true if interval is fully covered by parameter
  def is_fully_covered(self, interval):
    if((self.end <= interval.end) and (self.begin >= interval.begin)):
      return True
    else: 
      return False

class GenomicMap:
  """A Class to represent gapped genomic mappings"""
  def __init__(self, interval_a, interval_b, strand):
    self.source = interval_a
    self.destination = interval_b
    self.strand = strand
    
    
def check_junctions(readlist, output):
  junctions = []
  rlen=90 #pass this as parameter
  #Can only have an event per read (readlist are alignment blocks per read)
  max_size=0
  for n in range(0,len(readlist)):	
    read = readlist[n].source
    if(max_size<read.size()): max_size=read.size()
  for a in range(0,len(readlist)-1):
    reada = readlist[a].source
    desta = readlist[a].destination
    # One alignment begins with the first base of the read.
    if(reada.begin==1):
      for b in range(1,len(readlist)):
	#check junction
	readb = readlist[b].source
	destb = readlist[b].destination
	overlap = 0
	intersection = reada.intersect(readb)
	if(intersection != None): overlap = intersection.size()
	# Both alignments contain at least 5 read bases that do not overlap the other.
	reada_uniq = reada.size() - overlap	
	readb_uniq = readb.size() - overlap
	if(reada_uniq<5): break
	if(readb_uniq<5): break	
	# Both alignments together cover a number of bases in the read that is more than 2 bases longer than the length covered by any other single alignment.
	align_cover=reada_uniq+readb_uniq+overlap	
	if(align_cover<(max_size+2)): break		
	# One alignment contains at least 10 read bases that do not overlap the other.
	if((reada_uniq<10) and (readb_uniq<10)): break
	# There are at most 20 bp unique to the read between matches to the reference.
	if((rlen-align_cover)>20): break
	#If it survived all these tests, then declare a "winner"
	junctions.append(Interval(desta.end,destb.begin))
	#output.write(str(junction.begin)+"\t"+str(junction.end)+"\n")
	#print "Junction: "+str(desta.begin)+"-"+str(desta.end)+" : "+str(destb.begin)+"-"+str(destb.end)
  return junctions
  
def estimate_junctions(map_file,outfile,min_cov):
   
    output=open(outfile,"w")
    
    events = {}
    read=''
    begin=0
    end=0
    strand=''
    ref_begin=0
    ref_end=0
    possible=True
    readlist=[]
    
    junction_evidence = []
    
    map=open(map_file,"r")
    #First line is header
    map.readline()
    for nline in map:
      line=nline.strip()   
      flist=line.split("\t")
      
      #custom format
      #new_read=flist[3]
      #new_begin=int(flist[4])
      #new_end=int(flist[5])
      #new_strand=flist[6]
      #new_ref_begin=int(flist[1])
      #new_ref_end=int(flist[2])
      
      new_read=flist[3]
      new_begin=int(flist[4])
      new_end=new_begin+int(flist[2])
      new_strand=flist[5]
      new_ref_begin=int(flist[1])
      new_ref_end=new_ref_begin+int(flist[2])
      
      if(read==new_read):
	readlist.append(GenomicMap(Interval(new_begin,new_end),Interval(new_ref_begin,new_ref_end),new_strand))
	#if parts of the read are overlapping, then there is no junction
	#if((new_begin<=end) and (new_end>=begin)): possible = False
	#It may be just split alignments because of bad quality reads...
	#read_diff=new_begin-new_end
	#ref_diff=new_ref_begin-new_ref_end
	#if(read_diff == ref_diff): 
	#  possible=False
	#
	#if(possible):
	#  output.write(str(begin)+"-"+str(end)+"->"+str(ref_begin)+":"+str(ref_end)+"\t"+str(new_begin)+":"+str(new_end)+"->"+str(new_ref_begin)+":"+str(new_ref_end)+"\t"+new_strand+"\n")
      else:
	
	#check combinations of reads for possible events...
	for junction in check_junctions(readlist,output):
	  junction_evidence.append(junction)
	
	readlist=[]
	readlist.append(GenomicMap(Interval(new_begin,new_end),Interval(new_ref_begin,new_ref_end), new_strand))
	#set variables
	read=new_read
	begin=new_begin
	end=new_end
	strand=new_strand
	ref_begin=new_ref_begin
	ref_end=new_ref_end
	possible=True
    
    junction_counts={}
    #Then we just need to check the last one    
    for junction in junction_evidence:
      #Allow some variation around the point...
      #start = junction.begin / 10
      #end = junction.end / 10
      start = junction.begin
      end = junction.end
      if start not in junction_counts:
	junction_counts[start] = {}
	junction_counts[start][end] = 1
      elif end not in junction_counts[start]:
	junction_counts[start][end] = 1
      else:
	junction_counts[start][end] = junction_counts[start][end] + 1      
    
    for start, end_counts in sorted(junction_counts.iteritems()):
      for end, count in end_counts.iteritems():
	if(count>10):
	  output.write(str(start)+"<->"+str(end)+"\t"+str(count)+"\n")
    
    output.close()
    return
  

if __name__=="__main__":
    if len(sys.argv) == 3:
        estimate_junctions(sys.argv[1], sys.argv[2], 30)
    else:
        print "estimate_junctions.py map_file outfile"