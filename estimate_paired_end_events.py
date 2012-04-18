#!/usr/bin/env python
"""
This script estimates duplications from a Bam file
"""
import pysam
import sys
from Bio import SeqIO

def estimate_events(reference_file,bam_file,bin,min,max):
    samfile = pysam.Samfile(bam_file, "rb")
    
    reference = SeqIO.parse(reference_file, "fasta").next()
    reflen=len(reference.seq)
    refid=reference.id
    
    totals = {}
    events_ins = {}
    events_del = {}
    for n in range(0, reflen, bin):
	events_ins[n]=0
	events_del[n]=0	
	totals[n]=0
	
    #Binning (a little more than read length)
    for n in range(0, reflen, bin):
      for aligned_read in samfile.fetch(samfile.references[0],n,n+bin): #0-based
	num_events=0
	#just in case!
	if(aligned_read.is_unmapped): continue
	if aligned_read.is_paired:
	  #If we're doing pair number 2 then stop because we did it already!
	  if(aligned_read.is_read2): continue
	  #Only check if the mate is mapped
	  if(aligned_read.mate_is_unmapped):continue
	  totals[n] = totals[n] + 1
	  dist = abs(aligned_read.pnext - aligned_read.pos)+90
	  if(dist<min):
	    events_ins[n] = events_ins[n] + 1
	    #print str(aligned_read.pnext)+"\tIns/Dup\t"+str(min-dist)
	  if((abs(aligned_read.pnext - aligned_read.pos)+90)>max):
	    events_del[n] = events_del[n] + 1
	    #print str(aligned_read.pnext)+"\tDel/Jun\t"+str(dist-max)
	  #print str(aligned_read.pos)+"\t"+str(aligned_read.pnext)

    for n in range(0, reflen, bin):
      if(events_ins[n]>15): 
	if(float(events_ins[n])/float(totals[n])>0.75): print str(n)+"\tIns/Dup"
      if(events_del[n]>15): 
	if(float(events_del[n])/float(totals[n])>0.75): print str(n)+"\tDel/Jun"
      
    samfile.close()
    return
  

if __name__=="__main__":
    if len(sys.argv) == 3:
        estimate_events(sys.argv[1], sys.argv[2], 100, 450, 550)
    else:
        print "estimate_paired_end_envents.py reference.fa bam_file"