#!/usr/bin/env python
"""
This script estimates duplications from a Bam file
"""
import pysam
import sys
from Bio import SeqIO

def estimate_events(bam_file):
    samfile = pysam.Samfile(bam_file, "rb")
    
    #reference = SeqIO.parse(reference_file, "fasta").next()
    #reflen=len(reference.seq)
    #refid=reference.id

    for aligned_read in samfile.fetch(): #0-based
      num_events=0
      #just in case!
      if(aligned_read.is_unmapped): continue
      if aligned_read.is_paired:
	#If we're doing pair number 2 then stop because we did it already!
	if(aligned_read.is_read2): continue
	#Only check if the mate is mapped
	if(aligned_read.mate_is_unmapped):continue
	print str(aligned_read.pos)+"\t"+str(aligned_read.pnext)

	#If paired...
        #Check event...
        #if event: num_events++
          
     
    samfile.close()
    return
  

if __name__=="__main__":
    if len(sys.argv) == 2:
        # estimate snps and indels with more than 75% support  
        # and ignoring bases with coverage less than 30 (arbitrary)
        estimate_events(sys.argv[1], 400, 600)
    else:
        print "estimate_paired_end_envents.py bam_file"