#!/usr/bin/env python
"""
This script estimates duplications from a Bam file
"""
import pysam
import sys
from Bio import SeqIO

def estimate_duplications(reference_file,bam_file,min_fold):
    samfile = pysam.Samfile(bam_file, "rb")
    
    reference = SeqIO.parse(reference_file, "fasta").next()
    reflen=len(reference.seq)
    refid=reference.id
    
    #now usam pysam for fast indexed access (fasta must be indexed 'samtools faidx fasta_file')
    reference = pysam.Fastafile(reference_file)
    
    gcs = {}
    coverages = {}
    gc_coverages = {}
    gc_counts = {}
    for n in range(80, 161):
	gc_counts[n]=0
	gc_coverages[n]=0
	
    #binning...
    for n in range(0, reflen, 250):
      seq=reference.fetch(refid,n,n+250)
      gc = seq.count('G')+seq.count('C')
      #Bin extremities
      if(gc>160): 
	gc=160
      if(gc<80):
	gc=80
      gc_counts[gc] = gc_counts[gc]+1
      
      cov=0
      for read in samfile.fetch(refid, n, n+250):
	cov = cov + 1
      gc_coverages[gc] = gc_coverages[gc]+cov
      coverages[n]=cov
      gcs[n]=gc
      
    #Normalize by counts  
    for n in range(80, 161):
      gc_coverages[n] = float(gc_coverages[n]) / float(gc_counts[n])
      
    for n in range(0, reflen, 250):
      print str(n)+"\t"+str(float(coverages[n])/float(gc_coverages[gcs[n]]))
      #print str(n)+"\t"+str(coverages[n])
      
    #Carefull: cannot do this with eukaryotic/larger genome (or genome with more than one chromosome/contig)
    # e.g. 'gi|238899406|ref|NC_012759.1|'
    #for pileupcolumn in samfile.pileup(reference.id, 0, samfile.lengths[0]): #0-based

    #for i in sorted(gc_counts.iterkeys()):
    #  print str(i)+"\t"+str(gc_counts[i])+"\t"+str(gc_coverages[i])
     
    samfile.close()
    return
  

if __name__=="__main__":
    if len(sys.argv) == 3:
        # estimate snps and indels with more than 75% support  
        # and ignoring bases with coverage less than 30 (arbitrary)
        estimate_duplications(sys.argv[1],sys.argv[2],1.4)
    else:
        print instructions
