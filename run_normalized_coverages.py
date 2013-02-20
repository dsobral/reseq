#!/usr/bin/env python
"""
This script creates coverage file using a 250bp window from a Bam file
"""
import pysam
import sys
from Bio import SeqIO
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
                        FileTransferSpeed, FormatLabel, Percentage, \
                        ProgressBar, ReverseBar, RotatingMarker, \
                        SimpleProgress, Timer
                        
def run_normalized_coverages(reference_file,bam_file,outfile,window,group_min,group_max):
    samfile = pysam.Samfile(bam_file, "rb")
    output = open(outfile, "w")
      
    reference = SeqIO.parse(reference_file, "fasta").next()
    reflen=len(reference.seq)
    refid=reference.id
    
    widgets = ['Progress: ', Percentage(), ' ', Bar(marker=RotatingMarker()),' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=samfile.lengths[0]).start() 
    
    #now usam pysam for fast indexed access (fasta must be indexed 'samtools faidx fasta_file')
    reference = pysam.Fastafile(reference_file)
    
    gcs = {}
    coverages = {}
    gc_coverages = {}
    gc_counts = {}
    for n in range(group_min, group_max+1):
	gc_counts[n]=0
	gc_coverages[n]=0
	
    #binning...
    for n in range(0, reflen, window):
      pbar.update(n)
      seq=reference.fetch(refid,n,n+window)
      gc = seq.count('G')+seq.count('C')
      #Bin extremities
      if(gc>group_max): 
	gc=group_max
      if(gc<group_min):
	gc=group_min
      gc_counts[gc] = gc_counts[gc]+1
      
      cov=0
      #this does not work with split alignments
      for read in samfile.fetch(refid, n, n+window):
	cov = cov + 1
      gc_coverages[gc] = gc_coverages[gc]+cov
      coverages[n]=cov
      gcs[n]=gc
      
    #Normalize by counts  
    for n in range(group_min, group_max+1):
      gc_coverages[n] = float(gc_coverages[n]) / float(gc_counts[n])
      #output.write(str(n)+"\t"+str(gc_coverages[n])+"\n")
      
    
    for n in range(0, reflen, window):
      ##print str(n)+"\t"+str(float(coverages[n])/float(gc_coverages[gcs[n]]))
      output.write(str(n)+"\t"+str(float(coverages[n])/float(gc_coverages[gcs[n]]))+"\n")
      ##print str(n)+"\t"+str(coverages[n])
      
    #Carefull: cannot do this with eukaryotic/larger genome (or genome with more than one chromosome/contig)
    # e.g. 'gi|238899406|ref|NC_012759.1|'
    #for pileupcolumn in samfile.pileup(reference.id, 0, samfile.lengths[0]): #0-based

    #for i in sorted(gc_counts.iterkeys()):
    #  print str(i)+"\t"+str(gc_counts[i])+"\t"+str(gc_coverages[i])
     
    samfile.close()
    output.close()
    pbar.finish()
    return
  

if __name__=="__main__":
    if len(sys.argv) == 4:
        run_normalized_coverages(sys.argv[1],sys.argv[2],sys.argv[3],250,80,160)
        #estimate_duplications(sys.argv[1],sys.argv[2],sys.argv[3],100,30,70)
    else:
        print "run_normalized_coverages.py reference.fa bam_file outfile"
