#!/usr/bin/env python
"""
This is a simple script to estimate deletions based on absence of coverage
It needs low coverage and low quality mapping filtering
"""
import sys
import pysam
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
                        FileTransferSpeed, FormatLabel, Percentage, \
                        ProgressBar, ReverseBar, RotatingMarker, \
                        SimpleProgress, Timer
                        
def estimate_deletions(fname,outfile,read_len):
    samfile = pysam.Samfile(fname, "rb")
    output = open(outfile,"w")
    
    widgets = ['Progress: ', Percentage(), ' ', Bar(marker=RotatingMarker()),' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=samfile.lengths[0]).start() 
    
    deletion=False
    event_begin="0"
    
    for n in range(0,samfile.lengths[0],read_len):
      
      pbar.update(n)
      #print n
      #Should contain only one!
      pileups=samfile.count(samfile.references[0], n, n+read_len) #0-based
      
      if(pileups==0):
	if(not deletion):
	  deletion=True
	  event_begin=str(n)
      else:
	if(deletion):
	  output.write("deletion\t"+event_begin+"\t"+str(n)+"\n")
	  deletion=False
        
    samfile.close()
    output.close()
    pbar.finish()
    return

instructions = "e.g. python estimate_large_deletions.py lenient_alignments.bam outputfile.tab"

if __name__=="__main__":
    if len(sys.argv) == 3:
        estimate_deletions(sys.argv[1],sys.argv[2],45)
    else:
        print instructions
