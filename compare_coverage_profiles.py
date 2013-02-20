#!/usr/bin/env python
"""
This is a simple script to large genome events based on coverage files
Use run_coverage_window to obtain the coverage files
"""
import sys

def estimate_events(coverage_ref,coverage_mut,outfile,change):
    file_ref = open(coverage_ref, 'r')
    file_mut = open(coverage_mut, 'r')
    output=open(outfile,"w")
    deletion=False
    duplication=False
    event_begin="0"
    for line_ref in file_ref:
        line_ref.strip()
        fields_ref = line_ref.split("\t")
        point_ref=fields_ref[0]
        fold_ref=float(fields_ref[1])
        line_mut=file_mut.readline()
        line_mut.strip()
        fields_mut = line_mut.split("\t")
        point_mut=fields_mut[0]
        fold_mut=float(fields_mut[1])
        
        if(point_mut != point_ref): 
	  print "Error: coverage files do not seem compatible!"
	  sys.exit(1)
        
        if((fold_mut-fold_ref)>change):
	  #duplication
	  if(not duplication):
	    if(deletion):
	      output.write("deletion\t"+event_begin+"\t"+point_ref+"\n")
	      deletion=False
	    duplication=True
	    event_begin=point_ref
	elif((fold_ref-fold_mut)>change):
	  #deletion
	  if(not deletion):
	    if(duplication):
	      output.write("duplication\t"+event_begin+"\t"+point_ref+"\n")
	      duplication=False
	    deletion=True
	    event_begin=point_ref
	else:
	  if(duplication):
	    output.write("duplication\t"+event_begin+"\t"+point_ref+"\n")
	    duplication=False
	  if(deletion):
	    output.write("deletion\t"+event_begin+"\t"+point_ref+"\n")
	    deletion=False
	  
        
    file_ref.close()
    file_mut.close()
    output.close()
    return

instructions = "e.g. python estimate_large_deletions.py window_coverage_ref.cov window_coverage_mut.cov events_file.tab"

if __name__=="__main__":
    if len(sys.argv) == 4:
        estimate_events(sys.argv[1],sys.argv[2],sys.argv[3],0.5)
    else:
        print instructions