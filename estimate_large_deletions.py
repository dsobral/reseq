#!/usr/bin/env python
"""
This is a simple script to estimate deletions based on absence of coverage
It needs low coverage and low quality mapping filtering
"""
import sys

def estimate_deletions(fname):
    file = open(fname, 'r')
    cur = 0
    for line in file:
        line.rstrip()
        fields = line.split("\t")
        if(int(fields[0])>cur+1): 
            print "deletion: "+str(cur)+" - "+fields[0]
        cur = int(fields[0])
    file.close()
    return

instructions = "e.g. python estimate_large_deletions.py file.cov"

if __name__=="__main__":
    if len(sys.argv) == 2:
        estimate_deletions(sys.argv[1])
    else:
        print instructions
