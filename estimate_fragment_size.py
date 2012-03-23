#!/usr/bin/env python
"""
This is a simple script to estimate fragment length of paired reads
"""
import re
import sys

def print_frags(rlen, fname_a, fname_b):
    lim = 10000 #may pass as parameter
    count = 0
    #processes two sam files
    file_a = open(fname_a, 'r')
    file_b = open(fname_b, 'r')
    for line_a in file_a:
        if(count>lim): return
        line_b = file_b.readline()
        #ignore headers
        if(re.match('^@',line_a)): continue
        line_a.rstrip()
        line_b.rstrip()
        a_fields = line_a.split("\t")
        b_fields = line_b.split("\t")
        #MAPQ>=20
        if((int(a_fields[4])<20) or (int(b_fields[4])<20)): continue
        #Only perfect matches (ignores reads that are shorter than rlen
        if((a_fields[5]!=(str(rlen)+'M')) or (b_fields[5]!=(str(rlen)+'M'))): continue
        #Calculate insert size based on strand information: this also ignores when pairs when other bits are on.
        if(a_fields[1]=='16' and b_fields[1]=='0'):
            print str(int(a_fields[3])+rlen-int(b_fields[3]))
            count = count + 1
        if(a_fields[1]=='0' and b_fields[1]=='16'):
            print str(int(b_fields[3])+rlen-int(a_fields[3]))
            count = count + 1
    file_a.close()
    file_b.close()

    return

instructions = "e.g. python estimate_fragment_length.py 90 file_1.sam file_2.sam"

if __name__=="__main__":
    if len(sys.argv) == 4:
        print_frags(int(sys.argv[1]), sys.argv[2], sys.argv[3])
    else:
        print instructions
