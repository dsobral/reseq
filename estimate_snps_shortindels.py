#!/usr/bin/env python
"""
This script estimates SNPs and short indels from a Bam file
It requires a certain coverage to define an indel
"""
import pysam
import sys
from Bio import SeqIO

def estimate_snps_indels(reference,filename,min_cov,percent):
    samfile = pysam.Samfile(filename, "rb")
    reference = SeqIO.parse(reference, "fasta").next().seq

    #Carefull: cannot do this with eukaryotic/larger genome (or genome with more than one chromosome/contig)
    # e.g. 'gi|238899406|ref|NC_012759.1|'
    for pileupcolumn in samfile.pileup(samfile.references[0], 0, samfile.lengths[0]): #0-based
        
        #only do a test if we have enough coverage (number is arbitrary - can be passed as argument?)  
        #Also there can be strand bias... a true SNP should be detected by reads in both strands...
        if(pileupcolumn.n < min_cov): continue
        ref_base = reference[pileupcolumn.pos]
        total = pileupcolumn.n
        #total_himapqual = 0
        #total_hireadqual = 0
        
        mutant_count = 0
        #mutant_himapqual = 0
        #mutant_hireadqual = 0
        
        m_dic = { }
        indel_count = 0
        i_dic = { }
        for pileupread in pileupcolumn.pileups:
            #Ignore non-informative bases (careful with case sensitivity?): use BioPython?
            if(pileupread.alignment.seq[pileupread.qpos] == 'N'): continue
            #Ignore bases with deletions (as they should have been accounted for already)
            if(pileupread.is_del): continue
            #just in case!
            if(pileupread.alignment.is_unmapped): continue
            #Information about indels is in the previous base... (e.g. if it is an insertion)
            if(pileupread.indel != 0): 
                indel_count = indel_count + 1     
                if(pileupread.indel not in i_dic): i_dic[ pileupread.indel ] = 1
                else: i_dic[ pileupread.indel ] = i_dic[ pileupread.indel ] + 1
                
            if(ref_base != pileupread.alignment.seq[pileupread.qpos]): 
                mutant_count = mutant_count + 1
                note = ref_base + ">" + pileupread.alignment.seq[pileupread.qpos]
                if(note not in m_dic): m_dic[ note ] = 1
                else: m_dic[ note ] = m_dic[ note ] + 1
                    
        if((float(mutant_count) / float(total)) > percent): #we have a mutation
            note = 'SNP'
            #Is s a systematic mutation, or it changes?
            for snp in m_dic.keys():
                if(float(m_dic[snp]) / float(total) > percent):
                    note = note + " " + snp
                    break
            print '%s\t%s' % (pileupcolumn.pos + 1, note) #0-based
              
        if((float(indel_count) / float(total)) > percent): #we have an indel  
            note = 'indel'
            for len in i_dic.keys():
                if(float(i_dic[len]) / float(total) > percent):
                    if(len > 0): note = "insertion "+ str(len) + "bp"
                    else: note = "deletion\t" + str(len) + "bp"
                    break
                
            print '%s\t%s' % (pileupcolumn.pos + 2, note) #0-based + 1
      
    samfile.close()
    return
  

if __name__=="__main__":
    if len(sys.argv) == 3:
        # estimate snps and indels with more than 75% support  
        # and ignoring bases with coverage less than 30 (arbitrary)
        estimate_snps_indels(sys.argv[1],sys.argv[2],30,0.75)
    else:
        print instructions
