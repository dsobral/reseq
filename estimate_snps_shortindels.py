#!/usr/bin/env python
"""
This script estimates SNPs and short indels from a Bam file
It requires a certain coverage to define an indel
"""
import pysam
import sys
from Bio import SeqIO
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
                        FileTransferSpeed, FormatLabel, Percentage, \
                        ProgressBar, ReverseBar, RotatingMarker, \
                        SimpleProgress, Timer

def estimate_snps_indels(reference,filename,outfile,min_cov,min_qual_base,min_qual_map,percent,phred_add):
    samfile = pysam.Samfile(filename, "rb")
    reference = SeqIO.parse(reference, "fasta").next().seq
    output=open(outfile,"w")
    
    widgets = ['Progress: ', Percentage(), ' ', Bar(marker=RotatingMarker()),' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=samfile.lengths[0]).start() 
    
    #Carefull: cannot do this with eukaryotic/larger genome (or genome with more than one chromosome/contig)
    # e.g. 'gi|238899406|ref|NC_012759.1|'
    for pileupcolumn in samfile.pileup(samfile.references[0], 0, samfile.lengths[0]): #0-based
        
        #print pileupcolumn.pos
        pbar.update(pileupcolumn.pos)
        
        #only do a test if we have enough coverage (number is arbitrary - can be passed as argument?)  
        #Also there can be strand bias... a true SNP should be detected by reads in both strands...
        if(pileupcolumn.n < min_cov): continue
        ref_base = reference[pileupcolumn.pos]
        
        #print pileupcolumn.n
        
        total = 0
        #total = pileupcolumn.n
        
        mutant_count = 0
        
        m_dic = { }
        indel_count = 0
        i_dic = { }
        for pileupread in pileupcolumn.pileups:
            #Ignore non-informative bases (careful with case sensitivity?): use BioPython?
            if(pileupread.alignment.seq[pileupread.qpos] == 'N'): continue
            
            #Ignore reads with low mapping quality
            #print "Map Qual "+str(pileupread.alignment.mapq)
            if(pileupread.alignment.mapq < min_qual_map): continue
	    
	    # extra_conditions
	    if not pileupread.alignment.is_paired: continue	    	    
	    if(pileupread.alignment.mate_is_unmapped): continue
	  
            #Ignore positions of the reads with low quality
            read_qual = ord(pileupread.alignment.qual[pileupread.qpos])-phred_add
            #print "Read Qual "+str(read_qual)
            if(read_qual<min_qual_base): continue
            
            total = total + 1
                      
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
        
        #print total
        #sys.exit(0);
                    
        #Although pilecolumn.n may be more than min_cov, there may be other problems...
        if(total < min_cov): continue
        
        if((float(mutant_count) / float(total)) > percent): #we have a mutation
            note = 'SNP'
            #Is s a systematic mutation, or it changes?
            for snp in m_dic.keys():
                if(float(m_dic[snp]) / float(total) > percent):
                    note = note + " " + snp
                    break
            msg = '%s\t%s' % (pileupcolumn.pos + 1, note) #0-based
            output.write(msg+"\n")
            
        if((float(indel_count) / float(total)) > percent): #we have an indel  
            note = 'indel'
            for len in i_dic.keys():
                if(float(i_dic[len]) / float(total) > percent):
                    if(len > 0): note = "insertion "+ str(len) + "bp"
                    else: note = "deletion\t" + str(len) + "bp"
                    break
                
            msg = '%s\t%s' % (pileupcolumn.pos + 2, note) #0-based + 1
	    output.write(msg+"\n")
	    
    samfile.close()
    output.close()
    pbar.finish()
    return
  

if __name__=="__main__":
    if len(sys.argv) == 4:
        # estimate snps and indels with more than 75% support  
        # and ignoring bases with coverage less than 30 (arbitrary)
        # also ignores cases where the mapping quality < 20 (avoid repeats) 
        # and the read quality < 20 (avoid potential read error) 
        # Assume quality in solexa Phred+64 format 
        estimate_snps_indels(sys.argv[1],sys.argv[2],sys.argv[3],30,20,20,0.75,64)
    else:
        print "estimate_snps_shortindels.py ref_genome.fa alignment.bam output.tab"
