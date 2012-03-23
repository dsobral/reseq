import re
import pysam
import sys
import subprocess
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from sys import stderr
#import Bio.Align.Applications
#from Bio.Align.Applications import TCoffeeCommandline


#Get the reference (put in config file)
reference = SeqIO.parse("../../ecoli_NC012759.1_bw2952.fa", "fasta").next().seq

#consider using this
#pyreference = pysam.Fastafile("../../ecoli_NC012759.1_bw2952.fa")

#only need to do this once
#ssaha2Build -kmer 13 -skip 1 -save ecoli_NC012759.1_bw2952.fa_k13s1 ecoli_NC012759.1_bw2952.fa
# seeds 1 is important! (can also try -solexa params which are quite similar)
#ssaha2 -kmer 13 -skip 1 -seeds 1 -score 12 -cmatch 9 -ckmer 1 -save ecoli_NC012759.1_bw2952_k13s1.fa -pair 400,600 -output sam -outfile ecoli_NC012759.1_bw2952_Muc_k13s1_paired.sam ../raw_data/Muc_Cleandata_1.fq ../raw_data/Muc_Cleandata_2.fq
#Stringent mappings (pass as parameter)
samfile = pysam.Samfile("../alignments/ecoli_NC012759.1_bw2952_Anc_k13s1_paired.sorted.bam", "rb")

p = re.compile('[ACTG]')

rlen = 90 #pass as parameter

#Carefull: cannot do this with eukaryotic/larger genome (or genome with more than one chromosome/contig)
# 'gi|238899406|ref|NC_012759.1|'
#for pileupcolumn in samfile.pileup(samfile.references[0], 0, samfile.lengths[0]): #0-based
for pileupcolumn in samfile.pileup(samfile.references[0], 1289300, 1289350): #0-based
  #Print coverage
  #print '%s\t%s' % (pileupcolumn.pos , pileupcolumn.n)
  
  #only do a test if we have enough coverage (20 is somewhat arbitrary - can be passed as argument?)  
  #Also there can be strand bias... a true SNP should be detected by reads in both strands...
  if(pileupcolumn.n < 30): continue
  ref_base = reference[pileupcolumn.pos]
  total = pileupcolumn.n
  mutant = 0
  m_dic = { }
  indel = 0
  i_dic = { }
  for pileupread in pileupcolumn.pileups:
    #Ignore non-informative bases (careful with case sensitivity?)
    if(pileupread.alignment.seq[pileupread.qpos] == 'N'): continue
    #Ignore bases with deletions (as they should have been accounted for already)
    if(pileupread.is_del): continue
    #just in case!
    if(pileupread.alignment.is_unmapped): continue
    #Information about indels is in the previous base... (e.g. if it is an insertion)
    if(pileupread.indel != 0): 
        indel = indel + 1     
        if(pileupread.indel not in i_dic): i_dic[ pileupread.indel ] = 1
        else: i_dic[ pileupread.indel ] = i_dic[ pileupread.indel ] + 1
    if(ref_base != pileupread.alignment.seq[pileupread.qpos]): 
        mutant = mutant + 1
        note = ref_base+">"+pileupread.alignment.seq[pileupread.qpos]
        if(note not in m_dic): m_dic[ note ] = 1
        else: m_dic[ note ] = m_dic[ note ] + 1

  if((float(mutant) / float(total)) > 0.75): #we have a mutation
      note = 'SNP'
      for snp in m_dic.keys():
          if(float(m_dic[snp]) / float(total) > 0.75):
              note = note+" : "+ snp
              break
      print '%s\t%s' % (pileupcolumn.pos+1, note) #0-based

  if((float(indel) / float(total)) > 0.75): #we have an indel  
      note = 'indel' 
      for len in i_dic.keys():
          if(float(i_dic[len]) / float(total) > 0.75):
              if(len>0): note = "insertion: "+str(len)+"bp"
              else: note = "deletion: "+str(len)+"bp"
              break
      print '%s\t%s' % (pileupcolumn.pos+2, note) #0-based + 1
      
      #ref_start = pileupcolumn.pos-rlen+1
      #ref_end = pileupcolumn.pos+rlen+1
      #records = []
      #Get reference around the indel +- read length
      #records.append(SeqRecord(seq = reference[(ref_start):(ref_end)], id = "reference", description="reference"))
      #We need to run T-Coffee in this position
      #for pileupread in pileupcolumn.pileups:
      #    if(pileupread.alignment.is_reverse): continue
      #    #output sequences
      #    records.append(SeqRecord(seq = Seq(pileupread.alignment.seq,IUPAC.unambiguous_dna), id = pileupread.alignment.qname, description=note))
      
      #SeqIO.write(records, "test.fasta", "fasta")
      #command_line= ['t_coffee', 'test.fasta', '-type', 'dna', '-output','aln']
      #subprocess.check_call(command_line)
      #alignment = AlignIO.read("test.aln", "clustal")
      #for record in alignment: print "%s - %s" % (record.seq, record.id)

   
  #To use quality values and alignment score
  #print '%s\t%f' % (pileupcolumn.pos,float(mutant)/float(total))
  #  print '\t%s\t%s\t%s\t%s' % (reference[pileupcolumn.pos], pileupread.alignment.seq[pileupread.qpos], 
  #  pileupread.alignment.qual[pileupread.qpos], pileupread.alignment.opt('AS'))


samfile.close()

