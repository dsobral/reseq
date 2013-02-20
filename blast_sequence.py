#!/usr/bin/env python
"""
This is a simple script to extract a sequence from a fasta database
"""
from Bio import SeqIO
from Bio.Blast import NCBIWWW
import sys
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
                        FileTransferSpeed, FormatLabel, Percentage, \
                        ProgressBar, ReverseBar, RotatingMarker, \
                        SimpleProgress, Timer

def blast_sequence(fasta_file, outprefix):
  
  f_gseq=open(fasta_file,"rU")
  records=list(SeqIO.parse(f_gseq, "fasta"))
  num_records=len(records)
  
  widgets = ['Progress: ', Percentage(), ' ', Bar(marker=RotatingMarker()),' ', ETA()]
  pbar = ProgressBar(widgets=widgets, maxval=len(records)).start() 
  count_record=0  
  for record in records:
    #do not send very small sequences
    if(len(record.seq)<30): continue
    blast_handle = NCBIWWW.qblast('blastn', 'nr', record.seq, format_type='HTML', megablast=True)
    count_record=count_record+1
    pbar.update(count_record)
    blast_file = open(outprefix+"_"+record.id.replace('/','_')+'_'+record.description.replace(record.id,"").replace(" ","_")+'_blast-output.html', 'w')
    blast_file.write(blast_handle.read())
    blast_file.close() 
  f_gseq.close()
  
  pbar.finish()
      
  return


instructions = "e.g. python blast_sequence.py fastafile outprefix"

if __name__=="__main__":
    if len(sys.argv) == 3:
        blast_sequence(sys.argv[1],sys.argv[2])
    else:
        print instructions