#!/usr/bin/env python
"""
This script runs the complete pipeline for bacterial 
resequencing, detection of SNPs and indels
"""
from Bio import SeqIO
import sys, copy, subprocess
import pysam
import yaml

def run_experiment(genome, home, experiment):
  print "Running experiment "+experiment
  
  #run alignment using paired_end information
  file_1 = home+"/raw_data/"+experiment+"_Cleandata_1.fq"
  file_2 = home+"/raw_data/"+experiment+"_Cleandata_2.fq"

  bam_out=home+"/results/"+experiment+"_paired.bam"
  #using bwa  
  bwa_index=home+"/bwa/bwa_index/"+genome
  bwa_tmp_1 = home+"/tmp/"+experiment+"_1.aln"
  bwa_tmp_2 = home+"/tmp/"+experiment+"_2.aln"
  bwa_out=home+"/results/"+experiment+"_paired_bwa.sam"

  #check if file already exists to skip this step
  try:
    open(bwa_out)
    print "File "+bwa_out+".sam already exists... skipping alignment step"
  except IOError:
    print "Run alignment using paired-end information"  
    subprocess.check_call(["bwa","aln", "-f", bwa_tmp_1, bwa_index,file_1])
    subprocess.check_call(["bwa","aln", "-f", bwa_tmp_2, bwa_index,file_2])
    subprocess.check_call(["bwa","sampe", "-a", "600", "-f", bwa_out, bwa_index, bwa_tmp_1, bwa_tmp_2, file_1, file_2 ])
   
  bam_out=home+"/results/"+experiment+"_paired_bwa.bam"
  bam_out_pref=home+"/results/"+experiment+"_paired_bwa"
  tmp_bam=home+"/tmp/"+experiment+"_paired_bwa.tmp.bam"
  try:
    open(bam_out)
    print "File "+bam_out+" already exists... skipping sam to bam filtering step"
  except IOError:
    print "Run sam to bam filtering"  
    #Ignore unmapped reads, and only take reads with properly mapped pairs
    subprocess.check_call(["samtools","view","-S", "-b","-F","4", "-f", "2", bwa_out,"-o",tmp_bam])
    subprocess.check_call(["samtools","sort",tmp_bam,bam_out_pref])
    subprocess.check_call(["samtools","index",bam_out])
  
  #using paired ssaha2
  #ssaha_out=home+"/results/"+experiment+"_paired_ssaha2.sam"
  #subprocess.check_call(["ssaha2 -kmer 13 -skip 1 -seeds 1 -score 12 -cmatch 9 -ckmer 1 -best 1"+
  #" -save "+home+"/ssaha2_build/"+genome+".ssaha2_k13s1 -output sam"+
  #" -outfile "+ssaha_out+" -pair 400,600 "+file_1+" "+file_2"], shell=True)
  #need to add sam header in beginning because ssaha does not put header in sam
  
  #Infer SNPs and short indels
  snp_indels=home+"/results/"+experiment+"_indels.txt"
  try:
    open(snp_indels)
    print "File "+snp_indels+" already exists... skipping inference of SNPs and short Indels"
  except IOError:
    print "Run inference of SNPs and short indels" 
    subprocess.check_call(["python "+home+"/scripts/estimate_snps_shortindels.py "+home+"/genome/"+genome+" "+bam_out+" > "+snp_indels],shell=True)
  
  #confirm indels with multiple alignments around the site of interest(?)
  
  #run denovo assembly (using edena?)  
  # check genomic events from results
  
  # use relaxed mappings to check deletions
  
  # use unmapped reads and "broken" pairs to detect insertions?

def run_pipeline(config_file):
  
  with open(config_file) as in_handle:
    config = yaml.load(in_handle)
  
  genome = config['genome']
  home_dir = config['home']
  for experiment in config['experiments']:
    run_experiment(genome, home_dir, experiment['name'])
  #subprocess.check_call(["ls", "-l"])

instructions = "e.g. python run_pipeline.py config_file.yaml"

if __name__=="__main__":
    if len(sys.argv) == 2:
        run_pipeline(sys.argv[1])
    else:
        print instructions