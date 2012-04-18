#!/usr/bin/env python
"""
This script runs the complete pipeline for bacterial 
resequencing, detection of SNPs and indels
"""
from Bio import SeqIO
import sys, copy, subprocess
import pysam
import yaml

def run_experiment(genome, home, experiment, control):
  print "Running experiment "+experiment
  
  #run alignment using paired_end information
  file_1 = home+"/raw_data/"+experiment+"_Cleandata_1.fq"
  file_2 = home+"/raw_data/"+experiment+"_Cleandata_2.fq"

  bam_out=home+"/results/"+experiment+"_paired.bam"
  
  #using bwa for single end
  bwa_index=home+"/bwa/bwa_index/"+genome
  bwa_tmp_1 = home+"/tmp/"+experiment+"_1.aln"
  bwa_tmp_2 = home+"/tmp/"+experiment+"_2.aln"
  bwa_out_1 = home+"/tmp/"+experiment+"_1.bwa.sam"
  bwa_out_2 = home+"/tmp/"+experiment+"_2.bwa.sam"
  #bwa_out=home+"/results/"+experiment+"_paired_bwa.sam"

  #check if file already exists to skip this step
  try:
    open(bwa_out_1)
    print "File "+bwa_out_1+".sam already exists... skipping alignment step"
  except IOError:
    print "Run basic alignment (single-end)"  
    subprocess.check_call(["bwa","aln", "-f", bwa_tmp_1, bwa_index,file_1])
    subprocess.check_call(["bwa","aln", "-f", bwa_tmp_2, bwa_index,file_2])
    subprocess.check_call(["bwa","samse","-f", bwa_out_1, bwa_index, bwa_tmp_1, file_1])
    subprocess.check_call(["bwa","samse","-f", bwa_out_2, bwa_index, bwa_tmp_2, file_2])
    #for paired-end
    #subprocess.check_call(["bwa","sampe", "-a", "600", "-f", bwa_out, bwa_index, bwa_tmp_1, bwa_tmp_2, file_1, file_2 ])
   
  bam_out=home+"/results/"+experiment+"_bwa_merged.bam"
  bam_out_pref_1=home+"/tmp/"+experiment+"_bwa_1_sorted"
  bam_out_pref_2=home+"/tmp/"+experiment+"_bwa_2_sorted"
  tmp_bam_1=home+"/tmp/"+experiment+"_bwa_1.tmp.bam"
  tmp_bam_2=home+"/tmp/"+experiment+"_bwa_2.tmp.bam"
  try:
    open(bam_out)
    print "File "+bam_out+" already exists... skipping sam to bam filtering step"
  except IOError:
    print "Run sam to bam filtering"  
    subprocess.check_call(["samtools","view","-S", "-b","-F","4", bwa_out_1,"-o",tmp_bam_1])
    subprocess.check_call(["samtools","sort",tmp_bam_1,bam_out_pref_1])
    subprocess.check_call(["samtools","view","-S", "-b","-F","4", bwa_out_2,"-o",tmp_bam_2])
    subprocess.check_call(["samtools","sort",tmp_bam_2,bam_out_pref_2])
    subprocess.check_call(["samtools","merge",bam_out,bam_out_pref_1+".bam",bam_out_pref_2+".bam"])
    subprocess.check_call(["samtools","index",bam_out])
    
  #using paired ssaha2
  ssaha_out=home+"/tmp/"+experiment+"_paired_ssaha2.sam"
  ssaha_bam_tmp=home+"/tmp/"+experiment+"_paired_ssaha2.bam.tmp"
  ssaha_bam_prefix=home+"/results/"+experiment+"_paired_ssaha2"
  ssaha_bam_out=home+"/results/"+experiment+"_paired_ssaha2.bam"
  try:
    open(ssaha_bam_out)
    print "File "+ssaha_bam_out+" already exists... skipping ssaha paired mappings"
  except IOError:
    print "Run paired-end mapping using ssaha2"  
    subprocess.check_call(["ssaha2 -kmer 13 -skip 1 -seeds 1 -score 12 -cmatch 9 -ckmer 1 -best 1"+
    " -save "+home+"/ssaha2_build/"+genome+".ssaha2_k13s1 -output sam"+
    " -outfile "+ssaha_out+".tmp -pair 400,600 "+file_1+" "+file_2], shell=True)
    #need to add sam header in beginning because ssaha does not put header in sam
    sam_index = home+"/samtools/"+genome+".sam_index"
    subprocess.check_call(["cat "+sam_index+" "+ssaha_out+".tmp > "+ssaha_out], shell=True)
    #Ignore unmapped reads, and only take reads with properly mapped pairs
    subprocess.check_call(["samtools","view","-S", "-b","-F","4", "-f", "2", ssaha_out,"-o",ssaha_bam_tmp])
    subprocess.check_call(["samtools","sort",ssaha_bam_tmp,ssaha_bam_prefix])
    subprocess.check_call(["samtools","index",ssaha_bam_out])
    
  #Infer SNPs and short indels
  snp_indels=home+"/results/"+experiment+"_snps_30.tab"
  try:
    open(snp_indels)
    print "File "+snp_indels+" already exists... skipping inference of SNPs and short Indels"
  except IOError:
    print "Run inference of SNPs and short indels using paired-end mappings" 
    subprocess.check_call(["python "+home+"/scripts/estimate_snps_shortindels.py "+home+"/genome/"+genome+" "+ssaha_bam_out+" "+snp_indels],shell=True)

  #TODO confirm indels with multiple alignments around the site of interest(?)

  #Infer deletions from coverage of lennient mappings
  deletions_cov=home+"/results/"+experiment+"_deletions_cov.tab"
  try:
    open(deletions_cov)
    print "File "+deletions_cov+" already exists... skipping inference of deletions using coverage"
  except IOError:
    print "Run inference of deletions by coverage using lenient mappings" 
    subprocess.check_call(["python "+home+"/scripts/estimate_large_deletions.py "+bam_out+" "+deletions_cov],shell=True)
    
  coverage=home+"/results/"+experiment+"_paired_250bp.cov" 
  try:
    open(coverage)
    print "File "+coverage+" already exists... skipping building normalized coverage"
  except IOError:
    print "Build GC normalized coverage profile" 
    subprocess.check_call(["python "+home+"/scripts/run_normalized_coverages.py "+genome+" "+ssaha_bam_out+" "+coverage],shell=True)

  #Infer deletions or duplications using control
  if(control != experiment):
    events_control=home+"/results/"+experiment+"_events_"+control+".tab"
    try:
      open(events_control)
      print "File "+events_control+" already exists... skipping inference of genomic events using comparative coverage"
    except IOError:
      print "Run inference of deletions using comparative coverage" 
      control_coverage=home+"/results/"+control+"_paired_250bp.cov"
      try:
	open(control_coverage)
	subprocess.check_call(["python "+home+"/scripts/compare_coverage_profiles.py "+coverage+" "+control_coverage+" "+events_control],shell=True)
      except IOError:
	print "Warning: File "+control_coverage+" does not exist. Skipping inference of events by comparative coverage"
    	
  #run denovo assembly (using edena?)  
  # check genomic events from results
  
  # use "broken" pairs to detect insertions and deletions?
  
  # use unmapped reads to try to find junctions (e.g. using lastz?)

def run_pipeline(config_file):
  
  with open(config_file) as in_handle:
    config = yaml.load(in_handle)
  
  genome = config['genome']
  home_dir = config['home']
  for experiment in config['experiments']:
    run_experiment(genome, home_dir, experiment['name'],experiment['control'])
  #subprocess.check_call(["ls", "-l"])

instructions = "e.g. python run_pipeline.py config_file.yaml"

if __name__=="__main__":
    if len(sys.argv) == 2:
        run_pipeline(sys.argv[1])
    else:
        print instructions