RESEQUENCING PIPELINE

Author: 
  Daniel Sobral (dsobral@igc.gulbenkian.pt)
  Bioinformatics Unit, IGC

INTRODUCTION

  The purpose of this pipeline is to identify a series of genomic events in a bacterial strain, given a reference sequence.
  It is to be used only as a complement to breseq (http://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing)

  It can identify:
  - SNPs and short indels
  - large deletions
  - potential duplications (when comparing two strains e.g. evolved vs parental)
  - potential novel sequence (with a denovo assembly analysis)
  - potential junctions

  Basic description of the whole process (see more details later in the document):
    1 - Align individual reads to reference genome using bwa
    2 - Align reads using pair-end information with ssaha2 to reference genome
    3 - Calculate snps and short indels using stringent paired-end mappings (from step 2):
    4 - Calculate large deletions based on lack of coverage, using lenient single-end mappings (from step 1)
    5 - Calculate GC Normalized Coverage profiles in windows of 250bp (using alignments from step 1)
    6 - Calculate possible duplications / deletions by comparison of coverage of an evolved and an original strain 
      (if there is no original strain or coverage of original strain has not been calculated, this step is skipped)
    # Experimental steps
    7 - Perform a denovo assembly using edena
    8 - Map the denovo assembly against the reference (requiring 100% accuracy) and check for possible genomic events
    9 - Separate sequences of the denovo assembly that didn't match the reference and blast them against NCBI nr (uses internet connection)
    10 - Take unmapped reads (from step 1) and try to remap them against the reference allowing for split reads, to identify possible junction events

REQUIREMENTS

  You'll need a Linux environment (tested with Ubuntu 10.04.4 LTS).
    - It was also tested successfully in the Elephant server (Ubuntu 12.04.1 LTS)
    
  You need Python installed (tested with Python 2.6.5 and Python 2.7.3 in Elephant)
  You'll need the following Python libraries
    - pysam (to manipulate SAM and BAM alignment files)
    - pyyaml (only to use the config files with run_pipeline.py)
    - BioPython (for certain parts of the denovo sequencing analysis)
    - progressbar (to have a nice progress bar!)

  You'll also need the following programs:
    - bwa (tested with version 0.6.1-r104)
      http://bio-bwa.sourceforge.net/
    - ssaha2 (tested with version 2.5.5)
      http://www.sanger.ac.uk/resources/software/ssaha2/
    - samtools (tested with version 0.1.18 (r982:295))
      http://samtools.sourceforge.net/

  The denovo part of the pipeline also requires:
    - edena (tested with v3.121122)
      http://www.genomic.ch/edena.php
    - lastz (tested with version 1.02.00 built 20100112)
      http://www.bx.psu.edu/~rsharris/lastz/
    - bedtools (tested with v2.16.2)
      http://code.google.com/p/bedtools/

  All the tools in the versions used in the pipeline can also be located in pipeline/virtualenv/bin

VIRTUAL_ENV

  Although not required, I recommend using a python virtual environment (you need to install virtualenv to do so)

  You can recover the existing virtual env doing:
    pipeline> virtualenv virtualenv
    pipeline> source virtualenv/bin/activate

  You may need to reinstall the python libraries (due to changed python versions):
    (virtualenv)pipeline> easy_install pysam
    (virtualenv)pipeline> easy_install numpy
    (virtualenv)pipeline> easy_install biopython
    (virtualenv)pipeline> easy_install pyyaml
    (virtualenv)pipeline> easy_install progressbar

  run "deactivate" when you're done with running the pipeline

  You can also create a new virtual environment:
    pipeline> virtualenv newvirtualenv
    pipeline> source newvirtualenv/bin/activate
  Then you need to install all the python libraries (as above), as well as add the programs to newvirtualenv/bin (or $PATH)
  

DIRECTORY STRUCTURE AND REQUIRED DATA

  The pipeline requires a certain directory structure:
  - pipeline
    --- bwa
    ------ bwa_index (to store indexes for bwa)
    --- genome (to store reference sequences)
    --- raw_data (to store the raw reads)
    --- results (where the results of the pipeline will be stored
    --- samtools (to store indexes and other information required by samtools)
    --- scripts (where all the scripts to run the pipeline are located)
    --- ssaha2_build (to store indexes for ssaha2)
    --- tmp (required to store intermediate steps of the pipeline - you should delete its contents from time to time)
    --- virtualenv (where the python virtual environment is located)
    ------ bin (where the binaries for the virtual environment are stored, including programs needed to run the pipeline)
  
  The reference sequences are stored in pipeline/genome
    E.g. pipeline/genome/ecoli_NC_000913.2_MG1655.fa (MG1655 strain)
    E.g. pipeline/genome/ecoli_NC012759.1_bw2952.fa (BW2952 strain)

  Then you need indexes for alignments
  
  Indexes for bwa are located in pipeline/bwa/bwa_index
  To generate these indexes, you need to download the reference sequence (as fasta), and do the following:
    pipeline> bwa index genome/reference.fa
  Then move the results of this process to pipeline/bwa/bwa_index
    pipeline> mv genome/reference.fa.* bwa/bwa_index

  Indexes for ssaha are located in pipeline/ssaha2_build
  To generate these indexes, you need to download the reference sequence (as fasta) and do the following
    pipeline> ssaha2Build -k 13 -s 1 -save ssaha2_build/reference.fa.ssaha2_k13s1 genome/reference.fa
    Note: -s 1 is important for speed; -k 13 from ssaha manual for solexa, and also used in Breseq (v0.16)

  The alignments from bwa and ssaha2 are in SAM (text) format. To optimize further steps of the pipeline, 
  we need to convert from SAM to an optimized indexed binary format: BAM. For this we need samtools.

  Due to the way ssaha2 outputs sam, we need an external header to be added to the ssaha2 sam output before 
  running samtools. This header is located in pipeline/samtools
  For bacteria, it will basically consist of a single line with the identifier of the reference sequence and its size
    E.g. samtools/ecoli_NC_000913.2_MG1655.fa.sam_index 
    @SQ	SN:gi|49175990|ref|NC_000913.2|	LN:4639675
    E.g. samtools/ecoli_NC012759.1_bw2952.fa.sam_index
    @SQ	SN:gi|238899406|ref|NC_012759.1|	LN:4578159

  For steps 8, 9 and 10 of the pipeline, we access the reference fasta file. The same way that the SAM file can be 
  optimized using an indexed format, we can also accelerate fasta access using indexes.

  This is done using samtools:
    pipeline> samtools faidx genome/reference.fa
  This wil generate the following file: genome/reference.fa.fai
    E.g. genome/ecoli_NC012759.1_bw2952.fa.fai
    E.g. genome/ecoli_NC_000913.2_MG1655.fa.fai

  For parts of the pipeline involving the denovo assembly we require bedtools. This needs a header file with 
  information of genome size. This information is in genome/reference.fa.sizes.bed
    E.g. ecoli_NC_000913.2_MG1655.fa.sizes.bed
    gi	1	4639675
    (because of the way the id is read in these steps, only the part of the reference id until the first '|' is considered)


RUNNING THE PIPELINE

  The pipeline is run using scripts/run_pipeline.py

  It takes as imput a configuration file in yaml format.
  This configuration file has 3 main parameters:
  - home (indicating the location of the pipeline where we have the proper directory structure, as seen above)
  E.g. home: /home/dsobral/reseq/pipeline
  - genome (indicating the reference genome, located in home/genome)
  E.g. genome: ecoli_NC012759.1_bw2952.fa
  - experiments: (a set of experiments that are to be run). 
    It has 2 subparameters: name and control. 
    The name indicates the strain that we're currently analyzing;
    The control indicates the original strain from which name is derived (or name again if no original strain exists)
    There can be as many experiments as needed

  E.g. pipeline/config.yaml: 
    home: /home/dsobral/reseq/pipeline
    genome: ecoli_NC012759.1_bw2952.fa
    experiments:
      - name: Anc
        control: Anc
      - name: M1
        control: Anc

  For each different name, there should be 2 files with paired end reads in pipeline/raw_data:
    E.g. raw_data/Anc_Cleandata_1.fq raw_data/Anc_Cleandata_2.fq
    E.g. raw_data/M1_Cleandata_1.fq raw_data/M1_Cleandata_2.fq

    Note: The suffixes *_Cleandata_1.fq and *_Cleandata_2.fq are currently required for run_pipeline.py 
  
    The original strain should appear always before the evolved strains (or else step 6 of the pipeline is skipped)
    If the original strain was already done in a different run of the pipeline, and results are in pipeline/results 
      then you don't need to run it again

  Example of running the pipeline:
  (virtualenv)pipeline> python scripts/run_pipeline.py config.yaml


OUTPUT OF THE PIPELINE

  For each sample, the pipeline outputs a series of files (mostly text) in pipeline/results
  
  - Single end Mappings (from step 1): results/sample_name_bwa_merged.bam results/sample_name_bwa_merged.bam.bai
    This contains the binary BAM file for the alignments, as well as its companion (always necessary) index .bai file
    E.g. results/M1_bwa_merged.bam results/M1_bwa_merged.bam.bai

    Note: These mappings can be visualized in a browser such as IGV, together with the reference sequence

  - Paired end Mappings (from step 2): results/sample_name_paired_ssaha2.bam results/sample_name_paired_ssaha2.bam.bai
    This contains the binary BAM file for the alignments, as well as its companion (always necessary) index .bai file
    E.g. results/M1_paired_ssaha2.bam results/M1_paired_ssaha2.bam.bai

    Note: These mappings can be visualized in a browser such as IGV, together with the reference sequence

  - SNPs and small indels (from step 3): results/sample_name_snps30.tab
    This contains one line per event, with the position of the event, and a description of the event
    E.g. results/M1_snps_30.tab
      784643  SNP T>G
      2614573 insertion 2bp

  - large deletions (from step 4): results/sample_name_deletions_cov.tab    
    This contains genomic intervals where there are no reads aligned
    E.g. results/M1_deletions_cov.tab
      deletion        1286145 1288440

  - GC Normalized Coverage profiles in windows of 250bp (from step 5): results/sample_name_relaxed_250bp.cov
    This contains, for each non-overlapping 250bp window, the window's genomic position and the normalized coverage
    E.g. results/M1_relaxed_250bp.cov
      0	0.803635637235
      250	0.827299032849

    Note: This coverage file can be read into R and processed there to detect genomic events

  - Putative large indels based on comparative coverage (from step 6): results/sample_name_events_original_strain_relaxed.tab
    This contains possible genomic deletions or duplications of an evolved strain relative to an original strain
    E.g. results/M1_events_Anc_relaxed.tab
      deletion        2873500 2873750

    Note: in older versions of the pipeline there was a bug in that deletions were named duplications and vice-versa
    Some of the results of previous runs, particularly of older strains, may have this bug. 
    E.g. 13YFP_events_OYFP_relaxed.tab
      deletion        4512750 4513250 (this is actually a duplication)
    The new version is corrected.

   # The following are experimental steps, may produce useful information but may also contain a lot of noise

  - Denovo assembly (step 7): results/sample_name_edena_contigs.fasta results/sample_name_edena_*
    This contains the fasta file (as well as other files) of the edena denovo genome assembly using paired reads
    E.g. results/M1_edena_contigs.fasta

  - A map of the denovo assembly against the reference (step 8): results/sample_name_edena_genome_reference.fa.map
      This contains the mapping (using lastz) of the denovo contigs against the reference at 100% identity
      Each line of the file contains a set of tabular columns: 
      contig_name contig_start contig_end contig_strand	reference_name reference_start reference_end reference_strand align_cigarx
      E.g. results/M1_edena_ecoli_NC012759.1_bw2952.fa.map
	#name2  start2  end2    strand2 name1   start1  end1    strand1 cigarx
	results/M1_edena_1      1       548     +       gi      1285564 1286111 +       548=
	results/M1_edena_1      1       327573  -       gi      1541554 1869126 +       327573=

  - Possible genomic events based on mappings of denovo assembly against reference (step 8): results/sample_name_genome_events.tab
      This contains a list of possible events based on an interpretation of the mappings file
      Given that it is prone to false positives, it also provides extra evidence to help filter those cases
      E.g. results/M1_genome_events.tab
	...
	4103050-4103051 variation       (1bp)
	4131201-4180366 deletion        (49165bp)
	  region 4180360-4409444 is covered by results/M1_edena_10
	  Caution: event may be false

  - Sequences of the denovo assembly that didn't match the reference (step 9): results/sample_name_edena_novel_reference_genome.fa.fasta
    This contains sequences of the denovo assembly that could not be mapped by lastz against the reference
    It will also contain very small sequences corresponding to SNPs and small insertions (since we map at 100% similarity)
    E.g. results/M1_edena_novel_ecoli_NC012759.1_bw2952.fa.fasta
      >results/M1_edena_9 99990-102215
      CGGCATACCGCTGTTTATCCATGTGATCCTGCCTGGTGCCCTGCCCTCAATTATGGTCGG
      CGTGCGTTTTGCGTTGGGCCTGATGTGGCTGACGCTGATTGTTGCCGAAACCATTTCTGC
      CAATTCAGGCATTGGTTATCTGGCGATGAATGCGCGGGAGTTTTTGCAAACGGACGTGGT
      ...
      >results/M1_edena_10 2600-2602
      C

  - Html results of the blast of the unmapped contigs against NCBI nr (step 9): results/sample_name_results_sample_name_edena_*_blast-output.html
    For each of the unmapped sequences with more than 30bp (to avoid SNPs and small indels), there is a file with the blast result against NCBI nr
    The output is in html just for visualization purposes.
    E.g. results/M1_results_M1_edena_9__99990-102215_blast-output.html

  - Possible junctions based on remapping of unmapped reads allowing splits (step 10): results/sample_name_junctions.txt
    This file reports splits that occurr more than 10 times between the same coordinates (e.g. a control is the end and beginning of the chromosome) 
    For each possible junctions, it also reports the number of reads supporting that junction
    E.g. results/M1_junctions.txt
      1<->4578160     34
      951969<->3411611        25
      3311974<->3831376       25
  
  Finally, the pipeline also leaves temp files in pipeline/tmp that can be useful on their own, or otherwise deleted


DESCRIPTION OF INDIVIDUAL TOOLS 

  Each tool can be used independently:

  - Estimate SNPs and short indels
    (virtualenv) pipeline> python scripts/estimate_snps_shortindels.py reference.fa alignments.bam output_events.txt

    The result is a file with the the position of the event and the events description (output_events.txt)
    EVENT_POSITION	EVENT_DESCRIPTION

  - Estimate Large deletions based on a mapping file (usually relaxed alignments)
    (virtualenv) pipeline> python estimate_large_deletions.py lenient_alignments.bam output_deletions.tab

    The result  is a file with genomic intervals where no reads were detected (output_deletions.tab)
    deletion	DELETION_START	DELETION_END

  - Run GC Normalized coverage (only for 250bp windows)
    (virtualenv) pipeline> python scripts/run_normalized_coverages.py reference.fa alignments.bam output_coverage.cov

    The result is a file with normalized coverage for each 250bp window (output_coverage.cov)
    WINDOW_START	GC_NORMALIZED_COVERAGE

    This output_coverage.cov file can be read into R and processed there to detect genomic events

  - Estimate Large Indels based on the comparison between two normalized coverage profiles (usually relaxed alignments)  
    (virtualenv) pipeline> python (virtualenv) python compare_coverage_profiles.py coverage_ref.cov coverage_mut.cov output_events.tab

    The result is a file with intervals where the normalized coverage differs of more than 0.5 between the original and evolved 
    strains (output_events.tab). Event type can be either deletion or duplication
    EVENT_TYPE	EVENT_START	EVENT_END
        
   # Other scripts are more experimental and can be used at the user's risk. 
    Documentation in those scripts is likely to be scarse/incomplete/outdated


OBTAINING THE LATEST VERSION OF THE SCRIPTS

  The latest version of the scripts can be obtained from git:
    pipeline> git clone git@github.com:dsobral/reseq.git scripts
  To update:
    pipeline/scripts> git pull

  To use run_pipeline.py, scripts need to be in pipeline/scripts
  If you want to use individual scripts, they can be anywhere

DETAILS AND NOTES ON EACH OF THE STEPS OF THE PIPELINE:

  1 - Align individual reads to reference genome using bwa
  
    To align reads with bwa we use default parameters. The process is simple:
      pipeline> bwa aln raw_data\reads.fq > tmp/reads.aln
      pipeline> bwa samse tmp/reads.aln raw_data\reads.fq > tmp/alignments.sam

    Note: In bwa, if a read maps to multiple places, it will place it in one of the possibilities randomly
    These cases can be filtered a-posteriori, as the mapQ value will be 0
    If using other software this behavior may be different (e.g. for bowtie you need to specify parameter -M)

    I align the first pair and the second pair, independently, and then transform both to bam and merge them into a single bam
  
  2 - Align reads using pair-end information with ssaha2 to reference genome

    Parameters passed to ssaha2 were selected according to what the ssaha manual suggests for solexa reads.
    We also selected some parameters as they were used in breseq v0.16
    Parameters: -kmer 13 -skip 1 -seeds 1 -score 12 -cmatch 9 -ckmer 1 -pair 400,600

    Note: We need to pass ssaha a bound [lower, upper] for the fragment. I give a manual bound of [400, 600], so this will 
      not work if the fragment size is very different from 500bp (what's being used now)
      For the general case, we would need to add an extra step to estimate fragment length bounds to pass to ssaha
      bwa, in the paired end mode, does its own estimate of fragment lengths based on unique mappings.
      So basically we could replace ssaha2 with bwa for paired end reads too...

  3 - Calculate snps and short indels using stringent paired-end mappings (from step 2):
    A SNP or short indel is called in a position p of the genome if:
      a) The number of reads with a high quality (Q>20) base on position p is greater than 30
      b) The same event (e.g. SNP A>T or deletion 1bp) happens in >75% of those bases

  4 - Calculate large deletions based on lack of coverage, using lenient single-end mappings (from step 1)
    Returns intervals of consecutive bases of the genome for which there is no read mapped

    Because of spurious mappings, some large deletions may be broken in several smaller intervervals
    TODO: Need to introduce a noise component (e.g. convolution), where I'll accept that a few reads may still be there...

    Because of mappability issues, some deletions may be falsely called, but:
      - Using single-end reads and relaxed mapping should compensate for this effect, but will not eliminate it
      - Comparing an evolved strain against an original should compensate for that effect (as the problem should occurr in both) 

  5 - Calculate GC Normalized Coverage profiles in windows of 250bp
    I first traverse the genome in fixed, non-overlapping windows of 250bp.
    I bin these 250bp windows of 250bp based on the GC content of the reference for that window.
    Then I do a mean coverage per 250bp for each bin, independently
      Note: There is a nice linear correlation between coverage and higher GC 
      Note: This fixed binning may not work if the coverage is low, since some bins may have few or no reads
    Finally, I run the genome with non-overlapping 250bp windows, and for each I divide the absolute coverage by the mean coverage
      of a 250bp window within the same GC bin.
    
    TODO: This only takes GC into account. To more accurately detect duplications in individual strains we need to handle 
    positional bias (e.g. higher coverage near the orig). For this, we can build a model (e.g. a polynomial regression model)
    of coverage change along the genome and normalize based on that model. There are some R packages to do this type of models
    I wasn't particularly convinced with the ones I tried but maybe I didn't used them well

    Note: Comparing coverage between two strains removes (partially) the GC and positional biases.
      The positional bias may still be a problem here, for example, if the two strains have a different division 
      rate (where one strain has a greater positional bias than the other)

    Note: I tried PRISM (http://compbio.cs.toronto.edu/prism) and results were promising
    PRISM may also be useful to detect transpositions etc...


  6 - Calculate possible duplications / deletions, if there is an evolved and an original strain AND if the GC normalized 
  coverage of the original strain has been run previously

    Returns intervals of consecutive 250bp windows of the genome for which there is a difference of more than 0.5 between 
      the normalized coverage of the evolved strain versus the normalized coverage of original strain in the same window

    It is sensitive to noise, and it returns false positives, typically of small 250bp or 500bp intervals.
    Similarly, it may break larger (real) intervals into smaller, almost contiguous intervals. 
    E.g. duplication in 13YFP_events_OYFP_relaxed.tab

    TODO: Build a convolution model to improve signal/noise ratio. There are already some models for this in the literature.
    I wasn't particularly convinced with the ones I tried but maybe I didn't use them well

    TODO: Build a better statistical model to avoid the use of a perfectly arbitrary value (0.5 in this case). This will depend
    on the granularity (window size) and how we can filter out noise using convolution

  # All the pipeline is experimental, but these last steps are particularly experimental and conclusions should be taken with great care

  7 - Perform a denovo assembly using edena

    Edena is very simple to use and works well with bacterial genomes, particularly if coverage is high, or even moderate (~20-30x)
    It does not require many parameters, as it uses an read overlap model
    
    edena -p /tmp/sample_name_edena -paired raw_data/sample_name_Cleandata_1.fq raw_data/sample_name_Cleandata_2.fq
    edena -e /tmp/sample_name_edena.ovl -p /tmp/sample_name_edena -peHorizon 750
    
    The only parameter here is an expected horizon maximum fragment size.
    Edena is not very sensitive to that parameter, so 750 should work for most fragment sizes (200-500).

  8 - Map the denovo assembly against the reference (requiring 100% accuracy) and check for possible genomic events

    I use lastz to do the mappings. Sometimes lastz seem to miss one or two large sequences (still do not know exactly why). 
    These can be filtered when matching against NCBI nr in step 9. Nonetheless, these (few, usually one-two) large sequences 
    can cause false events to be detected

    Moreover, transposons are big confounders in this step, as they may suggest many false large events (because they cross-match). 
    Checking for contigs that completely cover a large event helps eliminating these false positive. To worsen this, usually
    the areas not covered in the denovo assembly are also transposon / repeat areas... so sometimes cannot reliably eliminate the event
  
    Nonetheless, sometimes there are events that are only partially covered by contigs. These are more likely to be true, 
    particularly if the event is only covered in a few bases on the border. 
  
  9 - Separate sequences of the denovo assembly that didn't match the reference and blast them against NCBI nr (uses internet connection)

    It successfully detected the construct used to insert the LacZ
    TODO: this step may also be modified done to work with a local database

  10 - Take unmapped reads (from step 1) and try to remap them against the reference allowing for split reads, to identify possible junction events
  
    One obvious control is end of chromosome against beginning of chromosome, which is a junction that should always appear.
    The number of reads supporting that junction should also be taken as an estimate of what is expected in a "real" junction...

    Transposons are big confounders in this step. E.g. many junctions will appear to happen when a transposon is jumping simply because 
    there are multiple copies of that transposon elsewhere in the genome. Nonetheless, when these multiple events happen, this is an indication
    that something may be happening...
