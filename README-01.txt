This repository currently contains code and scripts for alignment.
When we are done with that, we will use the results in this same repository
for variant calling and other downstream analyses, etc.  Keep this repository
around in your scratch space!

In this exercise, you are going to align reads from multiple FASTQ files
into SAM files, convert those to the binary BAM format, run them through
the samtools fixmate utility to compute insert lengths appropriately, etc,
sort those BAM files into coordinate order, and then merge all the BAM files
from a single individual and library prep and 
and then run the resulting BAM file through samtools markdup to identify
PCR duplicates. Finally, you will use samtools to index that BAM file.

The commands that I used to accomplish these tasks are in three following scripts.
You will need to modify the command lines within them to work on your own system.
I recommend working through the first two interactively to make sure that they
are working, before running script 3 as a batch job.

1. obtain-and-prep-data.sh:   In this script I download the Chinook salmon
genome and its BWA index needed for bwa to align the sequences, and then I also
download the FASTQ files to be used in the exercise.

2. map-N-files-from-K.sh:  This is a shell script that uses the information
in chinook-fastq-meta-data.tsv to cycle over N fastq files (starting from the K-th
file index in chinook-fastq-meta-data.tsv) and for each file it maps it, adds some 
mate information, sorts it in coordinate order.  At the end it merges BAMs from the
same individual and then marks PCR duplicates.  You should always make sure that,
when this is run, all the lanes of any individuals being mapped are included in the
run.  Basically this means that N should always be a multiple of 8 and K should always
be a multiple of 8 (including 8 * 0) plus 1.

3. run-single-job.sh:   A shell script with imbedded SBATCH directives that
runs map-N-files-from-K.sh on the first 128 pairs of FASTQ files listed in
chinook-fastq-meta-data.tsv. (Note that there are 1280 total pairs of FASTQ
files.  We will map the remaining 90% next week.)

Your mission is to review those scripts, understand how they work, and
modify them, if need be, for your own HPCC accounts, then use them to
map the reads in the first 128 pairs of FASTQ files.  (You
will align the remaining FASTQ files after we have
talked about SLURM job arrays.)

The data themselves are reads from 160 Chinook salmon that should (for the most part)
map to Chinook salmon chromosome 32.  This chromosome is about 13.5 Mb long.
Since we are focusing on reads that originally mapped to chromosome 32,
this is only about 1/200 of the total amount of data from the particular Chinooks
salmon sequencing project that generated the original data.
It should be enough to learn a lot about bioinformatics, without
overwhelming the compute resources we have available.

**chinook-fastq-meta-data.tsv**: We have mentioned a few
times a file called chinook-fastq-meta-data.tsv. This
is a TAB-delimited file in the top level of this repository
that contains information about the provenance of the
reads in each of the files.  Briefly, 160 salmon were prepared for sequencing 
(96 on Plate 1 and 64 on Plate 2).  All of these fish were prepared in a single
library prep (named Lib-1) and then DNA from each fish was split across 8 lanes
of a single flowcell on an Illumina Hi-Seq Machine.  The reads from each individual
from each lane are in two different FASTQ files (one for read 1 and the other for
read 2).  Thus, there are 160 x 8 x 2 = 2560 different gzipped fastq files.
The file chinook-fastq-meta-data.tsv holds 1280 rows (plus the header row).
Each of the columns carries some information, and we will use the ID, PU, PL,
SM, and LB columns to populate the read group information in each SAM file.

The files are named like this:

DPCh_plate2_E08_S136_L8_R1.fq.gz

Sample: DPCh_plate2_E08
Sample_index: 136
Lane: 8
Read: 1

But, you don't have to worry about parsing these names, because that is
effectively done for you in chinook-fastq-meta-data.tsv.

Parting words: this is going to call upon everything that we have
learned in the course so far, and I have set this up so that you will
have to look hard to understand what is being done at each step, and you
might have to modify the code to get it to work on your system.  

I strongly encourage everyone to run through each line of code in an
interactive session to make sure that it is working, before trying to
run run-single-job.sh as a BATCH job.  (Test it all first, so
that you don't make it a BOTCH job (to quote Nathan!)).

Please get started early, and ask your peers for help, if needed.
Those that find this relatively straightforward
will benefit from explaining it to others. 

** HANDING IN THE ASSIGNMENT **

This is due Thursday, April 7, before the beginning of class.
You are expected to have modified the three scripts above to
run on your system.  You will commit those three and push them
back up to GitHub.  

Also, you must print the file sizes of the bam files produced into a file called bam_file_sizes.txt
and the SHA1 hashes of the bam files produced into a file call bam_file_sha1s.txt,
and commit both of those files (bam_file_sizes.txt and bam_file_sha1s.txt) and
push them back to GitHub to get credit.

The commands to print the file sizes and MD5 hashes that you use must be these
(executed from the top level of this repository):

du -h mkdup/* > bam_file_sizes.txt
shasum mkdup/* > bam_file_sha1s.txt







