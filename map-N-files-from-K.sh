

# This is a script that takes two positional parameters N and K.
# It assumes that it is being launched from a directory that contains
# a TAB-delimited file called chinook-fastq-meta-data.tsv that gives a
# number, a file name prefix, and read group information for each of the
# pairs of FASTQ files.


# Check to make sure there are two appropriate positional parameters,
# and print some usage information
if [ $# -ne 2 ]; then
  echo "Wrong number of parameters"
  echo "Syntax "
  echo "   map-N-files-from-K.sh  N  K"
  echo "where:"
  echo "  N is the number of file pairs to align"
  echo "  K is the file pair number to start on"
  echo "This assumes there is a file named chinook-fastq-meta-data.tsv"
  echo "that holds information about each file pair, in the current"
  echo "working directory."
  echo "Example:
    map-N-files-from-K.sh  32  65
would align file pairs starting at 65 and going up to 96
    "
    
    exit 1;
fi


# set up your bioinf conda environemnt
source ~/.bashrc
conda activate bioinf


# Now, get the first file index and the last file index
START=$2
STOP=$(($START + $1 - 1))

# loop over those file indexes and do the alignment steps on all of them.
# NOTE: If you want to test the lines inside the loop just set Idx=5 (for example)
# on the command line and then run through the lines within the for loop.
for((Idx=$START; Idx<=$STOP; Idx++)); do

  echo "Starting work on file index $Idx at $(date)"
  
  # This awk script extracts the appropriate line from the chinook-fastq-meta-data.tsv
  # file and makes a command line that sets shell variables named after the column
  # headers in the file to their appropriate values:
  ASSIGNMENTS=$(awk -v LINE=$Idx '
    $1 == "index" {for(i=1; i<=NF; i++) vars[i]=$i; next}
    $1 == LINE {for(i=1; i<=NF; i++) printf("%s=%s; ", vars[i], $i)}
  ' chinook-fastq-meta-data.tsv)
  
  # then we have to eval the assignments in that command line:
  eval $ASSIGNMENTS
  
  # make some directories for output (if not there already)
  # the -p says "don't bark a warning if the directory already exists."
  mkdir -p bam
  mkdir -p stderr
  
  # Define variables for file names
  READ2=fastqs/${file_prefix}2.fq.gz
  READ1=fastqs/${file_prefix}1.fq.gz
  GENOME=genome/GCA_002872995.1_Otsh_v1.0_genomic.fna.gz
  FIXEDMATES=bam/${file_prefix}_FIXED.bam
  MARKED=bam/${file_prefix}_mkdup.bam
  BWA_STDERR=stderr/bwa_stderr_$file_prefix
  SAM_VIEW_STDERR=stderr/samtools_view_stderr_$file_prefix
  SAM_FIX_STDERR=stderr/samtools_fixmate_stderr_$file_prefix
  SAM_SORT_STDERR=stderr/samtools_sort_stderr_$file_prefix
  SAM_MARK_STDERR=stderr/samtools_markdup_stderr_$file_prefix
  
  # make a variable for the read-group string
  RGString="@RG\tID:$ID\tSM:$SM\tLB:$LB\tPU:$PU\tPL:$PL"
  
  # now, map, convert to bam, fix mates, sort in coordinate order, mark dupes, and
  # finally remove the intermediate fixed-mate file
  (bwa mem -R $RGString $GENOME $READ1 $READ2 2> $BWA_STDERR | \
    samtools view -b -1 - 2> $SAM_VIEW_STDERR | \
    samtools fixmate -m -O BAM - $FIXEDMATES 2> $SAM_FIX_STDERR) && \
  (samtools sort $FIXEDMATES 2> $SAM_SORT_STDERR | \
    samtools markdup - $MARKED) &&
  rm -f $FIXEDMATES &&
  echo "Done working on file index $Idx at $(date)" # this only gets printed if it was successful
  
done


