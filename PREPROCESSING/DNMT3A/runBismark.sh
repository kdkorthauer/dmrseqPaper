#!/bin/bash
#SBATCH -n 8 #Number of cores
#SBATCH -N 1 #Number of nodes
#SBATCH -t 500 #Runtime in minutes
#SBATCH -p irizarry,serial_requeue #Partition to submit to
#SBATCH --mem=125000   #125000 # MB
#SBATCH --mail-type=ALL      #Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=kdkorthauer@gmail.com  #Email to which notifications will be sent

# assume that the variable f houses the SRR id

# comment these out or (if necessary) add paths to where these tools are installed on your system
module load sratoolkit
module load bismark
module load cutadapt

# get fastq from sra
cd $DATDIR/fastq

# check if files already exist and don't recreate them if they do
file1=$f\_1.fastq
file2=$f\_2.fastq
if [[ ! -e $file1 ]] && [[ ! -e $file2 ]];
then
  fastq-dump -I --split-files ../SRA/$f
fi
 
# We use the fastqc command to get an idea about the sequencing 
# quality of the .fastq files. Using fastqc can determine if the reads
# contain adapter sequences, overrepresented sequences/contaminants, 
# differences between paired-end reads and library complexity.
# check if output files exist and don't recreate them if they do

file1=$f\_1_fastqc.zip
file2=$f\_2_fastqc.zip

if [[ ! -e $file1 ]];
then
  fastqc $f\_1.fastq
fi

if [[ ! -e $file2 ]];
then
  fastqc $f\_2.fastq
fi

# Look at QC report produced by fastqc

# firefox SRR949215_1_val_1_fastqc/fastqc_report.html &

# Adapter and quality trimming
# We use trim_galore to trim the reads to improve mapping efficiency 
# and reduce the chance of misalignments. Use the --paired parameter 
# to specify that this was a paired-end experiment. This will remove 
# base calls with a Phred score of 20 or lower, removes adapter sequences, 
# removes sequence pairs where either read became too short as a result 
# of trimming (less than 20 bp) and keeps paired-end sequence pairs in sync.

file1=$f\_1_val_1.fq
file2=$f\_2_val_2.fq

if [[ ! -e $file1 ]] && [[ ! -e $file2 ]];
then
  ~/trim_galore_zip/trim_galore --paired $f\_1.fastq $f\_2.fastq
fi

# this will create trimmed fastq files with names
# $f\_1_val_1.fq
# $f\_2_val_2.fq
# (check if they are already present and don't rerun trim_galore if so

#Bismark read aligner

# The paired-end reads for each sample are aligned with Bowtie 2 because it 
# scaled better for longer reads (i.e. from Illumina Hiseq 2000). 
# The parameter --multicore is used to speed up the mapping. 
# Here we use 12 cores. The parameter --bam writes the output file to a 
# BAM format instead of an uncompressed SAM file. 
# The folder /n/irizarryfs01_backed_up/kkorthauer/ReferenceGenomes/mouse/mm10/ 
# contains the bisulfite genome indexes for Bowtie 2 (--bowtie2). 
# The -1 and -2 parameters correspond to the paired-end .fastq files that
# have been trimmed.

file=$f\_1_val_1.fq_bismark_bt2_pe.bam

if [[ ! -e $file ]];
then
  bismark --multicore 10 --bowtie2 --bam $REFDIR -1 $f\_1_val_1.fq -2 $f\_2_val_2.fq
fi

# Each sample will take many hours (~10-20 hrs depending on the file size). 
# To assess the quality of the alignment, look at the Bismark mapping reports. 

# For example, the report for the SRR949215 can be viewed using cat.

# cat SRR949215_1_val_1.fq_bismark_bt2_PE_report.txt

# have to do combine the bam files from the same biosamples across different runs #
# samtools merge out.bam in1.bam in2.bam in3.bam

cd $DATDIR
readarray runs < $SCRIPTDIR/SraRunTable.txt

# loop through array to get experiment id

for i in "${runs[@]}"
do
   if grep -q "$f" <<< "$i" ;
   then
     exp=${i:13:10}
   fi
   # do whatever on $i
done

bams=()
# loop through array to get all SRR ids corresponding to that experiment
for i in "${runs[@]}"
do
   if grep -q "$exp" <<< "$i" ;
   then
     bams=("${bams[@]}" "${i:34:10}")
   fi
   # do whatever on $i
done

echo "${bams[@]}"

cd $DATDIR/fastq
ALLTHERE=TRUE
for i in "${bams[@]}"
do
   thisfile=$i\_1_val_1.fq_bismark_bt2_pe.bam
   if [[ ! -e $thisfile ]];
   then
     ALLTHERE=FALSE
     echo $thisfile
     echo $ALLTHERE
   fi
done
	
#Assessing the alignment

#The module bismark2report will generate an HTML report from the Bismark alignment,
# deduplication (if applicable), methylation extraction and M-bias reports. 
# The purpose is to get an idea of the experiment worked. The command attempts 
# to detect all relevant files in the current working directory 
# (but all files can be optionally listed here).
file=$f\_1_val_1.fq_bismark_bt2_pe.bam
reportfile=$f\_1_val_1.fq_bismark_bt2_PE_report.html

cd $DATDIR/fastq
if [[ -e $file ]] & [[ ! -e $reportfile ]];
then
  bismark2report --alignment_report $f\_1_val_1.fq_bismark_bt2_PE_report.txt
fi	
	
# Extracting methylation calls

# The cytosine methylation calls are extracted using the bismark_methylation_extractor 
# command for every single C analysed. Several output files will be produced including 
# strand and context specific cytosine output files (CpG, CHG or CHH). 
# The --bedGraph parameter generates coverage files which are useful as input 
# for the R/Bioconductor packge bsseq later on. In addition, an M-bias report 
# will be generated and an overall count report. 
# The command bismark_methylation_extractor will detect if the input BAM file
# is for single-end or paired-end (and will set the -p --no_overlap parameter 
# for paired-end in our case). The parameter --gzip returns gzipped context 
# specific cytosine output files.

# coverage (M + U) files $f\_1_val_1.fq_bismark_bt2_pe.bismark.cov
#srr=$f\_1_val_1.fq_bismark_bt2_pe.bismark.cov
#export srr

# need to add the generation of the cytosine report, since this preserves the strand information
# which is necessary to collapse CpGs on opposite strands (if just using the CpG coverage files
# CpGs on neighboring strands will be treated as independent CpG sites -> artificially 
# lowers coverage and increases CpG density

# only carry out this last step if this is the representative SRR (only one per experiment)
REP=FALSE
if [ $f = SRR3457812 ] || [ $f = SRR3457817 ] || [ $f = SRR3457822 ] || [ $f = SRR3457829 ] || [ $f = SRR3457830 ] || [ $f = SRR3457835 ];
then
  REP=TRUE
  echo $REP
fi 

if [ $ALLTHERE = TRUE ] && [ $REP = TRUE ];
then
  file=$exp\.bam
  if [[ ! -e $file ]] & [[ ! -e CpG_context_$exp\.txt.gz ]];
  then
    samtools cat -o $exp\.bam ${bams[@]/%/\_1_val_1.fq_bismark_bt2_pe.bam}
  fi	
  cd $DATDIR/fastq
  if [[ ! -e CpG_context_$exp\.txt.gz  ]];
  then
    bismark_methylation_extractor -p --no_overlap --comprehensive --multicore 8 --buffer_size 50G --counts --bedGraph --gzip $exp\.bam 
  fi
  if [[ ! -e $exp\.cytosine ]];
  then
    coverage2cytosine --merge_CpG --genome_folder $REFDIR -o $exp\.cytosine $exp\.bismark.cov 
  fi
  # convert the coverage file to a bsseq RData object and save it in one-up level directory
  #
  echo 
  export exp
  R CMD BATCH "cytosineToBSseq.R" $exp\.Rout
fi

