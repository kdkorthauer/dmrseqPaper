#!/bin/bash

# assume that the variable f houses the SRR id

# remove these if your system doesn't use modules
# specify path to these tools (if necessary)
module load sratoolkit
module load bismark
module load cutadapt
module load fastqc

# get fastq from sra
cd $DATDIR/fastq

# check if files already exist and don't recreate them if they do
#### UPDATE 8/15/2016:: REDO fastq dump since metadata has been corrected
# and adding the --skip-technical flag to ignore barcode reads from pooling
file1=$f\_1.fastq

fastqs=0

# count how many fastq files there are
if [[ -e $file1 ]]; then
  (( fastqs++ ))
fi

# there should be two.  If not, redo the fastq-dump step.
if [[ $fastqs -ne 1 ]];
then
  fastq-dump -I --split-files --skip-technical ../sra/$f.sra
fi

# if there are any empty files sitting in the current directory, delete them
#find . -type f -size 0 -delete 

# We use the fastqc command to get an idea about the sequencing 
# quality of the .fastq files. Using fastqc can determine if the reads
# contain adapter sequences, overrepresented sequences/contaminants, 
# differences between paired-end reads and library complexity.
# check if output files exist and don't recreate them if they do

file1=$f\_1_fastqc.zip

if [[ ! -e $file1 ]];
then
  fastqc $f\_1.fastq
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

file1=$f\_1_trimmed.fq # double check the name of these generated files
echo $file1

if [[ ! -e $file1 ]];
then
  echo $file1 " missing; beginning trimming" 
  ~/trim_galore_zip/trim_galore --fastqc $f\_1.fastq
  echo "trimming ended" 
else
  echo $file1 " is there; not trimming"  
fi

# this will create trimmed fastq files with names
# $f\_1_val_1.fq
# $f\_2_val_2.fq
# (check if they are already present and don't rerun trim_galore if so
#!/bin/sh

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
# don't carry out this step if the final bam file exists or if there are 
# temporary bam files in the directory (this means the previous attempt was stopped
# short 

# load perl module so can get the M-bias plots
module load perl-modules

file=$DATDIR/fastq/$f\_1_trimmed_bismark_bt2.bam
outdir=$DATDIR/fastq
tempdir=$DATDIR/temp

cd $DATDIR/fastq

# if the bam file exists, then check for temporary files. If present, remove them and the
# main bam file, and then perform the mapping step.
if [[ -e $file ]]; then
  if ls $f\_1_trimmed.fq.temp* 1> /dev/null 2>&1; then
  	echo "Temporary Bismark files exist.  Performing mapping step..."
  	rm $f*temp*    # remove any temporary files
  	rm $f\_1_trimmed_bismark_bt2.bam  # remove the bam file that is likely incomplete
  	$BISMARKPATH/bismark --bowtie2 --multicore 4 --bam -o $outdir --temp_dir $tempdir $REFDIR --se $DATDIR/fastq/$f\_1_trimmed.fq 
  fi  #removed --multicore 4 option for faster testing
else # if the main bam file does not exist, then need to carry it out
    echo "No bam file found.  Performing mapping step..."
    rm $f*temp*   # flush of temporary bam files... 
	$BISMARKPATH/bismark --bowtie2 --multicore 4 --bam -o $outdir --temp_dir $tempdir $REFDIR --se $DATDIR/fastq/$f\_1_trimmed.fq 
fi # removed --multicore 4 option for faster testing

# Each sample will take many hours (~10-20 hrs depending on the file size). 
# To assess the quality of the alignment, look at the Bismark mapping reports. 

# For example, the report for the SRR949215 can be viewed using cat.

# cat SRR949215_1_val_1.fq_bismark_bt2_PE_report.txt

# have to do combine the bam files from the same biosamples across different runs #
# samtools merge out.bam in1.bam in2.bam in3.bam

cd $DATDIR
readarray runs < $SCRIPTDIR/ID.table_public.txt

# loop through array to get biosample id

for i in "${runs[@]}"
do
   if grep -q "$f" <<< "$i" ;
   then
     exp=${i:0:12}  # format is index:startpos:length
     exp=$(echo $exp)
     echo "Biosample" $exp
   fi
   # do whatever on $i
done

bams=()
# loop through array to get all SRR ids corresponding to that experiment
for i in "${runs[@]}"
do
   if grep -q "$exp" <<< "$i" ;
   then
     bams=("${bams[@]}" "$(echo ${i:13:10})")
   fi
   # do whatever on $i
done

echo "${bams[@]}"

cd $DATDIR/fastq
ALLTHERE=TRUE
for i in "${bams[@]}"
do
   thisfile=$DATDIR/fastq/$i\_1_trimmed_bismark_bt2.bam
   if [[ ! -e $thisfile ]];
   then
     ALLTHERE=FALSE
     echo $thisfile
     echo $ALLTHERE
   fi
done
	
# prepend and append the bams array with file prefix and file suffix
# uncombined bams are in the fastq directory
bams=( "${bams[@]/%/_1_trimmed_bismark_bt2.bam}" )
bams=( "${bams[@]/#/$DATDIR/fastq/}" )
	
#Assessing the alignment

#The module bismark2report will generate an HTML report from the Bismark alignment,
# deduplication (if applicable), methylation extraction and M-bias reports. 
# The purpose is to get an idea of the experiment worked. The command attempts 
# to detect all relevant files in the current working directory 
# (but all files can be optionally listed here).
#file=$f\_1_val_1.fq_bismark_bt2_pe.bam
#reportfile=$f\_1_val_1.fq_bismark_bt2_PE_report.html


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

# coverage (M + U) files $f\_1_val_1.fq_bismark_bt2_pe.bismark.cov.gz

if [ $ALLTHERE = TRUE ];
then
  file=$DATDIR/bam/$exp\.bam
  if [[ ! -e $file ]];  # check for presence of combined bam file
  then
  	echo "Combining individual bam files..."
    samtools cat -o $file ${bams[@]}
    ls -l $file
  fi	
  cd $DATDIR/fastq
  if [[ ! -e $DATDIR/fastq/CpG_context_$exp\.txt.gz  ]];
  then
    echo "Starting extraction step..."
    rm $DATDIR/fastq/CpG_context_$exp\*.methXtractor.temp
    $BISMARKPATH/bismark_methylation_extractor -s --no_overlap --comprehensive --multicore 6 --buffer_size 50G --counts --bedGraph --gzip $file

  fi  # generates $exp\.bismark.cov which is used in the following step to create the input for read.bismark
  if [[ ! -e $DATDIR/fastq/$exp\.bismark.cov.gz ]];
  then
    echo "Starting bedGraph step..."
    rm $DATDIR/fastq/CpG_context_$exp\*.methXtractor.temp
    $BISMARKPATH/bismark2bedGraph --buffer_size 60G -o $exp\.bedGraph.gz $DATDIR/fastq/CpG_context_$exp\.txt.gz
  fi  # generates $exp\.bismark.cov which is used in the following step to create the input for read.bismark

  file2=$DATDIR/fastq/$exp\.cytosine
  if [[ ! -e $file2 ]]; then
    echo "Starting cytosine step..."
  	$BISMARKPATH/coverage2cytosine --merge_CpG --genome_folder $REFDIR -o $file2 $DATDIR/fastq/$exp\.bismark.cov.gz 
  fi # update 1/19/17 - use patched version of Bismark from Felix to handle scaffolds robustly
  # run the R script to read in the .cytosine files and convert into Bsseq R objects,
  # and save them in the project DATA directory for further use
  # export the environmental variable that houses the experiment ID so that the R 
  # script can access it  
  export exp
  cd $SCRIPTDIR
  R CMD BATCH "cytosineToBSseq_public.R" $exp\.Rout
fi