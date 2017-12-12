#!/bin/bash      

######################################################
### parameters to change to run on your own system ###
######################################################
# change the following to the root directory you'd like to download the wgbs data on your system
DATDIR='/n/irizarryfs01/kkorthauer/WGBS/DNMT3A/'
# change the following to where you'd like to save the reference genome assembly on your system
REFDIR='/n/irizarryfs01/kkorthauer/ReferenceGenomes/mouse/mm10/'
######################################################
###         end of parameters to change            ###
######################################################

# use SRA Toolkit to download SRA files from the
# list of SRRs contained in the file 'SRR_Acc_List.txt'

# comment these out or (if necessary) add paths to where these tools are installed on your system
module load sratoolkit
module load bismark

# save current directory (to access url list)
SCRIPTDIR=$PWD

cd $DATDIR
while read p; do
  cd SRA
  echo $p
  export p
  wget ${p}
  cd ..
done <$SCRIPTDIR/SRR_Acc_List_Url.txt

# get mm10 reference genome & unzip
cd $REFDIR
wget ftp://ftp.ensembl.org/pub/release-84/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

# Building the bisulfite genome indexes using Bowtie 2

# Next, we build the bisulfite genome indexes using Bowtie 2 for the human 
# genome (ENSEMBL GRCh38 build). This process can take several hours.

# The bismark_genome_preparation command needs to only be once for the 
# genome of interest for bisulfite alignments. This will create two 
# folders within /n/irizarryfs01_backed_up/kkorthauer/ReferenceGenomes/mouse/mm10/. 
# Within this folder needs to be a .fa or .fasta extension (single or 
# multiple entries per file). The two individual folders within this 
# directory are one for a C->T converted genome and the other one for 
# the G->A converted genome. After creating C->T and G->A versions of
# the genome they will be indexed in parallel using the indexer bowtie-build
# (or bowtie2-build). The --bowtie2 parameter will create bisulfite indexes
# for Bowtie 2 (default is Bowtie 1).
# Bowtie2 is optimized for newer NG sequencers such as Illumina HiSeq 2000
# (what this data was generated using)

bismark_genome_preparation --bowtie2 /n/irizarryfs01_backed_up/kkorthauer/ReferenceGenomes/mouse/mm10/



