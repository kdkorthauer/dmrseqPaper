#!/bin/bash
#SBATCH -N 1 #Number of nodes
#SBATCH --mail-type=ALL      #Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=kdkorthauer@gmail.com  #Email to which notifications will be sent

#!/bin/bash      

module load bismark

# get human grch38 reference genome and unzip
cd /n/irizarryfs01_backed_up/kkorthauer/ReferenceGenomes/human_bt1/
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

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
# Will use bowtie 1 since these sequences are all single end, and rather short
# (i.e. 50 to 100 basepairs long) and since bowtie 1 is faster

bismark_genome_preparation --bowtie1 /n/irizarryfs01_backed_up/kkorthauer/ReferenceGenomes/human/



