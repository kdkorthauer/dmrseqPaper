#!/bin/bash
######################################################
### parameters to change to run on your own system ###
######################################################
# this script calls the runBismark.sh script using the slurm scheduler by default
# If you do not have slurm, you can adapt this script to suit a different batch processing system,
# or change the following line and change to SLURM=F. This will run through the SRAs in serial
SLURM=T
# change the following to the root directory you'd like to download the wgbs data on your system
export DATDIR='/n/irizarryfs01/kkorthauer/WGBS/DNMT3A/'
# change the following to where you'd like to save the reference genome assembly on your system
export REFDIR='/n/irizarryfs01_backed_up/kkorthauer/ReferenceGenomes/mouse/mm10/'
######################################################
###         end of parameters to change            ###
######################################################

export SCRIPTDIR=$PWD
# for each accession number, run the script runBismark:
cd $DATDIR
for f in `cat SRR_Acc_List.txt`; do
	echo $f
  	export f
  	  # submit a batch job 
  	if [[ "$SLURM" =~ "T" ]]; then
  	  sbatch --job-name=${f} \
		-e /n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/DATA/SRA/SlurmErr/${f}.err.txt \
		-o /n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/DATA/SRA/SlurmOut/${f}.out.txt \
 		runBismark.sh;
    else 
      sh runBismark.sh;
    fi
done;