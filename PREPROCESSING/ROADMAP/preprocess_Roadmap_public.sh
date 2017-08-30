#!/bin/bash

# add in all slurm options to sbatch submission command instead of having them 
# buried inside the calling script to make sure that correct partition is used
# and the constraint to use the holylfs file system is satisfied

# query how far along the pipeline is by looking for the presence of key
# files along the way.  This assumes that all fastqs have been obtained through
# fastq-dump (if not, will still re-obtain them from the SRAs, but will do 
# so using a multi-core job and thus will not be efficient since that step
# of the pipeline only uses one core.  
cd /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/kkorthauer/fastq
noBam=()
noContext=()
for f in `cat /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SRR_Acc_List_public_heart_intestine.txt`; do
  cd /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/kkorthauer/fastq
  export f
  if [ $f != "Run_s" ]; then # skip the first line (just a header)
  	# check whether the SRR bam file from mapping exists
  echo $f
  if [[ ! -e $f\_1_trimmed.fq ]]; then
    noBam=("${noBam[@]}" "${f}")
    # if the bam file doesn't exist, check to see if it might be currently running
    isRunning=`showq-slurm -u kkorthauer | grep $f`
    if [[ -z $isRunning ]]; then
      # if it isn't currently running, then generate a job to perform the mapping step
      echo $f " needs fastq-dump or trimming; generating a new job to perform this step"
      sbatch --job-name=${f} \
            -e /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SlurmErr/${f}.err.txt \
            -o /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SlurmOut/${f}.out.txt \
            -t 1000 --mem=4000 -n 2 -N 1 -p general --constraint=holyib \
            --mail-type=ALL --mail-user=kdkorthauer@gmail.com \
            /n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/fastqToBismarkExtractCytosine_public.sh;
    fi	
  elif [[ ! -e $f\_1_trimmed_bismark_bt2.bam ]]; then
  	noBam=("${noBam[@]}" "${f}")
    # if the bam file doesn't exist, check to see if it might be currently running
    isRunning=`showq-slurm -u kkorthauer | grep $f`
    if [[ -z $isRunning ]]; then
      # if it isn't currently running, then generate a job to perform the mapping step
      echo $f " needs mapping; generating a new job to perform mapping step"
      sbatch --job-name=${f} \
            -e /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SlurmErr/${f}.err.txt \
            -o /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SlurmOut/${f}.out.txt \
            -t 6000 --mem=65000 -n 20 -N 1 -p general --constraint=holyib \
            --mail-type=ALL --mail-user=kdkorthauer@gmail.com \
            /n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/fastqToBismarkExtractCytosine_public.sh;
    fi	
  else # get the experiment ID for the current SRR id - loops through each SRA ID 
    readarray runs < /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/ID.table_public.txt
    for i in "${runs[@]}"
    do
      if grep -q "$f" <<< "$i" ; then
        exp=${i:0:12}  # format is index:startpos:length
        exp=$(echo $exp)
        echo "Biosample" $exp
      fi
    done
    # check if the combined bam file, CpG context extraction, and CpG coverage bed file
    # for that biosample id all exist (in the bam subdirectory).  If at least one of these
    # files is missing, then generate an extraction step job. 
    # files from the extraction will be placed in the bam subdirectory to distinguish 
    # from those generated in earlier (possibly truncated steps).
    cd /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/kkorthauer/fastq
  # uncomment the following line when not running the Rscript to get BSseq objects
    if [[ ! -e CpG_context_$exp\.txt.gz ]] || [[ ! -e $exp\.bismark.cov.gz ]] || [[ ! -e /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/kkorthauer/bam/$exp\.bam ]] || [[ ! -e $exp\.cytosine ]]; then  
      # check if CpG context file and coverage files are present
      # otherwise if the bam file is present, look for the CpG_context / coverage file 
      # that is output from the bismark_methylation_extractor step
      noContext=("${noContext[@]}" "${f}")
      echo $exp "Needs extraction"
      # if the CpG context file doesn't exist, check to see if it might be 
      # currently running under this, or one of the other SRR ids that are associated
      # with this experiment ID
      # first grab the other SRR IDs connected to this experiment ID
      # then check running jobs for each one
      isRunning=0
      for i in "${runs[@]}"; do
        if grep -q "$exp" <<< "$i" ; then
          SRRid=${i:13:10}
          SRRid=$(echo $SRRid)		#remove trailing whitespace by reassigning output of echo
          Running=`showq-slurm -u kkorthauer | grep $SRRid`
          if [[ ! -z $Running ]]; then
    	    ((++isRunning))
    	  fi
        fi
      done
      # if there are no running jobs associated with this experiment ID, 
      # then generate a job to carry out the meth extraction step
      if [[ $isRunning -eq 0 ]]; then
        if [[ ! -e CpG_context_$exp\.txt.gz ]] || [[ ! -e /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/kkorthauer/bam/$exp\.bam ]]; then 
          echo $f "job missing; generating a new job to perform extraction step(s)"
          sbatch --job-name=${f} \
            -e /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SlurmErr/${f}.err.txt \
            -o /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SlurmOut/${f}.out.txt \
            -t 6000 --mem=65000 -n 20 -N 1 -p general --constraint=holyib \
            --mail-type=ALL --mail-user=kdkorthauer@gmail.com \
            /n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/fastqToBismarkExtractCytosine_public.sh;
        elif [[ ! -e /n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/kkorthauer/fastq/$exp\.bismark.cov.gz ]]; then 
          echo $f "job missing; generating a new job to perform bismark2bedGraph step(s)"
          sbatch --job-name=${f} \
            -e /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SlurmErr/${f}.err.txt \
            -o /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SlurmOut/${f}.out.txt \
            -t 6000 --mem=65000 -p general --constraint=holyib \
            --mail-type=ALL --mail-user=kdkorthauer@gmail.com \
            /n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/fastqToBismarkExtractCytosine_public.sh;
        elif [[ ! -e $exp\.cytosine ]]; then
          echo $f "job missing; generating a new job to perform cytosine report step"
          sbatch --job-name=${f} \
            -e /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SlurmErr/${f}.err.txt \
            -o /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SlurmOut/${f}.out.txt \
            -t 600 --mem=65000 -n 20 -N 1 -p general --constraint=holyib \
            --mail-type=ALL --mail-user=kdkorthauer@gmail.com \
            /n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/fastqToBismarkExtractCytosine_public.sh;
        else
          echo $f "job missing; generating a new job to generate bsseq RData objects"
          # export the experiment ID variable so it can be accessed by the R script
          export exp
          sbatch --job-name=${f} \
            -e /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SlurmErr/${f}.err.txt \
            -o /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SlurmOut/${f}.out.txt \
            -t 60 --mem=65000 -n 1 \
            R CMD BATCH --quiet --no-restore --no-save /n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/cytosineToBSseq_public.R $exp\.Rout
        fi 
       fi	
    fi
  fi
  fi
  sleep 1
done;