######################################################
### parameters to change to run on your own system ###
######################################################
# change the following to where the data files are stored (see PREPROCESSING folder)
export DATDIR="/n/irizarryfs01/kkorthauer/WGBS/DNMT3A/"
# change the following to where you would like to save results output 
export OUTDIR0="/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/RESULTS/"
# change the following to the name of the file containing the RNA-seq counts (see PREPROCESSING folder)
export EXPDAT="/n/irizarryfs01/kkorthauer/WGBS/DNMT3A/RNAseq/GSE61969_cufflinks.gene.counts.txt"
# change the following to where metilene is installed on your system
export METPATH="/n/irizarryfs01_backed_up/kkorthauer/softwareTools/metilene_v0.2-7/metilene"
# change the following to TRUE if BSmooth, DSS and metilene should be run in default mode
export default=FALSE
# change the following to "area" or "avg" to use the respective naive test statistics
export STAT="stat" 
# change the following to F if not using slurm batch processing system
SLURM=T
# uncomment the following to specify number of cores if not using the slurm batch processing system 
#export SLURM_NTASKS=6
######################################################
###         end of parameters to change            ###
######################################################
JOB=0
for ARGS in 'time="dnmt3a" tiss1="KO_FLT3" tiss2="WT_FLT3" num.dmrs=0 sampleSize=2' \
'time="dnmt3a" tiss1="KO_FLT3" tiss2="WT_WTFL" num.dmrs=0 sampleSize=2' \
'time="dnmt3a" tiss1="WT_FLT3" tiss2="WT_WTFL" num.dmrs=0 sampleSize=2' \
; do
 ((JOB++))
 # add the name of each method to run here (dmrseq, BSmooth, DSS, metilene)
 for METHOD in dmrseq;
 	do
 	if [[ "$METHOD" =~ "dmrseq" ]]; then
 	  OUTDIR=$OUTDIR0$METHOD\_pkg
 	  if [[ "$STAT" =~ "avg"|"area" ]]; then
 	    OUTDIR=$OUTDIR\_naive_$STAT
 	  fi
 	else 
 	  OUTDIR=$OUTDIR0$METHOD
 	fi  
 	if [[ "$default" =~ "TRUE" ]] && [[ "$METHOD" =~ "metilene"|"DSS"|"BSmooth" ]]; then
 	  OUTDIR=$OUTDIR\_default
    fi
    mkdir -p ${OUTDIR}/Slurm
 	#if [[ "$JOB" -eq 2 ]]; then
        echo ${METHOD} ${JOB} ${ARGS} ${OUTDIR}
        export METHOD JOB ARGS OUTDIR
        if [[ "$SLURM" =~ "T" ]]; then
 		  sbatch -e ${OUTDIR}/Slurm/${METHOD}.DNMTe${JOB}.txt \
 		  -o ${OUTDIR}/Slurm/DNMTo${JOB}.${METHOD}.txt \
 		  --job-name=${METHOD}.Leuk.${JOB} -n 1 -N 1 -t 300 \
 		  -p irizarry,serial_requeue --mem=130000 --mail-type=ALL \
 		  --mail-user=kdkorthauer@gmail.com \
 		  R CMD BATCH --quiet --no-restore --no-save "--args $ARGS" dnmt3a_pkg.R ${OUTDIR}/Slurm/DNMTo${JOB}.${METHOD}.txt
 		else
 		  R CMD BATCH --quiet --no-restore --no-save "--args $ARGS" dnmt3a_pkg.R ${OUTDIR}/Slurm/DNMTo${JOB}.${METHOD}.txt
 		fi 
  	sleep 1 # pause to be kind to the scheduler
  	#fi
  	done
done