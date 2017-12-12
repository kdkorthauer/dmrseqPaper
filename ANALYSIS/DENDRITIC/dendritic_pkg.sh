######################################################
### parameters to change to run on your own system ###
######################################################
# change the following to where the raw data files have been downloaded (see PREPROCESSING folder)
export RAWDATDIR="/n/irizarryfs01/kkorthauer/WGBS/DENDRITIC/DATA/" 
# change the following to where you would like to save results output 
export OUTDIR0="/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/"
# change the following to where you would like to save simulated data 
export SIMDIR="/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/binSim"
# change the following to where metilene is installed on your system
export METPATH="/n/irizarryfs01_backed_up/kkorthauer/softwareTools/metilene_v0.2-7/metilene"
# change the following to TRUE for BSmooth, DSS, and metilene to be run in default mode
export default=TRUE
# change the following to "area" or "avg" to use the respective naive test statistics
export STAT="stat" 
# change the following to F if not using slurm batch processing system
SLURM=T
# uncomment the following to specify number of cores if not using the slurm batch processing system 
#export SLURM_NTASKS=6
######################################################
###         end of parameters to change            ###
######################################################
export JOB=1
for ARGS in 'sampleSize=2 num.dmrs=3000 cond="control"' \
   'sampleSize=2 num.dmrs=0 cond="control"' \
   'sampleSize=2 num.dmrs=3000 cond="infected"' \
   'sampleSize=2 num.dmrs=0 cond="infected"' \
   'sampleSize=3 num.dmrs=3000 cond="control"' \
   'sampleSize=3 num.dmrs=0 cond="control"' \
   'sampleSize=3 num.dmrs=3000 cond="infected"' \
   'sampleSize=3 num.dmrs=0 cond="infected"' \
   'sampleSize=6 num.dmrs=0 cond="infected"' \
   ;
   do
   #if [[  "$JOB" =~ "1"|"5"|"2"|"6" ]]; then
   if [[  "$JOB" =~ "1"|"5" ]]; then
   # add the name of each method to run here (dmrseq, BSmooth, DSS, metilene)
   for METHOD in BSmooth; 
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
    echo ${METHOD} ${JOB} ${ARGS} ${OUTDIR}
    export METHOD JOB ARGS OUTDIR
 	if [ $JOB -ge 5 ]; then
 		echo "Using larger job settings..."
 		if [[ "$SLURM" =~ "T" ]]; then
		  sbatch --job-name=${METHOD}.DCN.${JOB} \
		  -e ${OUTDIR}/Slurm/denN.${JOB}.${METHOD}.err.txt \
		  -o ${OUTDIR}/Slurm/denN.${JOB}.${METHOD}.out.txt \
		  -t 300 --mem=200000 -n 1 -N 1 -p serial_requeue,irizarry \
		  --mail-type=ALL --mail-user=kdkorthauer@gmail.com \
		  R CMD BATCH --quiet --no-restore --no-save "--args $ARGS" dendritic_pkg.R ${OUTDIR}/Slurm/denN.${JOB}.${METHOD}.out.txt
  		  sleep 1 # pause to be kind to the scheduler
  		else
  		  R CMD BATCH --quiet --no-restore --no-save "--args $ARGS" dendritic_pkg.R ${OUTDIR}/Slurm/denN.${JOB}.${METHOD}.out.txt
  		fi
  	else
  		echo "Using smaller job settings..."
  		if [[ "$SLURM" =~ "T" ]]; then
		  sbatch --job-name=${METHOD}.DCN.${JOB}. \
		  -e ${OUTDIR}/Slurm/denN.${JOB}.${METHOD}.err.txt \
		  -o ${OUTDIR}/Slurm/denN.${JOB}.${METHOD}.out.txt \
		  -t 200 --mem=180000 -n 1 -N 1 -p serial_requeue,irizarry \
		  --mail-type=ALL --mail-user=kdkorthauer@gmail.com \
		  R CMD BATCH --quiet --no-restore --no-save "--args $ARGS" dendritic_pkg.R ${OUTDIR}/Slurm/denN.${JOB}.${METHOD}.out.txt
  		  sleep 1 # pause to be kind to the scheduler
  		else
  		  R CMD BATCH --quiet --no-restore --no-save "--args $ARGS" dendritic_pkg.R ${OUTDIR}/Slurm/denN.${JOB}.${METHOD}.out.txt
  		fi
    fi
  done
  fi
  ((JOB++))
done