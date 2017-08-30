JOB=1
echo "${JOB}"
export JOB 
OUTDIR0="/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/"
echo "${OUTDIR0}"
export OUTDIR0
DATDIR="/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/binSim_preferential_subset/" 
echo "${DATDIR}"
export DATDIR
default=FALSE
export default
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
 echo ${ARGS}
 export ARGS 
 echo ${JOB}
 export JOB 
 #
 for METHOD in dmrseq; 
 #for METHOD in BSmooth DSS;
 #for METHOD in BSmooth;
 #for METHOD in DSS;
 #for METHOD in metilene 
 	do
 	OUTDIR=$OUTDIR0$METHOD\_pkg_naive_avg
 	#OUTDIR=$OUTDIR0$METHOD\_default
 	#OUTDIR=$OUTDIR0$METHOD #\_default
 	export OUTDIR
 	mkdir -p ${OUTDIR}/Slurm
 	if [[  "$JOB" =~ "1"|"5"|"2"|"6" ]]; then
 	#if [[  "$JOB" =~ "5"|"6" ]]; then
 	echo ${METHOD}
  	export METHOD
 	if [ $JOB -ge 5 ]; then
 		echo "Using larger job settings..."
		sbatch --job-name=${METHOD}.DCN.${JOB} \
		-e ${OUTDIR}/Slurm/denN.${JOB}.${METHOD}.err.txt \
		-o ${OUTDIR}/Slurm/denN.${JOB}.${METHOD}.out.txt \
		-t 1000 --mem=200000 -n 8 -N 1 -p irizarry,serial_requeue \
		--mail-type=ALL --mail-user=kdkorthauer@gmail.com \
		R CMD BATCH --quiet --no-restore --no-save "--args $ARGS" dendritic_pkg_avg.R ${OUTDIR}/Slurm/denN.${JOB}.${METHOD}.out.txt
  		sleep 1 # pause to be kind to the scheduler
  	else
  		echo "Using smaller job settings..."
		sbatch --job-name=dmrseq.DCN.${JOB}${METHOD}. \
		-e ${OUTDIR}/Slurm/denN.${JOB}.${METHOD}.err.txt \
		-o ${OUTDIR}/Slurm/denN.${JOB}.${METHOD}.out.txt \
		-t 400 --mem=180000 -n 10 -N 1 -p irizarry,serial_requeue \
		--mail-type=ALL --mail-user=kdkorthauer@gmail.com \
		R CMD BATCH --quiet --no-restore --no-save "--args $ARGS" dendritic_pkg_avg.R ${OUTDIR}/Slurm/denN.${JOB}.${METHOD}.out.txt
  		sleep 1 # pause to be kind to the scheduler
  	fi
  	fi
  done
  ((JOB++))
done