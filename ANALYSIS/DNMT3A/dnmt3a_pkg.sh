JOB=0
OUTDIR0="/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/RESULTS/"
export OUTDIR0
DATDIR="/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/DATA/"
echo "${DATDIR}"
export DATDIR
echo "${JOB}"
export JOB 
default=TRUE
export default
for ARGS in 'time="dnmt3a" tiss1="KO_FLT3" tiss2="WT_FLT3" num.dmrs=0 sampleSize=2' \
'time="dnmt3a" tiss1="KO_FLT3" tiss2="WT_WTFL" num.dmrs=0 sampleSize=2' \
'time="dnmt3a" tiss1="WT_FLT3" tiss2="WT_WTFL" num.dmrs=0 sampleSize=2' \
; do
 echo ${ARGS}
 export ARGS 
 ((JOB++))
 echo ${JOB}
 export JOB 
 for METHOD in dmrseq;
 #for METHOD in DSS BSmooth;
 #for METHOD in DSS;
 #for METHOD in metilene;
 	do
 	OUTDIR=$OUTDIR0$METHOD\_pkg
 	#OUTDIR=$OUTDIR0$METHOD\_default/
    echo "${OUTDIR}"
    export OUTDIR
    mkdir -p ${OUTDIR}/Slurm
 	#
 	if [[ "$JOB" -eq 1 ]];
 	then
 	 	echo ${METHOD}
  		export METHOD
 		sbatch -e ${OUTDIR}/Slurm/${METHOD}.DNMTe${JOB}.txt \
 		-o ${OUTDIR}/Slurm/DNMTo${JOB}.${METHOD}.txt \
 		--job-name=${METHOD}.Leuk.${JOB} -n 8 -N 1 -t 75 \
 		-p irizarry,serial_requeue --mem=130000 --mail-type=ALL \
 		--mail-user=kdkorthauer@gmail.com \
 		R CMD BATCH --quiet --no-restore --no-save "--args $ARGS" dnmt3a_pkg.R ${OUTDIR}/Slurm/DNMTo${JOB}.${METHOD}.txt
 		#dnmt3a_comparison.slurm
  	sleep 1 # pause to be kind to the scheduler
  	fi
  	done
done