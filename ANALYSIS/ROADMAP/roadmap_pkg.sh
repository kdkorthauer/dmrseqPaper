JOB=0
OUTDIR0="/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/RESULTS/"
export OUTDIR0
DATDIR="/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/"
echo "${DATDIR}"
export DATDIR
echo "${JOB}"
export JOB 
default=TRUE
export default
for ARGS in 'time="roadmap" tiss1="Cerebellum_tissue" tiss2="Ganglionic_eminence_neurosphere_culture" num.dmrs=0 sampleSize=2' \
'time="roadmap" tiss1="Myoepithelial_cell" tiss2="Luminal_epithelial_cell" num.dmrs=0 sampleSize=2' \
'time="roadmap" tiss1="Cerebellum_tissue" tiss2="Myoepithelial_cell" num.dmrs=0 sampleSize=2' \
'time="roadmap" tiss1="Ganglionic_eminence_neurosphere_culture" tiss2="Luminal_epithelial_cell" num.dmrs=0 sampleSize=2' \
'time="roadmap" tiss1="Right_Ventricle" tiss2="Small_Intestine" num.dmrs=0 sampleSize=2' \
'time="roadmap" tiss1="Right_Ventricle" tiss2="Sigmoid_Colon" num.dmrs=0 sampleSize=2' \
'time="roadmap" tiss1="Sigmoid_Colon" tiss2="Small_Intestine" num.dmrs=0 sampleSize=2' \
'time="roadmap" tiss1="Left_Ventricle" tiss2="Right_Ventricle" num.dmrs=0 sampleSize=2' \
'time="roadmap" tiss1="Left_Ventricle" tiss2="Small_Intestine" num.dmrs=0 sampleSize=2' \
'time="roadmap" tiss1="Left_Ventricle" tiss2="Sigmoid_Colon" num.dmrs=0 sampleSize=2' \
; do
 echo ${ARGS}
 export ARGS 
 ((JOB++))
 echo ${JOB}
 export JOB 
 for METHOD in dmrseq; 
 #for METHOD in BSmooth;
 #for METHOD in DSS;
 #for METHOD in BSmooth DSS;
 #for METHOD in metilene
 	do
 	OUTDIR=$OUTDIR0$METHOD\_pkg
 	#OUTDIR=$OUTDIR0$METHOD\_default/
    echo "${OUTDIR}"
    export OUTDIR
    mkdir -p ${OUTDIR}/Slurm
 	#
 	#if [[ "$JOB" -eq 1 ]];
 	#if [[ "$JOB" =~ "5"|"6"|"7"|"8"|"9"|"10" ]];
 	if [[ "$JOB" =~ "5"|"6"|"7"|"9"|"10" ]];
 	#if [[ "$JOB" =~ "5" ]];
 	then
 	 	echo ${METHOD}
  		export METHOD 
 		sbatch -e ${OUTDIR}/Slurm/${METHOD}.RM${JOB}.txt \
 		-o ${OUTDIR}/Slurm/RMo${JOB}.${METHOD}.txt \
 		--job-name=${METHOD}.RM${JOB} \
 		-n 6 -N 1 -t 120 -p irizarry,serial_requeue \
 		--mem=130000 --mail-type=ALL \
 		--mail-user=kdkorthauer@gmail.com \
 		R CMD BATCH --quiet --no-restore --no-save "--args $ARGS" roadmap_pkg.R ${OUTDIR}/Slurm/RMo${JOB}.${METHOD}.txt 
 		#roadmap_comparison.slurm 
  	sleep 1 # pause to be kind to the scheduler
  	fi
  	done
done