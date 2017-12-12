######################################################
### parameters to change to run on your own system ###
######################################################
# change the following to where the data files are stored (see PREPROCESSING folder)
export DATDIR="/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/"
# change the following to where you would like to save results output 
export OUTDIR0="/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/RESULTS/"
# change the following to the name of the file containing the RNA-seq counts (see PREPROCESSING folder)
export EXPDAT="/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/RNAseq/UCSD.mRNA-Seq.counts.txt"
# change the following to where metilene is installed on your system
export METPATH="/n/irizarryfs01_backed_up/kkorthauer/softwareTools/metilene_v0.2-7/metilene"
# change the following to TRUE if BSmooth, DSS and metilene should be run in default mode
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
JOB=0
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
 ((JOB++))
 export JOB 
  # add the name of each method to run here (dmrseq, BSmooth, DSS, metilene)
 for METHOD in BSmooth DSS metilene;
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
 	#if [[ "$JOB" -eq 5 ]];
 	if [[ "$JOB" =~ "5"|"6"|"7"|"8"|"9"|"10" ]];
 	then
 	    echo ${METHOD} ${JOB} ${ARGS}
        export METHOD JOB ARGS OUTDIR
        if [[ "$SLURM" =~ "T" ]]; then
 		  sbatch -e ${OUTDIR}/Slurm/${METHOD}.RM${JOB}.txt \
 		  -o ${OUTDIR}/Slurm/RMo${JOB}.${METHOD}.txt \
 		  --job-name=${METHOD}.RM${JOB} \
 		  -n 1 -N 1 -t 500 -p shared,irizarry,serial_requeue \
 		  --mem=130000 --mail-type=ALL \
 		  --mail-user=kdkorthauer@gmail.com \
 		  R CMD BATCH --quiet --no-restore --no-save "--args $ARGS" roadmap_pkg.R ${OUTDIR}/Slurm/RMo${JOB}.${METHOD}.txt 
 		else 
 		  R CMD BATCH --quiet --no-restore --no-save "--args $ARGS" roadmap_pkg.R ${OUTDIR}/Slurm/RMo${JOB}.${METHOD}.txt
 		fi
  	sleep 1 # pause to be kind to the scheduler
  	fi
  	done
done