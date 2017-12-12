#!/bin/bash

######################################################
### parameters to change to run on your own system ###
######################################################
# change the following to the root directory you'd like to download the wgbs data on your system
export DATDIR='/n/holylfs/EXTERNAL_REPOS/NCBI/RESTRICTED/kkorthauer/'
# change the following to the path where sratoolkit is installed on your system
export SRATOOLPATH='/n/home06/kkorthauer/sratoolkit.2.8.0-ubuntu64/'
######################################################
###         end of parameters to change            ###
######################################################


# read in the srr run ids
list=`cat SRR_Acc_List_public_heart_intestine.txt`
njobs=5

# use prefetch to download the sra files
for i in $list ; do
if [ $i != "Run_s" ];  # skip the first line (just a header)
	then
	f=$DATDIR/sra/$i.sra
	if [[ ! -e $f ]];
	then
 		echo -n "--> STARTING ASSESSION $i DOWNLOAD @ " ; date
 			$SRATOOLPATH/bin/prefetch --max-size 500000000 $i ./ > $i.log &
 			running=`jobs -r | grep -c Running`
 		while [ "$running" -gt "$njobs" ] ; do
        	sleep 1
        	running=`jobs -r | grep -c Running`
        	done
    fi
fi
done
