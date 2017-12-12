#!/bin/bash

######################################################
### parameters to change to run on your own system ###
######################################################
# change the following to the root directory you'd like to download the RNAseq data on your system
export DATDIR='/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/RNAseq/'
# change the following to where you'd like to save the reference genome assembly on your system
export REFDIR='/n/irizarryfs01/kkorthauer/ReferenceGenomes/human/'
# change the following the number of cores to use
export CORES=8
######################################################
###         end of parameters to change            ###
######################################################

SCRIPTDIR=$PWD
cd $DATDIR

# first download the bed file data from the selected tissue types from SRP000941(ftp address ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/roadmapepigenomics/) 

wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/roadmapepigenomics/by_experiment/mRNA-Seq/heart_left_ventricle/*.bed.gz'
wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/roadmapepigenomics/by_experiment/mRNA-Seq/heart_right_ventricle/*.bed.gz'
wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/roadmapepigenomics/by_experiment/mRNA-Seq/small_intestine/*.bed.gz'
wget 'ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/roadmapepigenomics/by_experiment/mRNA-Seq/sigmoid_colon/*.bed.gz'

# unzip the files
gunzip *gz 

# download the genome annotation files
cd $REFDIR
wget 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz'

# unzip the file
gunzip *gz 

# remove the following if you don't use a module system and (if necessary) specify the path of bedtools on your system
module load bedtools

cd $DATDIR
for f in *.bed.gz
do    
  outfile=$f\.bam 
  if [[ ! -e $outfile ]];
  then
	echo "converting file - $f to bam"
	bedToBam -i $f -g hg19.txt > $outfile
  else
  	echo "bam file $outfile already exists"
  fi
done

module load subread

# loop through bam files found and count features with featureCounts
for f in *.bam
do    
  fout=${f%%???????????}
  outfile=$fout\.counts.txt
  if [[ ! -e $outfile ]];
  then
	echo "counting features of - $fout"
	featureCounts -T $CORES -g gene_id -a $REFDIR/gencode.v19.chr_patch_hapl_scaff.annotation.gtf -o $outfile $f
  else
  	echo "count file $outfile already exists"
  fi
done


# R analysis
cat <<EOT >> combineCounts.R
dat.dir <- Sys.getenv("DATDIR")
setwd(paste0(dat.dir));
# list all files that have "UCSD" in the name and end in ".counts.txt";
# but do not contain ".summary" ;
count.files <- list.files(path = ".", pattern = "UCSD");
count.files <- count.files[count.files %in% list.files(path = ".", pattern = ".counts.txt")];
count.files <- count.files[!(count.files %in% list.files(path = ".", pattern = "summary"))];

# build count table (combine all samples into one table);
counts <- NULL;
it <- 1;
for (f in count.files){;
  tab <- read.table(file=f, header = TRUE, stringsAsFactors=FALSE);
    if (it > 1){;
      if (all.equal(counts$Geneid, tab$Geneid)){;
        counts <- cbind(counts, tab[,ncol(tab)]);
        colnames(counts)[ncol(counts)] <- colnames(tab)[ncol(tab)];
      }else{;
        stop("error: features not identical between count matrices");
      };
     }else{;
       counts <- tab;
     };
    it <- it + 1;
    rm(tab);
};

write.table(counts, file="UCSD.mRNA-Seq.counts.txt", quote=FALSE, 
   row.names=FALSE, sep="\t")
EOT

R CMD BATCH combineCounts.R