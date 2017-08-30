#!/bin/bash
#SBATCH -N 1 #Number of nodes
#SBATCH --mail-type=ALL      #Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=kdkorthauer@gmail.com  #Email to which notifications will be sent
#SBATCH --job-name=bedToBam
#SBATCH -e /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SlurmErr/bedToBam.err.txt
#SBATCH -o /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SlurmOut/bedToBam.out.txt
#SBATCH -t 1000 
#SBATCH --mem 150000
#SBATCH -p irizarry
#SBATCH -n 8

module load bedtools

cd /n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/RNAseq/

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
	featureCounts -T 8 -g gene_id -a /n/irizarryfs01_backed_up/kkorthauer/ReferenceGenomes/human/gencode.v19.chr_patch_hapl_scaff.annotation.gtf -o $outfile $f
  else
  	echo "count file $outfile already exists"
  fi
done


# R analysis
cat <<EOT >> combineCounts.R
setwd("/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/RNAseq/");
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