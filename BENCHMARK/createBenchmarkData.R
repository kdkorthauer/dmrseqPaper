######################################################
### parameters to change to run on your own system ###
######################################################
# change the following to where you would like to download the raw dendritic cell data files from GEO
raw.file.prefix <- data.file.prefix <- "/n/irizarryfs01/kkorthauer/WGBS/DENDRITIC/DATA/RAW/" 
# change the following to where you would like to save the resulting simulated
sim.file.prefix <- "/n/irizarryfs01/kkorthauer/WGBS/DENDRITIC/DATA/BENCHMARK/"
# change the following to the number of DMRs to add
num.dmrs <- 3000
######################################################
###         end of parameters to change            ###
######################################################

library(dmrseq)
library(readr)
library(rtracklayer)

# this will download the Dendritic data from GEO (GSE64177)
system(paste0("wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE64nnn/GSE64177/suppl/GSE64177_RAW.tar' -P ", 
              raw.file.prefix))

# untar and unzip the files
system(paste0("tar -xf ", raw.file.prefix, "GSE64177_RAW.tar -C ", raw.file.prefix))
system(paste0("for f in ", raw.file.prefix, "*.gz; do\n",
  'STEM=$(basename "${f}" .gz)\n',
  'gunzip -c "${f}" > ', raw.file.prefix,'"${STEM}"\n',
  'done\n'))
system(paste0("rm ", raw.file.prefix, "*.gz"))

 # functions to read in txt files to BSseq objects
 read.DC <- function(files, sampleNames, rmZeroCov = FALSE, verbose = TRUE){
	  ## Argument checking
	  if (anyDuplicated(files)){
		  stop("duplicate entries in 'files'")
	  }
	  if (length(sampleNames) != length(files) | anyDuplicated(sampleNames)){
		  stop("argument 'sampleNames' has to have the same length as argument 'files', without duplicate entries")
	  }
	  ## Process each file
	  idxes <- seq_along(files)
	  names(idxes) <- sampleNames
	  allOut <- lapply(idxes, function(ii){
		  if (verbose) {
			  cat(sprintf("[read.DC] Reading file '%s' ... ", files[ii]))
		  }
		  if (file.exists(files[ii])){
		  ptime1 <- proc.time()
   
		  raw <- read.DCFileRaw(thisfile = files[ii])

		  BS <- BSseq(pos = raw$start, chr = raw$chr, M = as.matrix(raw$M, ncol = 1),
							  Cov = as.matrix(raw$Cov, ncol = 1))     
		  sampleNames(BS) <- sampleNames[ii]
		  ptime2 <- proc.time()
		  stime <- (ptime2 - ptime1)[3]
		  if (verbose) {
			  cat(sprintf("done in %.1f secs\n", stime))  
		  }
		  print(paste0("class is: ", class(BS)))
		  return(BS) 
		  }else{
			 cat(sprintf("File missing; moving on to next..."))
		  } 
	  })
	  if (verbose) {
		  cat(sprintf("[read.DC] Joining samples ... "))
	  }
	  ptime1 <- proc.time()
	  OutC <- allOut[[1]]
	  show(OutC)
	  pData(OutC)
	  if (length(allOut)>1){
		for (l in 2:length(allOut)){
		  OutC <- bsseq::combine(OutC, allOut[[l]])
	  }}

	  ptime2 <- proc.time()
	  stime <- (ptime2 - ptime1)[3]
	  if (verbose) {
		  cat(sprintf("done in %.1f secs\n", stime))
	  }
	  OutC
  }

 read.DCFileRaw <- function(thisfile, verbose = TRUE){
   columnHeaders <- c("chr", "start", "M", "Cov")
   what0 <- replicate(length(columnHeaders), character(0))
   names(what0) <- columnHeaders
   int <- which(columnHeaders %in% c("start", "M", "Cov"))
   what0[int] <- replicate(length(int), integer(0))


   ## Read in the file
   if (grepl("\\.gz$", thisfile)) {
	   con <- gzfile(thisfile)
   }else{ 
	   con <- file(thisfile, open = "r")
   }	  
   out <- scan(file = con, what = what0, quote = "", sep="\t",
					 na.strings = "NA", quiet = FALSE)
   close(con)

   # return data.frame
   DataFrame(out)
 }


# This process will be carried out twice to generate two benchmark datasets:
# one with 2 samples per condition and one with 3 samples per condition 
# (since the dendritic cell data has 6 control samples total)
for (sampleSize in c(2,3)){

   # set a different seed for each sample size
   set.seed(848)
   if (sampleSize ==3){set.seed(848*5)}

   # list files & parse sample information from filenames 
   # select a subset of control ("NI" = not infected) files
   files <- list.files(path = raw.file.prefix, pattern = "*NI_5mC.txt")
   sampleNames <- sapply(strsplit(files, "_"), function(x) paste(x[[3]], x[[2]], sep="_"))
   Cell <- sapply(strsplit(files, "_"), function(x) x[[2]])
   selected <- seq_along(Cell)

   # exclude cell DC89 unless using all cells (coverage outlier)
   if (sampleSize == 2 ){
	  rmv <- grep("DC89", Cell)
	  selected <- selected[!(selected %in% rmv)]
   }

   # reorder odd/even in case of sample ordering effects
   selected <- selected[c(which(seq_along(selected)%%2==1), which(seq_along(selected)%%2==0))] 
   selected <- selected[1:(2*sampleSize)]
   files <- files[selected]
   sampleNames <- sampleNames[selected]
   Cell <- Cell[selected]

   # read in text files to a bsseq object and add sample metadata
   DC <- read.DC(paste0(raw.file.prefix, files), sampleNames, verbose = TRUE)
   pData(DC) <- data.frame(Sample=sampleNames, Cell=Cell) 

   # exclude sex chromosomes.
   DC <- chrSelectBSseq(DC, seqnames = paste0("chr", c(seq(1,22))), order = TRUE)

   # filter loci by coverage
   DC <- filterLoci(DC, minCoverage = 1, numSamples="all")
   show(DC)
 
   # simulate DMRS - see simDMRs documentation for details on output format
   simDC <- simDMRs(bs=DC, num.dmrs=num.dmrs)
   rm(DC)
   
   # save it as rds - read in later with readRDS(file)
   saveRDS(simDC, file = paste0(sim.file.prefix, "/sim.data.n", sampleSize, ".rds"))

   # also save it as one compressed txt file per sample
   ct <- 1
   for (samp in sampleNames){
	 gr <- simDC$bs[,sampleNames==samp]
	 df <- data.frame(seqnames=seqnames(gr),
	   starts=start(gr),
	   M=getCoverage(gr, type="M"),
	   U=getCoverage(gr, type="Cov"))

	 fz <- gzfile(paste0(sim.file.prefix, "/sim.data.n", sampleSize, "_", 
						samp, ".sample", ct, ".txt.gz"))
	 write_tsv(df, path=fz, col_names=FALSE)
   }

   # bed file of true DMRs
   export.bed(simDC$gr.dmrs, paste0(sim.file.prefix, "/trueDMRs.n", sampleSize, ".bed"))
   rm(simDC)
   ct <- ct + 1
}