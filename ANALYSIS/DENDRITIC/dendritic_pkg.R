source("../summarizeResultsFunctions.R")
raw.file.prefix <- data.file.prefix <- Sys.getenv("RAWDATDIR")
result.file.prefix <- Sys.getenv("OUTDIR")
sim.file.prefix <- Sys.getenv("SIMDIR")
metilene.path <- Sys.getenv("METPATH")
print(result.file.prefix) # make sure this is getting set correctly
print(raw.file.prefix)
print(sim.file.prefix)
ifelse(!dir.exists(result.file.prefix), dir.create(file.path(result.file.prefix)), FALSE) 
job <- as.numeric(Sys.getenv("JOB"))
print(job)
set.seed(848*job)

library(dmrseq)
STAT = Sys.getenv("STAT")

if (!file.exists(paste0(raw.file.prefix, "DC_BSSeq.RData"))){
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
	 # hack the code in 'read.bismarkFileRaw' function of bsseq to 
	 # read in ENCODE's pipeline output
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


	# list files
	files <- list.files(path = raw.file.prefix, pattern = "*.txt")
	if (length(grep("df", files))>0){ files <- files[-grep("df", files)] }
	sampleNames <- sapply(strsplit(files, "_"), function(x) paste(x[[3]], x[[2]], sep="_"))

	Infected <- sapply(strsplit(files, "_"), function(x) x[[3]])
	Cell <- sapply(strsplit(files, "_"), function(x) x[[2]])


	# read in one by one 
	# 4 columns: chr pos M M+U
	# need to keep only positions in common across samples
	DC <- read.DC(paste0(raw.file.prefix, files), sampleNames, verbose = TRUE)
	pData(DC) <- data.frame(Sample=sampleNames, Infected=Infected, Cell=Cell) 
	save(DC, file=paste0(raw.file.prefix, "DC_BSSeq.RData")) 

}

# First read in the arguments listed at the command line
args=(commandArgs(TRUE))

# args is now a list of character vectors
# First check to see if arguments are passed.
# Then cycle through each element of the list and evaluate the expressions.
  if (length (args) == 0) {
    stop("No arguments supplied.")
  } else {
    for (i in 1:length(args)) {
 	   eval (parse (text = args[[i]] ))
    }
  }
print(sampleSize)
print(num.dmrs)
print(cond)
workers <- as.numeric(Sys.getenv("SLURM_NTASKS"))
# register cores with #workers
library("BiocParallel")
register(MulticoreParam(workers))

print(workers)
time <- "dendritic"
time2 <- NULL

METHOD <- Sys.getenv("METHOD")
print(METHOD)

default <- as.logical(Sys.getenv("default"))
print(default)

 
# list files
files <- list.files(path = raw.file.prefix, pattern = "*5mC.txt")
sampleNames <- sapply(strsplit(files, "_"), function(x) paste(x[[3]], x[[2]], sep="_"))

Infected <- sapply(strsplit(files, "_"), function(x) x[[3]])
Cell <- sapply(strsplit(files, "_"), function(x) x[[2]])

# read in bsseq object
load(paste0(raw.file.prefix, "DC_BSSeq.RData"))

# exclude sex chromosomes.
DC <- chrSelectBSseq(DC, seqnames = paste0("chr", c(seq(1,22))), order = TRUE)

pData(DC) <- data.frame(Sample=sampleNames, Infected=Infected, Cell=Cell) 
print(show(DC))

 meth.levels.raw = getMeth(DC, type = "raw")
 nlocs <- nrow(meth.levels.raw)
 print(nlocs)
 
 tiss1 <- "NI" 
 tiss2 <- "MTB"
 tissue1 = grep(tiss1, pData(DC)$Infected) 
 tissue2 = grep(tiss2, pData(DC)$Infected)
 
 # exclude DC89 unless using all cells
 if (sampleSize == 2 | (cond == "infected" & sampleSize == 3)){
 	rmv <- grep("DC89", pData(DC)$Cell)
 	tissue1 <- tissue1[!(tissue1 %in% rmv)]
 	tissue2 <- tissue2[!(tissue2 %in% rmv)]
 }
 
 # compare two cells v two cells (both not infected)
 if (cond=="infected"){
 	tiss <- c(tissue1[1:sampleSize], 
  			  tissue2[(length(tissue2)-sampleSize+1):(length(tissue2))]) 
 	
 }else{
 	tiss <- tissue1
 	tiss <- tiss[c(1,3,5,2,4,6)]
 }
 tissue1 <- tiss[1:sampleSize]
 tissue2 <- tiss[(sampleSize+1):(2*sampleSize)]

 tissue = as.factor(c(rep(tiss1, sampleSize), rep(tiss2, sampleSize)))
 design = model.matrix(~tissue)
 
 print(nrow(meth.levels.raw))
 no.hits = which(is.na(rowMeans(meth.levels.raw[, c(tissue1, tissue2)])) == TRUE)
 str(no.hits)
 nrow(meth.levels.raw) - length(no.hits)

 nlocs = min(nlocs, (length(DC) - length(no.hits))) # nlocs orig specified in commandarg
  
 bs = DC[-no.hits][1:nlocs, c(tissue1, tissue2)]
 rm(DC); gc()
 meth.levels = meth.levels.raw[-no.hits, c(tissue1, tissue2)][1:nlocs, ]
 rm(meth.levels.raw); gc()
 metadata <- pData(bs)

# simulate DMRS
dmrs.true <- NULL
sim.file <- paste0(sim.file.prefix, "/sim.data.n", sampleSize, 
 					  ".", cond, ".all.rda")
if (num.dmrs>0){  
 if(!file.exists(sim.file)){
 	simDMRs(sim.file, bs, num.dmrs)
 }
 load(sim.file)
 bs <- sim.dat.red$bs
 pData(bs)$Infected <- design[,2]
 dmr.mncov <- sim.dat.red$dmr.mncov
 dmr.L <- sim.dat.red$dmr.L
 dmrs.true = sim.dat.red$gr.dmrs
 rm(sim.dat.red)
}

 
num.to.plot <- 500
genomeName <- "hg19"
pval.thresh <- 0.05
testCovariate="Infected"
minNumRegion=5
bpSpan=1000
maxGap=1000
maxGapSmooth=2500
minInSpan=30

if (METHOD=="dmrseq"){
  if (!file.exists(paste0(result.file.prefix, "/regions_", cond, "_", sampleSize, "_", num.dmrs, "DMRs.RData"))){
  regions <- dmrseq(bs, testCovariate=testCovariate, 
                   cutoff = 0.10, minNumRegion=minNumRegion,
                   smooth = TRUE,
                   bpSpan=bpSpan, minInSpan=minInSpan, 
                   maxGapSmooth=maxGapSmooth, maxGap = maxGap,   
                   verbose = TRUE,
                   parallel=TRUE,
                   stat=STAT)
                   
   if(STAT != "stat"){
     statCol <- which(colnames(regions)==STAT)
     regions$stat <- regions[,statCol]
   }
   
   save(regions, file=paste0(result.file.prefix, "/regions_", cond, "_", sampleSize, "_", num.dmrs, "DMRs.RData"))
   }else{
     load(paste0(result.file.prefix, "/regions_", cond, "_", sampleSize, "_", num.dmrs, "DMRs.RData"))
   }
   
   summaryPlotting(regions, bs, testCovariate=testCovariate, sampleSize,
   	 	num.dmrs, cond=cond, 
		pval.thresh, num.to.plot, genomeName, 
		result.file.prefix,
		dmrs.true)
	

  tryDifferentCutoffs(LOCI=regions, minNumRegion,
			  cutoffsQ=c(0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.75, 0.95, 1),
			  result.file.prefix, sampleSize, cond,
			  METHOD="dmrseq", sim.file)	

 # get subsetted power/fdr results
  subsets <- c("low.density", "high.density", "low.coverage", "high.coverage", "low.effsize", "high.effsize")
  for (sub in subsets){
    tryDifferentCutoffs(LOCI=regions, minNumRegion,
			  cutoffsQ=c(0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.75, 0.95, 1),
			  result.file.prefix, sampleSize, cond,
			  METHOD="dmrseq", sim.file, subset=sub)	
  }

  }else{
  if(default){
    source("../BSmooth_DSS_comparison_DEFAULT.R")
	system("cat ../BSmooth_DSS_comparison_DEFAULT.R")
  }else{
    source("../BSmooth_DSS_comparison.R")
    system("cat ../BSmooth_DSS_comparison.R")
  }
  
  if (METHOD == "metilene"){
   		dmrs <- dmrs[dmrs$qval < pval.thresh,]
  }		
  if (nrow(dmrs) == 0){
		message("No significant DMRs found")
  }
}

sessionInfo()




