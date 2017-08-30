# generate bsseq objects from BED files using the preprocessing script:
# "/n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/DNMT3A/preprocess_DNMT3A_chunks.R"

setwd("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/")
source("/n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/summarizeResultsFunctions.R")
data.file.prefix <- Sys.getenv("DATDIR")
result.file.prefix <- Sys.getenv("OUTDIR")
ifelse(!dir.exists(result.file.prefix), dir.create(file.path(result.file.prefix)), FALSE)

library(dmrseq)

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
time2 <- NULL
print(time) 
print(tiss1)  # eg "heart_e16.5"
print(tiss2)  # eg "heart_e11.5"
cond <- paste0(tiss1, tiss2, collapse="")
print(num.dmrs)
print(sampleSize)

METHOD <- Sys.getenv("METHOD")
workers <- as.numeric(Sys.getenv("SLURM_NTASKS"))
print(workers)
# register cores with #workers
library("BiocParallel")
register(MulticoreParam(workers))

job <- as.numeric(Sys.getenv("JOB"))
print(job)
set.seed(3935757)

default <- as.logical(Sys.getenv("default"))
print(default)

load(paste0(data.file.prefix, time, "_all_bsseq.RData"))

print(pData(bs))

if(is.null(sampleNames(bs))){
	sampleNames(bs) <- pData(bs)$SampleNames
}
 
meth.levels.raw = getMeth(bs, type = "raw")
 print(nrow(meth.levels.raw))
nlocs <- nrow(meth.levels.raw)
tissue1 = grep(tiss1, pData(bs)$Cond) 
tissue2 = grep(tiss2, pData(bs)$Cond)
if (length(tissue1) < sampleSize | length(tissue2) < sampleSize){
 	stop("Error: Fewer samples present than specified in 'sampleSize'.  Please specify a valid number of samples to use in each group.")
 }
 tissue1 <- tissue1[1:sampleSize]
 tissue2 <- tissue2[1:sampleSize]
 tissue = as.factor(c(rep(tiss1, length(tissue1)), rep(tiss2, length(tissue2))))
 design = model.matrix(~tissue)
 
 ##Subsetting by the number of CpG sites with a positive coverage 
 ##Generating data for first 100000 such CpG Locations
 no.hits = which(is.na(rowMeans(meth.levels.raw[, c(tissue1, tissue2)])) == TRUE)
 str(no.hits)
 nlocs = min(nlocs, (length(bs) - length(no.hits))) # nlocs orig specified in commandarg
 bs = bs[-no.hits][1:nlocs,c(tissue1, tissue2)]

 show(bs)
 rm(meth.levels.raw)
 gc()


 meth.mat = getCoverage(bs, type = "M")
 unmeth.mat = getCoverage(bs, type = "Cov") - meth.mat
 chr = as.character(seqnames(bs))
 pos = start(bs)
 colnames(meth.mat) <- colnames(unmeth.mat) <- pData(bs)$SampleNames
 chr <- paste0("chr", chr, sep="")
 bs <- BSseq(chr = chr, pos = pos,
               M = meth.mat, Cov = meth.mat+unmeth.mat, sampleNames = sampleNames(bs),
               pData=pData(bs))
 rm(meth.mat)
 rm(unmeth.mat)
 rm(chr)
 rm(pos)
 gc();

# peek at the workspace
print(ls())

pval.thresh <- 0.10
num.to.plot <- 500
genomeName <- "mm10"
testCovariate="Cond"
minNumRegion=5
bpSpan=1000
maxGap=1000
maxGapSmooth=2500
minInSpan=30

if (METHOD=="dmrseq"){

  if (!file.exists(paste0(result.file.prefix, "/regions_", cond, "_", 
  					sampleSize, "_", num.dmrs, "DMRs.RData"))){

   regions <- dmrseq(bs=bs, testCovariate=testCovariate,
		minInSpan=minInSpan, maxGapSmooth=maxGapSmooth, 
		minNumRegion=minNumRegion,
        cutoff = 0.10, maxGap = maxGap, 
        smooth = TRUE, bpSpan=bpSpan,  
        verbose = TRUE,
        parallel=TRUE)
   
   # plot summary of results using the regions results object
   save(regions, file=paste0(result.file.prefix, "/regions_", 
        cond, "_", sampleSize, "_", num.dmrs, "DMRs.RData"))

 }else{
     load(paste0(result.file.prefix, "/regions_", cond, "_", 
                 sampleSize, "_", num.dmrs, "DMRs.RData"))
 }
  

   summaryPlotting(regions, bs, testCovariate=testCovariate, sampleSize,
   	 	num.dmrs, cond=cond, 
		pval.thresh, num.to.plot, genomeName, 
		result.file.prefix,
		dmrs.true)


   dmrs <- regions[regions$qval < pval.thresh & !is.na(regions$qval),]
   # if no dmrs found, take top 2000 as significant
   if (nrow(dmrs) == 0){
	   message(paste0("No significant DMRs found by OLS; correlating ",
			"expression with top 2000 regions"))
	   dmrs <- dmrs[1:min(2000, nrow(dmrs)),]
   }
   default <- FALSE
		
}else{
  if(default){
    source("/n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/BSmooth_DSS_comparison_DEFAULT.R")
	system("cat /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/BSmooth_DSS_comparison_DEFAULT.R")
  }else{
    source("/n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/BSmooth_DSS_comparison.R")
    system("cat /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/BSmooth_DSS_comparison.R")
  }
	
}

source("/n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/DNMT3A/correlateDMRwithExpr_tx.R")
system("cat /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/DNMT3A/correlateDMRwithExpr_tx.R")

 if (METHOD=="dmrseq"){
    load(paste0(result.file.prefix, "/regions_", cond, "_", 
                 sampleSize, "_", num.dmrs, "DMRs.RData"))
    
	dmrs <- regions
	rm(regions)
	dmrs$beta <- dmrseq:::meanDiff(bs, dmrs, testCovariate, adjustCovariate=NULL)

 	corrWithExpression(dmrs, time=time, time2=time2, tiss1=tiss1, tiss2=tiss2, 
			   METHOD=METHOD, result.file.prefix, distance=2000, 
			   pval.thresh=0.10, fdr=TRUE)
    dmrs <- dmrs[dmrs$qval < pval.thresh & !is.na(dmrs$qval),]	
    dmrs$beta <- dmrseq:::meanDiff(bs, dmrs, testCovariate, adjustCovariate=NULL)
 }			    

  if (METHOD == "metilene"){
  		corrWithExpression(dmrs, time=time, time2=time2, tiss1=tiss1, tiss2=tiss2, 
			   METHOD=METHOD, result.file.prefix, distance=2000, 
			   pval.thresh=0.10, fdr=TRUE)
   		dmrs <- dmrs[dmrs$qval < pval.thresh,]
  }	
  
if (nrow(dmrs) > 100){
	corrWithExpression(dmrs, time=time, time2=time2, tiss1=tiss1, tiss2=tiss2, 
						METHOD=METHOD, result.file.prefix, distance=2000,
						pval.thresh=0.10, default=default, effectSize=TRUE)
	corrWithExpression(dmrs, time=time, time2=time2, tiss1=tiss1, tiss2=tiss2, 
						METHOD=METHOD, result.file.prefix, distance=2000,
						pval.thresh=0.10, default=default, effectSize=FALSE)
						
}

if (METHOD != "dmrseq"){
	message(nrow(dmrs), " DMRs found at FDR threshold level ", pval.thresh)
}else{
	message(nrow(dmrs), " DMRs found at FDR threshold level ", pval.thresh)
}




if (METHOD == "dmrseq"){
	# plotting of top 10 exclusive regions for each method
	dmrseq <- dmrs

	load(paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/RESULTS/BSmooth/",
	  "dmrs.BSmooth.n", sampleSize, ".", cond, ".Rdata"))
	bsmooth <- dmrs
	bsmooth$stat <- bsmooth$areaStat #/bsmooth$n
  	bsmooth$beta <- bsmooth$meanDiff
	bsmooth$indexStart <- bsmooth$idxStart
	bsmooth$indexEnd <- bsmooth$idxEnd

	load(paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/RESULTS/DSS/",
	  "dmrs.DSS.n", sampleSize, ".", cond, ".Rdata"))
	rm(dmlTest.sm)
	dss <- dmrs
	dss$stat <- dss$areaStat #/dss$nCG
  	dss$beta <- dss$diff.Methy
	dss$indexStart <- dss$start
	dss$indexEnd <- dss$end

	rm(dmrs)

	# compare dmrSeq with BSmooth
	compareOverlap(dmr1=dmrseq, dmr2=bsmooth, num.to.plot=500,
				   bs=bs, LIST1="dmrseq", LIST2="BSmooth", 
				   result.file.prefix=result.file.prefix, 
				   sampleSize=sampleSize, cond=cond,
				   testCovariate=testCovariate, genomeName=genomeName)
						
	# compare dmrSeq with DSS
	compareOverlap(dmr1=dmrseq, dmr2=dss, num.to.plot=500,
				   bs=bs, LIST1="dmrseq", LIST2="DSS", 
				   result.file.prefix=result.file.prefix, 
				   sampleSize=sampleSize, cond=cond,
				   testCovariate=testCovariate, genomeName=genomeName)
}
