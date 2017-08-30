# generate bsseq objects from BED files using the preprocessing script:
# /n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/cytosineToBSseq_public.R

setwd("/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA")
source("/n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/summarizeResultsFunctions.R")
data.file.prefix <- Sys.getenv("DATDIR")
result.file.prefix <- Sys.getenv("OUTDIR")
ifelse(!dir.exists(result.file.prefix), 
	dir.create(file.path(result.file.prefix)), FALSE)
	
library(dmrseq)

args=commandArgs(TRUE)
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

# replace "_" in tiss1, tiss2 with spaces
# (had to put "_" in tiss names since ARGS doesn't interpret spaces properly)
# these didn't work properly - don't change tiss1 or tiss2
tiss1 <- gsub("_", " ", tiss1); print(tiss1)
tiss2 <- gsub("_", " ", tiss2); print(tiss2)

cond <- paste0(tiss1, "_", tiss2)
print(num.dmrs)
print(sampleSize)

chromosome <- "all"
print(chromosome)
METHOD <- Sys.getenv("METHOD")

default <- as.logical(Sys.getenv("default"))
print(default)

workers <- as.numeric(Sys.getenv("SLURM_NTASKS"))
print(workers)
# register cores with #workers
library("BiocParallel")
register(MulticoreParam(workers))

job <- as.numeric(Sys.getenv("JOB"))
print(job)
set.seed(3935757)


if (!chromosome == "all"){
  chromosome <- as.numeric(chromosome)
  chr.set <- paste0("chr", c(seq(1,22), "X", "Y"))
  chromosome <- chr.set[chromosome]
}
	
load(paste0(data.file.prefix, time, "_", tiss1, "_all_bsseq_id.RData")); show(bs)
bs_tiss1 <- bs
rm(bs); gc()
meta <- pData(bs_tiss1)
rownames(meta) <- meta$SampleNames
pData(bs_tiss1) <- meta
load(paste0(data.file.prefix, time, "_", tiss2, "_all_bsseq_id.RData")); show(bs)
meta <- pData(bs)
rownames(meta) <- meta$SampleNames
pData(bs) <- meta
bs <- bsseq::combine(bs, bs_tiss1)
rm(bs_tiss1); gc()

print(pData(bs)$SampleNames)
print(pData(bs)$Cond)
print(pData(bs)$Rep)
print(pData(bs)$BioSamp)
print(pData(bs)$Age)

sampleNames(bs) <- pData(bs)$SampleNames

meth.levels.raw = getMeth(bs, type = "raw")
 print(nrow(meth.levels.raw))
nlocs <- nrow(meth.levels.raw)
tissue1 = grep(tiss1, pData(bs)$Cond) 
tissue2 = grep(tiss2, pData(bs)$Cond)
if (length(tissue1) < sampleSize | length(tissue2) < sampleSize){
 	stop("Error: Fewer samples present than specified in 'sampleSize'.  Please specify a valid number of samples to use in each group.")
 }
 
 # skip rep 2 of Small Intestine - is from a different sample than the other tissues
 if (tiss1 == "Small Intestine"){
 	tissue1 <- tissue1[pData(bs)$Age[tissue1] %in% c(3,34)]
 }else{
 	tissue1 <- tissue1[1:sampleSize]
 }
 
 if (tiss2 == "Small Intestine"){
 	tissue2 <- tissue2[pData(bs)$Age[tissue2] %in% c(3,34)]
 }else{
 	tissue2 <- tissue2[1:sampleSize]
 }
 tissue = as.factor(c(rep(tiss1, length(tissue1)), rep(tiss2, length(tissue2))))
 design = model.matrix(~tissue)
 
 ##Subsetting by the number of CpG sites with a positive coverage 
 ##Generating data for first 100000 such CpG Locations
 no.hits = which(is.na(rowMeans(meth.levels.raw[, c(tissue1, tissue2)])) == TRUE)
 str(no.hits)
 nlocs = min(nlocs, (length(bs) - length(no.hits))) # nlocs orig specified in commandarg
 bs = bs[-no.hits][1:nlocs,c(tissue1, tissue2)]
 meth.mat = getCoverage(bs, type = "M")
 unmeth.mat = getCoverage(bs, type = "Cov") - meth.mat
 
 show(bs)
 rm(meth.levels.raw)
 gc()
 chr = as.character(seqnames(bs))
 pos = start(bs)
 print(unique(chr))
 colnames(meth.mat) <- colnames(unmeth.mat) <- pData(bs)$SampleNames

# add "chr" prefix to chromosome names
chr <- paste0("chr", chr, sep="")
bs <- BSseq(chr = chr, pos = pos,
               M = meth.mat, Cov = meth.mat+unmeth.mat, sampleNames = sampleNames(bs),
               pData=pData(bs))
rm(meth.mat)
rm(unmeth.mat)
rm(pos)
rm(chr)
gc()

# peek at the workspace
print(ls())        

pval.thresh <- 0.05
num.to.plot <- 500
genomeName <- "hg38"
testCovariate="Cond"
minNumRegion=5
bpSpan=1000
maxGap=1000
maxGapSmooth=2500
minInSpan=30

if (METHOD=="dmrseq"){

  if (!file.exists(paste0(result.file.prefix, "/regions_", cond, "_", 
  					sampleSize, "_", num.dmrs, "DMRs.RData"))){
   regions <- dmrseq(bs, testCovariate= testCovariate, 
				   matchCovariate="Age",
                   cutoff = 0.10, minNumRegion=minNumRegion,
                   smooth = TRUE,
                   bpSpan=bpSpan, minInSpan=minInSpan, 
                   maxGapSmooth=maxGapSmooth, maxGap = maxGap,   
                   verbose = TRUE,
                   parallel=TRUE, workers=workers)
   
   save(regions, file=paste0(result.file.prefix, "/regions_", 
        cond, "_", sampleSize, "_", num.dmrs, "DMRs.RData"))

   }else{
     load(paste0(result.file.prefix, "/regions_", cond, "_", 
                 sampleSize, "_", num.dmrs, "DMRs.RData"))
   }

   summaryPlotting(regions, bs, testCovariate=testCovariate, sampleSize,
   	 	num.dmrs, cond=cond, 
		pval.thresh, num.to.plot, genomeName, 
		result.file.prefix)
   
   dmrs <- regions[regions$qval < pval.thresh & !is.na(regions$qval),]	
   
   message(paste0(nrow(dmrs), " DMRs found at the ", pval.thresh, " level"))
   # if no dmrs found, take top 1000 as significant
   if (nrow(dmrs) < 100){
	  message(paste0("No significant DMRs found; correlating ",
		   "expression with top 1000 regions"))
	  dmrs <- regions[1:min(15000, nrow(regions)),]
   }
   
default <- FALSE

}else{

  if(default){
    system("cat /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/BSmooth_DSS_comparison_DEFAULT.R")
	source("/n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/BSmooth_DSS_comparison_DEFAULT.R")
  }else{
    system("cat /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/BSmooth_DSS_comparison.R")
    source("/n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/BSmooth_DSS_comparison.R")
  }

}

if (METHOD == "dmrseq"){
	# plotting of top 10 exclusive regions for each method
	dmrseq <- dmrs

	load(paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/RESULTS/BSmooth_default/",
	  "dmrs.BSmooth.n", sampleSize, ".", cond, ".DEFAULT.Rdata"))
	rm(bs.tstat)
	bsmooth <- dmrs
	bsmooth$stat <- bsmooth$areaStat #/bsmooth$n
  	bsmooth$beta <- bsmooth$meanDiff
	bsmooth$indexStart <- bsmooth$idxStart
	bsmooth$indexEnd <- bsmooth$idxEnd
	bsmooth$L <- bsmooth$n

	load(paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/RESULTS/DSS_default/",
	  "dmrs.DSS.n", sampleSize, ".", cond, ".DEFAULT.Rdata"))
	rm(dmlTest.sm)
	dss <- dmrs
	dss$stat <- dss$areaStat #/dss$nCG
    dss$beta <- dss$diff.Methy 
	dss$indexStart <- dss$start
	dss$indexEnd <- dss$end
	dss$L <- dss$nCG
	rm(dmrs)
		
	met <- read.table(paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/RESULTS/metilene_default/",
	  "metilene_output_", gsub(" ", "_", cond), "_n", sampleSize, "_0DMRs.txt"), stringsAsFactors=FALSE)
	colnames(met) <- c("chr", "start", "end", "qval", "beta", "L", "pval1", 
	 				"pval2", "mean1", "mean2")
	met <- met[met$qval < pval.thresh,]
	met$stat <- met$beta * met$L
	met$indexStart <- match(paste0(met$chr, "_", met$start+1), 
	 					paste0(as.character(seqnames(bs)), "_", start(bs)))
    met$indexEnd <- match(paste0(met$chr, "_", met$end), 
	 					paste0(as.character(seqnames(bs)), "_", start(bs)))
	
    comp.plots <- paste0(result.file.prefix, "/compPlots/")
    ifelse(!dir.exists(comp.plots), dir.create(file.path(comp.plots)), FALSE)

    # compare dmrSeq with metilene
	compareOverlap(dmr1=dmrseq, dmr2=met, num.to.plot=25,
				   bs=bs, LIST1="dmrseq", LIST2="metilene", 
				   result.file.prefix=comp.plots, 
				   sampleSize=sampleSize, cond=cond,
				   testCovariate=testCovariate, genomeName=genomeName)

	# compare dmrSeq with BSmooth
	compareOverlap(dmr1=dmrseq, dmr2=bsmooth, num.to.plot=25,
				   bs=bs, LIST1="dmrseq", LIST2="BSmooth", 
				   result.file.prefix=comp.plots, 
				   sampleSize=sampleSize, cond=cond,
				   testCovariate=testCovariate, genomeName=genomeName)
						
	# compare dmrSeq with DSS
	compareOverlap(dmr1=dmrseq, dmr2=dss, num.to.plot=25,
				   bs=bs, LIST1="dmrseq", LIST2="DSS", 
				   result.file.prefix=comp.plots, 
				   sampleSize=sampleSize, cond=cond,
				   testCovariate=testCovariate, genomeName=genomeName)

    if (METHOD=="dmrseq"){
      dmrs <- dmrseq
    }else if(METHOD=="BSmooth"){
      dmrs <- bsmooth
    }else if(METHOD=="DSS"){
      dmrs <- dss
    }else if(METHOD=="metilene"){
      dmrs <- met
    }
    rm(dmrseq)
    rm(dss)
    rm(bsmooth)
    rm(met)				   
}

# correlate methylation difference with expression from RNA seq

if (nrow(dmrs) > 0){
 source("/n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/correlateDMRwithExpr_tx.R")
 system("cat /n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/correlateDMRwithExpr_tx.R")
 
 if (METHOD=="dmrseq"){
    load(paste0(result.file.prefix, "/regions_", cond, "_", 
                 sampleSize, "_", num.dmrs, "DMRs.RData"))
    dmrs <- regions
    rm(regions)
    dmrs$beta <- dmrseq:::meanDiff(bs, dmrs, testCovariate, adjustCovariate=NULL)

 	corrWithExpression(dmrs, time=time, time2=time2, tiss1=tiss1, tiss2=tiss2, 
			   METHOD=METHOD, result.file.prefix, distance=2000, 
			   pval.thresh=0.05, fdr=TRUE)
    dmrs <- dmrs[dmrs$qval < pval.thresh & !is.na(dmrs$qval),]	
    dmrs$beta <- dmrseq:::meanDiff(bs, dmrs, testCovariate, adjustCovariate=NULL)
 }			    
 
 if (METHOD == "metilene"){
	   corrWithExpression(dmrs, time=time, time2=time2, tiss1=tiss1, tiss2=tiss2, 
		  METHOD=METHOD, result.file.prefix, distance=2000, 
		  pval.thresh=0.05, fdr=TRUE)
	   dmrs <- dmrs[dmrs$qval < pval.thresh,]
 }	
 
 corrWithExpression(dmrs, time=time, time2=time2, tiss1=tiss1, tiss2=tiss2, 
			   METHOD=METHOD, result.file.prefix, distance=2000, 
			   pval.thresh=0.05, default=default, effectSize=TRUE)
  
 corrWithExpression(dmrs, time=time, time2=time2, tiss1=tiss1, tiss2=tiss2, 
			   METHOD=METHOD, result.file.prefix, distance=2000, 
			   pval.thresh=0.05, default=default, effectSize=FALSE)

}else{
	message("no DMRs found")
}

