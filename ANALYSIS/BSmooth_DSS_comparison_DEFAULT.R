## code to run comparison methods BSmooth and DSS
## Variable "METHOD" determines which of these will be run
## Results will be placed in a subfolder that is named after the method
##
library(data.table)

  # convert covariates to column numbers if characters
  if (is.character(testCovariate)){
    tc <- testCovariate
    testCovariate <- which(colnames(pData(bs)) == testCovariate)
    if (length(testCovariate)==0){
      stop(paste0("testCovariate named ", tc, 
                  " not found in pData(). ", 
                  "Please specify a valid testCovariate"))
    }
    rm(tc)
  }
  
message("Running BSmooth/DSS with Default settings")
 adjustCovariate=NULL
 coeff <- 2:(2+length(testCovariate)-1)

 testCov <- pData(bs)[,testCovariate]

 if (length(unique(testCov))==1){
 	message(paste0("Warning: only one unique value of the specified ",
 				"covariate of interest.  Assuming null comparison and ",
 				"splitting sample group into two equal groups"))
 	testCov <- rep(1, length(testCov))
    testCov[1:round(length(testCov)/2)] <- 0
 }
 if (!is.null(adjustCovariate)){
 	adjustCov <- pData(bs)[,adjustCovariate]
 	design <- model.matrix( ~ testCov + adjustCov)
 	colnames(design)[coeff] <- colnames(pData(bs))[testCovariate]
 	colnames(design)[,(max(coeff)+1):ncol(design)] <- colnames(pData(bs))[adjustCovariate]
 }else{
 	design <- model.matrix(~ testCov)
 	colnames(design)[coeff] <- colnames(pData(bs))[testCovariate]
 }

  # load simulated data if necessary
if (num.dmrs > 0){
   load(sim.file)			  
   bs <- sim.dat.red$bs
   pData(bs) <- metadata 
   rm(sim.dat.red); gc()
}
# rename column names to match sample names
sampleNames(bs) <- pData(bs)$Sample
M.mat <- assay(bs, "M")
Cov.mat <- assay(bs, "Cov")
colnames(M.mat) <- colnames(Cov.mat) <- sampleNames(bs)
assays(bs) <- SimpleList("M"=M.mat, "Cov"=Cov.mat)

# Run Metilene 
if (METHOD=="Metilene" | METHOD=="metilene"){
    # first generate input file
	Prop.mat <- M.mat / Cov.mat
	colnames(Prop.mat) <- paste0(c("g2_", "g1_")[design[,coeff]+1],
								colnames(Prop.mat))
	met.file <- paste0(data.file.prefix, "/metilene_input_", cond,
		"_n", sampleSize, "_", num.dmrs, "DMRs.txt")
	met.file <- gsub(" ", "_", met.file)
	
	Prop.mat <- data.frame(chr=as.character(seqnames(bs)), pos=start(bs), Prop.mat)
	if (!file.exists(met.file)){
		fwrite(Prop.mat, file = met.file, append = FALSE, quote = FALSE, sep = "\t")
  	}
  	met.file.out <- paste0("/metilene_output_", cond,
		"_n", sampleSize, "_", num.dmrs, "DMRs.txt")
	met.file.out <- gsub(" ", "_", met.file.out)
	
	if (!file.exists(paste0(result.file.prefix, met.file.out)) | 
		file.size(paste0(result.file.prefix, met.file.out)) == 0){
  		system(paste0("/n/irizarryfs01_backed_up/kkorthauer/softwareTools/",
  			"metilene_v0.2-7/metilene -a g1 -b g2 -t ", workers, " ",
  			 met.file, " > ", result.file.prefix, met.file.out))
  	}
  	
  	dmrs <- read.table(paste0(result.file.prefix, met.file.out),
  				stringsAsFactors=FALSE, sep="\t", header=FALSE)
  	colnames(dmrs) <- c("chr", "start", "end", "qval", "beta", "L", "pMU", 
  					    "pKS", "mu1", "mu2")
  	dmrs <- dmrs[order(dmrs$qval),]
  	message(paste0(sum(dmrs$qval < pval.thresh), " dmrs found at the ", 
  		pval.thresh, " level for ", METHOD))
  	dmrs$stat <- 1-dmrs$qval
  	pos <- as.numeric(start(bs))
  	dmrs$indexStart <- match(dmrs$start+1, pos)
    dmrs$indexEnd <- match(dmrs$end, pos)
    rm(pos)
    
     if (num.dmrs == 0){
	# summary plot and plots of individual regions
	summarizeResults(bs, OBS=dmrs, METHOD, which.sig=which(dmrs$qval < pval.thresh), 
				testCovariate=testCovariate, num.to.plot=500, 
				result.file.prefix=result.file.prefix, 
				sampleSize=sampleSize, cond=cond, 
				genomeName=genomeName)
	}else{
	  evaluateSimulation(OBS=dmrs, result.file.prefix, num.dmrs, sampleSize, cond,
				which.sig=which(dmrs$qval < pval.thresh), METHOD, num.to.plot=500, 
				sim.file, design)
			 
	  tryDifferentCutoffs(LOCI=dmrs, minNum,
			cutoffsQ=c(0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.75, 0.95, 1),
			result.file.prefix, num.dmrs, sampleSize, cond,
			METHOD, sim.file, maxGap=maxGap)	
   
	}
 
	if (nrow(dmrs) > 0){
	 message(paste0("Median n loci: ", median(dmrs$L)))
	 paste0("Mean n loci: ", mean(dmrs$L))
	 paste0("Range n loci: ", range(dmrs$L))
    }
}

rm(M.mat); rm(Cov.mat); gc()

if (METHOD=="BSmooth" | METHOD=="bsmooth"){


  library(bsseq)
  
  ##### Code to run BSmooth

   if (num.dmrs == 0){
	   smooth.file.name=paste0(result.file.prefix, "/smooth.BSmooth.n", 
		   sampleSize, ".", cond, ".DEFAULT.RData")
	   dmr.file.name=paste0(result.file.prefix, "/dmrs.BSmooth.n", 
	   		sampleSize, ".", cond, ".DEFAULT.Rdata")

   }else{
	   smooth.file.name=paste0(result.file.prefix, "/smooth.BSmooth.n", 
		   sampleSize, ".", cond, ".", num.dmrs, "DMRS.DEFAULT.RData")
	   dmr.file.name=paste0(result.file.prefix, "/dmrs.BSmooth.n", 
	   	   sampleSize, ".", cond, ".", num.dmrs, "DMRS.DEFAULT.Rdata")
   }


  # smooth & save results object
  if(!file.exists(smooth.file.name) ){
  	# run BSmooth with default parameters
  	# with specialized smoothing parameters as in Perm, there are a small
  	# number of smoothing errors out of locfit; BSmooth doesn't handle these well
  	# so leave params as default
  	# Bsmooth uses more memory than other methods; use fewer cores to stay under constraints
	 bs.SMOOTH <- BSmooth(bs, mc.cores = ceiling(workers/2), verbose = TRUE,
	 						#ns = minInSpan, 
	 						h = bpSpan, 
	 						#maxGap = maxGapSmooth,
	 						parallelBy="chromosome")
	 save(bs.SMOOTH, file=smooth.file.name)
  }else{
	 load(smooth.file.name)  
  }

  # print out how many loci are missing smoothing estimates
  
  
  # get group names
  # make sure they are in the same order as the design mat - 
  # such that the direction will be tissue with x=1 minus tissue with x=0
  g1 <- sampleNames(bs)[which(design[,coeff] == 1)]
  g2 <- sampleNames(bs)[which(design[,coeff] == 0)]
    

  if (!file.exists(dmr.file.name) ){
	 # compute t-statistics
	 # local.correct = FALSE so that analysis is focused on the 
	 # local scale only and no low-frequency signal is subtracted
	 
	 bs.tstat <- BSmooth.tstat(bs.SMOOTH,
								 group1 = g1,
								 group2 = g2,
								 local.correct = FALSE,
								 verbose = TRUE)

	 # DMR calling
	 dmrs0 <- dmrFinder(bs.tstat, maxGap=maxGap,
	 					stat="tstat")
	 dmrs <- subset(dmrs0, n >= minNum & abs(meanDiff) >= 0)
	 nrow(dmrs)
	 head(dmrs)	
	 dmrs <- dmrs[order(-abs(dmrs$areaStat)/dmrs$n),]

	 save(dmrs, file=dmr.file.name)
  }else{
  	 load(dmr.file.name)
  }
  
  # plotting
  # set up colors to represent the covariate of interest
  colors <- c("red", "blue")
  colors <- cbind( c(0,1), colors)

  # map each row of the design matrix to a color
  z <- colors[,2][match(design[,coeff], colors[,1])]
  pData(bs)$col <- as.character(z)
  
  # search for tissue/condition names in relevant design matrix column
  # doesn't assume pData already contains the condition labels
  t1 <- grep(tiss1, colnames(design)[coeff])
  t2 <- grep(tiss2, colnames(design)[coeff])

  # 0 = red, 1 = blue
  if (length(t1) > 0){
	  pData(bs)$label <- tiss2
	  pData(bs)$label[pData(bs)$col == "blue"] <- tiss1
  }else if(length(t2) > 0){
	  pData(bs)$label <- tiss1
	  pData(bs)$label[pData(bs)$col == "blue"] <- tiss2
  }

  
  dmrs$stat <- dmrs$areaStat/dmrs$n
  dmrs$beta <- dmrs$meanDiff
  dmrs$chr <- as.character(dmrs$chr)
  dmrs$indexStart <- dmrs$idxStart
  dmrs$indexEnd <- dmrs$idxEnd
  if (length(grep("chr", dmrs$chr)) == 0 ){
   		dmrs$chr <- paste0("chr", dmrs$chr)
  }	
  if (num.dmrs == 0){
	# summary plot and plots of individual regions
	summarizeResults(bs, OBS=dmrs, METHOD, which.sig=1:nrow(dmrs), 
				testCovariate=testCovariate, num.to.plot=500, 
				result.file.prefix=result.file.prefix, 
				sampleSize=sampleSize, cond=cond, 
				genomeName=genomeName)
  }else{
    evaluateSimulation(OBS=dmrs, result.file.prefix, num.dmrs, sampleSize, cond,
			  which.sig=1:nrow(dmrs), METHOD, num.to.plot=500, 
			  sim.file, design)
	# try different cutoffs and save #Detected and FDR
		tryDifferentCutoffs(LOCI=bs.tstat, minNum,
			  cutoffsQ=c(0.00001,0.0001,0.00025,0.0005,0.001,0.002,0.005,0.01,0.025,0.05),
			  result.file.prefix, num.dmrs, sampleSize, cond,
			  METHOD, sim.file, maxGap=maxGap)
 }
 
  if (nrow(dmrs) > 0){
   message(paste0("Median n loci: ", median(dmrs$n)))
   paste0("Mean n loci: ", mean(dmrs$n))
   paste0("Range n loci: ", range(dmrs$n))
 }
}


if (METHOD=="DSS" | METHOD=="dss"){
   ##### Code to run DSS
   library(DSS)

	#  results file name will include the number of dmrs
	# if this is not zero.
	
   if (num.dmrs == 0){
	   smooth.file.name=paste0(result.file.prefix, "/smooth.DSS.n", 
		   sampleSize, ".", cond, ".DEFAULT.RData")
	   dmr.file.name=paste0(result.file.prefix, "/dmrs.DSS.n", 
	   		sampleSize, ".", cond, ".DEFAULT.Rdata")

   }else{
	   smooth.file.name=paste0(result.file.prefix, "/smooth.DSS.n", 
		   sampleSize, ".", cond, ".", num.dmrs, "DMRS.DEFAULT.RData")
	   dmr.file.name=paste0(result.file.prefix, "/dmrs.DSS.n", 
	   	   sampleSize, ".", cond, ".", num.dmrs, "DMRS.DEFAULT.Rdata")
   }

   if(!file.exists(smooth.file.name) ){

   		# smooth methylation levels and test each site for DML
   		dmlTest.sm <- DMLtest(bs, 
   							group1=sampleNames(bs)[which(design[,coeff] == 1)],
							group2=sampleNames(bs)[which(design[,coeff] == 0)],
							smoothing=TRUE)

		save(dmlTest.sm, file=smooth.file.name)
	}else{
		load(smooth.file.name)
	}
	
   				
   # determine which DMLS are significant according some pvalue threshold		
   # how to pick the pvalue threshold? is this already adjusted for multiple comparisons?
 
  
   if (!file.exists(dmr.file.name) ){
	  # merge nearby DMLs to construct DMRs (add hoc) - use default threshold settings
	  # default p.threshold is too stringent when not merging nearby regions
	  # the way the callDMR function interprets the pct.sig number does not allow for
	  # specifying pct.sig=1 (since it will require that the proportion significant in
	  # the region is strictly greater than pct.sig.  Code at line 140 should be changed 
	  # from if(mean(pp) > pct.sig) to if(mean(pp) >= pct.sig) so that the function 
	  # behaves as expected.  Instead, need to put in a number just less than 1 but greater
	  # than the number of loci in the region minus 1 divided by the total number of loci
	  
	  dmrs <- callDMR(dmlTest.sm, minCG=minNumRegion-1)

	  # order by area stat
	  dmrs <- dmrs[order(-abs(dmrs$areaStat)/dmrs$nCG),]

	  # save dmr result object
	  # file name will contain number of dmrs if not zero
 	  save(dmrs, file=dmr.file.name)
   }else{
   	  load(dmr.file.name)
   }
   
   # plotting
   # set up colors to represent the covariate of interest
   colors <- c("red", "blue")
   colors <- cbind( c(0,1), colors)

   # map each row of the design matrix to a color
   z <- colors[,2][match(design[,coeff], colors[,1])]
   pData(bs)$col <- as.character(z)

   # search for tissue/condition names in relevant design matrix column
   # doesn't assume pData already contains the condition labels
   t1 <- grep(tiss1, colnames(design)[coeff])
   t2 <- grep(tiss2, colnames(design)[coeff])

   # 0 = red, 1 = blue
   if (length(t1) > 0){
	   pData(bs)$label <- tiss2
	   pData(bs)$label[pData(bs)$col == "blue"] <- tiss1
   }else if(length(t2) > 0){
	   pData(bs)$label <- tiss1
	   pData(bs)$label[pData(bs)$col == "blue"] <- tiss2
   }
 	
   dmrs$stat <- dmrs$areaStat/dmrs$nCG
   dmrs$beta <- dmrs$diff.Methy 
   dmrs$chr <- as.character(dmrs$chr)
   dmrs$indexStart <- dmrs$start
   dmrs$indexEnd <- dmrs$end
   if (length(grep("chr", dmrs$chr)) == 0 ){
   		dmrs$chr <- paste0("chr", dmrs$chr)
   }

   
   if(num.dmrs==0){
		# summary plot and plots of individual regions
		summarizeResults(bs, OBS=dmrs, METHOD, which.sig=1:nrow(dmrs), 
				testCovariate=testCovariate, num.to.plot=500, 
				result.file.prefix=result.file.prefix, 
				sampleSize=sampleSize, cond=cond, 
				genomeName=genomeName)
    }else{
    	evaluateSimulation(OBS=dmrs, result.file.prefix, num.dmrs, sampleSize, cond,
				which.sig=1:nrow(dmrs), METHOD, num.to.plot=500,
				sim.file, design)
		# try different cutoffs and save #Detected and FDR
			tryDifferentCutoffs(LOCI=dmlTest.sm, minNum,
			  cutoffsQ=c(1e-6, 1e-5, 1e-4, 0.00025, 0.001, 0.002, 0.01, 0.025, 0.05, 0.1),
			  result.file.prefix, num.dmrs, sampleSize, cond,
			  METHOD, sim.file, maxGap=maxGap)
    }
   
   if (nrow(dmrs) > 0){
	 message(paste0("Median n loci: ", median(dmrs$nCG)))
	 paste0("Mean n loci: ", mean(dmrs$nCG))
	 paste0("Range n loci: ", range(dmrs$nCG))
   }
}
