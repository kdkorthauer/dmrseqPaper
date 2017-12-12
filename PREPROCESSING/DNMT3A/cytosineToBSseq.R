# R script to convert individual sample .cov file from Bismark to bsseq object (RData obj)

# get env variable with the name of the .cov file for this particular sample 
exp <- Sys.getenv("exp")


# load bsseq package with readin function
library(bsseq)

if (file.exists(paste0(exp, ".cytosine"))){
  
    # read in the SRA run table metadata
    # to extract the condition name and replicate number for the given experiment
    meta <- read.table("../SraRunTable.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
	meta$cond <- NA
	meta$cond[meta$genotype_s == "Mx1-cre+;Dnmt3af/f;MSCV_FLT3-ITD-GFP"] <- "KO_FLT3"
	meta$cond[meta$genotype_s == "Mx1-cre+;Dnmt3a+/+;MSCV_FLT3-ITD-GFP"] <- "WT_FLT3"
	meta$cond[meta$genotype_s == "Mx1-cre+;Dnmt3a+/+;MSCV_GFP"] <- "WT_WTFL"
	
	meta$rep <- NA
	for (r in 1:nrow(meta)){
	    biosamps <- unique(meta$BioSample_s[meta$cond ==  meta$cond[r]])
		meta$rep[r] <- which(meta$BioSample_s[r] == biosamps)
	}
	
	#exp <- meta$Experiment_s[which(meta$Run_s=="SRR3457815")[1]]
	#srr <- meta$Run_s[meta$Experiment == exp]
    #fils <- paste0(srr, "_1_val_1.fq_bismark_bt2_pe.bismark.cov")
    #bs <- read.bismark(fils, srr, verbose=TRUE)
  
	Cond <- meta$cond[meta$Experiment_s == exp][1]
	Rep <- meta$rep[meta$Experiment_s == exp][1]
	Samp <- paste(Cond, Rep, sep="_rep")
	Tissue <- meta$source_name_s[meta$Experiment_s == exp][1]
	
	bs <- read.bismark(paste0(exp, ".cytosine"),
				 Samp, fileType="cytosineReport",
				 strandCollapse=TRUE,
				 rmZeroCov=FALSE,
				 verbose = TRUE)

	# remove chromosomes that aren't in 1:23,X,Y
	bs <- chrSelectBSseq(bs, seqnames = c(seq(1,23),"X","Y"), order = TRUE)
    
    message(paste0("writing to file ", "../dnmt3a_", Cond, "_rep", Rep, "_bsseq.RData"))
	
	save(bs, file=paste0("../dnmt3a_", Cond, "_rep", Rep, "_bsseq.RData"))
	
	allSamps <- unique(paste(meta$cond, "rep", meta$rep, meta$source_name_s, 
						sep="_"))
	processed <- rep(FALSE, length(allSamps))
	for (s in 1:length(allSamps)){
	    allSamps2 <- paste0(c(substr(allSamps[s],1,11),
	    					  substr(allSamps[s],13,13)), collapse="")
	    processed[s] <- file.exists(paste0("../dnmt3a_",
		 allSamps2,
		 "_bsseq.RData"))
	}

    # if all there, join together
	if (sum(processed)==length(processed)){
		message("Merging samples into one BSseq object")
		allSamps2 <- paste0(c(substr(allSamps[1],1,11),
					  substr(allSamps[1],13,13)), collapse="")
		load(paste0("../dnmt3a_", allSamps2, "_bsseq.RData"))
		bs_all <- bs
		sampleNames(bs) <- allSamps2
		
		message(paste0("Sample 1 done"))
		for (i in 2:length(processed)){
			allSamps2 <- paste0(c(substr(allSamps[i],1,11),
						  substr(allSamps[i],13,13)), collapse="")
			load(paste0("../dnmt3a_", allSamps2, "_bsseq.RData"))
			sampleNames(bs) <- allSamps2
			bs_all <- combine(bs_all, bs)
			message(paste0("Sample ", i, " done"))
		}
		
		bs <- bs_all
		rm(bs_all)
		# add nicely formatted metadata
		pData(bs) <- DataFrame(SampleNames=allSamps, 
							   Cond=substr(allSamps, 1, 7),
							   Rep=substr(allSamps, 13, 13),
							   Tissue=substring(allSamps, first=15))
							   
		save(bs, file="../dnmt3a_all_bsseq.RData")
	}

}else {
	message(paste0("Coverage file for experiment ", exp, " doesn't exist"))
}

