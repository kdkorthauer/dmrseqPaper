# R script to convert individual sample .cytosine file from Bismark to bsseq object (RData obj)

# get env variable with the name of the .cytosine file for this particular sample
# the experiment id
exp <- Sys.getenv("exp")

# list all merged cytosine coverage files
system("ls -l *cytosine.merged_CpG_evidence.cov")

# load bsseq package with readin function
library(bsseq)

if (file.exists(paste0(exp, ".cytosine.merged_CpG_evidence.cov"))){
  
    # read in the SRA run table metadata
    # to extract the condition name and replicate number for the given experiment
    meta <- read.table("/n/irizarryfs01_backed_up/kkorthauer/WGBS/CODE/ROADMAP/SraRunTable_public_heart_intestine.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
	meta$cond <- meta$TISSUE_TYPE_s
		
	meta$rep <- NA
	for (r in 1:nrow(meta)){
	    biosamps <- unique(meta$BioSample_s[meta$cond ==  meta$cond[r]])
		meta$rep[r] <- which(meta$BioSample_s[r] == biosamps)
	}
  
	Cond <- meta$cond[meta$BioSample_s == exp][1]
	Rep <- meta$rep[meta$BioSample_s == exp][1]
	Samp <- paste(Cond, Rep, sep="_rep")
	BioSampleID <- paste(exp)
	Age <- meta$DONOR_AGE_s[meta$BioSample_s == exp][1]
	
	if (!file.exists(paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/",
    			 "roadmap_", Cond, "_rep", Rep, "_bsseq.RData"))){
			bs <- read.bismark(paste0(exp, ".cytosine.merged_CpG_evidence.cov"),
						 Samp, fileType="cov",
						 strandCollapse=TRUE,
						 rmZeroCov=FALSE,
						 verbose = TRUE)

			# remove chromosomes that aren't in 1:23,X,Y
			bs <- chrSelectBSseq(bs, seqnames = c(seq(1,23),"X","Y"), order = TRUE)

			message(paste0("writing individual replicate to file /n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/",
						 "roadmap_", Cond, "_rep", Rep, "_bsseq.RData"))
			
			save(bs, file=paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/",
						 "roadmap_", Cond, "_rep", Rep, "_bsseq.RData"))
	}else{
		load(paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/",
    			 "roadmap_", Cond, "_rep", Rep, "_bsseq.RData"))
    	pData(bs) <- data.frame(SampleNames=Samp, Cond=Cond, Rep=Rep, 
    							BioSamp=BioSampleID, Age=Age)
    	sampleNames(bs) <- Samp
    	save(bs, 			 
    		file=paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/",
						 "roadmap_", Cond, "_rep", Rep, "_bsseq.RData"))
	}
		
	if (!file.exists(paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/",
						  "roadmap_", Cond, "_all_bsseq_id.RData"))){
	  # vector of all the reps that exist for that cond
	  allreps <- unique(meta$rep[meta$cond==Cond])
	  allSamps <- paste(Cond, allreps, sep="_rep")
	  processed <- rep(FALSE, length(allSamps))
	  for (s in 1:length(allSamps)){
		  # check if that condition/rep combination has been processed and file exists
		  processed[s] <- file.exists(paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/",
				   "roadmap_", allSamps[s], "_bsseq.RData"))
	  }


	  # if all there, join together (only if there is more than one rep)
	  if (sum(processed)==length(processed) & length(processed) > 1){
		  message(paste0("Merging all reps from ", Cond, " into one BSseq object"))
		  message(paste0("writing all replicates for condition ", Cond, 
				   " to file /n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/",
				   "roadmap_", Cond, "_all_bsseq.RData"))
		  load(paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/",
				   "roadmap_", allSamps[1], "_bsseq.RData"))
		  bs_all <- bs
		  #sampleNames(bs) <- allSamps[1]
		
		  message(paste0("Rep 1 done"))
		  for (i in 2:length(processed)){
			  load(paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/",
				   "roadmap_", allSamps[i], "_bsseq.RData"))
			  #sampleNames(bs) <- allSamps[i]
			  bs_all <- combine(bs_all, bs)
			  message(paste0("Rep ", i, " done"))
		  }
		
		  bs <- bs_all
		  rm(bs_all)
							   
		  save(bs, file=paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/",
						  "roadmap_", Cond, "_all_bsseq_id.RData"))
		}
	}
}else {
	message(paste0("Coverage file for experiment ", exp, " doesn't exist"))
}

