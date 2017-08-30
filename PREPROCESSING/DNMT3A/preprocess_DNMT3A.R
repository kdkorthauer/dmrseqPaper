## 6/2016
## Read in mouse DNMT3a WGBS files (from GEO) in bed format
## Details of format are listed in README file in the WGBS/DNMT3A directory
## Use two functions modified from bsseq package to wrangle them into the BSseq 
##  object format - includes collapsing information over strand
## Once data is in BSseq format, save the resulting .RData file so it can be 
##  loaded into another R session for analysis

## first run the PREPROCESSING/DNMT3A/filterCpGs.sh script to extract only
## the Cytosines in CpG context

 
 library(GenomicRanges)
 library(bsseq)

 read.dnmt <- function(files, sampleNames, rmZeroCov = FALSE, verbose = TRUE){
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
             cat(sprintf("[read.encode] Reading file '%s' ... ", files[ii]))
         }
         if (file.exists(files[ii])){
         ptime1 <- proc.time()
         
         raw <- read.dnmtFileRaw(thisfile = files[ii])
                  
		 ## Now we need to handle that the data has separate lines for each strand
		 ## We join these
		 if (length(unique(raw$strand)) == 2){
			 tmp <- raw[raw$strand == "+",]
			 BS.forward <- BSseq(pos = tmp$pos, chr = tmp$chr, 
								 M = as.matrix(tmp$C_count, ncol = 1),
								 Cov = as.matrix(tmp$CT_count, ncol=1),
									 sampleNames = "forward")
			 tmp <- raw[raw$strand == "-",]
			 BS.reverse <- BSseq(pos = tmp$pos, chr = tmp$chr, 
								 M = as.matrix(tmp$C_count, ncol = 1),
								 Cov = as.matrix(tmp$CT_count, ncol=1),
									 sampleNames = "reverse")
			 BS <- combine(BS.forward, BS.reverse)
			 cn <- rep(sampleNames[ii], 2)
			 names(cn) <- NULL
			 BS <- collapseBSseq(BS, columns = cn)
			 rm(raw)
			 rm(tmp)
		 }else{
		 	BS <- BSseq(pos = raw$pos, chr = raw$chr, 
								 M = as.matrix(raw$C_count, ncol = 1),
								 Cov = as.matrix(raw$CT_count, ncol=1),
									 sampleNames = sampleNames)
		 	rm(raw)
		 }
		 
         ptime2 <- proc.time()
         stime <- (ptime2 - ptime1)[3]
         if (verbose) {
             cat(sprintf("done in %.1f secs\n", stime))  
         }
         print(paste0("class is:", class(BS)))
         return(BS) 
         }else{
         	cat(sprintf("File missing; moving on to next..."))
         } 
     })
     if (verbose) {
         cat(sprintf("[read.encode] Joining samples ... "))
     }
     ptime1 <- proc.time()
     print( sapply(allOut, function(x) class(x) != "BSseq") )
     nullentry <- sapply(allOut, is.null)
     nullentry <- sapply(allOut, function(x) class(x) != "BSseq")
     allOut <- allOut[!nullentry]
     allOut <- combineList(allOut)
     ptime2 <- proc.time()
     stime <- (ptime2 - ptime1)[3]
     if (verbose) {
         cat(sprintf("done in %.1f secs\n", stime))
     }
     allOut
 }
 
 
 
 read.dnmtFileRaw <- function(thisfile, verbose = TRUE){
 # hack the code in 'read.bismarkFileRaw' function of bsseq to 
 # read in ENCODE's pipeline output
  columnHeaders <- c("chr", "pos", "strand", "context", "ratio", "eff_CT_count",
 		"C_count", "CT_count", "rev_G_count", "rev_GA_count", "CI_lower", "CI_upper")
  what0 <- replicate(length(columnHeaders), character(0))
  names(what0) <- columnHeaders
  int <- which(columnHeaders %in% 
  	c("pos","C_count", "CT_count" ))
  num <- which(columnHeaders %in% 
    c("ratio", "eff_CT_count", "CI_lower", "CI_upper"))
  what0[int] <- replicate(length(int), integer(0))
  what0[num] <- replicate(length(num), numeric(0))
 

  ## Read in the file
  if (grepl("\\.gz$", thisfile)) {
      con <- gzfile(thisfile)
  }else{ 
      con <- file(thisfile, open = "r")
  }	  
  # set skip=1 since the files contain a header
  out <- scan(file = con, what = what0, quote = "", sep="\t",
  					na.strings = "NA", quiet = FALSE, skip=1)
  close(con)

  # return data.frame
  DataFrame(out)
}



job <- as.numeric(Sys.getenv("JOB")); print(job)
print(job) # indicates which file will be read in and saved

# get list of files 
setwd("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/DATA/CpG/")
filesToGet <- list.files(".", pattern="*.txt")[job]
Type <- substr(filesToGet, 1, 7)
repl <- as.numeric(substr(filesToGet, 9, 9))
sampNames <-  substr(filesToGet, 1, 9)
table(Type)

print(filesToGet)
print(Type)
print(repl)
print(sampNames)

if (!file.exists(paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/dnmt3a_",
						Type, "_rep", repl, "_bsseq.RData"))){

	dnmt3a <- read.dnmt(filesToGet, sampleNames=sampNames)
	pData <- pData(dnmt3a)
	pData$Sample <- sampNames
	pData$Type <- Type
	pData$Rep <- repl
	pData(dnmt3a) <- pData

	save(dnmt3a, 
		file=paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/dnmt3a_",
							Type, "_rep", repl, "_bsseq.RData"))
	rm(dnmt3a)
}

# if they are all there, combine them
setwd("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/")
Rdats <- list.files(".", pattern="*.RData")
if(length(Rdats)>=6){
  filesToGet <- list.files("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/DATA/CpG/",
  						   pattern="*.txt")
  Type <- substr(filesToGet, 1, 7)
  repl <- as.numeric(substr(filesToGet, 9, 9))
  sampNames <-  substr(filesToGet, 1, 9)
  
  for (j in 1:length(filesToGet)){
  		message(paste0("Joining sample ", j, "..."))
  		if (file.exists(paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/dnmt3a_",
						Type[j], "_rep", repl[j], "_bsseq.RData"))){
		   load(paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/dnmt3a_",
						   Type[j], "_rep", repl[j], "_bsseq.RData"))
		   # combine with previous object (if there is one)
		   if (j > 1){
			   dnmt3a <- combine(dnmt3a_old, dnmt3a)
		   }
		   dnmt3a_old <- dnmt3a
		   rm(dnmt3a)
		}else{
			stop("not all samples processed into their own BSseq object")
		}
  }

dnmt3a <- dnmt3a_old
rm(dnmt3a_old)
sampleNames(dnmt3a) <- sampNames

save(dnmt3a, 
	file=paste0("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/dnmt3a_all_bsseq.RData"))

}
