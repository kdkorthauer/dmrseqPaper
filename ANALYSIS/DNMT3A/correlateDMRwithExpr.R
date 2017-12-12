library(ggplot2)
library(annotatr)
library(rtracklayer)
library(biomaRt)
library(DESeq2)
library(fmsb) 
library(readr)   
   
corrWithExpression <- function(dmrs, time, time2=NULL, tiss1, tiss2, METHOD, 
							result.file.prefix, distance=2000, 
							pval.thresh=0.10,
							min.length = NULL, max.length = NULL,
							default=FALSE, effectSize=FALSE,  fdr=FALSE,
							expr.file=NULL, Window="promoter"){
							
  message("Assessing overlap of DMRs with ", Window)
  if(effectSize){
  	 dmrs <- dmrs[order(-abs(dmrs$beta)),]
   }else{
   	 dmrs <- dmrs[order(-abs(dmrs$stat)),]
   }
    
   # start looking at expression data to explore correlation between methylation differences
   # and expression changes

   # if min.length or max.length are non-null, restrict the set of dmrs to 
   # associate by these parameters
   if (METHOD=="BSmooth"){
   	  	dmrs$L <- dmrs$n
   }else if (METHOD == "DSS"){
   	 	dmrs$L <- dmrs$nCG
   }
    
    
   if(default & (METHOD=="BSmooth" | METHOD=="DSS")){
   	   METHOD <- paste0(METHOD, "_default")
   } 
    
    subfolder <- "odds"
   if (effectSize){
   	 subfolder <- paste0(subfolder, "_effectSize")
   }else if(fdr){
   	 subfolder <- paste0(subfolder, "_fdr_Cumulative")
   }else{
   	 subfolder <- paste0(subfolder, "_statistic")
   }
   if (default & (grepl("BSmooth", METHOD) | grepl("DSS", METHOD))){
   	 subfolder <- paste0(subfolder, "_default")
   }
    
   if(!is.null(min.length)){
	  if(min.length > 100){  # interpret as #bps not nCpGs
	   dmrs$L2 <- dmrs$end - dmrs$start + 1
	   subfolder <- "odds_ln"
	  }else{
	   dmrs$L2 <- dmrs$L
	   subfolder <- "odds_nCpG"
	  }
   }
   
   if(!is.null(max.length)){
	  if(max.length > 100){  # interpret as #bps not nCpGs
	    dmrs$L2 <- dmrs$end - dmrs$start + 1
	    subfolder <- "odds_ln"
	  }else{
	   dmrs$L2 <- dmrs$L
	   subfolder <- "odds_nCpG"
	  }
   }
   
    ifelse(!dir.exists(paste0(result.file.prefix, "/", subfolder)), 
		dir.create(file.path(paste0(result.file.prefix, "/", subfolder))), FALSE)	

   
   if (!is.null(min.length)){
   		if (!is.null(max.length)){
   			if (max.length > min.length){
   			    if (sum(dmrs$L2 > min.length) > 0){
   					dmrs <- dmrs[dmrs$L2 > min.length & dmrs$L2 <= min.length,,drop=FALSE]
   				}else{
   					stop("Error: no DMRs of specified length")
   				}
   			}else{
   				stop("Error: max.length must be greater than min.length")
   			}
   		}else{
   		 	if (sum(dmrs$L2 > min.length) > 0){
   				dmrs <- dmrs[dmrs$L2 > min.length,,drop=FALSE]
   			}else{
   				stop("Error: no DMRs of specified length")
   			}
   		}
   }else if (!is.null(max.length)){
     dmrs <- dmrs[dmrs$L2 <= max.length,,drop=FALSE]
   }
   
   if (nrow(dmrs) < 500){
   	 stop(paste0("Error: not enough candidates with min.length ", min.length, 
   	 	"and max.length ", max.length))
   }
   
   message(paste0(nrow(dmrs), " dmrs being evaluated"))
 
   
   # start looking at expression data to explore correlation between methylation differences
   # and expression changes

   # use the text file downloaded from GEO (contains the counts of reads mapping
   # to each transcript obtained from cufflinks
   # file: /n/irizarryfs01/kkorthauer/WGBS/DNMT3A/RNAseq/GSE61969_cufflinks.gene.counts.txt
   # this file has genes in rows, with the gene name in the first column, and each subsequent
   # column is a sample (with a header containing the sample names of each column in the first
   # line.  Some of the samples represent the combination of multiple technical replicates.

   expr <- read.table(expr.file,
   				header=TRUE, stringsAsFactors=FALSE)
   	
   rownames(expr) <- expr$Gene
   expr <- expr[,-1]
   	
   # colIDs for files selected time points and tissue types
   if (!exists("time2")){
	   time2 <- time
   }else if (is.null(time2)){
	   time2 <- time
   }

	# remove KI (knock-in) samples (don't have methylation data for these)
   expr <- expr[, -grep("KI", colnames(expr))]	
   
   # remove the KO + WT FLT2 samples (don't have methylation data for these)
   expr <- expr[, -grep("KO3a_GFP", colnames(expr))]

   colnames(expr)[grep("GFP", colnames(expr))] <- 
   				paste0(colnames(expr)[grep("GFP", colnames(expr))], "_WTFL")
    	
   c1a <- grep(substr(tiss1,1,2), colnames(expr)) 
   c1b <- grep(substr(tiss1,4,7), colnames(expr))
   c1 <- c1a[c1a %in% c1b]
   
   c2a <- grep(substr(tiss2,1,2), colnames(expr)) 
   c2b <- grep(substr(tiss2,4,7), colnames(expr))
   c2 <- c2a[c2a %in% c2b]
 
# there has to be at least one sample of each type in order to go on
if (length(c1)>0 & length(c2)>0){
   tiss <- c(rep(tiss1, length(c1)),
			 rep(tiss2, length(c2)))
		  
   expr <- expr[,c(c1,c2)]
   colnames(expr) <- tiss

   # calculate a fold change from tiss1 to tiss2 for each of the genes in the 
   # expr matrix

   expr.ct <- ceiling(expr)
   sampleNames <- paste0(colnames(expr.ct), c(1:length(c1), 1:length(c2)))
   colData <- data.frame(samp=sampleNames, rep=factor(c(1:length(c1), 1:length(c2))),
						 condition=factor(colnames(expr.ct)))
   rownames(colData) <- sampleNames
   colnames(expr.ct) <- sampleNames
   condition <- colData$cond
   dds <- DESeqDataSetFromMatrix(countData = expr.ct,
								 colData = colData,
								 design = ~ condition)
   dds <- DESeq(dds)
   log2FC <- results(dds)$log2FoldChange
   if (sum(is.na(log2FC)) > 0){
   		log2FC[is.na(log2FC)] <- 0
   }
   p.FC <- results(dds)$padj
   if (sum(p.FC < pval.thresh, na.rm=TRUE)==0){
	   return(NULL)
   }   
   names(log2FC) <- rownames(expr)
   
   message(paste0("Number of genes total: ", sum(!is.na(p.FC))))
   message(paste0("Number of DE genes total: ", sum(p.FC < pval.thresh, na.rm=TRUE)))
    
   # associate each row of dmr object with log2FC of genes within
   # 2kb of TSS, if any
   dmrs$chr <- as.character(dmrs$chr)
   if (length(grep("chr", dmrs$chr)) == 0 ){
   		dmrs$chr <- paste0("chr", dmrs$chr)
   }

   # Get CpG and promoter / TSS / gene body annotation  
   annot_CpG = c('mm10_cpgs')
    for (attempt in 1:5){
	  annotations = try(build_annotations(genome = 'mm10', annotations = annot_CpG), 
				  silent=TRUE)  
	  if(!is(annotations, 'try-error')){
	  	  message("Download of CpG annotation successful!")
		  break;
	  }else{
		  message(paste0("Download of CpG annotation failed. ",
		          5-attempt, " attempts left"))
	  }
   }
  
   # adjust annot based on window type
   if (Window=="genebody"){
      for (attempt in 1:5){
	     annot = try(build_annotations(genome = 'mm10', 
                   annotations = c('mm10_genes_5UTRs','mm10_genes_cds','mm10_genes_3UTRs')), 
				 silent=TRUE)  
		 annot <- unlist(reduce(split(annot, annot$symbol),min.gapwidth=1e8))
		 annot$symbol <- names(annot)
		 
	   if(!is(annot, 'try-error')){
		 message("Download of gene body annotation successful!")
		 break;
	   }else{
		 message(paste0("Download of gene body annotation failed. ",
				 5-attempt, " attempts left"))
	   }
     }
     
   }else if(Window=="islandshore" | Window=="promoter"){
    for (attempt in 1:5){
	 annot = try(build_annotations(genome = 'mm10', 
                   annotations = 'mm10_genes_promoters'), 
				 silent=TRUE)  
	 if(!is(annot, 'try-error')){
		 message("Download of Promoter annotation successful!")
		 break;
	 }else{
		 message(paste0("Download of Promoter annotation failed. ",
				 5-attempt, " attempts left"))
	 }
    }
         
     if(Window=="islandshore"){
       annot <- trim(resize(annot,  width = distance, fix = 'end'))
     }
    
     annot <- unlist(reduce(split(annot, annot$symbol)))
     annot$symbol <- names(annot)

   }else{
     stop("please specify a valid window (from genebody, islandshore or promoter)")
   }
   
   x <- match(annot$symbol,rownames(expr))   
   annot$log2FC <- log2FC[x]
   annot$p.FC <- p.FC[x]
   annot <- annot[!is.na(annot$p.FC),]
   
   centers <- makeGRangesFromDataFrame(dmrs)

   dm_annotated = annotate_regions(regions = centers,
   	annotations = annotations, # CpG annotations from annotatr
   	ignore.strand = TRUE,
   	quiet = FALSE)
   dm_annotated$type <- dm_annotated@elementMetadata@listData$annot$type
   
   m <- findOverlaps(centers, dm_annotated)
   tmp <- dm_annotated$type[m@to]
   tmp <- split(tmp, m@from)
   tmp <- sapply(tmp, function(x) paste(x,collapse=""))
   
   dmrs$CpGcat <- NA
   dmrs$CpGcat[grepl("inter", tmp)] <- "Open Sea"
   dmrs$CpGcat[grepl("shelves", tmp)] <- "Shelf"
   dmrs$CpGcat[grepl("shores", tmp)] <- "Shore"
   dmrs$CpGcat[grepl("islands", tmp)] <- "Island"
   rm(tmp)
 
 topRanks <- c(500,1000,2000,5000,10000,15000,25000,50000,75000,100000,
 				125000,150000,175000,200000)
 topRanks <- c(topRanks[topRanks<nrow(dmrs)], nrow(dmrs))
 
 if (fdr){
 	topRanks <- c(1000,2000,5000,10000,25000,50000,100000,
 				150000,200000,250000,300000,400000)
    if (cond %in% c("KO_FLT3WT_FLT3")){
 		topRanks <- c(2000,4000,8000,25000,50000,100000,
 				150000,200000,250000,300000,400000)
 	}
    smallest <- min(topRanks)
    
    if(sum(dmrs$qval==1, na.rm=TRUE) > 0){
    	topRanks <- topRanks[topRanks < min(which(dmrs$qval==1))]
    }
     
    topRanks <- topRanks[topRanks < nrow(dmrs)]
        
    if (length(topRanks) == 0){
    	topRanks <- round(smallest/2)
    }
	
	if (nrow(dmrs)-topRanks[length(topRanks)] < 10000 & length(topRanks) > 1){
    	topRanks <- topRanks[-length(topRanks)]
    }

    topRanks <- c(topRanks, nrow(dmrs))
 	fdrRanks <- dmrs$qval[topRanks]
 	fdrDiff <- diff(fdrRanks)
 	if(sum(fdrDiff==0, na.rm=T) > 0){
 		topRanks <- topRanks[-(which(fdrDiff==0)+1)]
 		fdrRanks <- fdrRanks[-(which(fdrDiff==0)+1)]
 	}
    
 	fdrRanks[length(fdrRanks)] <- 1
 }
 

 odd.rank <- odd.rank.sig <- CI.low <- CI.hi <- propDE <- rep(NA, length(topRanks))
 pIslands <- pShores <- pShelves <- pOpen <- rep(NA, length(topRanks))
 deIslands <- deShores <- deShelves <- deOpen <- nDML <- lDMR <- rep(NA, length(topRanks))
 ct <- 1
 for(rnk in topRanks){
  if (fdr & ct > 1 & 2 < 1){
   		st <- topRanks[ct-1]
   		if (fdrRanks[ct] > 0.10){
   			st <- min(topRanks[fdrRanks>=0.10])
   		}
   		if (st == topRanks[ct]){
   			st <- topRanks[ct-1]
   		}
   }else{
   		st <- 1
   }
   
   centers <- makeGRangesFromDataFrame(dmrs[st:rnk,])
   overlaps <- findOverlaps(centers, annot)
    
   nDML[ct] <- sum(dmrs[1:rnk,]$L)
   
   if(ct==1){
   		lDMR[ct] <- mean(dmrs[1:rnk,]$L)
   }else{
   		lDMR[ct] <- mean(dmrs[(topRanks[ct-1]):rnk,]$L)
   }
   
    # Intersect the regions we read in with the annotations
   	topDMRs <- dmrs[st:rnk,]

   # Split by type
   
   regionTab <- table(topDMRs$CpGcat)/nrow(topDMRs)
     	
   if (length(grep("Open Sea", names(regionTab)))>0){
   	pOpen[ct] <- regionTab[grep("Open Sea", names(regionTab))]
   }else{
   	pOpen[ct] <- 0
   }
   
   if (length(grep("Shelf", names(regionTab)))>0){
   	pShelves[ct] <- regionTab[grep("Shelf", names(regionTab))]
   }else{
   	pShelves[ct] <- 0
   }
   
   if (length(grep("Shore", names(regionTab)))>0){
   	pShores[ct] <- regionTab[grep("Shore", names(regionTab))]
   }else{
    pShores[ct] <- 0
   }
   
   if (length(grep("Island", names(regionTab)))>0){
   	pIslands[ct] <- regionTab[grep("Island", names(regionTab))]
   }else{
	pIslands[ct] <- 0
   }

  if(length(overlaps) > 0){  
   
   meth.diff <- topDMRs$beta[overlaps@from]
   expr.diff <- annot$log2FC[overlaps@to]
   
   pvals <- annot$p.FC[overlaps@to]
   pvals[is.na(pvals)] <- 1
   
   corr.me <- round(cor(meth.diff, expr.diff, use="pairwise", method="spearman"), 4)
   p.me <- cor.test(meth.diff, expr.diff, method="spearman")$p.value
   
   a <- sum(expr.diff > 0 & meth.diff < 0) + 1
   b <- sum(expr.diff > 0 & meth.diff > 0) + 1
   c <- sum(expr.diff < 0 & meth.diff < 0) + 1
   d <- sum(expr.diff < 0 & meth.diff > 0) + 1
  
   or.me <- oddsratio(a,b,c,d)[["estimate"]]
   ci.me <- oddsratio(a,b,c,d)[["conf.int"]]
   
   odds.me <- (a+d) / (b+c) 
   odds.me.CI <- c( exp(log(odds.me) - 1.96*sqrt((odds.me + 1)^2/(odds.me*sum(a,b,c,d)))),
				   exp(log(odds.me) + 1.96*sqrt((odds.me + 1)^2/(odds.me*sum(a,b,c,d)))) )
   
   
   
 if (sum(pvals < pval.thresh, na.rm=TRUE) > 1){
    
    sigGenes <- annot[annot$p.FC < pval.thresh,]

    propDE[ct] <- sum(countOverlaps(centers,
    			sigGenes)>0) / length(st:rnk)
     dmrs.DE <- topDMRs[countOverlaps(makeGRangesFromDataFrame(topDMRs), 
    			sigGenes)>0,]
   
    regionALLTab <- regionTab*nrow(topDMRs)
    regionDETab <-  table(dmrs.DE$CpGcat)
    
    tmp <- rep(0, length(regionALLTab))
    names(tmp) <- names(regionALLTab)
    for(k in 1:length(regionDETab)){
   		tmp[names(tmp) == names(regionDETab)[k]] <- regionDETab[k]
    }
    regionDETab <- tmp
    rm(tmp)		
   		
    regionDETab <- regionDETab / regionALLTab
  
   if (length(grep("Open Sea", names(regionDETab)))>0){
   	deOpen[ct] <- (regionDETab[grep("Open Sea", names(regionDETab))])
   }else{
   	deOpen[ct] <- 0
   }
   
   if (length(grep("Shelf", names(regionDETab)))>0){
   	deShelves[ct] <- (regionDETab[grep("Shelf", names(regionDETab))])
   }else{
   	deShelves[ct] <- 0
   }
   
   if (length(grep("Shore", names(regionDETab)))>0){
   	deShores[ct] <- (regionDETab[grep("Shore", names(regionDETab))])
   }else{
    deShores[ct] <- 0
   }
   
   if (length(grep("Island", names(regionDETab)))>0){
   	deIslands[ct] <- (regionDETab[grep("Island", names(regionDETab))])
   }else{
	deIslands[ct] <- 0
   }
   
   
   	corr.me.de <- round(cor(meth.diff[pvals < pval.thresh], 
   						   expr.diff[pvals < pval.thresh], use="pairwise", method="spearman"), 4)
   	p.me.de <- cor.test(meth.diff[pvals < pval.thresh], 
   							expr.diff[pvals < pval.thresh], method="spearman")$p.value
   
    a <- sum(expr.diff[pvals < pval.thresh] > 0 & meth.diff[pvals < pval.thresh] < 0) + 1
    b <- sum(expr.diff[pvals < pval.thresh] > 0 & meth.diff[pvals < pval.thresh] > 0) + 1
    c <- sum(expr.diff[pvals < pval.thresh] < 0 & meth.diff[pvals < pval.thresh] < 0) + 1
    d <- sum(expr.diff[pvals < pval.thresh] < 0 & meth.diff[pvals < pval.thresh] > 0) + 1
   
   	or.me.de <- oddsratio(a,b,c,d)[["estimate"]]
   	ci.me.de <- oddsratio(a,b,c,d)[["conf.int"]]
   	
   	odds.me.de <- (a+d) / (b+c) 
   	odds.me.de.CI <- c( exp(log(odds.me.de) - 1.96*sqrt((odds.me.de + 1)^2/(sum(a,b,c,d)*odds.me.de))),
   					exp(log(odds.me.de) + 1.96*sqrt((odds.me.de + 1)^2/(sum(a,b,c,d)*odds.me.de))) )
 
   }else{
   	corr.me.de <- NA
   	p.me.de <- 1
   	or.me.de <- NA
   	por.me.de <- NA
   	odds.me <- odds.me.de <- NA
   	odds.me.de.CI <- NA

   }
   
   
   print(paste0("Correlation of methylation and expression: ", corr.me))
   print(paste0("p-value : ", p.me))
   print(paste0("Correlation of methylation and expression (DE): ", corr.me.de))
   print(paste0("p-value : ", p.me.de))
   print(paste0("Number of genes neighboring DMRs: ", length(expr.diff)))
   print(paste0("Number of DE genes neighboring DMRs: ", sum(pvals < pval.thresh)))

   print(paste0("Odds (CI): ", round(odds.me,3), 
   				"(", round(odds.me.CI[1], 3), "-", round(odds.me.CI[2], 3), ")" ))
   print(paste0("DE Odds (CI): ", round(odds.me.de,3), 
   				"(", round(odds.me.de.CI[1], 3), "-", round(odds.me.de.CI[2], 3), ")" ))
   
   print(design)
   
   odd.rank[ct] <- odds.me.de
   odd.rank.sig[ct] <- sum(odds.me.de.CI > 1)==2
   CI.low[ct] <- odds.me.de.CI[1] 
   CI.hi[ct] <- odds.me.de.CI[2]
   
   xlimits <- range(meth.diff)
      
   save(meth.diff, expr.diff, file=paste0(result.file.prefix, "/", subfolder, "/corrWithExpression.", METHOD, ".n", 
		   sampleSize, ".", cond, ".", num.dmrs, "DMRS.top", rnk, "_min",
		   min.length, "_max", max.length, ".", Window, ".RData"))
   # correlate log2 foldchanges with DMR test statistic
   diff.df <- data.frame(meth.diff=meth.diff, expr.diff=expr.diff,
   						 DE = pvals < pval.thresh)
   diff.df <- diff.df[order(diff.df$DE),]
   diff.df$DE <- as.factor(diff.df$DE)
   levels(diff.df$DE) <- list(True="TRUE", False="FALSE")
   
   pdf(paste0(result.file.prefix, "/", subfolder, "/corrWithExpressionPlot.", METHOD, ".n", 
		   sampleSize, ".", cond, ".", num.dmrs, "DMRS.top", rnk, "_min",
		   min.length, "_max", max.length, ".", Window, ".pdf"),
		   height=5, width=6)
		gp <- ggplot(aes(x=meth.diff, y=expr.diff, fill=DE, colour=DE), 
							  data=diff.df) + 
			geom_vline(xintercept=0, colour="grey", alpha=0.8, lwd=1.5) + 
   			geom_hline(yintercept=0, colour="grey", alpha=0.8, lwd=1.5) +
			geom_point(alpha=0.3) +
   			xlab("Methylation difference (statistic)") +
   			ylab("Expression difference (log2 FC estimate)") +
   			scale_fill_manual(values=c("red", "black")) +
   			scale_colour_manual(values=c("red", "black")) +
   			ggtitle(paste0("Corr (DE only): ", round(corr.me,3), "(", round(corr.me.de,3), rep("*", as.numeric(p.me.de < 0.05)), ")",
   			               ", DE Odds: ", round(odds.me.de,3), "(", round(odds.me.de.CI[1],3), "-", round(odds.me.de.CI[2],3), ")" )) 	+ 
   			theme_classic()
   		print(gp)
   dev.off()
   
      save(gp, file=paste0(result.file.prefix, "/", subfolder, "/corrWithExpressionPlot.", 
    		METHOD, ".n", sampleSize, ".", cond, ".", num.dmrs, "DMRS.top", rnk, ".", Window, ".RData"))
 }else{
    message("No overlaps between DMRs and TSS within ", distance, " bps.")
 }
    ct <- ct + 1   
 }
 }
 
 
if (fdr){
 	topRanks <- fdrRanks
} 
 
rank.df <- data.frame(ranks=topRanks, odds = odd.rank, 
							  Significant = odd.rank.sig,
							  CI.low = CI.low,
							  CI.hi = CI.hi,
							  propDE = propDE,
							  Open = pOpen,
							  Islands = pIslands,
							  Shores = pShores,
							  Shelves = pShelves,
							  OpenDE = deOpen,
							  IslandsDE = deIslands,
							  ShoresDE = deShores,
							  ShelvesDE = deShelves,
							  nDML = nDML,
							  lDMR = lDMR)    
save(rank.df, file=paste0(result.file.prefix, "/", subfolder, "/oddsByRanking.", METHOD, ".n", 
		   sampleSize, ".", cond, ".", num.dmrs, "DMRS.top", rnk, "_min",
		   min.length, "_max", max.length, ".", Window,  ".RData"))
	rank.df <- rank.df[rank.df$ranks <= nrow(dmrs),]
    
	pdf(paste0(result.file.prefix, "/", subfolder, "/oddsByRankingPlot.", METHOD, ".n", 
		   sampleSize, ".", cond, ".", num.dmrs, "DMRS.", "_min",
		   min.length, "_max", max.length, ".", Window, ".pdf"),
		   height=4, width=7)
		gp <- ggplot(aes(x=ranks, y=odds), 
							  data=rank.df) + 
   			geom_hline(yintercept=1, colour="grey", alpha=0.8, lwd=1.5, lty=2) +
   			geom_line(color="grey") +
			geom_point(alpha=0.9, aes(color=Significant)) +
   			xlab("Number of top-ranked DMRs") +
   			ylab("Odds of DE genes in discordant direction") +
   			scale_fill_manual(values=c("black", "red")) +
   			scale_colour_manual(values=c("black", "red")) +
   			ggtitle("Odds of DE discordant with Methylation by number of top DMRs") + 
   			theme_classic()
   		print(gp)
   dev.off()
   
   # for plotting when calling BSmooth or DSS - need to hard specify which 
   # version of the dmrseq results to use for plotting. 
   if (grepl("DSS", METHOD) | grepl("BSmooth", METHOD) | grepl("metilene", METHOD)){
   		pfile1 <- list.files(paste0(result.root.dir,"/dmrseq_pkg/", subfolder, "/"),
   						pattern=paste0("oddsByRanking.dmrseq.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
   		pfile2 <- list.files(paste0(result.root.dir,"/dmrseq_pkg/", subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", Window, ".RData")) 		  
   
		if (length(pfile1[pfile1 %in% pfile2]) == 1){	   					 					
		 pfile <- paste0(paste0(result.root.dir,"/dmrseq_pkg/", subfolder, "/"), 
				   pfile1[pfile1 %in% pfile2])
		}else if (length(pfile1[pfile1 %in% pfile2]) > 1){	
			thisOne <- grep(nrow(dmrs), pfile1[pfile1 %in% pfile2])
			pfile <- paste0(result.file.prefix, "/", subfolder, "/", pfile1[pfile1 %in% pfile2][thisOne])
		}else{
		   pfile <- "null"
		}
		rank.df$Method <- METHOD
   }else{
   # for plotting when calling dmrseq
   		pfile1 <- list.files(path=paste0(result.file.prefix, "/", subfolder, "/"),
   						pattern=paste0("oddsByRanking.dmrseq.n", 
		   				sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
  	 	pfile2 <- list.files(path=paste0(result.file.prefix, "/", subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", Window, ".RData")) 
   		if (length(pfile1[pfile1 %in% pfile2]) == 1){	   					 					
   			pfile <- paste0(result.file.prefix, "/", subfolder, "/", 
   		    pfile1[pfile1 %in% pfile2])
		}else if (length(pfile1[pfile1 %in% pfile2]) > 1){	  #if > 1, pick most recently modified
			details <- file.info(pfile1[pfile1 %in% pfile2])
			details <- details[with(details, order(as.POSIXct(mtime), decreasing=TRUE)), ]
			pfile <- paste0(result.file.prefix, "/", subfolder, "/", rownames(details)[1])
			
			thisOne <- grep(nrow(dmrs), pfile1[pfile1 %in% pfile2])
			pfile <- paste0(result.file.prefix, "/", subfolder, "/", pfile1[pfile1 %in% pfile2][thisOne])
		}else{
		   pfile <- "null"
		}
		load(pfile)
		rank.df$Method <- "dmrseq"		
		ranks.last <- rank.df[nrow(rank.df),,drop=FALSE]
		ranks.last$Method <- "dmrseq"  
   }
   
   if(fdr){
       bfile1 <- list.files(path=paste0(result.root.dir,"/BSmooth_default/", "odds_statistic_default", "/"),
							pattern=paste0("oddsByRanking.BSmooth_default.n", 
								sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
	   bfile2 <- list.files(path=paste0(result.root.dir,"/BSmooth_default/", "odds_statistic_default", "/"),
							pattern=paste0("_min", min.length, "_max", max.length, ".", Window, ".RData")) 		  
	   if (length(bfile1[bfile1 %in% bfile2]) == 1){	   					 					
	   bfile <- paste0(result.root.dir,"/BSmooth_default/", "odds_statistic_default", "/", 
				  bfile1[bfile1 %in% bfile2])
	   }else{
		  bfile <- "null"
	   }

	   dfile1 <- list.files(path=paste0(result.root.dir,"/DSS_default/", "odds_statistic_default", "/"),
							pattern=paste0("oddsByRanking.DSS_default.n", 
								sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
	   dfile2 <- list.files(path=paste0(result.root.dir,"/DSS_default/", "odds_statistic_default", "/"),
							pattern=paste0("_min", min.length, "_max", max.length, ".", Window,  ".RData")) 		  
	   if (length(dfile1[dfile1 %in% dfile2]) == 1){	   					 					
	   dfile <- paste0(result.root.dir,"/DSS_default/", "odds_statistic_default", "/", 
				  dfile1[dfile1 %in% dfile2])
	   }else{
		  dfile <- "null"
	   }
   }else{		   
	   bfile1 <- list.files(path=paste0(result.root.dir,"/BSmooth/", subfolder, "/"),
							pattern=paste0("oddsByRanking.BSmooth.n", 
								sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
	   bfile2 <- list.files(path=paste0(result.root.dir,"/BSmooth/", subfolder, "/"),
							pattern=paste0("_min", min.length, "_max", max.length, ".", Window,  ".RData")) 		  
	   if (length(bfile1[bfile1 %in% bfile2]) == 1){	   					 					
	   bfile <- paste0(result.root.dir,"/BSmooth/", subfolder, "/", 
				  bfile1[bfile1 %in% bfile2])
	   }else{
		  bfile <- "null"
	   }

	   dfile1 <- list.files(path=paste0(result.root.dir,"/DSS/", subfolder, "/"),
							pattern=paste0("oddsByRanking.DSS.n", 
								sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
	   dfile2 <- list.files(path=paste0(result.root.dir,"/DSS/", subfolder, "/"),
							pattern=paste0("_min", min.length, "_max", max.length,".", Window,  ".RData")) 		  
	   if (length(dfile1[dfile1 %in% dfile2]) == 1){	   					 					
	   dfile <- paste0(result.root.dir,"/DSS/", subfolder, "/", 
				  dfile1[dfile1 %in% dfile2])
	   }else{
		  dfile <- "null"
	   }
   }
   
   pdf(paste0(result.file.prefix, "/", subfolder, "/lengthDistribution.", METHOD, ".n", 
		   sampleSize, ".", cond, ".", num.dmrs, "DMRS.", "_min",
		   min.length, "_max", max.length, ".", Window, ".pdf"),
		   height=4, width=7)
		gp <- ggplot(data=dmrs, aes(log(L))) +
				geom_density() +
	   		    xlab("Log dmr length (nCpGs)") +
   			    ggtitle("Region Length Distribution") + 
   				theme_classic()
   		print(gp)
   dev.off()
   
 if (METHOD == "dmrseq"){
   
   if(file.exists(pfile) | file.exists(bfile) | file.exists(dfile)){
   		message("plotting OR comparisons")
   		
   		if(file.exists(pfile)){
		  load(pfile)
		  if (fdr){
 			rank.df$ranks <- fdrRanks
		  } 
 
		  ranks.all <- rank.df
		  ranks.last <- rank.df[nrow(rank.df),,drop=FALSE]
		  ranks.all$Method <- "dmrseq"
		  ranks.last$Method <- "dmrseq"
		}
				
   		if(file.exists(bfile)){
		  load(bfile)
		  coldiff <- ncol(ranks.all)-ncol(rank.df)
		  if (coldiff>1){
		  	for (m in 1:(coldiff-1)){
		  		rank.df <- cbind(rank.df, rep(NA, nrow(rank.df)))
		  	}
		  }
		  rank.df$Method <- "BSmooth"
		  last <- rank.df[nrow(rank.df),,drop=FALSE]
		  last$Method <- "BSmooth"
		  
		  if (fdr){
		  	rank.df <- last
		  }	
		  
		  if (!exists("ranks.all")){
		  	ranks.all <- rank.df
		  }else{
		  	ranks.all <- rbind(ranks.all, rank.df)
		  	ranks.last <- rbind(ranks.last, last)
		  }
		}	   
   		
   		if(file.exists(dfile)){
		  load(dfile)
		  coldiff <- ncol(ranks.all)-ncol(rank.df)
		  if (coldiff>1){
		  	for (m in 1:(coldiff-1)){
		  		rank.df <- cbind(rank.df, rep(NA, nrow(rank.df)))
		  	}
		  }
		  
		  rank.df$Method <- "DSS"
		  last <- rank.df[nrow(rank.df),,drop=FALSE]
		  last$Method <- "DSS"
		  
		  if (fdr){
		  	rank.df <- last
		  }	
		  if (!exists("ranks.all")){
		  	ranks.all <- rank.df
		  }else{
		  	ranks.all <- rbind(ranks.all, rank.df)
		  	ranks.last <- rbind(ranks.last, last)
		  }
		}
		
   		ranks.all <- ranks.all[order(ranks.all$Method),]
   		ranks.last <- ranks.last[order(ranks.last$Method),]
   		
   	pdf(paste0(result.file.prefix, "/", subfolder, "/oddsByRankingPlot.COMPARISON.n", 
		   sampleSize, ".", cond, ".", num.dmrs, "DMRS.",  "_min",
		   min.length, "_max", max.length, ".", Window, ".pdf"),
		   height=4, width=7)
    
    if(!fdr){ 
		gp <- ggplot(aes(x=ranks, y=log2(odds)), 
							  data=ranks.all) + 
   			geom_hline(yintercept=0, colour="grey", alpha=0.6, lwd=1.5, lty=2) +
   			geom_line(aes(group=Method, color=Method), lwd=1.1) +
			geom_point(alpha=0.6, aes(color=Method)) +
			geom_errorbar(aes(ymin=log2(CI.low), ymax=log2(CI.hi), colour=Method),
                 alpha=0.7, size=1.1, width=0) +
   			xlab("Number of top-ranked DMRs") +
   			ylab("Log2 Odds of DE genes in discordant direction") +
   			#scale_fill_manual(values=c("black", "red")) +
   			#scale_colour_manual(values=c("black", "red")) +
   			ggtitle("Odds of DE discordant with Methylation by number of top DMRs") + 
   			theme_classic()
   		gpl <- ggplot(aes(x=log2(ranks), y=log2(odds)), 
							  data=ranks.all) + 
   			geom_hline(yintercept=0, colour="grey", alpha=0.6, lwd=1.5, lty=2) +
   			geom_line(aes(group=Method, color=Method), lwd=1.1) +
			geom_point(alpha=0.6, aes(color=Method)) +
			geom_errorbar(aes(ymin=log2(CI.low), ymax=log2(CI.hi), colour=Method),
                 alpha=0.7, size=1.1, width=0) +
   			xlab("Log2 Number of top-ranked DMRs") +
   			ylab("Log2 Odds of DE genes in discordant direction") +
   			#scale_fill_manual(values=c("black", "red")) +
   			#scale_colour_manual(values=c("black", "red")) +
   			ggtitle("Odds of DE discordant with Methylation by number of top DMRs") + 
   			theme_classic()	
   	    gpl2 <- ggplot(aes(x=log2(ranks), y=log2(odds)), 
							  data=ranks.all) + 
			geom_ribbon(aes(ymin=log2(CI.low), ymax=log2(CI.hi), fill=Method),
                 alpha=0.3) +
   			geom_hline(yintercept=0, colour="grey", alpha=0.6, lwd=1.5, lty=2) +
   			geom_line(aes(group=Method, color=Method), lwd=1.1) +
			#geom_point(alpha=0.6, aes(color=Method)) +
   			xlab("Log2 Number of top-ranked DMRs") +
   			ylab("Log2 Odds of DE genes in discordant direction") +
   			#scale_fill_manual(values=c("black", "red")) +
   			#scale_colour_manual(values=c("black", "red")) +
   			ggtitle("Odds of DE discordant with Methylation by number of top DMRs") + 
   			theme_classic()			
   				
   		gC <- ggplot(aes(x=log2(ranks), y=log2(nDML)), 
							  data=ranks.all) + 
   			geom_line(aes(group=Method, color=Method), lwd=1.1) +
			geom_point(alpha=0.6, aes(color=Method)) +
   			xlab("Log2 Number of top-ranked DMRs") +
   			ylab("Log2 Number of CpG loci covered in top-ranked DMRs") +
   			ggtitle("Number of DMRs by number of CpG Loci covered") + 
   			theme_classic()	
   			
   		gL <- ggplot(aes(x=log2(ranks), y=lDMR), 
							  data=ranks.all) + 
   			geom_line(aes(group=Method, color=Method), lwd=1.1) +
			geom_point(alpha=0.6, aes(color=Method)) +
   			xlab("Log2 Number of top-ranked DMRs") +
   			ylab("Mean Number of CpG loci in top-ranked DMRs") +
   			ggtitle("Mean length of DMRs") + 
   			theme_classic()		
   			
   		print(gp)
   		print(gpl)
   		print(gpl2)
   		print(gC)
   		print(gL)	
   		
   		}else if (fdr){
   		    ranks.all$Method <- as.factor(ranks.all$Method)
   		    ranks.all$ranks[ranks.all$ranks > 1] <- min(fdrRanks)
   		    add2 <- ranks.all[ranks.all$Method != "dmrseq",]
   		    if(nrow(add2)>0){
   		    	add2$ranks <- 1
   		    	ranks.all <- rbind(ranks.all, add2)
   		 	}
   		    ranks.all$fdr <- 1
   		    ranks.all$fdr[ranks.all$Method == "dmrseq"] <- 0
   		    ranks.all$fdr <- as.factor(ranks.all$fdr)
   		    
   			pD <- ggplot(aes(x=ranks, y=log2(odds)), 
							  data=ranks.all) + 
				geom_ribbon(aes(ymin=log2(CI.low), ymax=log2(CI.hi), fill=Method),
					 alpha=0.3) +
				#geom_hline(yintercept=0, colour="grey", alpha=0.6, lwd=1.5, lty=1) +
				#geom_vline(xintercept=0.10, colour="grey", alpha=0.6, lwd=1.2, lty=2) +
				geom_line(aes(group=Method, color=Method, linetype=fdr), lwd=1.1) +
				#geom_point(alpha=0.6, aes(color=Method)) +
				xlab("FDR threshold") +
				ylab("log2 Odds of DE genes in discordant direction") +
				#scale_fill_manual(values=c("black", "red")) +
				#scale_colour_manual(values=c("black", "red")) +
				ggtitle("Odds of DE discordant with Methylation by number of top DMRs") + 
				theme_classic() +
				guides(linetype=FALSE)
   		
   			print(pD)
   			
   			pDl <- ggplot(aes(x=log2(ranks), y=log2(odds)), 
							  data=ranks.all) + 
				geom_ribbon(aes(ymin=log2(CI.low), ymax=log2(CI.hi), fill=Method),
					 alpha=0.3) +
				#geom_hline(yintercept=0, colour="grey", alpha=0.6, lwd=1.5, lty=1) +
				#geom_vline(xintercept=0.10, colour="grey", alpha=0.6, lwd=1.2, lty=2) +
				geom_line(aes(group=Method, color=Method, linetype=fdr), lwd=1.1) +
				#geom_point(alpha=0.6, aes(color=Method)) +
				xlab("log2 FDR threshold") +
				ylab("log2 Odds of DE genes in discordant direction") +
				#scale_fill_manual(values=c("black", "red")) +
				#scale_colour_manual(values=c("black", "red")) +
				ggtitle("Odds of DE discordant with Methylation by number of top DMRs") + 
				theme_classic() +
				guides(linetype=FALSE)
   		
   			print(pDl)
   		}
   		   		
   		
   		if ("propDE" %in% colnames(ranks.all)){
   		
   			
   			p4 <-  ggplot(aes(x=Method, y=propDE), data=ranks.last) + 
   			geom_bar(stat="identity", aes(fill=Method)) +
   			ylab("Proportion of DMRs near a DE gene") +
   			ggtitle("Proportion of DMRs overlapping the promoter region of a DE gene") + 
   			theme_classic()
   			
   			CpG.cats <- reshape2::melt(ranks.last[,c(7:10,17)])
   			colnames(CpG.cats)[2] <- "CpG_Category"
   			
   			p5 <-  ggplot(aes(x=Method, y=value), data=CpG.cats) + 
   			geom_bar(stat="identity", aes(fill=CpG_Category)) +
   			ylab("Proportion of CpGs in DMRs") +
   			ggtitle("CpGs in DMRs annotated by CpG category") + 
   			theme_classic()
   			
   			DE.cats <- reshape2::melt(ranks.all[,c(11:14, 17)])
   			colnames(DE.cats)[2] <- "CpG_Category"
   			nchars <- sapply(as.character(DE.cats$CpG_Category), nchar)
   			DE.cats$CpG_Category <- substr(DE.cats$CpG_Category, 1, nchars-2)
   			
   			p6 <-  ggplot(aes(x=CpG_Category, y=value, group=Method), data=DE.cats) + 
   			geom_bar(stat="identity", aes(fill=Method), position="dodge") +
   			ylab("Proportion of DMRs near a DE gene") +
   			ggtitle("Proportion of DMRs near a DE gene by CpG category and Method") + 
   			theme_classic()
   			
   			DE.cats <- reshape2::melt(ranks.all[,c(1,11:14,17)], id=c("ranks", "Method"))
   			colnames(DE.cats)[3] <- "CpG_Category"
   			nchars <- sapply(as.character(DE.cats$CpG_Category), nchar)
   			DE.cats$CpG_Category <- substr(DE.cats$CpG_Category, 1, nchars-2)
   			
   			p7 <-  ggplot(aes(x=log2(ranks), y=value, colour=Method), 
   				data=DE.cats) + 
   			geom_line(aes(group=interaction(CpG_Category,Method)), 
   				alpha=0.6, lwd=1.1) +
   			facet_wrap( ~ CpG_Category, ncol=2) +
			geom_point(alpha=0.6, size=3) +
   			ylab("Proportion of covered CpGs near a DE gene") +
   			xlab("Log2 Number of Top-ranked DMRs") +
   			ggtitle("Proportion of DMR-covered CpGs near a DE gene by CpG category") + 
   			theme_classic()
   			
   			# add 95% confidence bounds
   		    DE.cats$CI.low <- DE.cats$value - 
   		    	1.96 * sqrt(DE.cats$value*(1-DE.cats$value) / DE.cats$ranks)
   		    DE.cats$CI.hi <-  DE.cats$value + 
   		    	1.96 * sqrt(DE.cats$value*(1-DE.cats$value) / DE.cats$ranks)
   		    	
   			p7b <-  ggplot(aes(x=log2(ranks), y=value), 
   				data=DE.cats) + 
   			geom_ribbon(aes(ymin=(CI.low), ymax=(CI.hi), fill=Method),
                 alpha=0.3) +
   			geom_line(aes(group=interaction(CpG_Category,Method), colour=Method), 
   				 lwd=1.1) +
   			facet_wrap( ~ CpG_Category, ncol=2) +
   			ylab("Proportion of covered CpGs near a DE gene") +
   			xlab("Log2 Number of Top-ranked DMRs") +
   			ggtitle("Proportion of DMR-covered CpGs near a DE gene by CpG category") + 
   			theme_classic()
   			
   			# add 95% confidence bounds
   		    ranks.all$CI.pDE.low <- ranks.all$propDE - 
   		    	1.96 * sqrt(ranks.all$propDE*(1-ranks.all$propDE) / ranks.all$ranks)
   		    ranks.all$CI.pDE.hi <-  ranks.all$propDE + 
   		    	1.96 * sqrt(ranks.all$propDE*(1-ranks.all$propDE) / ranks.all$ranks)
   		    
   			p3 <- ggplot(aes(x=log2(ranks), y=propDE), 
							  data=ranks.all) + 
   			geom_line(aes(group=Method, color=Method), lwd=1.1) +
			geom_point(alpha=0.6, aes(color=Method)) +
   			xlab("Log2 Number of top-ranked DMRs") +
   			ylab("Proportion of DMRs near a DE gene") +
   			ggtitle("Proportion of DMRs overlapping the promoter region of a DE gene") + 
   			geom_errorbar(aes(ymin=CI.pDE.low, ymax=CI.pDE.hi, colour=Method),
                 alpha=0.7, size=1.1, width=0) +
   			theme_classic()
   			
   			print(p3)
   			print(p4)
   			print(p5)
   			print(p6)
   			print(p7)
   			print(p7b)
   		}
   		
   dev.off()
   
   }}
   
   return(NULL)
}

