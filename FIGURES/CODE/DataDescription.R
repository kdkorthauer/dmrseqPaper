# generate summary figures and data for summary table
# to give main summary/data description
library(bsseq)
library(ggplot2)
library(reshape2)
library(cowplot)
result.file.prefix <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/dmrseqPaper/FIGURES/out/"

plotBetaDistribution <- function(bs, main=""){

	 meth.mat = getCoverage(bs, type = "M")
	 unmeth.mat = getCoverage(bs, type = "Cov") - meth.mat
     
     if("Cell" %in% colnames(pData(bs))){
     	colnames(meth.mat) <- colnames(unmeth.mat) <- pData(bs)$Cell
     }
     
     if("Condition" %in% colnames(pData(bs))){
     	colnames(meth.mat) <- colnames(unmeth.mat) <- paste0(pData(bs)$Condition,
     				"_", pData(bs)$Rep)
     }
      
	 #plot empirical distribution of meth levels (betas)
	 # random sample of 5M rows
	 rows.plot <- sample(1:nrow((meth.mat)), min(nrow((meth.mat)),5e6), replace=FALSE)
	 meth.levelsm <- data.frame(meth.mat[rows.plot,] / (meth.mat + unmeth.mat)[rows.plot,])
	 cov.matm <- data.frame((meth.mat + unmeth.mat)[rows.plot,])
	 
	 coverageBySample <- colSums(cov.matm)
	 weight.mat <- melt(t(apply(as.matrix(cov.matm), 1, function(x) x/coverageBySample)))
	 	 
	 meth.levelsm <- melt(meth.levelsm)
	 cov.matm <- melt(cov.matm)
     
     if (grepl("_", cov.matm$variable[1])){
	   meth.levelsm$Tissue <- unlist(sapply(strsplit(as.character(meth.levelsm$variable), "_"),
								  function(x) x[[1]]))
	   meth.levelsm$Replicate <- unlist(sapply(strsplit(as.character(meth.levelsm$variable), "_"),
								  function(x) x[[2]]))
	   meth.levelsm$Replicate[meth.levelsm$Replicate == "rep3"] <- "rep2"
	   meth.levelsm$Replicate[meth.levelsm$Replicate == "rep1"] <- "1"
	   meth.levelsm$Replicate[meth.levelsm$Replicate == "rep2"] <- "2"
	 
	   cov.matm$Tissue <- unlist(sapply(strsplit(as.character(cov.matm$variable), "_"),
								  function(x) x[[1]]))
	   cov.matm$Replicate <- unlist(sapply(strsplit(as.character(cov.matm$variable), "_"),
								  function(x) x[[2]]))
	   cov.matm$Replicate[cov.matm$Replicate == "rep3"] <- "rep2"
	   cov.matm$Replicate[cov.matm$Replicate == "rep1"] <- "1"
	   cov.matm$Replicate[cov.matm$Replicate == "rep2"] <- "2"

	   meth.levelsm$Tissue <- gsub("\\.", "_", as.character(meth.levelsm$Tissue))
	   cov.matm$Tissue <- gsub("\\.", "_", as.character(cov.matm$Tissue))
	   
     }else{
       meth.levelsm$Cell <- as.character(meth.levelsm$variable)
	   cov.matm$Cell <- as.character(cov.matm$variable)
	   meth.levelsm$Replicate <- meth.levelsm$Cell
	   cov.matm$Replicate <- cov.matm$Cell
     }
     
     meth.levelsm$Coverage <- weight.mat$value
     
     if (length(unique(meth.levelsm$Tissue))>1 & !("Cell" %in% colnames(meth.levelsm))){
		 p0 <- ggplot(meth.levelsm, aes(value, colour=Tissue, 
					  group=interaction(Tissue, Replicate))) +
					geom_density(aes(linetype=Replicate, weight=Coverage), 
							adjust=1.75, 
							alpha=0.6, size=1)+
					xlab("Methylation Proportion") + 
					ggtitle(paste0(main)) +
					#geom_hline(yintercept=0, colour="white", size=1.6) +
					#geom_vline(xintercept=0, colour="white", size=1.6) + 
					#geom_vline(xintercept=1, colour="white", size=1.6) +
					theme_classic() + 
					theme(plot.title = element_text(face="bold"))# colour=condition
		 p1 <- ggplot(cov.matm, aes(log(value), colour=Tissue, 
					  group=interaction(Tissue, Replicate))) +
					geom_line(aes(linetype=Replicate), adjust=2.5, 
							alpha=0.6, stat="density", size=1)+
					xlab("Coverage (log scale)") + 
					ggtitle(paste0(main)) +
					theme_classic() + 
					theme(plot.title = element_text(face="bold"))# colour=condition
	    if("Condition" %in% colnames(pData(bs))){	
		    meth.levelsm$Condition <- meth.levelsm$Tissue
		    cov.matm$Condition <- cov.matm$Tissue
	       
				 p0 <- ggplot(meth.levelsm, aes(value, colour=Condition, 
					  group=interaction(Condition, Replicate))) +
					geom_density(aes(linetype=Replicate, 
					weight=Coverage), adjust=4.5, 
							alpha=0.6, size=1)+
					xlab("Methylation Proportion") + 
					#geom_hline(yintercept=0, colour="white", size=1.6) +
					#geom_vline(xintercept=0, colour="white", size=1.6) + 
					#geom_vline(xintercept=1, colour="white", size=1.6) +
					ggtitle(paste0(main)) +
					theme_classic() + 
					theme(plot.title = element_text(face="bold"))# colour=condition
		 		p1 <- ggplot(cov.matm, aes(log(value), colour=Condition, 
					  group=interaction(Condition, Replicate))) +
					geom_line(aes(linetype=Replicate), adjust=2.5, 
							alpha=0.6, stat="density", size=1)+
					xlab("Coverage (log scale)") + 
					ggtitle(paste0(main)) +
					theme_classic() + 
					theme(plot.title = element_text(face="bold"))# colour=condition
		 }
	   		
	   }else{
	     meth.levelsm$Cell <- meth.levelsm$Replicate
	     cov.matm$Cell <- cov.matm$Replicate
	     
		 p0 <- ggplot(meth.levelsm, aes(value, colour=Cell, 
					   group=Cell)) +
					 geom_density(adjust=4, aes(weight=Coverage),
							 alpha=0.6, size=1)+
					 xlab("Methylation Proportion") + 
					 ggtitle(paste0(main)) +
					 #geom_hline(yintercept=0, colour="white", size=1.6) +
					 #geom_vline(xintercept=0, colour="white", size=1.6) + 
					 #geom_vline(xintercept=1, colour="white", size=1.6) +
					 theme_classic()+ 
					 theme(plot.title = element_text(face="bold")) # colour=condition
					 
		  p1 <- ggplot(cov.matm, aes(log(value), colour=Cell, 
					   group=Cell)) +
					 geom_line(adjust=2.5, 
							 alpha=0.6, stat="density", size=1)+
					 xlab("Coverage (log scale)") + 
					 ggtitle(paste0(main)) +
					 theme_classic() + 
					 theme(plot.title = element_text(face="bold"))# colour=condition

	   }		  
	   
	rm(meth.levelsm) 
	rm(rows.plot)
	return(list(p0, p1))
}

getCoverageStats <- function(bs){
	cov.mat <- getCoverage(bs, type = "Cov")
	med.cov <- apply(cov.mat, 2, median)
	max.cov <- apply(cov.mat, 2, max)
	return(list(med.cov, max.cov))
}


			
summary.table <- vector("list", 6)

# Roadmap data - 4 tissues
data.file.prefix <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/"

	 
# Left Ventricle
load(paste0(data.file.prefix, "roadmap_Left Ventricle_all_bsseq_id.RData")); 
meth.levels.raw = getMeth(bs, type = "raw")
no.hits = which(is.na(rowMeans(meth.levels.raw)) == TRUE)
bs = bs[-no.hits]
LV <- bs 
rm(bs); rm(meth.levels.raw); gc();

# Right Ventricle
load(paste0(data.file.prefix, "roadmap_Right Ventricle_all_bsseq_id.RData")); 
meth.levels.raw = getMeth(bs, type = "raw")
no.hits = which(is.na(rowMeans(meth.levels.raw)) == TRUE)
bs = bs[-no.hits]
RV <- bs 
rm(bs); rm(meth.levels.raw); gc();

# Sigmoid Colon
load(paste0(data.file.prefix, "roadmap_Sigmoid Colon_all_bsseq_id.RData")); 
meth.levels.raw = getMeth(bs, type = "raw")
no.hits = which(is.na(rowMeans(meth.levels.raw)) == TRUE)
bs = bs[-no.hits]
SC <- bs 
rm(bs); rm(meth.levels.raw); gc();

# Left Ventricle
load(paste0(data.file.prefix, "roadmap_Small Intestine_all_bsseq_id.RData")); 
bs <- bs[,-which(pData(bs)$Age==30)]
meth.levels.raw = getMeth(bs, type = "raw")
no.hits = which(is.na(rowMeans(meth.levels.raw)) == TRUE)
bs = bs[-no.hits]
SI <- bs 
rm(bs); rm(meth.levels.raw); gc();

# summary table

CpGs <- c(length(LV), length(RV), length(SC), length(SI))
LV.stats <- getCoverageStats(LV)
RV.stats <- getCoverageStats(RV)
SC.stats <- getCoverageStats(SC)
SI.stats <- getCoverageStats(SI)
med.cov <- (rbind(unlist(LV.stats[[1]]), unlist(RV.stats[[1]]), 
					unlist(SC.stats[[1]]), unlist(SI.stats[[1]])))
max.cov <- (rbind(unlist(LV.stats[[2]]), unlist(RV.stats[[2]]), 
					unlist(SC.stats[[2]]), unlist(SI.stats[[2]])))

summary.table[[1]] <- cbind(CpGs, med.cov, max.cov)
rownames(summary.table[[1]]) <- c("LV", "RV", "SC", "SI")
colnames(summary.table[[1]]) <- c("CpGs", "Med1", "Med2", "Max1", "Max2")

bs <- bsseq::combine(LV, RV, SC, SI)
rm(LV)
rm(RV)
rm(SC)
rm(SI)

# add "chr" prefix to chromosome names
chr <- paste0("chr", as.character(seqnames(bs)), sep="")
bs <- BSseq(chr = chr, pos = start(bs),
               M = getCoverage(bs, type="M"), 
               Cov = getCoverage(bs, type="Cov"), 
               sampleNames = sampleNames(bs),
               pData=pData(bs))

roadmap <- plotBetaDistribution(bs, "Human Tissues")
rm(bs)

### dendritic data
data.file.prefix <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/DATA/"

# 3 vs 3 controls
load(paste0(data.file.prefix, "DC_BSSeq.RData"))
DC <- chrSelectBSseq(DC, seqnames = paste0("chr", c(seq(1,22))), order = TRUE)
DC <- DC[,pData(DC)$Infected == "NI"]

meth.levels.raw = getMeth(DC, type = "raw")
no.hits = which(is.na(rowMeans(meth.levels.raw)) == TRUE)
DC = DC[-no.hits]

CpGs <- length(DC)
DC.stats <- getCoverageStats(DC)

med.cov <- unlist(DC.stats[[1]])
max.cov <- unlist(DC.stats[[2]])

newtab <- cbind(CpGs, med.cov, max.cov)
summary.table[[2]] <- newtab

dendritic.control3 <- plotBetaDistribution(DC, "Dendritic Controls")
meta3 <- pData(DC)

# 3vs3 Simulated
load("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/binSim_preferential/sim.data.n3.control.all.rda")
DC <- sim.dat.red$bs
rm(sim.dat.red)
pData(DC) <- meta3
DC <- DC[,c(1,4,2,5,3,6)]
pData(DC) <- meta3

CpGs <- length(DC)
DC.stats <- getCoverageStats(DC)

med.cov <- unlist(DC.stats[[1]])
max.cov <- unlist(DC.stats[[2]])

newtab <- cbind(CpGs, med.cov, max.cov)
rownames(newtab) <- paste0("Dendritic_Spikein_3", rep(1:nrow(newtab)))
summary.table[[4]] <- newtab

dendritic.sim3 <- plotBetaDistribution(DC, "Dendritic Simulated")
rm(DC)

## DNMT3a
load("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/DATA/dnmt3a_all_bsseq.RData")
pData(bs)$Condition <- pData(bs)$Cond
pData(bs)$Condition[pData(bs)$Condition=="KO_FLT3"] <- "ALL"
pData(bs)$Condition[pData(bs)$Condition=="WT_FLT3"] <- "AML"
pData(bs)$Condition[pData(bs)$Condition=="WT_WTFL"] <- "Control"

ALL <- bs[,pData(bs)$Condition == "ALL"]
AML <- bs[,pData(bs)$Condition == "AML"]
Control <- bs[,pData(bs)$Condition == "Control"]

rm(bs)
meth.levels.raw = getMeth(ALL, type = "raw")
no.hits = which(is.na(rowMeans(meth.levels.raw)) == TRUE)
ALL = ALL[-no.hits]

meth.levels.raw = getMeth(AML, type = "raw")
no.hits = which(is.na(rowMeans(meth.levels.raw)) == TRUE)
AML = AML[-no.hits]

meth.levels.raw = getMeth(Control, type = "raw")
no.hits = which(is.na(rowMeans(meth.levels.raw)) == TRUE)
Control = Control[-no.hits]

CpGs <- c(length(ALL), length(AML), length(Control))
ALL.stats <- getCoverageStats(ALL)
AML.stats <- getCoverageStats(AML)
Control.stats <- getCoverageStats(Control)
med.cov <- (rbind(unlist(ALL.stats[[1]]), unlist(AML.stats[[1]]), 
					unlist(Control.stats[[1]])))
max.cov <- (rbind(unlist(ALL.stats[[2]]), unlist(AML.stats[[2]]), 
					unlist(Control.stats[[2]])))

summary.table[[6]] <- cbind(CpGs, med.cov, max.cov)
rownames(summary.table[[6]]) <- c("ALL", "AML", "Control")
colnames(summary.table[[6]]) <- c("CpGs", "Med1", "Med2", "Max1", "Max2")

rm(ALL)
rm(AML)
rm(Control)

load("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/DATA/dnmt3a_all_bsseq.RData")
pData(bs)$Condition <- pData(bs)$Cond
pData(bs)$Condition[pData(bs)$Condition=="KO_FLT3"] <- "ALL"
pData(bs)$Condition[pData(bs)$Condition=="WT_FLT3"] <- "AML"
pData(bs)$Condition[pData(bs)$Condition=="WT_WTFL"] <- "Control"

meth.levels.raw = getMeth(bs, type = "raw")
no.hits = which(is.na(rowMeans(meth.levels.raw)) == TRUE)
bs = bs[-no.hits]

mouse <- plotBetaDistribution(bs, "Murine Leukemia")
rm(bs)

print(summary.table)


# initialize plot file
plotFile <- paste0(result.file.prefix, 
						"/supp_fig2.pdf")
 pdf(file=plotFile, height=6, width=10)	
  
  print(plot_grid(dendritic.control3[[1]], 
  		   dendritic.sim3[[1]], 
  		   roadmap[[1]],
 		   mouse[[1]], 
 		   labels = c("(A)", "(B)", "(C)", "(D)"),
 		   align="hv"))
 		   
 dev.off()
 