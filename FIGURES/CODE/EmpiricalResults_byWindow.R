# Code to generate Figure 3, Supplementary Figures S6-S9

# source("/n/irizarryfs01_backed_up/kkorthauer/WGBS/PAPER/FIGURES/CODE/EmpiricalResults_byWindow.R")

# plots associations of expression with methylation broken up by type of window used to
# overlap the DMRs with: promoter, islandshores (extended promoter regions), and 
# gene body (cds plus 3 prime UTR plus 5 prime UTR)

######################################################
### parameters to change to run on your own system ###
######################################################
# change the following the dmrseq results directory - assumes dmrseq results 
# are in a folder called 'dmrseq_pkg', and the directories for other methods 
# are called 'BSmooth', 'DSS', and 'metilene'. 
rm.result.dir <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/RESULTS/"
# change the following to the dnmt3a results directory 
dn.result.dir <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/RESULTS/"
# change the following to where you'd like to save the figure output
result.file.prefix <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/PAPER/FIGURES/out/"
######################################################
###         end of parameters to change            ###
######################################################

##### Roadmap Results Plots
library(ggplot2)
library(reshape2)
library(cowplot)
library(ggthemes)

for(win in c("promoter", "genebody", "islandshore")){

subfolder <- "odds_effectSize"
sampleSize <- 2
num.dmrs <- 0
min.length <- max.length <- NULL
allConditions <- c("Right Ventricle_Small Intestine",
				   "Right Ventricle_Sigmoid Colon",
				   "Left Ventricle_Small Intestine",	
				   "Left Ventricle_Sigmoid Colon",	   
				   "Sigmoid Colon_Small Intestine")
ct <- 1
odds <- covg <- near <- fdr <- tit <- vector("list", length(allConditions))
grob.dmrseq <- grob.bsmooth <- grob.dss <- grob.met <- vector("list", length(allConditions))
labs <- c("A", "B", "C", "D", "E")

for(cond in allConditions){

   # for plotting when calling dmrseq
   pfile1 <- list.files(path=paste0(rm.result.dir, "/dmrseq_pkg/", subfolder, "/"),
				   pattern=paste0("oddsByRanking.dmrseq.n", 
					   sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
   pfile2 <- list.files(path=paste0(rm.result.dir, "/dmrseq_pkg/", subfolder, "/"),
				   pattern=paste0("_min", min.length, "_max", max.length, ".", win, ".RData")) 
   if (length(pfile1[pfile1 %in% pfile2]) == 1){	   					 					
	   pfile <- paste0(rm.result.dir, "/dmrseq_pkg/", subfolder, "/", 
		 pfile1[pfile1 %in% pfile2])
   }else if (length(pfile1[pfile1 %in% pfile2]) > 1){	  #if > 1, pick most recently modified
	   details <- file.info(pfile1[pfile1 %in% pfile2])
	   details <- details[with(details, order(as.POSIXct(mtime), decreasing=TRUE)), ]
	   pfile <- paste0(rm.result.dir, "/dmrseq_pkg/", subfolder, "/", rownames(details)[1])
	   
	   thisOne <- grep(nrow(dmrs), pfile1[pfile1 %in% pfile2])
	   pfile <- paste0(rm.result.dir, "/dmrseq_pkg/", subfolder, "/", pfile1[pfile1 %in% pfile2][thisOne])
   }else{
	  pfile <- "null"
   }
   load(pfile)
   rank.df$Method <- "dmrseq"		
   ranks.last <- rank.df[nrow(rank.df),,drop=FALSE]
   ranks.last$Method <- "dmrseq"  

 
   bfile1 <- list.files(path=paste0(rm.result.dir, "/BSmooth/", subfolder, "/"),
   						pattern=paste0("oddsByRanking.BSmooth.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
   bfile2 <- list.files(path=paste0(rm.result.dir, "/BSmooth/", subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", win, ".RData")) 		  
   if (length(bfile1[bfile1 %in% bfile2]) == 1){	   					 					
   bfile <- paste0(rm.result.dir, "/BSmooth/", subfolder, "/", 
   		      bfile1[bfile1 %in% bfile2])
   }else{
   	  bfile <- "null"
   }

   dfile1 <- list.files(path=paste0(rm.result.dir,"/DSS/", subfolder, "/"),
   						pattern=paste0("oddsByRanking.DSS.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
   dfile2 <- list.files(path=paste0(rm.result.dir, "/DSS/", subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", win,  ".RData")) 		  
   if (length(dfile1[dfile1 %in% dfile2]) == 1){	   					 					
   dfile <- paste0(rm.result.dir, "/DSS/", subfolder, "/", 
   		      dfile1[dfile1 %in% dfile2])
   }else{
   	  dfile <- "null"
   }
  
   mfile1 <- list.files(path=paste0(rm.result.dir, "/metilene/", subfolder, "/"),
   						pattern=paste0("oddsByRanking.metilene.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
   mfile2 <- list.files(path=paste0(rm.result.dir, "/metilene/", subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", win, ".RData")) 		  
   if (length(mfile1[mfile1 %in% mfile2]) == 1){	   					 					
   mfile <- paste0(rm.result.dir, "/metilene/", subfolder, "/", 
   		      mfile1[mfile1 %in% mfile2])
   }else{
   	  mfile <- "null"
   }
  
   if(file.exists(pfile)){
	 load(pfile)
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
	 if (!exists("ranks.all")){
	   ranks.all <- rank.df
	 }else{
	   ranks.all <- rbind(ranks.all, rank.df)
	   ranks.last <- rbind(ranks.last, last)
	 }
   }
   
  if(file.exists(mfile)){
	 load(mfile)
	 coldiff <- ncol(ranks.all)-ncol(rank.df)
	 if (coldiff>1){
	   for (m in 1:(coldiff-1)){
		   rank.df <- cbind(rank.df, rep(NA, nrow(rank.df)))
	   }
	 }
	 rank.df$Method <- "metilene"
	 last <- rank.df[nrow(rank.df),,drop=FALSE]
	 last$Method <- "metilene"
	 if (!exists("ranks.all")){
	   ranks.all <- rank.df
	 }else{
	   ranks.all <- rbind(ranks.all, rank.df)
	   ranks.last <- rbind(ranks.last, last)
	 }
   }
  
   ranks.all <- ranks.all[order(ranks.all$Method),]
   
   odds[[ct]] <- ggplot(aes(x=ranks, y=odds), 
						 data=ranks.all) + 
	   geom_ribbon(aes(ymin=CI.low, ymax=CI.hi, fill=Method),
			alpha=0.3) +
	   geom_hline(yintercept=0, colour="grey", alpha=0.6, lwd=1.5, lty=2) +
	   geom_line(aes(group=Method, color=Method), lwd=1.1) +
	   xlab("Number of top-ranked DMRs") +
	   ylab("Odds") +
	   ggtitle("Odds of inverse DE association with Methylation") + 
	   theme_classic()	+
     scale_colour_colorblind()+
     scale_fill_colorblind() +
     scale_x_continuous(trans="log2") +
     scale_y_continuous(trans="log2")
  
   covg[[ct]] <- ggplot(aes(x=ranks, y=nDML), 
						 data=ranks.all) + 
	   geom_line(aes(group=Method, color=Method), lwd=1.1) +
	   geom_point(aes(color=Method), size=2) +
	   xlab("Number of top-ranked DMRs") +
	   ylab("Number of CpG loci") +
	   ggtitle("Number CpG Loci in top DMRs") + 
	   geom_ribbon(aes(ymin=nDML, ymax=nDML, fill=Method),
				   alpha=0.3)+
	   theme_classic()+
     scale_colour_colorblind()+
     scale_fill_colorblind() +
     scale_x_continuous(trans="log2") +
     scale_y_continuous(trans="log2")
  
	  
   DE.cats <- melt(ranks.all[,c(1,11:14,17)], id=c("ranks", "Method"))
   colnames(DE.cats)[3] <- "CpG_Category"
   nchars <- sapply(as.character(DE.cats$CpG_Category), nchar)
   DE.cats$CpG_Category <- substr(DE.cats$CpG_Category, 1, nchars-2)
  
   # add 95% confidence bounds
   DE.cats$CI.low <- DE.cats$value - 
	   1.96 * sqrt(DE.cats$value*(1-DE.cats$value) / DE.cats$ranks)
   DE.cats$CI.hi <-  DE.cats$value + 
	   1.96 * sqrt(DE.cats$value*(1-DE.cats$value) / DE.cats$ranks)
	  
   near[[ct]] <-  ggplot(aes(x=ranks, y=value), 
	   data=DE.cats) + 
	   geom_ribbon(aes(ymin=(CI.low), ymax=(CI.hi), fill=Method),
			alpha=0.3) +
	   geom_line(aes(color=Method,
					group=interaction(CpG_Category,Method)), 
					lwd=1.1) +
	   facet_wrap( ~ CpG_Category, ncol=2) +
	   ylab("Proportion DMRs near a DE gene") +
	   xlab("Number of Top-ranked DMRs") +
	   ggtitle("Proportion of DMRs near a DE gene") + 
	   theme_classic() +
     scale_colour_colorblind()+
     scale_fill_colorblind() +
     scale_x_continuous(trans="log2")
	
  labs[ct] <- paste0("(", labs[ct], ") ", gsub("_", " vs ", cond))
 	   
  ############################################################################ 
  # correlation (supplementary results)
  # use up to top 10000 DMRs (or max number)

  load(paste0(rm.result.dir, "/dmrseq_pkg/", subfolder, "/corrWithExpressionPlot.", 
		 "dmrseq.n", sampleSize, ".", cond, ".", num.dmrs, "DMRS.top", 
		  min(10000, max(ranks.all$ranks[ranks.all$Method=="dmrseq"])), ".", win, ".RData"))
  grob.dmrseq[[ct]] <- gp + xlab("Methylation difference (effect size)")
  
  
  load(paste0(rm.result.dir, "/BSmooth/",
  		  subfolder, "/corrWithExpressionPlot.", 
		 "BSmooth.n", sampleSize, ".", cond, ".", num.dmrs, "DMRS.top", 
		  min(10000, max(ranks.all$ranks[ranks.all$Method=="BSmooth"])), ".", win, ".RData"))
  grob.bsmooth[[ct]] <- gp + xlab("Methylation difference (effect size)")
		  
		
  load(paste0(rm.result.dir, "/DSS/", 
  		  subfolder, "/corrWithExpressionPlot.", 
		 "DSS.n", sampleSize, ".", cond, ".", num.dmrs, "DMRS.top", 
		  min(10000, max(ranks.all$ranks[ranks.all$Method=="DSS"])), ".", win, ".RData"))
  grob.dss[[ct]] <- gp + xlab("Methylation difference (effect size)")
  
  load(paste0(rm.result.dir, "/metilene/", 
  		  subfolder, "/corrWithExpressionPlot.", 
		 "metilene.n", sampleSize, ".", cond, ".", num.dmrs, "DMRS.top", 
		  min(10000, max(ranks.all$ranks[ranks.all$Method=="metilene"])), ".", win, ".RData"))
  grob.met[[ct]] <- gp + xlab("Methylation difference (effect size)")
  
    ############################################################################
  # fdr rankings / default bsmooth/dss
  subfolder <-  "odds_fdr_Cumulative"
  
  pfile1 <- list.files(path=paste0(rm.result.dir, "/dmrseq_pkg/", subfolder, "/"),
   						pattern=paste0("oddsByRanking.dmrseq.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
  pfile2 <- list.files(path=paste0(rm.result.dir, "/dmrseq_pkg/", subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", win, ".RData")) 
  if (length(pfile1[pfile1 %in% pfile2]) == 1){	   					 					
	  pfile <- paste0(rm.result.dir, "/dmrseq_pkg/", subfolder, "/", 
		pfile1[pfile1 %in% pfile2])
  }else if (length(pfile1[pfile1 %in% pfile2]) > 1){	  #if > 1, pick most recently modified
	  details <- file.info(pfile1[pfile1 %in% pfile2])
	  details <- details[with(details, order(as.POSIXct(mtime), decreasing=TRUE)), ]
	  pfile <- paste0(rm.result.dir, "/dmrseq_pkg/", subfolder, "/", rownames(details)[1])
	  
	  thisOne <- grep(nrow(dmrs), pfile1[pfile1 %in% pfile2])
	  pfile <- paste0(rm.result.dir, "/dmrseq_pkg/", subfolder, "/", pfile1[pfile1 %in% pfile2][thisOne])
  }else{
	 pfile <- "null"
  }
  load(pfile)
  rank.df$Method <- "dmrseq"		
  ranks.last <- rank.df[nrow(rank.df),,drop=FALSE]
  ranks.last$Method <- "dmrseq"  
  fdrRanks <- rank.df$ranks
  ranks.all <- rank.df
  
  mfile1 <- list.files(path=paste0(rm.result.dir, "/metilene_default/",
  						subfolder, "/"),
   						pattern=paste0("oddsByRanking.metilene.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
  mfile2 <- list.files(path=paste0(rm.result.dir, "/metilene_default/",
  				 		subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", win, ".RData")) 
  if (length(mfile1[mfile1 %in% mfile2]) == 1){	   					 					
	  mfile <- paste0(rm.result.dir, "/metilene_default/",
	    subfolder, "/", 
		mfile1[mfile1 %in% mfile2])
  }else if (length(mfile1[mfile1 %in% mfile2]) > 1){	  #if > 1, pick most recently modified
	  details <- file.info(mfile1[mfile1 %in% mfile2])
	  details <- details[with(details, order(as.POSIXct(mtime), decreasing=TRUE)), ]
	  mfile <- paste0(rm.result.dir, "/metilene_default/",
	  			 		subfolder, "/", rownames(details)[1])
	  
	  thisOne <- grep(nrow(dmrs), mfile1[mfile1 %in% mfile2])
	  mfile <- paste0(rm.result.dir, "/metilene_default/",
	   					subfolder, "/", mfile1[mfile1 %in% mfile2][thisOne])
  }else{
	 mfile <- "null"
  }
  load(mfile)
  rank.df$Method <- "metilene"		
  last <- rank.df[nrow(rank.df),,drop=FALSE]
  ranks.all <- rbind(ranks.all, rank.df)
  ranks.last <- rbind(ranks.last, last)
  fdrRanks <- ranks.all$ranks
  
  subfolder <-  "odds_statistic_default"
   
   bfile1 <- list.files(path=paste0(rm.result.dir, "/BSmooth_default/", subfolder, "/"),
   						pattern=paste0("oddsByRanking.BSmooth_default.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
   bfile2 <- list.files(path=paste0(rm.result.dir, "/BSmooth_default/", subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", win, ".RData")) 		  
   if (length(bfile1[bfile1 %in% bfile2]) == 1){	   					 					
   bfile <- paste0(rm.result.dir, "/BSmooth_default/", subfolder, "/", 
   		      bfile1[bfile1 %in% bfile2])
   }else{
   	  bfile <- "null"
   }

   dfile1 <- list.files(path=paste0(rm.result.dir, "/DSS_default/", subfolder, "/"),
   						pattern=paste0("oddsByRanking.DSS_default.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
   dfile2 <- list.files(path=paste0(rm.result.dir, "/DSS_default/", subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", win, ".RData")) 		  
   if (length(dfile1[dfile1 %in% dfile2]) == 1){	   					 					
   dfile <- paste0(rm.result.dir, "/DSS_default/", subfolder, "/", 
   		      dfile1[dfile1 %in% dfile2])
   }else{
   	  dfile <- "null"
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
	 rank.df <- last
  
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
	 rank.df <- last
	 
	 if (!exists("ranks.all")){
	   ranks.all <- rank.df
	 }else{
	   ranks.all <- rbind(ranks.all, rank.df)
	   ranks.last <- rbind(ranks.last, last)
	 }
   }
  
   ranks.all <- ranks.all[order(ranks.all$Method),]
   
   ranks.all$Method <- as.factor(ranks.all$Method)
   ranks.all$ranks[ranks.all$ranks > 1] <- min(fdrRanks)
   if (length(unique(ranks.all$ranks[ranks.all$Method=="metilene"])) > 1){
   		grad <- c("dmrseq", "metilene")
   }else{
   		grad <- c("dmrseq")
   		ranks.all <- ranks.all[-(which(ranks.all$Method=="metilene")[1]),]
   		ranks.all$ranks[ranks.all$Method=="metilene"] <- min(fdrRanks)
   }
   add2 <- ranks.all[!(ranks.all$Method %in% grad),]
   if(nrow(add2)>0){
	   add2$ranks <- 1
	   ranks.all <- rbind(ranks.all, add2)
   }
   
   ranks.all$fdr <- 1
   ranks.all$fdr[ranks.all$Method %in% grad] <- 0
   ranks.all$fdr <- as.factor(ranks.all$fdr)
   minx <- min(ranks.all$ranks)
      
   fdr[[ct]] <- ggplot(aes(x=ranks, y=odds), 
					 data=ranks.all) + 
	   geom_ribbon(aes(ymin=CI.low, ymax=CI.hi, fill=Method),
			alpha=0.3) +
	   geom_line(aes(group=Method, color=Method, linetype=fdr), lwd=1.1) +
	   scale_y_continuous(trans="log2", limits=c(2^-0.065, 2^3.2),
	                      breaks=c(1,2,4,8,16)) +
	   xlab("FDR threshold (square root scaled)") +
	   ylab("Odds") +
	   ggtitle(substr(labs[ct], 5, nchar(labs[ct]))) + 
	   theme_classic() +
	   scale_x_sqrt(breaks=c(0,0.01,0.1,0.25,0.5,0.75,1)) +
	   guides(linetype=FALSE)+
     scale_colour_colorblind() +
     scale_fill_colorblind() 
  
  ct <- ct + 1
  subfolder <- "odds_effectSize"
}   			

 p1 <- plot_grid( odds[[1]] + theme(legend.position="none", axis.title.x=element_blank()), 
 covg[[1]] + theme(legend.position="none", axis.title.x=element_blank()), 
 odds[[2]] + theme(legend.position="none", plot.title =element_blank(), 
 					axis.title.x=element_blank()), 
 covg[[2]] + theme(legend.position="none", plot.title =element_blank(), 
 					axis.title.x=element_blank()), 
 odds[[3]] + theme(legend.position="none", plot.title =element_blank(), 
 					axis.title.x=element_blank()), 
 covg[[3]] + theme(legend.position="none", plot.title =element_blank(),
 					axis.title.x=element_blank()), 
 odds[[4]] + theme(legend.position="none", plot.title =element_blank(), 
 					axis.title.x=element_blank()), 
 covg[[4]] + theme(legend.position="none", plot.title =element_blank(), 
 					axis.title.x=element_blank()), 
 odds[[5]] + theme(legend.position="none", plot.title =element_blank()), 
 covg[[5]] + theme(legend.position="none", plot.title =element_blank()), 
                    labels = c("", "", labs[2], "", labs[3], "", labs[4], "", labs[5], ""),
                    hjust = 0,
                    nrow = 5,
                    align = "hv",
                    rel_widths = c(4,2.5))

                   
p2 <- plot_grid(near[[1]] + theme(legend.position="none", axis.title.x=element_blank()),
					near[[2]] + theme(legend.position="none", plot.title =element_blank(), axis.title.x=element_blank()),
					near[[3]] + theme(legend.position="none", plot.title =element_blank(), axis.title.x=element_blank()),
					near[[4]] + theme(legend.position="none", plot.title =element_blank(), axis.title.x=element_blank()),
					near[[5]] + theme(legend.position="none", plot.title =element_blank()),
                   hjust = -1,
                   nrow = 5)
legend <- get_legend(covg[[1]] + 
                theme(legend.position="bottom", 
                	  legend.title=element_text(size=14), 
                	  legend.text=element_text(size=13)))
top.title <- ggdraw() + draw_label(paste0(labs[1]), x=0.003, hjust=0, 
								fontface="bold")

pdf(paste0(result.file.prefix, "/supp_fig6.", win, ".pdf"),
  height=16, width=16)
	p <- plot_grid( p1, p2, rel_widths = c(2,1.3), nrow=1)
	p <- plot_grid( top.title, p, legend, nrow=3, rel_heights = c(0.02,1,0.025), hjust=0)
	print(p)
dev.off()



p1 <- plot_grid(grob.dmrseq[[1]] + theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()), 
				grob.bsmooth[[1]] + theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()), 
				grob.dss[[1]]+ theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()),
				grob.met[[1]]+ theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()),
				grob.dmrseq[[2]]+ theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()), 
				grob.bsmooth[[2]]+ theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()), 
				grob.dss[[2]]+ theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()),
				grob.met[[2]]+ theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()),
				grob.dmrseq[[3]]+ theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()), 
				grob.bsmooth[[3]]+ theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()), 
				grob.dss[[3]]+ theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()),
				grob.met[[3]]+ theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()),
				grob.dmrseq[[4]]+ theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()), 
				grob.bsmooth[[4]]+ theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()), 
				grob.dss[[4]]+ theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()),
				grob.dss[[4]]+ theme(legend.position="none", plot.title =element_text(size=8.25), axis.title.x=element_blank()),
				grob.dmrseq[[5]]+ theme(legend.position="none", plot.title =element_text(size=8.25)), 
				grob.bsmooth[[5]]+ theme(legend.position="none", plot.title =element_text(size=8.25)), 
				grob.dss[[5]]+ theme(legend.position="none", plot.title =element_text(size=8.25)),
				grob.met[[5]]+ theme(legend.position="none", plot.title =element_text(size=8.25)),
				nrow=5, hjust=-0.05, vjust=-1.25,
				labels = c(NA, "", "", "", labs[2], "", "", "", labs[3], "", "",  "",	
						   labs[4], "", "", "", labs[5], "", "", ""),
				align="hv",
				scale=0.95)
				
legend <- suppressMessages(get_legend(grob.dmrseq[[1]] + 
 				scale_colour_manual(name="",
                            values=c("red", "black"),
                            labels=c("DE", "Not DE")) +
                scale_fill_discrete(name="",
                         breaks=c(TRUE, FALSE),
                         labels=c("DE", "Not DE"))+ 
                theme(legend.position="bottom", 
                	  legend.title=element_text(size=14), 
                	  legend.text=element_text(size=13))))
top.title <- ggdraw() + draw_label(paste0(labs[1]), 
							x=0.003, hjust=-0.027, vjust=0.2, fontface="bold")
				
pdf(paste0(result.file.prefix, "/supp_fig7.", win, ".pdf"), 
  height=18.333, width=17)
	p <- plot_grid( top.title, p1, nrow=2, rel_heights = c(0.036,1))
	p <- plot_grid( p, legend, rel_heights = c(3,0.06), nrow=2)
	p <- p  + draw_label("dmrseq", x = 0.16, y = 0.973,
            	vjust = 1, hjust = 1, size = 14, fontface = 'bold') +
            draw_label("BSmooth", x = 0.412, y = 0.973,
            	vjust = 1, hjust = 1, size = 14, fontface = 'bold') +
            draw_label("DSS", x = 0.65, y = 0.973,
            	vjust = 1, hjust = 1, size = 14, fontface = 'bold') +
            draw_label("metilene", x = 0.908, y = 0.973,
            	vjust = 1, hjust = 1, size = 14, fontface = 'bold')
	print(p)
dev.off()   


p1 <- plot_grid( fdr[[1]] + theme(legend.position="none", 
					plot.title = element_text(size = 12)), 
				 fdr[[2]]+ theme(legend.position="none", 
 					axis.title.y=element_blank(), 
 					plot.title = element_text(size = 12)), 
				 fdr[[3]] + theme(legend.position="none", 
 					axis.title.y=element_blank(), 
 					plot.title = element_text(size = 12)),
				 fdr[[4]] + theme(legend.position="none",
				  	plot.title = element_text(size = 12)), 
				 fdr[[5]] + theme(legend.position="none",
 					axis.title.y=element_blank(),
 					plot.title = element_text(size = 12)),
				 nrow=2 )
	
legend <- get_legend(fdr[[1]] + 
                theme(legend.position="bottom", 
                	  legend.title=element_text(size=14), 
                	  legend.text=element_text(size=13)))
top.title <- ggdraw() + draw_label("FDR versus odds of inverse DE association with Methylation", 
							x=0.003, hjust=-0.29, fontface="bold", vjust=0.8) 
top.title2 <- ggdraw() + draw_label("(A) Roadmap Tissue Comparisons", 
							x=0.003, hjust=0, fontface="bold", vjust=0.8)                   	  
p <- plot_grid(top.title, p1, legend, rel_heights = c(0.2,3,0.06), nrow=3) 


p1.fdr.rm <- plot_grid(top.title2, p1, rel_heights = c(0.2,3), nrow=2) 

}

############################################################################ 
############################################################################ 
############################################################################    
############################################################################ 
######################    Do the same for mouse  ###########################
############################################################################ 
############################################################################ 
############################################################################ 
############################################################################ 

for(win in c("promoter", "genebody", "islandshore")){

sampleSize <- 2
num.dmrs <- 0
min.length <- max.length <- NULL
allConditions <- c("KO_FLT3WT_WTFL",
				   "WT_FLT3WT_WTFL",
				   "KO_FLT3WT_FLT3")
				   

ct <- 1
odds <- covg <- near <- fdr <- tit <- vector("list", length(allConditions))
grob.dmrseq <- grob.bsmooth <- grob.dss <- grob.met <- vector("list", length(allConditions))
labs <- c("(A) AML vs Control", "(B) ALL vs Control", "(C) AML vs ALL")

for(cond in allConditions){
   subfolder <- "odds_effectSize"
   # for plotting when calling dmrseq
   		pfile1 <- list.files(path=paste0(dn.result.dir, "/dmrseq_pkg/", subfolder, "/"),
   						pattern=paste0("oddsByRanking.dmrseq.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
  	 	pfile2 <- list.files(path=paste0(dn.result.dir, "/dmrseq_pkg/", subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", win, ".RData")) 
   		if (length(pfile1[pfile1 %in% pfile2]) == 1){	   					 					
   			pfile <- paste0(dn.result.dir, "/dmrseq_pkg/", subfolder, "/", 
   		      pfile1[pfile1 %in% pfile2])
		}else if (length(pfile1[pfile1 %in% pfile2]) > 1){	  #if > 1, pick most recently modified
			details <- file.info(pfile1[pfile1 %in% pfile2])
			details <- details[with(details, order(as.POSIXct(mtime), decreasing=TRUE)), ]
			pfile <- paste0(dn.result.dir, "/dmrseq_pkg/", subfolder, "/", rownames(details)[1])
			
			thisOne <- grep(nrow(dmrs), pfile1[pfile1 %in% pfile2])
			pfile <- paste0(dn.result.dir, "/dmrseq_pkg/", subfolder, "/", pfile1[pfile1 %in% pfile2][thisOne])
		}else{
		   pfile <- "null"
		}
		load(pfile)
		rank.df$Method <- "dmrseq"		
		ranks.last <- rank.df[nrow(rank.df),,drop=FALSE]
		ranks.last$Method <- "dmrseq"  

 
   bfile1 <- list.files(path=paste0(dn.result.dir, "/BSmooth/", subfolder, "/"),
   						pattern=paste0("oddsByRanking.BSmooth.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
   bfile2 <- list.files(path=paste0(dn.result.dir, "/BSmooth/", subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", win, ".RData")) 		  
   if (length(bfile1[bfile1 %in% bfile2]) == 1){	   					 					
   bfile <- paste0(dn.result.dir, "/BSmooth/", subfolder, "/", 
   		      bfile1[bfile1 %in% bfile2])
   }else{
   	  bfile <- "null"
   }

   dfile1 <- list.files(path=paste0(dn.result.dir, "/DSS/", subfolder, "/"),
   						pattern=paste0("oddsByRanking.DSS.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
   dfile2 <- list.files(path=paste0(dn.result.dir, "/DSS/", subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", win, ".RData")) 		  
   if (length(dfile1[dfile1 %in% dfile2]) == 1){	   					 					
   dfile <- paste0(dn.result.dir, "/DSS/", subfolder, "/", 
   		      dfile1[dfile1 %in% dfile2])
   }else{
   	  dfile <- "null"
   }
   
   mfile1 <- list.files(path=paste0(dn.result.dir, "/metilene/", subfolder, "/"),
   						pattern=paste0("oddsByRanking.metilene.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
   mfile2 <- list.files(path=paste0(dn.result.dir, "/metilene/", subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", win, ".RData")) 		  
   if (length(mfile1[mfile1 %in% mfile2]) == 1){	   					 					
   mfile <- paste0(dn.result.dir, "/metilene/", subfolder, "/", 
   		      mfile1[mfile1 %in% mfile2])
   }else{
   	  mfile <- "null"
   }
  
  
   if(file.exists(pfile)){
	 load(pfile)
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
	 if (!exists("ranks.all")){
	   ranks.all <- rank.df
	 }else{
	   ranks.all <- rbind(ranks.all, rank.df)
	   ranks.last <- rbind(ranks.last, last)
	 }
   }
   
   if(file.exists(mfile)){
	 load(mfile)
	 coldiff <- ncol(ranks.all)-ncol(rank.df)
	 if (coldiff>1){
	   for (m in 1:(coldiff-1)){
		   rank.df <- cbind(rank.df, rep(NA, nrow(rank.df)))
	   }
	 }
	 rank.df$Method <- "metilene"
	 last <- rank.df[nrow(rank.df),,drop=FALSE]
	 last$Method <- "metilene"
	 if (!exists("ranks.all")){
	   ranks.all <- rank.df
	 }else{
	   ranks.all <- rbind(ranks.all, rank.df)
	   ranks.last <- rbind(ranks.last, last)
	 }
   }
  
  
   ranks.all <- ranks.all[order(ranks.all$Method),]
   
   odds[[ct]] <- ggplot(aes(x=ranks, y=odds), 
						 data=ranks.all) + 
	   geom_ribbon(aes(ymin=(CI.low), ymax=(CI.hi), fill=Method),
			alpha=0.3) +
	   geom_hline(yintercept=0, colour="grey", alpha=0.6, lwd=1.5, lty=2) +
	   geom_line(aes(group=Method, color=Method), lwd=1.1) +
	   scale_y_continuous(trans="log2") +
	   scale_x_continuous(trans="log2") +
	   xlab("Number of top-ranked DMRs") +
	   ylab("Odds") +
	   ggtitle("Odds of inverse DE association with Methylation") + 
	   theme_classic()		+
     scale_colour_colorblind()+
     scale_fill_colorblind()
  
   covg[[ct]] <- ggplot(aes(x=ranks, y=nDML), 
						 data=ranks.all) + 
	   geom_line(aes(group=Method, color=Method), lwd=1.1) +
	   geom_point(aes(color=Method), size=2) +
	   scale_y_continuous(trans="log2") +
	   scale_x_continuous(trans="log2") +
	   xlab("Number of top-ranked DMRs") +
	   ylab("Number of CpG loci") +
	   ggtitle("Number CpG Loci in top DMRs") + 
	   geom_ribbon(aes(ymin=nDML, ymax=nDML, fill=Method),
				   alpha=0.3)+
	   theme_classic()+
     scale_colour_colorblind()		+
     scale_fill_colorblind()
  
	  
   DE.cats <- melt(ranks.all[,c(1,11:14,17)], id=c("ranks", "Method"))
   colnames(DE.cats)[3] <- "CpG_Category"
   nchars <- sapply(as.character(DE.cats$CpG_Category), nchar)
   DE.cats$CpG_Category <- substr(DE.cats$CpG_Category, 1, nchars-2)
  
   # add 95% confidence bounds
   DE.cats$CI.low <- DE.cats$value - 
	   1.96 * sqrt(DE.cats$value*(1-DE.cats$value) / DE.cats$ranks)
   DE.cats$CI.hi <-  DE.cats$value + 
	   1.96 * sqrt(DE.cats$value*(1-DE.cats$value) / DE.cats$ranks)
	  
   near[[ct]] <-  ggplot(aes(x=ranks, y=value), 
	   data=DE.cats) + 
	   geom_ribbon(aes(ymin=CI.low, ymax=CI.hi, fill=Method),
			alpha=0.3) +
	   geom_line(aes(color=Method,
					group=interaction(CpG_Category,Method)), 
					lwd=1.1) +
	   facet_wrap( ~ CpG_Category, ncol=2) +
	   ylab("Proportion DMRs near a DE gene") +
	   xlab("Number of Top-ranked DMRs") +
	   ggtitle("Proportion of DMRs near a DE gene") + 
	   theme_classic()+
     scale_colour_colorblind()+
     scale_fill_colorblind() +
	 scale_x_continuous(trans="log2") 
	   
  ############################################################################ 
  # correlation (supplementary results)
  # use up to top 10000 DMRs (or max number)

  load(paste0(dn.result.dir, "/dmrseq_pkg/", subfolder, "/corrWithExpressionPlot.", 
		 "dmrseq.n", sampleSize, ".", cond, ".", num.dmrs, "DMRS.top", 
		  min(5000, max(ranks.all$ranks[ranks.all$Method=="dmrseq"])), ".", win, ".RData"))
  grob.dmrseq[[ct]] <- gp + xlab("Methylation difference (effect size)")
  
  
  load(paste0(dn.result.dir, "/BSmooth/",
  		  subfolder, "/corrWithExpressionPlot.", 
		 "BSmooth.n", sampleSize, ".", cond, ".", num.dmrs, "DMRS.top", 
		  min(5000, max(ranks.all$ranks[ranks.all$Method=="BSmooth"])), ".", win, ".RData"))
  grob.bsmooth[[ct]] <- gp + xlab("Methylation difference (effect size)")
		  
		
  load(paste0(dn.result.dir, "/DSS/", 
  		  subfolder, "/corrWithExpressionPlot.", 
		 "DSS.n", sampleSize, ".", cond, ".", num.dmrs, "DMRS.top", 
		  min(5000, max(ranks.all$ranks[ranks.all$Method=="DSS"])), ".", win, ".RData"))
  grob.dss[[ct]] <- gp + xlab("Methylation difference (effect size)")
  
  load(paste0(dn.result.dir, "/metilene/", 
  		  subfolder, "/corrWithExpressionPlot.", 
		 "metilene.n", sampleSize, ".", cond, ".", num.dmrs, "DMRS.top", 
		  min(5000, max(ranks.all$ranks[ranks.all$Method=="metilene"])), ".", win, ".RData"))
  grob.met[[ct]] <- gp + xlab("Methylation difference (effect size)")
  
  ############################################################################
  # fdr rankings / default bsmooth/dss
  subfolder <-  "odds_fdr_Cumulative"
  
  pfile1 <- list.files(path=paste0(dn.result.dir, "/dmrseq_pkg/", subfolder, "/"),
   						pattern=paste0("oddsByRanking.dmrseq.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
  pfile2 <- list.files(path=paste0(dn.result.dir, "/dmrseq_pkg/", subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", win, ".RData")) 
  if (length(pfile1[pfile1 %in% pfile2]) == 1){	   					 					
	  pfile <- paste0(dn.result.dir, "/dmrseq_pkg/", subfolder, "/", 
		pfile1[pfile1 %in% pfile2])
  }else if (length(pfile1[pfile1 %in% pfile2]) > 1){	  #if > 1, pick most recently modified
	  details <- file.info(pfile1[pfile1 %in% pfile2])
	  details <- details[with(details, order(as.POSIXct(mtime), decreasing=TRUE)), ]
	  pfile <- paste0(dn.result.dir, "/dmrseq_pkg/", subfolder, "/", rownames(details)[1])
	  
	  thisOne <- grep(nrow(dmrs), pfile1[pfile1 %in% pfile2])
	  pfile <- paste0(dn.result.dir, "/dmrseq_pkg/", subfolder, "/", pfile1[pfile1 %in% pfile2][thisOne])
  }else{
	 pfile <- "null"
  }
  load(pfile)
  rank.df$Method <- "dmrseq"		
  ranks.last <- rank.df[nrow(rank.df),,drop=FALSE]
  ranks.last$Method <- "dmrseq"  
  fdrRanks <- rank.df$ranks
  ranks.all <- rank.df
   
  mfile1 <- list.files(path=paste0(dn.result.dir, "/metilene_default/",
  						subfolder, "/"),
   						pattern=paste0("oddsByRanking.metilene.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
  mfile2 <- list.files(path=paste0(dn.result.dir, "/metilene_default/",
  				 		subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", win, ".RData")) 
  if (length(mfile1[mfile1 %in% mfile2]) == 1){	   					 					
	  mfile <- paste0(dn.result.dir, "/metilene_default/",
	    subfolder, "/", 
		mfile1[mfile1 %in% mfile2])
  }else if (length(mfile1[mfile1 %in% mfile2]) > 1){	  #if > 1, pick most recently modified
	  details <- file.info(mfile1[mfile1 %in% mfile2])
	  details <- details[with(details, order(as.POSIXct(mtime), decreasing=TRUE)), ]
	  mfile <- paste0(dn.result.dir, "/metilene_default/",
	  			 		subfolder, "/", rownames(details)[1])
	  
	  thisOne <- grep(nrow(dmrs), mfile1[mfile1 %in% mfile2])
	  mfile <- paste0(dn.result.dir, "/metilene_default/",
	   					subfolder, "/", mfile1[mfile1 %in% mfile2][thisOne])
  }else{
	 mfile <- "null"
  }
  load(mfile)
  rank.df$Method <- "metilene"		
  last <- rank.df[nrow(rank.df),,drop=FALSE]
  ranks.all <- rbind(ranks.all, rank.df)
  ranks.last <- rbind(ranks.last, last)
  fdrRanks <- ranks.all$ranks
  
  subfolder <-  "odds_statistic_default"
   
   bfile1 <- list.files(path=paste0(dn.result.dir, "/BSmooth_default/", subfolder, "/"),
   						pattern=paste0("oddsByRanking.BSmooth_default.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
   bfile2 <- list.files(path=paste0(dn.result.dir, "/BSmooth_default/", subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", win, ".RData")) 		  
   if (length(bfile1[bfile1 %in% bfile2]) == 1){	   					 					
   bfile <- paste0(dn.result.dir, "/BSmooth_default/", subfolder, "/", 
   		      bfile1[bfile1 %in% bfile2])
   }else{
   	  bfile <- "null"
   }

   dfile1 <- list.files(path=paste0(dn.result.dir, "/DSS_default/", subfolder, "/"),
   						pattern=paste0("oddsByRanking.DSS_default.n", 
		   					sampleSize, ".", cond, ".", num.dmrs, "DMRS.top")) 
   dfile2 <- list.files(path=paste0(dn.result.dir, "/DSS_default/", subfolder, "/"),
   						pattern=paste0("_min", min.length, "_max", max.length, ".", win, ".RData")) 		  
   if (length(dfile1[dfile1 %in% dfile2]) == 1){	   					 					
   dfile <- paste0(dn.result.dir, "/DSS_default/", subfolder, "/", 
   		      dfile1[dfile1 %in% dfile2])
   }else{
   	  dfile <- "null"
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
	 rank.df <- last
  
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
	 rank.df <- last
	 
	 if (!exists("ranks.all")){
	   ranks.all <- rank.df
	 }else{
	   ranks.all <- rbind(ranks.all, rank.df)
	   ranks.last <- rbind(ranks.last, last)
	 }
   }
  
   ranks.all <- ranks.all[order(ranks.all$Method),]
   
   ranks.all$Method <- as.factor(ranks.all$Method)
   ranks.all$ranks[ranks.all$ranks > 1] <- min(fdrRanks)
   if (length(unique(ranks.all$ranks[ranks.all$Method=="metilene"])) > 1){
   		grad <- c("dmrseq", "metilene")
   }else{
   		grad <- c("dmrseq")
   		ranks.all <- ranks.all[-(which(ranks.all$Method=="metilene")[1]),]
   		ranks.all$ranks[ranks.all$Method=="metilene"] <- min(fdrRanks)
   }
   add2 <- ranks.all[!(ranks.all$Method %in% grad),]
   if(nrow(add2)>0){
	   add2$ranks <- 1
	   ranks.all <- rbind(ranks.all, add2)
   }
   
   ranks.all$fdr <- 1
   ranks.all$fdr[ranks.all$Method %in% grad] <- 0
   ranks.all$fdr <- as.factor(ranks.all$fdr)
   minx <- min(ranks.all$ranks)
   
   ll <- 2^-1
   if (win=="genebody"){ ll <- 2^-0.5}   
   fdr[[ct]] <- ggplot(aes(x=ranks, y=odds), 
                       data=ranks.all) + 
     geom_ribbon(aes(ymin=CI.low, ymax=CI.hi, fill=Method),
                 alpha=0.3) +
     geom_line(aes(group=Method, color=Method, linetype=fdr), lwd=1.1) +
     xlab("FDR threshold (square root scaled)") +
     ylab("Odds") +
     ggtitle(substr(labs[ct], 5, nchar(labs[ct]))) + 
     theme_classic() +
     scale_x_sqrt(breaks=c(0,0.01,0.1,0.25,0.5,0.75,1)) +
     guides(linetype=FALSE)  +
     scale_colour_colorblind() +
     scale_fill_colorblind() +
     scale_y_continuous(trans="log2", limits=c(ll, 2^3),
	                      breaks=c(0.5, 1,2,4,8,16)) 
  	  
  ct <- ct + 1
}   			

 p1 <- plot_grid( odds[[1]] + theme(legend.position="none", axis.title.x=element_blank()), 
 covg[[1]] + theme(legend.position="none", axis.title.x=element_blank()), 
 odds[[2]] + theme(legend.position="none", plot.title =element_blank(), 
 					axis.title.x=element_blank()), 
 covg[[2]] + theme(legend.position="none", plot.title =element_blank(), 
 					axis.title.x=element_blank()), 
 odds[[3]] + theme(legend.position="none", plot.title =element_blank()), 
 covg[[3]] + theme(legend.position="none", plot.title =element_blank()), 
                    labels = c("", "", labs[2], "", labs[3], ""),
                    hjust = 0,
                    nrow = 3,
                    align = "hv",
                    rel_widths = c(4,2.5))

                   
p2 <- plot_grid(near[[1]] + theme(legend.position="none", axis.title.x=element_blank()),
					near[[2]] + theme(legend.position="none", plot.title =element_blank(), axis.title.x=element_blank()),
					near[[3]] + theme(legend.position="none", plot.title =element_blank()),
                   hjust = -1,
                   nrow = 3)
legend <- get_legend(covg[[1]] + 
                theme(legend.position="bottom", 
                	  legend.title=element_text(size=14), 
                	  legend.text=element_text(size=13)))
top.title <- ggdraw() + draw_label(paste0(labs[1]), x=0.003, hjust=0, fontface="bold")

pdf(paste0(result.file.prefix, "/supp_fig8.", win, ".pdf"),
  height=9.6, width=16)
	p <- plot_grid( p1, p2, rel_widths = c(2,1.3), nrow=1)
	p <- plot_grid( top.title, p, legend, nrow=3, rel_heights = c(0.04,1,0.04), hjust=0)
	print(p)
dev.off()
   

p1 <- plot_grid(grob.dmrseq[[1]] + theme(legend.position="none", plot.title =element_text(size=9), axis.title.x=element_blank()), 
				grob.bsmooth[[1]] + theme(legend.position="none", plot.title =element_text(size=9), axis.title.x=element_blank()), 
				grob.dss[[1]]+ theme(legend.position="none", plot.title =element_text(size=9), axis.title.x=element_blank()),
				grob.met[[1]]+ theme(legend.position="none", plot.title =element_text(size=9), axis.title.x=element_blank()),
				grob.dmrseq[[2]]+ theme(legend.position="none", plot.title =element_text(size=9), axis.title.x=element_blank()), 
				grob.bsmooth[[2]]+ theme(legend.position="none", plot.title =element_text(size=9), axis.title.x=element_blank()), 
				grob.dss[[2]]+ theme(legend.position="none", plot.title =element_text(size=9), axis.title.x=element_blank()),
				grob.met[[2]]+ theme(legend.position="none", plot.title =element_text(size=9), axis.title.x=element_blank()),
				grob.dmrseq[[3]]+ theme(legend.position="none", plot.title =element_text(size=9)), 
				grob.bsmooth[[3]]+ theme(legend.position="none", plot.title =element_text(size=9)), 
				grob.dss[[3]]+ theme(legend.position="none", plot.title =element_text(size=9)),
				grob.met[[3]]+ theme(legend.position="none", plot.title =element_text(size=9)),
				nrow=3, hjust=-0.05, vjust=-1.25,
				labels = c(NA, "", "", "", labs[2], "", "", "", labs[3], "", "", ""),
				align="hv",
				scale=0.95)
				
legend <- suppressMessages(get_legend(grob.dmrseq[[1]] + 
 				scale_colour_manual(name="",
                            values=c("red", "black"),
                            labels=c("DE", "Not DE")) +
                scale_fill_discrete(name="",
                         breaks=c(TRUE, FALSE),
                         labels=c("DE Gene", "Not DE")) + 
                theme(legend.position="bottom", 
                	  legend.title=element_text(size=14), 
                	  legend.text=element_text(size=13))))
top.title <- ggdraw() + draw_label(paste0(labs[1]), 
							x=0.003, hjust=-0.027, vjust=0.2, fontface="bold")
				
pdf(paste0(result.file.prefix, "/supp_fig9.", win, ".pdf"),
  height=11, width=17)
	p <- plot_grid( top.title, p1, nrow=2, rel_heights = c(0.05,1))
	p <- plot_grid( p, legend, rel_heights = c(3,0.1), nrow=2)
	p <- p  + draw_label("dmrseq", x = 0.151, y = 0.965,
            	vjust = 1, hjust = 1, size = 14, fontface = 'bold') +
            draw_label("BSmooth", x = 0.405, y = 0.965,
            	vjust = 1, hjust = 1, size = 14, fontface = 'bold') +
            draw_label("DSS", x = 0.640, y = 0.965,
            	vjust = 1, hjust = 1, size = 14, fontface = 'bold') +
            draw_label("metilene", x = 0.905, y = 0.965,
            	vjust = 1, hjust = 1, size = 14, fontface = 'bold')
	print(p)
dev.off()   


p1 <- plot_grid( fdr[[1]] + theme(legend.position="none"),
				 fdr[[2]]+ theme(legend.position="none", 
 					axis.title.y=element_blank()), 
				 fdr[[3]] + theme(legend.position="none", 
 					axis.title.y=element_blank()),
				 nrow=1 )
				 
legend <- get_legend(fdr[[1]] + 
                theme(legend.position="bottom", 
                	  legend.title=element_text(size=14), 
                	  legend.text=element_text(size=13)))
top.title <- ggdraw() + draw_label("FDR versus odds of inverse DE association with Methylation", 
							x=0.003, hjust=-0.37, fontface="bold", vjust=0.8)     
top.title2 <- ggdraw() + draw_label("(B) Murine Leukemia Models", 
							x=0.003, hjust=0, fontface="bold", vjust=0.8)            	  
p <- plot_grid(top.title, p1, legend, rel_heights = c(0.2,1.5,0.06), nrow=3) 

p.fdr.m <- plot_grid(top.title2, p1, legend, rel_heights = c(0.2,1.5,0.12), nrow=3) 


p1 <- plot_grid(p1.fdr.rm, p.fdr.m, nrow=2, rel_heights=c(7, 3.75))
	          
pdf(paste0(result.file.prefix, "/fig3.", win, ".pdf"),
  height=10.75, width=9)
	print(p1)
dev.off()   
   
}