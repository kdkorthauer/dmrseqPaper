# code to generate Supplementary Figures S10-12 (Fig 2 subsetted by various factors

# source("/n/irizarryfs01_backed_up/kkorthauer/WGBS/PAPER/FIGURES/CODE/SimulationResults_subset.R")

######################################################
### parameters to change to run on your own system ###
######################################################
# change the following the root results directory,
outdir <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/"
# change the following to represent the dmrseq results directory
result.file.prefix <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/dmrseq_pkg/"
# change the following to where you'd like to save the figure output
figure.file.prefix <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/PAPER/FIGURES/out/"
######################################################
###         end of parameters to change            ###
######################################################



METHOD <- "dmrseq"
num.dmrs <- 3000
min.length <- max.length <- NULL
cond <- "control"
allSampleSize <- c(2,3)
rocs <- rocs.d <- rocs.n <- fdr.n <- vector("list", length(allSampleSize))				   
				   
library(ggplot2) 
library(cowplot)
library(ggthemes)

default <- data.frame(Method = c("dmrseq", "dmrseq", "BSmooth", "BSmooth", 
						"DSS", "DSS", "metilene", "metilene"),
					  sampleSize = rep(c(2,3),4),
					  dmrs = rep(NA, 8),
					  fp = rep(NA, 8), 
					  tp = rep(NA, 8),
					  power = rep(NA, 8),
					  fdr = rep(NA, 8))

ct <- 1
fnum <- 10

subset.list <- vector("list", 2)

for (type in c("effsize", "density", "coverage")){
 
 t.title <- "Effect Size"
 if(type=="density"){ t.title <- "CpG density"}
 if(type=="coverage"){ t.title <- "Coverage"}
 
 d2 <- d3 <- vector("list", 2)
 names(d2) <- names(d3) <- c("low", "high")

 for (lev in c("low", "high")){ 
   l.title <- "Low "
   if(lev=="high"){ l.title <- "High "}
   
 for (sampleSize in allSampleSize){
 # defaults
 f.bsmooth.d <- list.files(path=paste0(outdir, "/BSmooth_default/"),
                           pattern=paste0("PowerFDRtable.n", sampleSize, ".",
			               cond, ".*", "DMRs.BSmooth.", lev, ".", type, ".sim.txt"),
			               full.names = TRUE)
 f.dss.d <- list.files(path=paste0(outdir, "/DSS_default/"),
                       pattern=paste0("PowerFDRtable.n", sampleSize, ".",
			                    cond, ".*", "DMRs.DSS.", lev, ".", type, ".sim.txt"),
			               full.names = TRUE)
 f.met.d <- list.files(path=paste0(outdir, "/metilene_default/"),
                       pattern=paste0("PowerFDRtable.n", sampleSize, ".",
			           cond, ".*", "DMRs.metilene.", lev, ".", type, ".sim.txt"),
			               full.names = TRUE)
 f.perm <- list.files(path=result.file.prefix, 
                      pattern=paste0("PowerFDRtable.n", sampleSize, ".",
			          cond, ".*", "DMRs.dmrseq.", lev, ".", type, ".sim.txt"),
			               full.names = TRUE)
			   
 tab.bsmooth.d <- read.table(f.bsmooth.d, stringsAsFactors=FALSE,
							 header=TRUE)
 tab.dss.d <- read.table(f.dss.d, stringsAsFactors=FALSE,
							 header=TRUE)
 tab.met.d <- read.table(f.met.d, stringsAsFactors=FALSE,
							 header=TRUE)
 tab.perm <- read.table(f.perm, stringsAsFactors=FALSE,
							 header=TRUE)		
							 						 							 	   
 default[default$sampleSize == sampleSize & default$Method=="BSmooth",]$dmrs <- 
 			tab.bsmooth.d$nDMR[rownames(tab.bsmooth.d)=="0.025"]
 default[default$sampleSize == sampleSize & default$Method=="BSmooth",]$tp <- 
 			round(tab.bsmooth.d$power[rownames(tab.bsmooth.d)=="0.025"]*num.dmrs)
 default[default$sampleSize == sampleSize & default$Method=="BSmooth",]$fp <- 
 			round(tab.bsmooth.d$fdr[rownames(tab.bsmooth.d)=="0.025"]*
 			tab.bsmooth.d$nDMR[rownames(tab.bsmooth.d)=="0.025"])
 default[default$sampleSize == sampleSize & default$Method=="BSmooth",]$power <-
 			tab.bsmooth.d$power[rownames(tab.bsmooth.d)=="0.025"]
 default[default$sampleSize == sampleSize & default$Method=="BSmooth",]$fdr <-
 			tab.bsmooth.d$fdr[rownames(tab.bsmooth.d)=="0.025"]
 			
 Qs <- as.character(c(1e-6, 1e-5, 1e-4, 0.00025, 0.001, 0.002, 0.01, 0.025, 0.05, 0.1))
 default[default$sampleSize == sampleSize & default$Method=="DSS",]$dmrs <- 
 			tab.dss.d$nDMR[!rownames(tab.dss.d) %in% Qs]
 default[default$sampleSize == sampleSize & default$Method=="DSS",]$tp <- 
 			round(tab.dss.d$power[!rownames(tab.dss.d) %in% Qs]*num.dmrs)
 default[default$sampleSize == sampleSize & default$Method=="DSS",]$fp <- 
 			round(tab.dss.d$fdr[!rownames(tab.dss.d) %in% Qs]*
 			tab.dss.d$nDMR[!rownames(tab.dss.d) %in% Qs])
 default[default$sampleSize == sampleSize & default$Method=="DSS",]$power <- 
 			tab.dss.d$power[!rownames(tab.dss.d) %in% Qs]
 default[default$sampleSize == sampleSize & default$Method=="DSS",]$fdr <- 
 			tab.dss.d$fdr[!rownames(tab.dss.d) %in% Qs]
 	
 			
 default[default$sampleSize == sampleSize & default$Method=="dmrseq",]$dmrs <- 
 			tab.perm$nDMR[rownames(tab.perm)=="0.05"]
 default[default$sampleSize == sampleSize & default$Method=="dmrseq",]$tp <- 
 			round(tab.perm$power[rownames(tab.perm)=="0.05"]*num.dmrs)
 default[default$sampleSize == sampleSize & default$Method=="dmrseq",]$fp <- 
 			round(tab.perm$fdr[rownames(tab.perm)=="0.05"]*
 			tab.perm$nDMR[rownames(tab.perm)=="0.05"])
 default[default$sampleSize == sampleSize & default$Method=="dmrseq",]$power <- 
 			tab.perm$power[rownames(tab.perm)=="0.05"]
 default[default$sampleSize == sampleSize & default$Method=="dmrseq",]$fdr <- 
 			tab.perm$fdr[rownames(tab.perm)=="0.05"]
 
 default[default$sampleSize == sampleSize & default$Method=="metilene",]$dmrs <- 
 			tab.met.d$nDMR[rownames(tab.met.d)=="0.05"]
 default[default$sampleSize == sampleSize & default$Method=="metilene",]$tp <- 
 			round(tab.met.d$power[rownames(tab.met.d)=="0.05"]*num.dmrs)
 default[default$sampleSize == sampleSize & default$Method=="metilene",]$fp <- 
 			round(tab.met.d$fdr[rownames(tab.met.d)=="0.05"]*
 			tab.met.d$nDMR[rownames(tab.met.d)=="0.05"])
 default[default$sampleSize == sampleSize & default$Method=="metilene",]$power <- 
 			tab.met.d$power[rownames(tab.met.d)=="0.05"]
 default[default$sampleSize == sampleSize & default$Method=="metilene",]$fdr <- 
 			tab.met.d$fdr[rownames(tab.met.d)=="0.05"]

 f.bsmooth <- list.files(path=paste0(outdir, "/BSmooth/"),
                           pattern=paste0("PowerFDRtable.n", sampleSize, ".",
			               cond, ".*", "DMRs.BSmooth.", lev, ".", type, ".sim.txt"),
			               full.names = TRUE)
 f.dss <- list.files(path=paste0(outdir, "/DSS/"),
                       pattern=paste0("PowerFDRtable.n", sampleSize, ".",
			                    cond, ".*", "DMRs.DSS.", lev, ".", type, ".sim.txt"),
			               full.names = TRUE)
 f.met <- list.files(path=paste0(outdir, "/metilene/"),
                       pattern=paste0("PowerFDRtable.n", sampleSize, ".",
			           cond, ".*", "DMRs.metilene.", lev, ".", type, ".sim.txt"),
			               full.names = TRUE)
			   
 tab.bsmooth <- read.table(f.bsmooth, stringsAsFactors=FALSE,
							 header=TRUE)
 tab.dss <- read.table(f.dss, stringsAsFactors=FALSE,
							 header=TRUE)
 tab.met <- read.table(f.met, stringsAsFactors=FALSE,
							 header=TRUE)
 
 tab.bsmooth <- tab.bsmooth[tab.bsmooth$nDMR > 0,]
 tab.bsmooth$Method <- "BSmooth"
 
 tab.dss <- tab.dss[tab.dss$nDMR > 0,]
 tab.dss$Method <- "DSS"
 
 tab.met <- tab.met[tab.met$nDMR > 0,]
 tab.met$Method <- "metilene"
 
 tab.perm <- tab.perm[tab.perm$nDMR > 0,]
 tab.perm$Method <- "dmrseq"
  
 tab.all <- rbind(tab.perm, tab.bsmooth, tab.dss, tab.met)
 
 tab.all$nFalse <- tab.all$fdr * tab.all$nDMR 
 tab.all$nTrue <- tab.all$power * num.dmrs
 
 tab.all$nHitsPerDMR <- tab.all$nDMR * (1-tab.all$fdr) / (tab.all$power * num.dmrs)
 
 tab.all$nHitsTrue <- tab.all$nDMR * (1-tab.all$fdr)
 tab.all <- tab.all[order(tab.all$Method),]
 
 default.current <- default[default$sampleSize == sampleSize,]
 default.current <- default.current[order(default.current$Method),]
 
 
 aorb <- "(A)"
 if (sampleSize==3){ aorb <- "(B)" }			  
 rocs[[ct]] <- ggplot(data=tab.all, aes(x=nFalse, y=nTrue, group=Method, 
 				fill=Method, colour=Method)) +
	 geom_line(alpha=I(0.65)) +
	 geom_point(alpha=I(0.65), stroke=0) + 
	 coord_cartesian(xlim=c(0, 5000), 
					 ylim=c(0, 2500)) +
	 ggtitle(paste0(aorb, " Simulation D", sampleSize, collapse="")) +
	 labs(x="False Positives", y="True Positives")+ 
	 theme_classic()+
   scale_colour_colorblind()
 
 if (sampleSize==2){
 	metTab <- tab.met[tab.met$nDMR >0,]
 	metTab$Method <- "metilene"
 	metTab$Sample.Size <- "D2"
 	metTab$level <- as.numeric(rownames(metTab))
 	
 	bsRows <- tab.bsmooth.d 
 	bsRows$Method <- "BSmooth"
    bsRows$Sample.Size <- "D2"
    
	# same for dss
	Qs <- as.character(c(1e-6, 1e-5, 1e-4, 0.00025, 0.001, 0.002, 0.01, 0.025, 0.05, 0.1))
	dsRows <- tab.dss.d
	dsRows$Method <- "DSS"
	dsRows$Sample.Size <- "D2"

 }else{
    tab.met$Method <- "metilene"
    tab.met$Sample.Size <- "D3"
    tab.met$level <- as.numeric(rownames(tab.met))
 	metTab <- rbind(metTab, tab.met[tab.met$nDMR >0,])
 	metTab$Sample.Size <- factor(metTab$Sample.Size)
 	
    tab.bsmooth.d$Method <- "BSmooth"
    tab.bsmooth.d$Sample.Size <- "D3"
    bsRows <- rbind(bsRows[rownames(bsRows)=="0.025",], tab.bsmooth.d[rownames(tab.bsmooth.d)=="0.025",])
    
	tab.dss.d$Method <- "DSS"
	tab.dss.d$Sample.Size <- "D3"
	
	bsRows <- rbind(bsRows, dsRows[!(rownames(dsRows) %in% Qs),], tab.dss.d[!(rownames(tab.dss.d) %in% Qs),])
 }
 
 if (sampleSize == 2){
 	alt.file <- list.files(path=result.file.prefix, 
                      pattern=paste0("PowerFDRtable.n3", ".",
			          cond, ".*", "DMRs.dmrseq.", lev, ".", type, ".sim.txt"),
			               full.names = TRUE)
 
	alt.tab.perm <- read.table(alt.file, stringsAsFactors=FALSE,
								header=TRUE)
	tab.perm <- tab.perm[,c(1:5)]
	tab.perm$Sample.Size <- sampleSize
	alt.tab.perm$Sample.Size <- 3
	tab.perm$level <- as.numeric(rownames(tab.perm))
	alt.tab.perm$level <- as.numeric(rownames(alt.tab.perm))
	tab.perm <- rbind(tab.perm, alt.tab.perm)
	tab.perm$Method <- "dmrseq"
	tab.perm$Sample.Size[tab.perm$Sample.Size == 2] <- "D2"
    tab.perm$Sample.Size[tab.perm$Sample.Size == 3] <- "D3"
	tab.perm$Sample.Size <- factor(tab.perm$Sample.Size)
	permTab <- tab.perm
 }
 
 ct <- ct + 1
}


## fdr vs power with defaults
permTab <- permTab[,-which(grepl("level", colnames(permTab)))]
permTab$Method <- "dmrseq"
metTab <- metTab[,-which(grepl("level", colnames(metTab)))]
metTab <- metTab[,c(1:5,7,6)]

permTab <- rbind(permTab, bsRows, metTab)

let <- c("(A) ", "(B) ")
if (lev == "high"){ let <- c("(C) ", "(D) ") }

d2[[lev]] <- ggplot(data=permTab[permTab$Sample.Size=="D2",], aes(x=fdr, y=power, 
		group=Method, colour=Method)) +
	geom_line(alpha=I(0.6)) +
	geom_point(alpha=I(0.6)) +
	ggtitle(paste0(let[1], l.title, t.title, ", Simulation D2")) +
	coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
	labs(x="FDR", y="Power")+ 
	theme_classic() +
	theme(plot.title = element_text(hjust = 1, face="bold"),
	      legend.position="bottom",
	      legend.text = element_text(size=12)) + 
	guides(colour = guide_legend(override.aes =
		 	list(linetype = c("blank", "solid", "blank", "solid")))) +
  scale_colour_colorblind()
	
d3[[lev]] <- ggplot(data=permTab[permTab$Sample.Size=="D3",], aes(x=fdr, y=power, 
		group=Method, colour=Method)) +
	geom_line(alpha=I(0.6), show.legend=FALSE) +
	geom_point(alpha=I(0.6)) +
	ggtitle(paste0(let[2], l.title, t.title, ", Simulation D3")) +
	coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
	labs(x="FDR", y="Power")+ 
	theme_classic() +
	theme(legend.position="none", axis.title.y=element_blank(),
		plot.title = element_text(hjust = 1, face="bold")) +
  scale_colour_colorblind()
} # end loop over level
			
# change to a four-panel figure with (A-B) low and (C-D) high
pdf(paste0(figure.file.prefix, "/supp_fig", fnum, ".pdf"),
  height=6.5, width=6.5)

  legend <- get_legend(d2[[2]])
  
  main_grid <- (plot_grid(d2[[1]] + theme(legend.position="none"), d3[[1]],
                  d2[[2]] + theme(legend.position="none"), d3[[2]],
                  nrow=2, ncol=2,
  				  rel_widths=c(1,1), rel_heights=c(1,1),
   				  hjust=-0.4, vjust=0.75,
   				  scale=0.98))
   				  
  print(plot_grid(main_grid,
                  legend, nrow=2, ncol=1,
  				  rel_heights=c(1,0.05),
   				  hjust=-0.4, vjust=0.75,
   				  scale=0.98))
dev.off()		

fnum<-fnum+1
}


 