# code to generate Figure 2, Figure 3, Supplementary Figures S3-S5

# source("/n/irizarryfs01_backed_up/kkorthauer/WGBS/dmrseqPaper/FIGURES/CODE/SimulationResults.R")

METHOD <- "dmrseq"
outdir <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/"
result.file.prefix <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/dmrseq_pkg/"
naive.avg <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/dmrseq_pkg_naive_avg/"
naive.sum <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/dmrseq_pkg_naive_sum/"

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
for (sampleSize in allSampleSize){

 # defaults
 f.bsmooth.d <- paste0(outdir, "/BSmooth_default/PowerFDRtable.n", sampleSize, ".",
			   cond, ".", num.dmrs, "DMRs.BSmooth.sim.txt")
 f.dss.d <- paste0(outdir, "/DSS_default/PowerFDRtable.n", sampleSize, ".",
			   cond, ".", num.dmrs, "DMRs.DSS.sim.txt")
 f.met.d <- paste0(outdir, "/metilene_default/PowerFDRtable.n", sampleSize, ".",
			   cond, ".", num.dmrs, "DMRs.metilene.sim.txt")
 f.perm <- paste0(result.file.prefix, "/PowerFDRtable.n", sampleSize, ".",
			   cond, ".", num.dmrs, "DMRs.dmrseq.sim.txt")
 f.perm.avg <- paste0(naive.avg, "/PowerFDRtable.n", sampleSize, ".",
			   cond, ".", num.dmrs, "DMRs.dmrseq.sim.txt")
 f.perm.sum <- paste0(naive.sum, "/PowerFDRtable.n", sampleSize, ".",
			   cond, ".", num.dmrs, "DMRs.dmrseq.sim.txt")
			   
 tab.bsmooth.d <- read.table(f.bsmooth.d, stringsAsFactors=FALSE,
							 header=TRUE)
 tab.dss.d <- read.table(f.dss.d, stringsAsFactors=FALSE,
							 header=TRUE)
 tab.met.d <- read.table(f.met.d, stringsAsFactors=FALSE,
							 header=TRUE)
 tab.perm <- read.table(f.perm, stringsAsFactors=FALSE,
							 header=TRUE)		
 tab.perm.avg <- read.table(f.perm.avg, stringsAsFactors=FALSE,
							 header=TRUE)		
 tab.perm.sum <- read.table(f.perm.sum, stringsAsFactors=FALSE,
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

 f.bsmooth <- paste0(outdir, "/BSmooth/PowerFDRtable.n", sampleSize, ".",
			   cond, ".", num.dmrs, "DMRs.BSmooth.sim.txt")
 f.dss <- paste0(outdir, "/DSS/PowerFDRtable.n", sampleSize, ".",
			   cond, ".", num.dmrs, "DMRs.DSS.sim.txt")
 f.met <- paste0(outdir, "/metilene/PowerFDRtable.n", sampleSize, ".",
			   cond, ".", num.dmrs, "DMRs.metilene.sim.txt")
			   
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
 
 tab.perm.avg <- tab.perm.avg[tab.perm.avg$nDMR > 0,]
 tab.perm.avg$Method <- "Mean Difference"

 tab.perm.sum <- tab.perm.sum[tab.perm.sum$nDMR > 0,]
 tab.perm.sum$Method <- "Area Statistic"
 
 tab.all <- rbind(tab.perm, tab.bsmooth, tab.dss, tab.met)
 
 tab.all$nFalse <- tab.all$fdr * tab.all$nDMR 
 tab.all$nTrue <- tab.all$power * num.dmrs
 
 tab.all$nHitsPerDMR <- tab.all$nDMR * (1-tab.all$fdr) / (tab.all$power * num.dmrs)
 
 tab.all$nHitsTrue <- tab.all$nDMR * (1-tab.all$fdr)
 tab.all <- tab.all[order(tab.all$Method),]
 
 default.current <- default[default$sampleSize == sampleSize,]
 default.current <- default.current[order(default.current$Method),]
 
 tab.naive <- rbind(tab.perm, tab.perm.avg, tab.perm.sum)
 tab.naive$nFalse <- tab.naive$fdr * tab.naive$nDMR 
 tab.naive$nTrue <- tab.naive$power * num.dmrs
 tab.naive$level <- as.numeric(rownames(tab.naive))
 tab.naive$level[tab.naive$level > 1] <- 1 
 
 aorb <- "(A)"
 if (ct==2){ aorb <- "(B)" }			  
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
 
 rocs.n[[ct]] <- ggplot(data=tab.naive, aes(x=fdr, y=power, group=Method, 
 				fill=Method, colour=Method)) +
	 geom_line(alpha=I(0.65)) +
	 geom_point(alpha=I(0.65), stroke=0) + 
	 coord_cartesian(xlim=c(0, 1), 
					 ylim=c(0, 1)) +
	 ggtitle(paste0(aorb, " Simulation D", sampleSize, collapse="")) +
	 labs(x="FDR", y="Power")+ 
	 theme_classic()+
   scale_colour_colorblind()
	 
 fdr.n[[ct]] <- ggplot(data=tab.naive, aes(x=level, y=fdr, group=Method, 
 				fill=Method, colour=Method)) +
 	 geom_abline(slope=1, intercept=0, colour="black", size=0.35, linetype=2, alpha=0.7) + 
	 geom_line(alpha=I(0.65)) +
	 geom_point(alpha=I(0.65), stroke=0) + 
	 coord_cartesian(xlim=c(0, 0.4), 
					 ylim=c(0, 0.4)) +
	 ggtitle(paste0(aorb, " Simulation D", sampleSize, collapse="")) +
	 labs(x="Specified FDR level", y="Observed FDR level") +
	 theme_classic()+
   scale_colour_colorblind()
 
  if (ct ==1){
 	metTab <- tab.met[tab.met$nDMR >0,]
 	metTab$Method <- "metilene"
 	metTab$Sample.Size <- "D2"
 	metTab$level <- as.numeric(rownames(metTab))
 }else{
    tab.met$Method <- "metilene"
    tab.met$Sample.Size <- "D3"
    tab.met$level <- as.numeric(rownames(tab.met))
 	metTab <- rbind(metTab, tab.met[tab.met$nDMR >0,])
 	metTab$Sample.Size <- factor(metTab$Sample.Size)
 	
 	p1m <- ggplot(data=metTab, aes(x=level, y=fdr, group=Sample.Size, colour=Sample.Size)) +
		geom_abline(slope=1, intercept=0, colour="black", size=0.35, linetype=2, alpha=0.7) + 
		geom_line(alpha=I(0.6)) +
		geom_point(alpha=I(0.6), stroke=0) +
		coord_cartesian(xlim=c(0,0.4), ylim=c(0,0.4)) +
		ggtitle("FDR control by Metilene") + 
		labs(x="Specified FDR level", y="Observed FDR level") +
		scale_color_manual(values=c("#E69F00", "#0072B2")) +
		labs(colour="Simulation") +
		theme_classic() 
 }
 
 
 if (sampleSize == 2){
 	alt.file <- paste0(result.file.prefix, "/PowerFDRtable.n3.",
			  cond, ".", num.dmrs, "DMRs.dmrseq.sim.txt")
 
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
    
	p1 <- ggplot(data=tab.perm, aes(x=level, y=fdr, group=Sample.Size, colour=Sample.Size)) +
		geom_abline(slope=1, intercept=0, colour="black", size=0.35, linetype=2, alpha=0.7) + 
		geom_line(alpha=I(0.6)) +
		geom_point(alpha=I(0.6), stroke=0) +
		coord_cartesian(xlim=c(0,0.4), ylim=c(0,0.4)) +
		ggtitle("FDR control by dmrseq") + 
		labs(x="Specified FDR level", y="Observed FDR level") +
		scale_color_manual(values=c("#E69F00", "#0072B2")) +
		labs(colour="Simulation") +
		theme_classic()
		
 }
 
 ct <- ct + 1
}


pdf("/n/irizarryfs01_backed_up/kkorthauer/WGBS/dmrseqPaper/FIGURES/out/supp_fig4.pdf",
  height=3, width=6.5)
  legend <- get_legend(rocs[[1]])
  print(plot_grid(rocs[[1]]+ theme(legend.position="none", 
  				  	plot.title = element_text(hjust = -0.85, face="bold")), 
   				  rocs[[2]]+ theme(legend.position="none", axis.title.y=element_blank(), plot.title = element_text(hjust =  -0.45, face="bold")),
   				  legend, nrow=1, rel_widths=c(1,1,0.3),
   				  hjust=0, vjust=0.75,
   				  scale=0.98))
 
dev.off()


pdf("/n/irizarryfs01_backed_up/kkorthauer/WGBS/dmrseqPaper/FIGURES/out/supp_fig5.pdf",
  height=2.75, width=6.5)
  legend <- get_legend(rocs.n[[1]])
  print(plot_grid(rocs.n[[1]]+ theme(legend.position="none", 
  				  	plot.title = element_text(hjust = -0.85, face="bold")), 
   				  rocs.n[[2]]+ theme(legend.position="none", axis.title.y=element_blank(), plot.title = element_text(hjust =  -0.45, face="bold")),
   				  legend, nrow=1, rel_widths=c(1,1,0.5),
   				  hjust=0, vjust=0.75,
   				  scale=0.98))
 
dev.off()

pdf("/n/irizarryfs01_backed_up/kkorthauer/WGBS/dmrseqPaper/FIGURES/out/fig2.pdf",
  height=3, width=3.75)

  legend <- get_legend(p1m)
   
  pfdr <- plot_grid(p1+ theme(legend.position="none", plot.title = element_text(hjust = -0.375, face="bold")), 
  					legend,
  				  nrow=1,
  				  rel_widths = c(1,0.3),
   				  hjust=0, vjust=0.75,
   				  scale=0.98)
  print(pfdr)
dev.off()		
	
pdf("/n/irizarryfs01_backed_up/kkorthauer/WGBS/dmrseqPaper/FIGURES/out/supp_fig3.pdf",
    height=3, width=3.75)

legend <- get_legend(p1m)

pfdr <- plot_grid(p1m+ theme(legend.position="none", 
                             plot.title = element_text(hjust = -0.375, face="bold")), 
                  legend,
                  nrow=1,
                  rel_widths = c(1,0.3),
                  hjust=0, vjust=0.75,
                  scale=0.98)
print(pfdr)
dev.off()		


# print numbers of dmrs found by each method in the NULL comparisons: 

message("dmrseq at 0.05 FDR:") 
load("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/dmrseq_pkg/regions_control_2_0DMRs.RData")
message(paste0("Sample Size 2: ", sum(regions$qval <= 0.05)))
rm(regions)

load("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/dmrseq_pkg/regions_control_3_0DMRs.RData")
message(paste0("Sample Size 3: ", sum(regions$qval <= 0.05)))
rm(regions)

message("BSmooth with defaults:")
load("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/BSmooth_default/dmrs.BSmooth.n2.control.DEFAULT.Rdata")
message(paste0("Sample Size 2: ", nrow(dmrs)))

load("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/BSmooth_default/dmrs.BSmooth.n3.control.DEFAULT.Rdata")
message(paste0("Sample Size 3: ", nrow(dmrs)))

message("DSS with defaults: ")
load("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/DSS_default/dmrs.DSS.n2.control.DEFAULT.Rdata")
message(paste0("Sample Size 2: ", nrow(dmrs)))

load("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/DSS_default/dmrs.DSS.n3.control.DEFAULT.Rdata")
message(paste0("Sample Size 3: ", nrow(dmrs)))

message("metilene at 0.05 FDR:") 
dmrs <- read.table("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/metilene_default/metilene_output_control_n2_0DMRs.txt", stringsAsFactors=FALSE)
message(paste0("Sample Size 2: ", sum(dmrs$V4 <= 0.05)))
rm(dmrs)

dmrs <- read.table("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/metilene_default/metilene_output_control_n3_0DMRs.txt", stringsAsFactors=FALSE)
message(paste0("Sample Size 3: ", sum(dmrs$V4 <= 0.05)))
rm(dmrs)

## fdr vs power with defaults
permTab <- permTab[,-which(grepl("level", colnames(permTab)))]
permTab$Method <- "dmrseq"
metTab <- metTab[,-which(grepl("level", colnames(metTab)))]
metTab <- metTab[,c(1:5,7,6)]

bsRows <- read.table("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/BSmooth_default/PowerFDRtable.n2.control.3000DMRs.BSmooth.sim.txt", stringsAsFactors=FALSE,
							 header=TRUE)
bsRows$Method <- "BSmooth"
bsRows$Sample.Size <- "D2"

tmp <- read.table("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/BSmooth_default/PowerFDRtable.n3.control.3000DMRs.BSmooth.sim.txt", stringsAsFactors=FALSE,
							 header=TRUE)
tmp$Method <- "BSmooth"
tmp$Sample.Size <- "D3"
bsRows <- rbind(bsRows[rownames(bsRows)=="0.025",], tmp[rownames(tmp)=="0.025",])


# same for dss
Qs <- as.character(c(1e-6, 1e-5, 1e-4, 0.00025, 0.001, 0.002, 0.01, 0.025, 0.05, 0.1))
dsRows <- read.table("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/DSS_default/PowerFDRtable.n2.control.3000DMRs.DSS.sim.txt", stringsAsFactors=FALSE,
							 header=TRUE)
dsRows$Method <- "DSS"
dsRows$Sample.Size <- "D2"


tmp <- read.table("/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/DSS_default/PowerFDRtable.n3.control.3000DMRs.DSS.sim.txt", stringsAsFactors=FALSE,
							 header=TRUE)
tmp$Method <- "DSS"
tmp$Sample.Size <- "D3"


bsRows <- rbind(bsRows, dsRows[!(rownames(dsRows) %in% Qs),], tmp[!(rownames(tmp) %in% Qs),])


permTab <- rbind(permTab, bsRows, metTab)

p2 <- ggplot(data=permTab, aes(x=fdr, y=power, 
		group=interaction(Method, Sample.Size), 
		colour=Sample.Size)) +
	geom_line(alpha=I(0.6)) +
	geom_point(alpha=I(0.6), aes(shape=Method)) +
	ggtitle("Power versus FDR") +
	coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
	labs(x="FDR", y="Power", colour="Simulation")+ 
	#scale_color_manual(values=c("#E69F00", "#0072B2")) +
	scale_shape_manual(values = c(17,16,2,1)) +
	theme_classic()+
  scale_colour_colorblind()	

d2 <- ggplot(data=permTab[permTab$Sample.Size=="D2",], aes(x=fdr, y=power, 
		group=Method, colour=Method)) +
	geom_line(alpha=I(0.6)) +
	geom_point(alpha=I(0.6)) +
	ggtitle("(A) Simulation D2") +
	coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
	labs(x="FDR", y="Power")+ 
	theme_classic() +
	theme(plot.title = element_text(hjust = -0.75, face="bold")) + 
	guides(colour = guide_legend(override.aes =
		 	list(linetype = c("blank", "solid", "blank", "solid")))) +
  scale_colour_colorblind()
	
d3 <- ggplot(data=permTab[permTab$Sample.Size=="D3",], aes(x=fdr, y=power, 
		group=Method, colour=Method)) +
	geom_line(alpha=I(0.6), show.legend=FALSE) +
	geom_point(alpha=I(0.6)) +
	ggtitle("(B) Simulation D3") +
	coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
	labs(x="FDR", y="Power")+ 
	theme_classic() +
	theme(legend.position="none", axis.title.y=element_blank(),
		plot.title = element_text(hjust = -0.4, face="bold")) +
  scale_colour_colorblind()
			
pdf("/n/irizarryfs01_backed_up/kkorthauer/WGBS/dmrseqPaper/FIGURES/out/fig3.pdf",
  height=3, width=6.5)

  legend <- get_legend(d2)
  print(plot_grid(d2 + theme(legend.position="none"), d3, legend, nrow=1,
  				  rel_widths=c(1,1,0.3),
   				  hjust=0, vjust=0.75,
   				  scale=0.98))
dev.off()		


print(default)



 