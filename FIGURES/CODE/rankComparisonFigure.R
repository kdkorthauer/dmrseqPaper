# Code to generate intermediate files used to make 
# Figure 5 (see 'fig4.Rmd' for code to produce Figure 4, which calls this script)
# source("rankComparisonFigure.R")


######################################################
### parameters to change to run on your own system ###
######################################################
# change the following to represent the main results directory for Roadmap
rm.result.dir <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/RESULTS/"
# change the following to represent the Roadmap data directory
rm.data.dir <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/ROADMAP/DATA/"
# change the following to represent the main results directory for DNMT3a
dn.result.dir <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/RESULTS/"
# change the following to represent the Roadmap data directory
dn.data.dir <- "/n/irizarryfs01/kkorthauer/WGBS/DNMT3A/"

# change the following to where you'd like to save the figure output
figure.file.prefix <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/PAPER/FIGURES/out/"
######################################################
###         end of parameters to change            ###
######################################################

##### Roadmap Results Plots
library(dmrseq)
library(GenomicRanges)
library(dplyr)
library(corrplot)
library(viridis)

sampleSize <- 2
num.dmrs <- 0
num.to.plot <- 1  #(how many top regions are plotted in the compPlots)
min.length <- max.length <- NULL
pval.thresh <- 0.05
allConditions <- c("Right Ventricle_Small Intestine",
                   "Right Ventricle_Sigmoid Colon",
                   "Left Ventricle_Small Intestine",
                   "Left Ventricle_Sigmoid Colon",	   
                   "Sigmoid Colon_Small Intestine")
ct <- 1
testCovariate <- "Cond"
adjustCovariate=NULL
genomeName="hg38"
time <- "roadmap"
ord <- "stat"


for(cond in allConditions){
  
  tiss1 <- strsplit(cond, "_")[[1]][1]
  tiss2 <- strsplit(cond, "_")[[1]][2]
  
  annot <- getAnnot(genomeName) 
  
  # construct bsseq object for this comparison
  load(paste0(rm.data.dir, time, "_", tiss1, "_all_bsseq_id.RData")); show(bs)
meta <- pData(bs)
bs <- BSseq(chr = as.character(seqnames(bs)), pos = start(bs),
            M =  as.matrix(getCoverage(bs, type = "M")), 
            Cov = as.matrix(getCoverage(bs, type = "Cov")), 
            sampleNames=sampleNames(bs))
bs_tiss1 <- bs
rm(bs); gc()
rownames(meta) <- meta$SampleNames
pData(bs_tiss1) <- meta
load(paste0(rm.data.dir, time, "_", tiss2, "_all_bsseq_id.RData")); show(bs)
meta <- pData(bs)
bs <- BSseq(chr = as.character(seqnames(bs)), pos = start(bs),
            M =  as.matrix(getCoverage(bs, type = "M")), 
            Cov = as.matrix(getCoverage(bs, type = "Cov")), 
            sampleNames=sampleNames(bs))
rownames(meta) <- meta$SampleNames
pData(bs) <- meta
bs <- bsseq::combine(bs, bs_tiss1)
rm(bs_tiss1); gc()

  
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
  
  # construct the design matrix using the pData of bs
  if (ncol(pData(bs)) < max(testCovariate, adjustCovariate)){
    stop(paste0("Error: pData(bs) has too few columns.  Please specify valid ",
                "covariates to use in the analysis"))
  }
  
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
    colnames(design)[,(max(coeff)+1):ncol(design)] <- 
      colnames(pData(bs))[adjustCovariate]
  }else{
    design <- model.matrix( ~ testCov)
    colnames(design)[coeff] <- colnames(pData(bs))[testCovariate]
  }
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  # set up colors and sample labels for plotting 
  # one unique color per unique value of the covariate of interest
  cov.unique <- unique(design[,coeff])
  colors <- gg_color_hue(length(cov.unique))
  if (length(cov.unique) == 2){
    colors <- c("mediumblue", "deeppink1")
  }
  colors <- cbind(cov.unique, colors[rank(as.numeric(cov.unique))])
  z <- colors[,2][match(design[,coeff], colors[,1])]
  pData(bs)$col <- as.character(z)
  pData(bs)$label <- paste0(pData(bs)[,testCovariate])
  
  # load DMR objects for each method
  load(paste0(rm.result.dir, "/dmrseq_pkg/regions_", cond, "_", 
              sampleSize, "_", num.dmrs, "DMRs.RData"))
  
  dmrseq <- regions[regions$qval < pval.thresh & !is.na(regions$qval),]
  dmrseq.all <- regions
  rm(regions)
  
  load(paste0(rm.result.dir, "/BSmooth_default/",
              "dmrs.BSmooth.n", sampleSize, ".", cond, ".DEFAULT.Rdata"))
  bsmooth <- dmrs
  bsmooth$stat <- bsmooth$areaStat 
  bsmooth$beta <- bsmooth$meanDiff
  bsmooth$indexStart <- bsmooth$idxStart
  bsmooth$indexEnd <- bsmooth$idxEnd
  bsmooth$L <- bsmooth$n
  
  load(paste0(rm.result.dir, "/DSS_default/",
              "dmrs.DSS.n", sampleSize, ".", cond, ".DEFAULT.Rdata"))
  rm(dmlTest.sm)
  dss <- dmrs
  dss$stat <- dss$areaStat 
  dss$beta <- dss$diff.Methy 
  dss$indexStart <- dss$start
  dss$indexEnd <- dss$end
  dss$L <- dss$nCG
  rm(dmrs)
  
  met <- read.table(paste0(rm.result.dir, "/metilene_default/",
                           "metilene_output_", gsub(" ", "_", cond), "_n", sampleSize, "_0DMRs.txt"), stringsAsFactors=FALSE)
  colnames(met) <- c("chr", "start", "end", "qval", "beta", "L", "pval1", 
                     "pval2", "mean1", "mean2")
  met.all <- met
  met.all$stat <- met.all$beta * met.all$L
  met.all$indexStart <- match(paste0(met$chr, "_", met$start+1), 
                          paste0(as.character(seqnames(bs)), "_", start(bs)))
  met.all$indexEnd <- match(paste0(met$chr, "_", met$end), 
                        paste0(as.character(seqnames(bs)), "_", start(bs)))
  met <- met.all[met.all$qval < pval.thresh,]
  
  ## examine overlaps
  o.bsmooth <- findOverlaps(makeGRangesFromDataFrame(dmrseq.all), 
                            makeGRangesFromDataFrame(bsmooth))
  o.dss <- findOverlaps(makeGRangesFromDataFrame(dmrseq.all), 
                            makeGRangesFromDataFrame(dss))
  o.met <- findOverlaps(makeGRangesFromDataFrame(dmrseq.all), 
                            makeGRangesFromDataFrame(met.all))
  i.all <- unique(o.bsmooth@from[(o.bsmooth@from %in% o.dss@from) & 
                          (o.bsmooth@from %in% o.met@from)])
  dmrseq.common <- dmrseq.all[i.all,]
  dmrseq.common$dmrseqRank <- rank(-abs(dmrseq.common$stat))/nrow(dmrseq.common)
  
  o.bsmooth <- findOverlaps(makeGRangesFromDataFrame(dmrseq.common), 
                            makeGRangesFromDataFrame(bsmooth))
  o.dss <- findOverlaps(makeGRangesFromDataFrame(dmrseq.common), 
                        makeGRangesFromDataFrame(dss))
  o.met <- findOverlaps(makeGRangesFromDataFrame(dmrseq.common), 
                        makeGRangesFromDataFrame(met.all))
  bsmooth.common <- bsmooth[o.bsmooth@to,]
  bsmooth.common$i.dmrseq <- o.bsmooth@from
  dss.common <- dss[o.dss@to,]
  dss.common$i.dmrseq <- o.dss@from
  met.common <- met.all[o.met@to,]
  met.common$i.dmrseq <- o.met@from
  
  bsmooth.common$area.rank <- rank(-abs(bsmooth.common$areaStat))/nrow(bsmooth.common)
  bsmooth.common$mean.rank <- rank(-abs(bsmooth.common$meanDiff))/nrow(bsmooth.common)
  dss.common$area.rank <- rank(-abs(dss.common$areaStat))/nrow(dss.common)
  dss.common$mean.rank <- rank(-abs(dss.common$diff.Methy))/nrow(dss.common)
  met.common$rank.qval <- rank(met.common$qval)/nrow(met.common)
  
  # add these ranks to common region table
  # average over all overlapping regions if there are more than 1
  dmrseq.common <- cbind(dmrseq.common,
    (bsmooth.common %>% group_by(i.dmrseq) %>% 
    summarize(BSmooth.areaRank=mean(area.rank),
              BSmooth.avgRank=mean(mean.rank)))[,-1],
    (dss.common %>% group_by(i.dmrseq) %>% 
       summarize(DSS.areaRank=mean(area.rank),
                 DSS.avgRank=mean(mean.rank)))[,-1],
    (met.common %>% group_by(i.dmrseq) %>% 
        summarize(Metilene.qvalRank=mean(rank.qval)))[,-1])
  
  col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                             "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))
  ranks <- dmrseq.common[,11:16]
 
  # find an enriched set of regions that are given high ranking by other 
  # method(s) but have high(er) FDR by dmrseq
  bsmooth.high <- dmrseq.common[(dmrseq.common$BSmooth.areaRank < 0.1 |
                                 dmrseq.common$BSmooth.avgRank < 0.1) &
                                 dmrseq.common$qval > 0.5, ]
  dss.high <- dmrseq.common[(dmrseq.common$DSS.areaRank < 0.1 |
                             dmrseq.common$DSS.avgRank < 0.1) &
                             dmrseq.common$qval > 0.50, ]
  met.high <- dmrseq.common[(dmrseq.common$Metilene.areaRank < 0.1 |
                             dmrseq.common$Metilene.avgRank < 0.1 |
                             dmrseq.common$Metilene.qvalRank < 0.05) &
                             dmrseq.common$qval > 0.50, ]
  all.high <- dmrseq.common[dmrseq.common$BSmooth.areaRank < 0.1 &
                              dmrseq.common$DSS.areaRank < 0.1 &
                              dmrseq.common$Metilene.qvalRank < 0.1 &
                              dmrseq.common$qval < 0.1, ]
  avg.high <- dmrseq.common[(dmrseq.common$BSmooth.avgRank < 0.2 &
                              dmrseq.common$DSS.avgRank < 0.2 & 
                              dmrseq.common$qval > 0.5), ]
  area.high <- dmrseq.common[(dmrseq.common$BSmooth.areaRank < 0.1 &
                               dmrseq.common$DSS.areaRank < 0.1 & 
                               dmrseq.common$qval > 0.5), ]
  
  
  # plot these guys
  if (nrow(bsmooth.high)>0){
    compareTrack <- GenomicRangesList(
      dmrseq=makeGRangesFromDataFrame(bsmooth.high, keep.extra.columns=TRUE),
      BSmooth=makeGRangesFromDataFrame(bsmooth.common[findOverlaps(makeGRangesFromDataFrame(bsmooth.common),
                                                                   makeGRangesFromDataFrame(bsmooth.high))@from,],
                                       keep.extra.columns=TRUE),
      DSS=makeGRangesFromDataFrame(dss.common[findOverlaps(makeGRangesFromDataFrame(dss.common), 
                                                           makeGRangesFromDataFrame(bsmooth.high))@from,],
                                   keep.extra.columns=TRUE),
      Metilene=makeGRangesFromDataFrame(met.common[findOverlaps(makeGRangesFromDataFrame(met.common), 
                                                                makeGRangesFromDataFrame(bsmooth.high))@from,],
                                        keep.extra.columns=TRUE))
    
    pdf(paste0(figure.file.prefix, "/rankPlots/",
               "DiscordantRanks_BSmooth_", cond, ".pdf"), width=6, height=3)
    plotDMRs(bs, regions=bsmooth.high, 
             extend=(bsmooth.high$end - bsmooth.high$start + 1)/2, 
             addRegions=bsmooth.high, regionCol=NULL,
             addPoints=TRUE, pointsMinCov=1, qval=FALSE,
             stat=FALSE, includeYlab=TRUE, compareTrack=compareTrack,
             labelCols=c("area.rank", "mean.rank",  "qval"))
    dev.off()
  }
  
  if (nrow(dss.high)>0){
    compareTrack <- GenomicRangesList(
      dmrseq=makeGRangesFromDataFrame(dss.high, keep.extra.columns=TRUE),
      BSmooth=makeGRangesFromDataFrame(bsmooth.common[findOverlaps(makeGRangesFromDataFrame(bsmooth.common),
                                                                   makeGRangesFromDataFrame(dss.high))@from,],
                                       keep.extra.columns=TRUE),
      DSS=makeGRangesFromDataFrame(dss.common[findOverlaps(makeGRangesFromDataFrame(dss.common), 
                                                           makeGRangesFromDataFrame(dss.high))@from,],
                                   keep.extra.columns=TRUE),
      Metilene=makeGRangesFromDataFrame(met.common[findOverlaps(makeGRangesFromDataFrame(met.common), 
                                                                makeGRangesFromDataFrame(dss.high))@from,],
                                        keep.extra.columns=TRUE))
    
    pdf(paste0(figure.file.prefix, "/rankPlots/",
               "DiscordantRanks_DSS_", cond, ".pdf"), width=6, height=3)
    plotDMRs(bs, regions=dss.high, 
             extend=(dss.high$end - dss.high$start + 1)/2, 
             addRegions=dss.high, regionCol=NULL,
             addPoints=TRUE, pointsMinCov=1,qval=FALSE,
             stat=FALSE, includeYlab=TRUE, compareTrack=compareTrack,
             labelCols=c("area.rank", "mean.rank",  "qval"))
    dev.off()
  }
  
  if (nrow(met.high)>0){
    compareTrack <- GenomicRangesList(
      dmrseq=makeGRangesFromDataFrame(met.high, keep.extra.columns=TRUE),
      BSmooth=makeGRangesFromDataFrame(bsmooth.common[findOverlaps(makeGRangesFromDataFrame(bsmooth.common),
                                                                   makeGRangesFromDataFrame(met.high))@from,],
                                       keep.extra.columns=TRUE),
      DSS=makeGRangesFromDataFrame(dss.common[findOverlaps(makeGRangesFromDataFrame(dss.common), 
                                                           makeGRangesFromDataFrame(met.high))@from,],
                                   keep.extra.columns=TRUE),
      Metilene=makeGRangesFromDataFrame(met.common[findOverlaps(makeGRangesFromDataFrame(met.common), 
                                                                makeGRangesFromDataFrame(met.high))@from,],
                                        keep.extra.columns=TRUE))
    
    pdf(paste0(figure.file.prefix, "/rankPlots/",
               "DiscordantRanks_Metilene_", cond, ".pdf"), width=6, height=3)
    plotDMRs(bs, regions=met.high, 
             extend=(met.high$end - met.high$start + 1)/2, 
             addRegions=met.high, regionCol=NULL,
             addPoints=TRUE, pointsMinCov=1,qval=FALSE,
             stat=FALSE, includeYlab=TRUE, compareTrack=compareTrack,
             labelCols=c("area.rank", "mean.rank",  "qval"))
    dev.off()
  }
  
  if (nrow(all.high)>0){
    all.high <- all.high[1:min(nrow(all.high), 25),]
    compareTrack <- GenomicRangesList(
      dmrseq=makeGRangesFromDataFrame(all.high, keep.extra.columns=TRUE),
      BSmooth=makeGRangesFromDataFrame(bsmooth.common[findOverlaps(makeGRangesFromDataFrame(bsmooth.common),
                                                                   makeGRangesFromDataFrame(all.high))@from,],
                                       keep.extra.columns=TRUE),
      DSS=makeGRangesFromDataFrame(dss.common[findOverlaps(makeGRangesFromDataFrame(dss.common), 
                                                           makeGRangesFromDataFrame(all.high))@from,],
                                   keep.extra.columns=TRUE),
      Metilene=makeGRangesFromDataFrame(met.common[findOverlaps(makeGRangesFromDataFrame(met.common), 
                                                                makeGRangesFromDataFrame(all.high))@from,],
                                        keep.extra.columns=TRUE))
    
    pdf(paste0(figure.file.prefix, "/rankPlots/",
               "ConcordantRanks_", cond, ".pdf"), width=6, height=3)
    plotDMRs(bs, regions=all.high, 
             extend=(all.high$end - all.high$start + 1)/2, 
             addRegions=all.high, regionCol=NULL,
             addPoints=TRUE, pointsMinCov=1,qval=FALSE,
             stat=FALSE, includeYlab=TRUE, compareTrack=compareTrack,
             labelCols=c("area.rank", "mean.rank",  "qval"))
    dev.off()
  }
  
  if (nrow(avg.high)>0){
    compareTrack <- GenomicRangesList(
      dmrseq=makeGRangesFromDataFrame(avg.high, keep.extra.columns=TRUE),
      BSmooth=makeGRangesFromDataFrame(bsmooth.common[findOverlaps(makeGRangesFromDataFrame(bsmooth.common),
                                                                   makeGRangesFromDataFrame(avg.high))@from,],
                                       keep.extra.columns=TRUE),
      DSS=makeGRangesFromDataFrame(dss.common[findOverlaps(makeGRangesFromDataFrame(dss.common), 
                                                           makeGRangesFromDataFrame(avg.high))@from,],
                                   keep.extra.columns=TRUE),
      Metilene=makeGRangesFromDataFrame(met.common[findOverlaps(makeGRangesFromDataFrame(met.common), 
                                                                makeGRangesFromDataFrame(avg.high))@from,],
                                        keep.extra.columns=TRUE))
    
    pdf(paste0(figure.file.prefix, "/rankPlots/",
               "DiscordantRanks_", "Avg_", cond, ".pdf"), width=6, height=3)
    plotDMRs(bs, regions=avg.high, 
             extend=(avg.high$end - avg.high$start + 1)/2, 
             addRegions=avg.high, regionCol=NULL,
             addPoints=TRUE, pointsMinCov=1,qval=FALSE,
             stat=FALSE, includeYlab=TRUE, compareTrack=compareTrack,
             labelCols=c("area.rank", "mean.rank",  "qval"))
    dev.off()
  }

if (nrow(area.high)>0){
  compareTrack <- GenomicRangesList(
    dmrseq=makeGRangesFromDataFrame(area.high, keep.extra.columns=TRUE),
    BSmooth=makeGRangesFromDataFrame(bsmooth.common[findOverlaps(makeGRangesFromDataFrame(bsmooth.common),
                                                                 makeGRangesFromDataFrame(area.high))@from,],
                                     keep.extra.columns=TRUE),
    DSS=makeGRangesFromDataFrame(dss.common[findOverlaps(makeGRangesFromDataFrame(dss.common), 
                                                         makeGRangesFromDataFrame(area.high))@from,],
                                 keep.extra.columns=TRUE),
    Metilene=makeGRangesFromDataFrame(met.common[findOverlaps(makeGRangesFromDataFrame(met.common), 
                                                              makeGRangesFromDataFrame(area.high))@from,],
                                      keep.extra.columns=TRUE))
  
  pdf(paste0(figure.file.prefix, "/rankPlots/",
             "DiscordantRanks_", "Area_", cond, ".pdf"), width=6, height=3)
  plotDMRs(bs, regions=area.high, 
           extend=(area.high$end - area.high$start + 1)/2, 
           addRegions=area.high, regionCol=NULL,
           addPoints=TRUE, pointsMinCov=1,qval=FALSE,
           stat=FALSE, includeYlab=TRUE, compareTrack=compareTrack,
           labelCols=c("area.rank", "mean.rank",  "qval"))
  dev.off()
}
}

############################################################################
############################################################################
################# do the same thing for the mouse experiment (DNMT3A)
############################################################################
############################################################################

sampleSize <- 2
num.dmrs <- 0
num.to.plot <- 1  #(how many top regions are plotted in the compPlots)
min.length <- max.length <- NULL
pval.thresh <- 0.10
allConditions <- c("KO_FLT3.WT_WTFL",
                   "WT_FLT3.WT_WTFL",
                   "KO_FLT3.WT_FLT3")

ct <- 1
testCovariate <- "Cond"
adjustCovariate=NULL
genomeName="mm10"
time <- "dnmt3a"

# location of bsseq data objects for each condition
annot <- getAnnot(genomeName) 

for(cond in allConditions){
  
  tiss1 <- strsplit(cond, "\\.")[[1]][1]
  tiss2 <- strsplit(cond, "\\.")[[1]][2]
  cond <- gsub("\\.", "", cond)
  
  load(paste0(dn.data.dir, time, "_all_bsseq.RData"))
  
  print(pData(bs))
  
  if(is.null(sampleNames(bs))){
    sampleNames(bs) <- pData(bs)$SampleNames
  }
  
  meth.levels.raw = getMeth(bs, type = "raw")
  print(nrow(meth.levels.raw))
  nlocs <- nrow(meth.levels.raw)
  tissue1 = grep(tiss1, pData(bs)$Cond) 
  tissue2 = grep(tiss2, pData(bs)$Cond)
  if (length(tissue1) < sampleSize | length(tissue2) < sampleSize){
    stop("Error: Fewer samples present than specified in 'sampleSize'.  Please specify a valid number of samples to use in each group.")
  }
  tissue1 <- tissue1[1:sampleSize]
  tissue2 <- tissue2[1:sampleSize]
  tissue = as.factor(c(rep(tiss1, length(tissue1)), rep(tiss2, length(tissue2))))
  design = model.matrix(~tissue)
  
  ##Subsetting by the number of CpG sites with a positive coverage 
  ##Generating data for first 100000 such CpG Locations
  no.hits = which(is.na(rowMeans(meth.levels.raw[, c(tissue1, tissue2)])) == TRUE)
  str(no.hits)
  nlocs = min(nlocs, (length(bs) - length(no.hits))) # nlocs orig specified in commandarg
  bs = bs[-no.hits][1:nlocs,c(tissue1, tissue2)]
  
  show(bs)
  rm(meth.levels.raw)
  gc()
  
  
  meth.mat = getCoverage(bs, type = "M")
  unmeth.mat = getCoverage(bs, type = "Cov") - meth.mat
  chr = as.character(seqnames(bs))
  pos = start(bs)
  colnames(meth.mat) <- colnames(unmeth.mat) <- pData(bs)$SampleNames
  chr <- paste0("chr", chr, sep="")
  bs <- BSseq(chr = chr, pos = pos,
              M = meth.mat, Cov = meth.mat+unmeth.mat, sampleNames = sampleNames(bs),
              pData=pData(bs))
  rm(meth.mat)
  rm(unmeth.mat)
  rm(chr)
  rm(pos)
  gc();
  
  
  # load DMR objects for each method
  load(paste0(dn.result.dir, "/dmrseq_pkg/regions_", cond, "_", 
              sampleSize, "_", num.dmrs, "DMRs.RData"))
  
  dmrseq <- regions[regions$qval < pval.thresh & !is.na(regions$qval),]
  dmrseq.all <- regions
  rm(regions)
  
  load(paste0(dn.result.dir, "/BSmooth_default/",
              "dmrs.BSmooth.n", sampleSize, ".", cond, ".DEFAULT.Rdata"))
  rm(bs.tstat)
  bsmooth <- dmrs
  bsmooth$stat <- bsmooth$areaStat #/bsmooth$n
  bsmooth$beta <- bsmooth$meanDiff
  bsmooth$indexStart <- bsmooth$idxStart
  bsmooth$indexEnd <- bsmooth$idxEnd
  bsmooth$L <- bsmooth$n
  
  load(paste0(dn.result.dir, "/RESULTS/DSS_default/",
              "dmrs.DSS.n", sampleSize, ".", cond, ".DEFAULT.Rdata"))
  rm(dmlTest.sm)
  dss <- dmrs
  dss$stat <- dss$areaStat #/dss$nCG
  dss$beta <- dss$diff.Methy 
  dss$indexStart <- dss$start
  dss$indexEnd <- dss$end
  dss$L <- dss$nCG
  rm(dmrs)
  
  met <- read.table(paste0(dn.result.dir, "/metilene_default/",
                           "metilene_output_", gsub(" ", "_", cond), "_n", sampleSize, "_0DMRs.txt"), stringsAsFactors=FALSE)
  colnames(met) <- c("chr", "start", "end", "qval", "beta", "L", "pval1", 
                     "pval2", "mean1", "mean2")
  met$stat <- met$beta * met$L
  met$indexStart <- match(paste0(met$chr, "_", met$start+1), 
                          paste0(as.character(seqnames(bs)), "_", start(bs)))
  met$indexEnd <- match(paste0(met$chr, "_", met$end), 
                        paste0(as.character(seqnames(bs)), "_", start(bs)))
  met.all <- met
  met <- met[met$qval < pval.thresh,]
  
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
  
  # construct the design matrix using the pData of bs
  if (ncol(pData(bs)) < max(testCovariate, adjustCovariate)){
    stop(paste0("Error: pData(bs) has too few columns.  Please specify valid ",
                "covariates to use in the analysis"))
  }
  
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
    colnames(design)[,(max(coeff)+1):ncol(design)] <- 
      colnames(pData(bs))[adjustCovariate]
  }else{
    design <- model.matrix( ~ testCov)
    colnames(design)[coeff] <- colnames(pData(bs))[testCovariate]
  }
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  # set up colors and sample labels for plotting 
  # one unique color per unique value of the covariate of interest
  cov.unique <- unique(design[,coeff])
  colors <- gg_color_hue(length(cov.unique))
  if (length(cov.unique) == 2){
    colors <- c("mediumblue", "deeppink1")
  }
  colors <- cbind(cov.unique, colors[rank(as.numeric(cov.unique))])
  z <- colors[,2][match(design[,coeff], colors[,1])]
  pData(bs)$col <- as.character(z)
  pData(bs)$label <- paste0(pData(bs)[,testCovariate])
  pData(bs)$label[pData(bs)$label=="KO_FLT3"] <- "AML"
  pData(bs)$label[pData(bs)$label=="WT_FLT3"] <- "ALL"
  pData(bs)$label[pData(bs)$label=="WT_WTFL"] <- "Control"
  
  ## examine overlaps
  o.bsmooth <- findOverlaps(makeGRangesFromDataFrame(dmrseq.all), 
                            makeGRangesFromDataFrame(bsmooth))
  o.dss <- findOverlaps(makeGRangesFromDataFrame(dmrseq.all), 
                        makeGRangesFromDataFrame(dss))
  o.met <- findOverlaps(makeGRangesFromDataFrame(dmrseq.all), 
                        makeGRangesFromDataFrame(met.all))
  i.all <- unique(o.bsmooth@from[(o.bsmooth@from %in% o.dss@from) & 
                                   (o.bsmooth@from %in% o.met@from)])
  dmrseq.common <- dmrseq.all[i.all,]
  dmrseq.common$dmrseqRank <- rank(-abs(dmrseq.common$stat))/nrow(dmrseq.common)
  
  o.bsmooth <- findOverlaps(makeGRangesFromDataFrame(dmrseq.common), 
                            makeGRangesFromDataFrame(bsmooth))
  o.dss <- findOverlaps(makeGRangesFromDataFrame(dmrseq.common), 
                        makeGRangesFromDataFrame(dss))
  o.met <- findOverlaps(makeGRangesFromDataFrame(dmrseq.common), 
                        makeGRangesFromDataFrame(met.all))
  bsmooth.common <- bsmooth[o.bsmooth@to,]
  bsmooth.common$i.dmrseq <- o.bsmooth@from
  dss.common <- dss[o.dss@to,]
  dss.common$i.dmrseq <- o.dss@from
  met.common <- met.all[o.met@to,]
  met.common$i.dmrseq <- o.met@from
  
  bsmooth.common$area.rank <- rank(-abs(bsmooth.common$areaStat))/nrow(bsmooth.common)
  bsmooth.common$mean.rank <- rank(-abs(bsmooth.common$meanDiff))/nrow(bsmooth.common)
  dss.common$area.rank <- rank(-abs(dss.common$areaStat))/nrow(dss.common)
  dss.common$mean.rank <- rank(-abs(dss.common$diff.Methy))/nrow(dss.common)
  met.common$rank.qval <- rank(met.common$qval)/nrow(met.common)
  
  # add these ranks to common region table
  # average over all overlapping regions if there are more than 1
  dmrseq.common <- cbind(dmrseq.common,
                         (bsmooth.common %>% group_by(i.dmrseq) %>% 
                            summarize(BSmooth.areaRank=mean(area.rank),
                                      BSmooth.avgRank=mean(mean.rank)))[,-1],
                         (dss.common %>% group_by(i.dmrseq) %>% 
                            summarize(DSS.areaRank=mean(area.rank),
                                      DSS.avgRank=mean(mean.rank)))[,-1],
                         (met.common %>% group_by(i.dmrseq) %>% 
                            summarize(Metilene.qvalRank=mean(rank.qval)))[,-1])
  
  col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                                 "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))
  ranks <- dmrseq.common[,11:16]
  
  # find an enriched set of regions that are given high ranking by other 
  # method(s) but have high(er) FDR by dmrseq
  bsmooth.high <- dmrseq.common[(dmrseq.common$BSmooth.areaRank < 0.1 |
                                   dmrseq.common$BSmooth.avgRank < 0.1) &
                                  dmrseq.common$qval > 0.5, ]
  dss.high <- dmrseq.common[(dmrseq.common$DSS.areaRank < 0.1 |
                               dmrseq.common$DSS.avgRank < 0.1) &
                              dmrseq.common$qval > 0.50, ]
  met.high <- dmrseq.common[dmrseq.common$Metilene.qvalRank < 0.05 &
                              dmrseq.common$qval > 0.50, ]
  all.high <- dmrseq.common[dmrseq.common$BSmooth.areaRank < 0.1 &
                              dmrseq.common$DSS.areaRank < 0.1 &
                              dmrseq.common$Metilene.qvalRank < 0.1 &
                              dmrseq.common$qval < 0.1, ]
  avg.high <- dmrseq.common[(dmrseq.common$BSmooth.avgRank < 0.2 &
                               dmrseq.common$DSS.avgRank < 0.2 & 
                               dmrseq.common$qval > 0.5), ]
  area.high <- dmrseq.common[(dmrseq.common$BSmooth.areaRank < 0.1 &
                                dmrseq.common$DSS.areaRank < 0.1 & 
                                dmrseq.common$qval > 0.5), ]
  
  # plot these guys
  if (nrow(bsmooth.high)>0){
    compareTrack <- GenomicRangesList(
      dmrseq=makeGRangesFromDataFrame(bsmooth.high, keep.extra.columns=TRUE),
      BSmooth=makeGRangesFromDataFrame(bsmooth.common[findOverlaps(makeGRangesFromDataFrame(bsmooth.common),
                                                                   makeGRangesFromDataFrame(bsmooth.high))@from,],
                                       keep.extra.columns=TRUE),
      DSS=makeGRangesFromDataFrame(dss.common[findOverlaps(makeGRangesFromDataFrame(dss.common), 
                                                           makeGRangesFromDataFrame(bsmooth.high))@from,],
                                   keep.extra.columns=TRUE),
      Metilene=makeGRangesFromDataFrame(met.common[findOverlaps(makeGRangesFromDataFrame(met.common), 
                                                                makeGRangesFromDataFrame(bsmooth.high))@from,],
                                        keep.extra.columns=TRUE))
    
    pdf(paste0(figure.file.prefix, "/rankPlots/",
               "DiscordantRanks_BSmooth_", cond, ".pdf"), width=6, height=3)
    plotDMRs(bs, regions=bsmooth.high, 
             extend=(bsmooth.high$end - bsmooth.high$start + 1)/2, 
             addRegions=bsmooth.high, regionCol=NULL,
             addPoints=TRUE, pointsMinCov=1,qval=FALSE,
             stat=FALSE, includeYlab=TRUE, compareTrack=compareTrack,
             labelCols=c("area.rank", "mean.rank",  "qval"))
    dev.off()
  }
  
  if (nrow(dss.high)>0){
    compareTrack <- GenomicRangesList(
      dmrseq=makeGRangesFromDataFrame(dss.high, keep.extra.columns=TRUE),
      BSmooth=makeGRangesFromDataFrame(bsmooth.common[findOverlaps(makeGRangesFromDataFrame(bsmooth.common),
                                                                   makeGRangesFromDataFrame(dss.high))@from,],
                                       keep.extra.columns=TRUE),
      DSS=makeGRangesFromDataFrame(dss.common[findOverlaps(makeGRangesFromDataFrame(dss.common), 
                                                           makeGRangesFromDataFrame(dss.high))@from,],
                                   keep.extra.columns=TRUE),
      Metilene=makeGRangesFromDataFrame(met.common[findOverlaps(makeGRangesFromDataFrame(met.common), 
                                                                makeGRangesFromDataFrame(dss.high))@from,],
                                        keep.extra.columns=TRUE))
    
    pdf(paste0(figure.file.prefix, "/rankPlots/",
               "DiscordantRanks_DSS_", cond, ".pdf"), width=6, height=3)
    plotDMRs(bs, regions=dss.high, 
             extend=(dss.high$end - dss.high$start + 1)/2, 
             addRegions=dss.high, regionCol=NULL,
             addPoints=TRUE, pointsMinCov=1,qval=FALSE,
             stat=FALSE, includeYlab=TRUE, compareTrack=compareTrack,
             labelCols=c("area.rank", "mean.rank",  "qval"))
    dev.off()
  }
  
  if (nrow(met.high)>0){
    compareTrack <- GenomicRangesList(
      dmrseq=makeGRangesFromDataFrame(met.high, keep.extra.columns=TRUE),
      BSmooth=makeGRangesFromDataFrame(bsmooth.common[findOverlaps(makeGRangesFromDataFrame(bsmooth.common),
                                                                   makeGRangesFromDataFrame(met.high))@from,],
                                       keep.extra.columns=TRUE),
      DSS=makeGRangesFromDataFrame(dss.common[findOverlaps(makeGRangesFromDataFrame(dss.common), 
                                                           makeGRangesFromDataFrame(met.high))@from,],
                                   keep.extra.columns=TRUE),
      Metilene=makeGRangesFromDataFrame(met.common[findOverlaps(makeGRangesFromDataFrame(met.common), 
                                                                makeGRangesFromDataFrame(met.high))@from,],
                                        keep.extra.columns=TRUE))
    
    pdf(paste0(figure.file.prefix, "/rankPlots/",
               "DiscordantRanks_Metilene_", cond, ".pdf"), width=6, height=3)
    plotDMRs(bs, regions=met.high, 
             extend=(met.high$end - met.high$start + 1)/2, 
             addRegions=met.high, regionCol=NULL,
             addPoints=TRUE, pointsMinCov=1,qval=FALSE,
             stat=FALSE, includeYlab=TRUE, compareTrack=compareTrack,
             labelCols=c("area.rank", "mean.rank",  "qval"))
    dev.off()
  }
  
  if (nrow(all.high)>0){
    all.high <- all.high[1:min(nrow(all.high), 25),]
    compareTrack <- GenomicRangesList(
      dmrseq=makeGRangesFromDataFrame(all.high, keep.extra.columns=TRUE),
      BSmooth=makeGRangesFromDataFrame(bsmooth.common[findOverlaps(makeGRangesFromDataFrame(bsmooth.common),
                                                                   makeGRangesFromDataFrame(all.high))@from,],
                                       keep.extra.columns=TRUE),
      DSS=makeGRangesFromDataFrame(dss.common[findOverlaps(makeGRangesFromDataFrame(dss.common), 
                                                           makeGRangesFromDataFrame(all.high))@from,],
                                   keep.extra.columns=TRUE),
      Metilene=makeGRangesFromDataFrame(met.common[findOverlaps(makeGRangesFromDataFrame(met.common), 
                                                                makeGRangesFromDataFrame(all.high))@from,],
                                        keep.extra.columns=TRUE))
    
    pdf(paste0(figure.file.prefix, "/rankPlots/",
               "ConcordantRanks_", cond, ".pdf"), width=6, height=3)
    plotDMRs(bs, regions=all.high, 
             extend=(all.high$end - all.high$start + 1)/2, 
             addRegions=all.high, regionCol=NULL,
             addPoints=TRUE, pointsMinCov=1,qval=FALSE,
             stat=FALSE, includeYlab=TRUE, compareTrack=compareTrack,
             labelCols=c("area.rank", "mean.rank",  "qval"))
    dev.off()
  }
  
  if (nrow(avg.high)>0){
    compareTrack <- GenomicRangesList(
      dmrseq=makeGRangesFromDataFrame(avg.high, keep.extra.columns=TRUE),
      BSmooth=makeGRangesFromDataFrame(bsmooth.common[findOverlaps(makeGRangesFromDataFrame(bsmooth.common),
                                                                   makeGRangesFromDataFrame(avg.high))@from,],
                                       keep.extra.columns=TRUE),
      DSS=makeGRangesFromDataFrame(dss.common[findOverlaps(makeGRangesFromDataFrame(dss.common), 
                                                           makeGRangesFromDataFrame(avg.high))@from,],
                                   keep.extra.columns=TRUE),
      Metilene=makeGRangesFromDataFrame(met.common[findOverlaps(makeGRangesFromDataFrame(met.common), 
                                                                makeGRangesFromDataFrame(avg.high))@from,],
                                        keep.extra.columns=TRUE))
    
    pdf(paste0(figure.file.prefix, "/rankPlots/",
               "DiscordantRanks_", "Avg_", cond, ".pdf"), width=6, height=3)
    plotDMRs(bs, regions=avg.high, 
             extend=(avg.high$end - avg.high$start + 1)/2, 
             addRegions=avg.high, regionCol=NULL,
             addPoints=TRUE, pointsMinCov=1,qval=FALSE,
             stat=FALSE, includeYlab=TRUE, compareTrack=compareTrack,
             labelCols=c("area.rank", "mean.rank",  "qval"))
    dev.off()
  }
  
  if (nrow(area.high)>0){
    compareTrack <- GenomicRangesList(
      dmrseq=makeGRangesFromDataFrame(area.high, keep.extra.columns=TRUE),
      BSmooth=makeGRangesFromDataFrame(bsmooth.common[findOverlaps(makeGRangesFromDataFrame(bsmooth.common),
                                                                   makeGRangesFromDataFrame(area.high))@from,],
                                       keep.extra.columns=TRUE),
      DSS=makeGRangesFromDataFrame(dss.common[findOverlaps(makeGRangesFromDataFrame(dss.common), 
                                                           makeGRangesFromDataFrame(area.high))@from,],
                                   keep.extra.columns=TRUE),
      Metilene=makeGRangesFromDataFrame(met.common[findOverlaps(makeGRangesFromDataFrame(met.common), 
                                                                makeGRangesFromDataFrame(area.high))@from,],
                                        keep.extra.columns=TRUE))
    
    pdf(paste0(figure.file.prefix, "/rankPlots/",
               "DiscordantRanks_", "Area_", cond, ".pdf"), width=6, height=3)
    plotDMRs(bs, regions=area.high, 
             extend=(area.high$end - area.high$start + 1)/2, 
             addRegions=area.high, regionCol=NULL,
             addPoints=TRUE, pointsMinCov=1,qval=FALSE,
             stat=FALSE, includeYlab=TRUE, compareTrack=compareTrack,
             labelCols=c("area.rank", "mean.rank",  "qval"))
    dev.off()
  }
}


