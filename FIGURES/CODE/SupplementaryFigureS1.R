# Code to generate Supplementary Figure S1
# To compare AR(1) with CAR(1)

# Simulation D2 (because then can evaluate sens & spec)
setwd("/n/irizarryfs01_backed_up/kkorthauer/WGBS/dmrseqPaper/FIGURES/out/")

# load necessary packages and region stat function (internal to dmrseq)
library(dmrseq)
library(nlme)
library(reshape2)
library(outliers)

asin.gls.cov <- function(ix, design, coeff,
                           correlation=corAR1(form = ~ 1 | sample),
                           correlationSmall=corCAR1(form = ~ L | sample),
                           weights=varPower(form=~1/MedCov, fixed=0.5)){
    sampleSize <- nrow(design)/2
    dat <- data.frame(
      g.fac=factor(as.vector(sapply(design[,coeff], 
                                    function(x) rep(x,length(ix))))),
      sample=factor(as.vector(sapply(1:(sampleSize*2), 
                                     function(x) rep(x,length(ix))))),
      meth=melt(meth.mat[ix,])$value,
      cov=melt(cov.mat[ix,])$value,
      L = as.vector(rep(pos[ix], nrow(design)))
    )
    
    # condition to remove regions with constant methylation / unmeth values
    if ( ! ((length(unique(dat$meth)) == 1 & dat$meth[1] == 0) |
            (length(unique(dat$cov-dat$meth)) == 1 & 
             (dat$cov-dat$meth)[1] == 0)) ){ 
      
      dat$pos <- as.numeric(factor(dat$L))   				
      X <- model.matrix( ~ dat$g.fac )
      colnames(X)[2] <- "grp"
      
      if(nrow(X) <= ncol(X))
        stop(paste0("Not enough degree of freedom to fit the linear ",
                    "model. Drop some terms in formula"))
      
      Y <- as.matrix(dat$meth)
      N <- as.matrix(dat$cov)
      
      dat$MedCov <- rep(as.numeric(by(dat$cov, dat$pos, median)), sampleSize*2)
      
      # tol is the convergence criteria (for difference Between two
      # iterations of the dispersion parameter estimate
      
      ## small constants to bound p and phi
      # pick this such that min value of z[m>0] 
      # is greater than all values of z[m==0]
      c0 = 0.05
      c1 = 0.001
      
      ## check to make sure data is complete
      ixn <- N > 0
      if(mean(ixn) < 1) { ## has missing entries
        X <- X[ixn,,drop=FALSE]
        Y <- Y[ixn]
        N <- N[ixn]
        ## check design
        if(nrow(X) < ncol(X) + 1) ## not enough df for regression
          return(NULL)
        ## design is not of full rank because of missing. Skip
        if(any(abs(svd(X)$d) <1e-8)) 
          return(NULL)
      }
      
      ## Transform the methylation levels. 
      #Add a small constant to bound away from 0/1.
      dat$Z = asin(2*(Y+c0)/(N+2*c0) - 1) 
      
      #Add a tiny amt of jitter to avoid numerically constant Z vals
      # across a sample over the entire region
      if (max(table(dat$Z, dat$sample)) >= length(ix)-1){
        dat$Z = asin(2*(Y+c0)/(N+2*c0) - 1) + 
          runif(length(Y), -c1, c1)
      }
      
      
      if (length(ix) >= 40){
        fit <- tryCatch({summary(gls(Z ~ g.fac + factor(L), 
                                     weights=weights,
                                     data=dat, 
                                     correlation=correlation))},
                        error=function(e) { return(NA)})
      }else{
        # check for presence of 1-2 coverage outliers that could end up driving
        # the difference between the groups
        if (length(unique(dat$MedCov[1:length(ix)])) > 1 & length(ix) <= 10 ){
          grubbs.one <- suppressWarnings(
              grubbs.test(dat$MedCov[1:length(ix)])$p.value)
          grubbs.two <- suppressWarnings(
              grubbs.test(dat$MedCov[1:length(ix)], type=20)$p.value)
        }else{
          grubbs.one <- grubbs.two <- 1
        }
        
        if (grubbs.one < 0.01 | grubbs.two < 0.01){
          weights = varIdent(form = ~ 1)
        }
        
        fit <- tryCatch({summary(gls(Z ~ g.fac + factor(L),
                                     weights=weights,
                                     data=dat,
                                     correlation=correlationSmall))},
                        error=function(e) { return(NA)})
        
        # error handling in case of false convergence (don't
        # include first variance weighting, and then corr str)
        if (sum(is.na(fit)) == length(fit)){
          fit <- tryCatch({summary(gls(Z ~ g.fac + factor(L),
                                       data=dat,
                                       correlation=correlationSmall))},
                          error=function(e) { return(NA)})
          if (sum(is.na(fit)) == length(fit)){
            fit <- tryCatch({summary(gls(Z ~ g.fac + factor(L),
                                         data=dat))},
                            error=function(e) { return(NA)})
          }
        }
      }
      
      if (!(sum(is.na(fit)) == length(fit))){
        stat <- fit$tTable[2,3]
        beta <- fit$tTable[2,1]
      }else{
        stat <- beta <- NA
      }
      
      return(data.frame(beta=beta,stat=stat))
    }else{
      return(NULL)
    }
}

data.file.prefix <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/binSim/" 

result.file.prefix <- "/n/irizarryfs01_backed_up/kkorthauer/WGBS/DENDRITIC/RESULTS/dmrseq_pkg"

sampleSize=2 
num.dmrs=3000 
cond="control"

workers <- as.numeric(Sys.getenv("SLURM_NTASKS"))
# register cores with #workers
library("BiocParallel")
register(MulticoreParam(workers))

time <- "dendritic"
time2 <- NULL
coeff <- 2
# load in simulated data
sim.file <- paste0(data.file.prefix, "/sim.data.n", sampleSize, 
 					  ".", cond, ".all.rda")
if (num.dmrs>0){  
 if(!file.exists(sim.file)){
 	simDMRs(sim.file, bs, num.dmrs)
 }
 load(sim.file)
 bs <- sim.dat.red$bs
 design <- model.matrix(~c(0,0,1,1))
 pData(bs)$Infected <- design[,2]
 dmrs.true = sim.dat.red$gr.dmrs
 rm(sim.dat.red)
}
 meth.mat <- as.matrix(getCoverage(bs, type="M"))
 cov.mat <- as.matrix(getCoverage(bs, type="Cov"))
 pos <- start(bs)
 chr <- as.character(seqnames(bs))
 meta <- pData(bs)
 rm(bs)

# load in candidate regions
load(paste0(result.file.prefix, "/regions_", cond, "_", 
                 sampleSize, "_", num.dmrs, "DMRs.RData"))
 
# construct Indices of candidate regions
ind <- 1:length(pos)
Indexes <- vector("list", nrow(regions))
for (r in 1:length(Indexes)){
	Indexes[[r]] <- regions$indexStart[r]:regions$indexEnd[r]
}   

if (!file.exists("corstructure_comparison.rda")){
  ar1 <- do.call("rbind", bplapply(Indexes, function(Index)
		asin.gls.cov(ix=ind[Index],design=design,coeff=coeff,
					 correlation=corAR1(form = ~ 1 | sample),
					 correlationSmall=corAR1(form = ~ 1 | sample))))

  car1 <- do.call("rbind", bplapply(Indexes, function(Index)
		asin.gls.cov(ix=ind[Index],design=design,coeff=coeff,
					 correlation=corCAR1(form = ~ L | sample),
					 correlationSmall=corCAR1(form = ~ L | sample))))
						  
  save(ar1, car1, file="corstructure_comparison.rda")
}else{
  load("corstructure_comparison.rda")
}

library(ggplot2)
library(cowplot)
regions$stat.AR1 <- ar1$stat
regions$stat.CAR1 <- car1$stat
regions$stat.diff <- regions$stat.AR1-regions$stat.CAR1
regions$width <- regions$end-regions$start + 1

pdf("supp_fig1.pdf")
p1 <- ggplot(regions, aes(x=L, y=stat.diff)) +
    xlim(c(5,100)) +
    ylim(c(-10,10)) +
    geom_point(shape=1) +    # Use hollow circles
    geom_smooth() +          # Add a loess smoothed fit curve with confidence region
	xlab("Number of CpG loci in the region") +
	ylab("AR(1) statistic - CAR(1) statistic") +
	ggtitle("") +
	theme_classic()
p2 <- ggplot(regions, aes(x=width, y=stat.diff)) +
    ylim(c(-10,10)) +
    geom_point(shape=1) +    # Use hollow circles
    geom_smooth() +          # Add a loess smoothed fit curve with confidence region
     ggtitle("") +
	xlab("Number of basepairs spanning the region") +
	ylab("AR(1) statistic - CAR(1) statistic") +
	theme_classic()
	
  print(plot_grid(p1, p2, nrow=2, hjust=0, 
  	    labels=c("(A) Statistics under AR(1) vs CAR(1) by number of CpG loci",
  	    		 "(B) Statistics under AR(1) vs CAR(1) by region width")))
dev.off()
          	          
	          	          
