######################
# SUMMARY STATISTICS #
######################

# Calculate Watterson theta for one replicate
# sapply(stats_neut$seg_sites, calculate.thetaW) to compute for all replicates
calculate.thetaW <- function(segsites){
        S <- ncol(segsites)
        n <- nrow(segsites)
        thetaW <- S/sum(1/(1:(n-1)))
        return(thetaW)
}


# Calculate PI ie Tajima's theta for one replicate
# sapply(stats_neut$seg_sites, calculate.thetaPi) to compute for all replicates
calculate.thetaPi <- function(segsites){
  S <- ncol(segsites)
  Nsam <- nrow(segsites)
  SFS <- table(apply(segsites,2,function(x){sum(as.numeric(x))}))
  DiffMoy <- as.numeric(unlist(lapply(as.numeric(names(SFS)),function(x){(x*(Nsam-x))/(Nsam*(Nsam-1)/2)})))
  thetaPi <- sum(DiffMoy*SFS)
return(thetaPi)
}


# Calculate Tajima's D
# the vector sapply(stats_neut$seg_sites, calculate.Dt) is equivalent to stats_neut$tajimas_d
calculate.Dt<-function(segsites){
  S <- ncol(segsites)
  n <- nrow(segsites)
  a1 = sum(1/(1:(n-1)))
  a2 = sum(1/(1:(n-1))^2)
  b1 = (n+1)/(3*(n-1))
  b2 = 2*(n^2+n+3)/(9*n*(n-1))
  c1 = b1-1/a1
  c2 = b2-(n+2)/(a1*n)+a2/(a1^2)
  e1 = c1/a1
  e2 = c2/(a1^2+a2)
  VarDt <- sqrt(e1*S+e2*S*(S-1))
  Dt <- (calculate.thetaPi(segsites) - calculate.thetaW(segsites)) / VarDt
  return(Dt)
}

# Calculate the variance of Tajima's D 
# (1) cut the segsite matrix into nwin consecutive windows with equal number of snps
# (2) compute Dt for each window 
# (3) return the variance 
# setting nwin to 10000 for the realdata allow to compute the variance of Dt for
# short loci as it was done for the simulations (in the simu locus length ~
# 2000bp)
calculate.Dt.var = function(segsites, nwin=1e4) {
  nsnps = ncol(segsites$snps)
  step = floor(nsnps/nwin)
  taj_d_win = sapply(0:(nwin-1), FUN=function(i) calculate.Dt(segsites$snps[,(i*step+1):((i+1)*step)]))
  return(var(taj_d_win))
}

# Draw the site frequency spectrum 
# sapply(stats_neut$seg_sites, draw.sfs) will draw the SFS for each replicate 
# (not very useful because it's too stochastic)
# better to compute the SFS based on ALL replicates
draw.sfs <- function(segsites){
    plot(table(apply(segsites,2,function(x){sum(as.numeric(x))})),type="h",xlab="Derived Allele Frequency",ylab="Count")
}

# Represent the coalescent tree
# sapply(stats_neut$seg_sites, plot.coalesc.tree) will draw the tree for each replicate 
library(ape)
plot.coalesc.tree <- function(segsites, ...){
  segsites <- as.matrix(segsites)
  Genotype <- rbind(segsites,rep(0,times=ncol(segsites))) # Add a line filled with zero for mimicking an outgroup with ancestral alleles
  #Genotype <- Genotype[,colSums(Genotype)>1 & colSums(Genotype)<(max(colSums(Genotype))-1)] #what was this line for?
  print(dim(Genotype))
  if (is.null(rownames(segsites))) {
    rownames(segsites) <- 1:nrow(segsites)
  }
  rownames(Genotype) <- c(rownames(segsites),"outgroup")
  Dist <- dist.gene(Genotype)
  Tree <- nj(Dist)
  Tree <- root(Tree,outgroup="outgroup")
  plot(Tree,...)
}


# calculate thetaL
# sapply(stats_neut$seg_sites, calculate.thetaL)
calculate.thetaL <- function(segsites){
    segsites <- as.matrix(segsites)
    n <- nrow(segsites)
    Epsilon <- colSums(segsites)
    thetaL <- table(Epsilon)
    thetaL <- sum(as.numeric(names(thetaL))*thetaL)/(n-1)
    return(thetaL)
}

#calculate Fay Wu H modified by Zeng
calculate.Hfw <- function(segsites){
    n <- nrow(segsites)
    S <- ncol(segsites)
    thetaL <- calculate.thetaL(segsites)
    thetaPi <- calculate.thetaPi(segsites)
    thetaW <- calculate.thetaW(segsites)
    an = sum(1/(1:(n-1)))
    bn = sum(1/((1:(n-1))^2))
    theta2 <- S*(S-1)/(an^2+bn)
    denom = (9*n*(n-1)*(n-1))
    VarTLW <- thetaW*(n-2)/(6*(n-1)) + theta2*(18*n*n*(3*n+2)*sum(1/(1:n)^2) - (88*n*n*n+9*n*n-13*n+6)) / denom
    return( (thetaL-thetaW) / VarTLW )
}


#' Calculate heterozygosity for each diploid individual 
#' @title Het
#' @param segsites segsites object (or a matrix representing snps)
#' @param sampleinfo info about the samples, i.e a matrix having the same number of lines than segsites
#' @return list(het=het, info=info) het = one value per diploid individual, info= code/population/region of ech diploid individual
#' @note diploid individuals haplotypes are adjacent in the table 
#' (eg line 1 = paternal chromosome of the first individual and line 2 = maternal chromosome of the first individual)
#' 
individual.het <-function(segsites, sampleinfo=NULL) {
  # 2 successive lines are considered as the 2 haplotypes of a same individual
  # output vector sizes = the nb of indivuals = nrow(segsites)/2
  segsites <- as.matrix(segsites)
  if (nrow(segsites)%%2!=0) 
    stop("The snp matrix should have an even number of lines to be able to create diploid individuals from 2 successive haplotypes")
  odd <- seq(1,nrow(segsites),by=2)
  het <- sapply(odd, FUN=function(i) mean(segsites[i,]!=segsites[i+1,]))
  info <- NULL
  if (!is.null(sampleinfo)) info <- sampleinfo[odd,]
  if (!is.null(rownames(segsites))) names(het) <- rownames(segsites)[odd]
  return(list(het=het, info=info))
}



#########################
# HANDLING REAL DATASET #
#########################


#' Create an object msobj from the human population Complete Genomics dataset
#' @title Read haplotype file and create msobj
#' @note additionnally to the attributes of a classic msobj created from ms simulation outputs
#' this sobj contain the absolute position along the cromosome (in base pairs)
#' and information about each haploid sample (individual ID, POP, REGION)
#' @param hap.input Name of file containing the genetic data
#' @param info.input Name of file containing the sampling information/details
#' @param skip number of lines to skip in the hap.input file
#' @param nrows total number of lines to read from hap.input file
#' @param ... all additional arguments are passed to the read.table funtion (you should not need that)
#' 
#' @return data  with attributes
#' \item{segsites} obsject of class segsites (from coala library) 
#' that has a snps and a position attributes
#'   - snps: contains binary genotypes (0=ancestral allele, 1=derived), rows=haploid individuals and columns=SNPs
#'   - position: relative position of the SNPs (between 0 and 1)
#'\item{region_bp_lim} start and end of the region in absolute base pairs (to transform back the relative snp positions to bp)
#' \item{info} matrix nindiv x 3 containing information about each haploid sample:
#'             individual code, population and region (ID, POP, REGION)
#'             
#'@example 
#' realdata <- read.hap.geno("subset.chr2.CG_54genomes_shapeit_phased.haps", info.input = "CG_54genomes_indiv.txt")
#'             
read.hap.geno <- function(hap.input, info.input=NULL, skip=0, nrows=-1, ...) {
  snps = read.table(hap.input, skip=skip, nrows=nrows, ...)
  absolute_pos = snps[,3]
  snps = t(as.matrix(snps[,-(1:5)])) 
  bp_start = min(absolute_pos)
  bp_end = max(absolute_pos)
  pos = (absolute_pos - bp_start)/(bp_end-bp_start)
  segsites = create_segsites(snps, pos, check = TRUE)

  sampleinfo = NULL
  if (!is.null(info.input)) {
    sampleinfo = (read.table(info.input,header=T,stringsAsFactors = F))
    rownames(segsites$snps) = apply(sampleinfo, 1, function(x) paste(x,collapse="_"))
  } 
  return(list(segsites=segsites, region_bp_lim=c(bp_start, bp_end), info=sampleinfo))
}


#' Crop genotype for individuals/populations
#' @title subsample genotype data at individuals selected from population or region label
#' Useful for calculating statistics for a specific population 
#' and for removing sites that are fixed in this population
select.label <- function(data, label){
  datacrop <- data
  selected <- (data$info$POP == label) | (data$info$REGION == label)
  snps <- data$segsites$snps[selected,]
  # Create segsite and filter out sites that are fixed for this subselection (via check=TRUE)
  segsites <- create_segsites(snps,  data$segsites$position, check=TRUE)
  datacrop$segsites <- segsites
  datacrop$info <- data$info[selected,]
  return(datacrop)
}

#' Crop genotype for snp
#' @title subsample genotype at given snps
#' Useful for testing on small subset of data
#' give either a nb of snp to keep, or the list/vector of snp indexes, or a minor allele frequency threshold
select.snps <- function(data, nsnps=NULL, keep=NULL, maf=NULL){
  datacrop <- data
  if (!is.null(maf)) {
    # keep snps with allele frequency higher than maf
    keep = which(colSums(data$segsites)/nrow(data$segsites) >= maf) 
  } else if (is.null(keep)) {
    # keep the nsnps first nsnps
    keep = 1:nsnps 
  }
  datacrop$segsites <- create_segsites(data$segsites$snps[,keep], data$segsites$position[keep])
  return(datacrop)
}


