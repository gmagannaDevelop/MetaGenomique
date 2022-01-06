#' ---
#' title: "ABC"
#' author:
#'   - Gustavo MAGAÑA LOPEZ
#'   - Théo RONCALLI
#' output:
#'   pdf_document:
#'     keep_tex: true
#' ---
#'

#+ echo=F, message=F, warning=F
# declare dependencies
.dependencies <- 
  c("here", "abc", "coala", 
    "magrittr", "ggplot2", "dplyr",
    "tibble", "foreach", "doRNG", "stringr"
    )
new.packages <- 
  .dependencies[!(.dependencies %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)
.. <- sapply(.dependencies, function(x){ library(x, character.only = T)})
# load utility functions
source(here("ABC/Fonctions_TP_demo_real.R")) # this script has slighly changed since TP1
# load data
hap.fname <- here("ABC/chr22.CG_54genomes_shapeit_phased.txt.haps")
info.fname <- here("ABC/CG_54genomes_indiv.txt")
realdata <- read.hap.geno(hap.input=hap.fname, info.input = info.fname)
set.seed(1234)


# sumstats
stat.1000g <- NULL
for (pop in unique(realdata$info$POP)) {
  # subset the dataset
  pop.subset <- select.label(realdata, label = pop)
  # compute estimates
  
  size.norm.factor <- pop.subset$region_bp_lim[2] - pop.subset$region_bp_lim[1]
  theta_watterson <- calculate.thetaW(pop.subset$segsites) / size.norm.factor
  theta_tajima <- calculate.thetaPi(pop.subset$segsites) / size.norm.factor
  tajimas_d <- calculate.Dt(pop.subset$segsites)
  tajimas_d_var <- calculate.Dt.var(pop.subset$segsites)
  # create the dataframe
  .pop.stats <- 
    data.frame(
      Theta.Tajima=theta_tajima, D.Tajima=tajimas_d, 
      Var.D.Tajima=tajimas_d_var,
      Theta.Watterson=theta_watterson
    )
  stat.1000g <- rbind(stat.1000g, .pop.stats)
}
rownames(stat.1000g) <- unique(realdata$info$POP)
stat.1000g.tb <- stat.1000g %>% rownames_to_column("Population") %>% tibble()

loaded <- load(file = here("ABC/simu.expeA.RData"))
loaded
# [1] "sumstats.expeA" "models.expeA"   "params.bott"
# sapply(loaded, function(x){ str(eval(as.symbol(x))) }, simplify = T)


#' Data analysis
dim(sumstats.expeA)
head(sumstats.expeA) # each line is a set of parameters... for a different model ?
table(models.expeA) # equally divided bott, const, exp -> 50_000 each
head(models.expeA)
dim(params.bott)  
head(params.bott)


# Visualise model parameters
{
  par(mfrow=c(1,3))
  boxplot(sumstats.expeA$pi~models.expeA)
  boxplot(sumstats.expeA$tajimas_d~models.expeA)
  boxplot(sumstats.expeA$tajimas_d_var~models.expeA)
  par(mfrow=c(1,1))
}


# Explore main ABC function
?postpr


# Explore
{
  random_idx <- 
    sample.int(length(models.expeA), size = 1)
  random_idx
  
  abc.mod.sel <- postpr(
    target = sumstats.expeA[random_idx,],   # pseudo-observed sumstats
    index = models.expeA[-random_idx],      # known models for reference db
    sumstat = sumstats.expeA[-random_idx,], # sumstats for reference db
    method = 'rejection',                   # type of ABC algorithm
    tol=0.05                                #
  )
  prob <- summary(abc.mod.sel)$Prob
  prob
}

# Gridsearch
  random_idx <- 
    sample.int(length(models.expeA), size = 1)
  random_idx
  
  abc.mod.sel <- postpr(
    target = sumstats.expeA[random_idx,],   # pseudo-observed sumstats
    index = models.expeA[-random_idx],      # known models for reference db
    sumstat = sumstats.expeA[-random_idx,], # sumstats for reference db
    method = 'rejection',                   # type of ABC algorithm
    tol=0.05                                #
  )
  prob <- summary(abc.mod.sel)$Prob
  prob
  

cv4postpr(
  index = models.expeA, 
  sumstat = sumstats.expeA, 
  nval=10, 
  tols=c(.0005,.1,.5), 
  method="rejection"
) -> cv.modsel.reject

parallel.cv.modsel.reject <- 
  parallel_cv4postpr(
    models.expeA, sumstats.expeA, nval=20, tols=seq(0.005, 0.95, length.out = 30)
  )

parallel.cv.modsel.reject2 <- 
  parallel_cv4postpr(
    models.expeA, sumstats.expeA, nval=20, 
    tols=seq(0.005, 0.95, length.out = 30),
    seed = 444
  )

parallel.cv.modsel.reject3 <- 
  parallel_cv4postpr(
    models.expeA, sumstats.expeA, nval=25, tols=seq(0.005, 0.95, length.out = 60)
  )

# plot(cv.modsel) # matrice de confusion pour chaque seuil de toleranc
compute_misclassif_from_cv_model <- function(cv.modsel.reject){
  if (!require(stringr, quietly = T)){ 
    stop("this function requires the package 'stringr'")
  }
  misclassif.reject <- list()
  for (tol in names(cv.modsel.reject$estim)) {
    misclassif.reject[[tol]] = mean(cv.modsel.reject$estim[[tol]] != cv.modsel.reject$true)
  }
  .tols <- as.numeric(str_remove(names(misclassif.reject), "tol"))
  data.frame(tolerance=.tols, perc.misclassified=as.numeric(misclassif.reject))
}

compute_misclassif_from_cv_model(cv.modsel.reject)  
rej <- compute_misclassif_from_cv_model(parallel.cv.modsel.reject)  
rej2 <- compute_misclassif_from_cv_model(parallel.cv.modsel.reject2)  
rej3 <- compute_misclassif_from_cv_model(parallel.cv.modsel.reject3)  
rej.full <- compute_misclassif_from_cv_model(parallel.cv.modsel.reject.full)  

rej.full %>% ggplot(aes(x=tolerance, y=perc.misclassified)) + geom_point() + geom_smooth()
  
# Parallel gridsearch
source(here("ABC/customfuncs.R"))




