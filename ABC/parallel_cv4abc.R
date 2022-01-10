#' Université Paris-Saclay
#' Faculté des Sciences d'Orsay
#' M2 AMI2B : Biologie Computationnelle
#' 
#' Métagénomique et génétique des populations
#' 
#' Gustavo Magaña López
#' Theó Roncalli
#' 

# Load libraries
shhh <- suppressPackageStartupMessages # It's a library, so shhh!
shhh(library(here))
shhh(library(glue))
shhh(library(ggplot2))

# Source our functions
source(here("ABC/cv4abc.R"))
source(here("ABC/customfuncs.R"))

# parse command line arguments
args <- commandArgs(trailingOnly=T)
if (length(args) != 6){
  print("usage: ")
  print("Rscript parallel_cv4abc.R nval tol.low tol.high n.tols threads method")
  stop("Incorrect number of arguments")
} else {
  nval <- as.numeric(args[1])
  tol.low <- as.numeric(args[2])
  tol.high <- as.numeric(args[3])
  n.tols <- as.numeric(args[4])
  threads <- as.numeric(args[5])
  method <- args[6]
  if (threads > n.tols){
    stop("You have more threads than tolerances to simulate!")
  }
}

# load simulated data
loaded <- load(file = here("ABC/simu.expeA.RData"))
# [1] "sumstats.expeA" "models.expeA"   "params.bott"
sumstats.bott <- subset(sumstats.expeA, subset=models.expeA=="bott")

# Perform parallel cross validation
parallel.cv.modsel <-
  parallel_cv4abc(
    param = params.bott, 
    sumstat = sumstats.bott, 
    nval=nval, 
    tols=seq(tol.low, tol.high, length.out = n.tols),
    nthreads=threads,
    method=method,
    abc.out = NULL
  )

# Compute misclassification percentages
rej.f.tol <- compute_rmse_from_cv_model(parallel.cv.modsel)  

# Automatically save results to the data directory 
write.csv(
  rej.f.tol, 
  file=here(
    glue::glue("ABC/data/cv4abc_{method}_nval_{nval}_tols_{tol.low}_{tol.high}_{n.tols}.csv")
  ), 
  quote = F, 
  row.names = F
)
