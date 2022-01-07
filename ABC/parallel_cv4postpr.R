library(here)
source(here("ABC/cv4postpr.R"))

args <- commandArgs(trailingOnly=T)
# basic input-output parameters
workdir <- args[1]
IPfile <- args[2]
controlFile <- args[3]


parallel.cv.modsel.reject.full <-
  parallel_cv4postpr(
    models.expeA, sumstats.expeA, nval=50, tols=seq(0.005, 0.95, length.out = 60)
  )

