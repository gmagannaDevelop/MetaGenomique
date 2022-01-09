#' Université Paris-Saclay
#' Faculté des Sciences d'Orsay
#' M2 AMI2B : Biologie Computationnelle
#' 
#' Métagénomique et génétique des populations
#' 
#' Gustavo Magaña López
#' Theó Roncalli
#' 

library(stringr)

compute_misclassif_from_cv_model <- function(cv.modsel.reject){
  #' Compute percentage of misclassified observations
  #' Takes a 'cv4postpr' or 'cv4postpr.parallel' object
  #' 
  misclassif.reject <- list()
  for (tol in names(cv.modsel.reject$estim)) {
    misclassif.reject[[tol]] = mean(cv.modsel.reject$estim[[tol]] != cv.modsel.reject$true)
  }
  .tols <- as.numeric(str_remove(names(misclassif.reject), "tol"))
  data.frame(tolerance=.tols, perc.misclassified=as.numeric(misclassif.reject))
}

my.read.table <- function(x){
  #' drop-in replacement for base read.table
  read.table(x, header = T, sep = ",")
}

plot_cv_rej <- function(rej.df){
  #' Scatterplot for cv4postpr
  #' Percentage misclassified ~ tolerance
  cv.plot <- rej.df %>% 
    ggplot(aes(x=tolerance, y=perc.misclassified)) + 
    geom_point()
  cv.plot
}

gof_for_population <- function (population, tol = 0.05){
  #' Calculate Goodness-of-fit for a given population
  #' This function is partially hardcoded and therefore unreliable
  #' Proceed with caution
  
  # We added a statistic which is not 
  # present in sumstats, that's why we have to take 
  # only the first three columns
  estim_x <- stat.1000g[population, ] %>% as.matrix 
  
  modsel_x <- postpr(target = estim_x, # pseudo-observed sumstats
                     index = models.expeA, # known models for reference db
                     sumstat = sumstats.expeA, # sumstats for reference db
                     method='rejection', # type of ABC algorithm
                     tol=tol) # proportion of kept simulations
  
  prob_x <- summary(modsel_x)$Prob
  gof_x <- gfitpca(target = estim_x, # pseudo-observed sumstats
                   index = models.expeA, # known models for reference db
                   sumstat = sumstats.expeA, # sumstats for reference db
                   main = glue::glue("ABC - goodness of fit - Pop : {population} "))
  
  list(
    postpr.summary = summary(modsel_x),
    posterior.probs = prob_x,
    gfitpca = gof_x
  )
}
