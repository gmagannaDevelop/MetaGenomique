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
library(tibble)

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

compute_rmse_from_cv_model <- function(abc.cv){
  #' Return a list containing the RMSE
  #' for each one of the original sumstats
  
  .yy <- abc.cv$true # ground truth
  .se.f <- function(y_hat){ (y_hat - .yy)^2 } # squared error function
  .root.mean.f <- function(xx){ # root of mean function
    apply(xx, 2, function(yy){ sqrt(mean(yy)) }) 
  }
  
  # get squared error for each tolerance in the abc.cv$estim list
  # (ground_truth_matrix - estimation_matrix)^2
  .diffs.by.tol.ls <- sapply(abc.cv$estim, .se.f, simplify = F, USE.NAMES = T)
  # build a list with the errors
  .rmse.ls <- sapply(.diffs.by.tol.ls, .root.mean.f, simplify = F, USE.NAMES = T)
  # Transform it into a data.frame
  .rmse.df <- do.call(rbind, .rmse.ls) %>% 
      as.data.frame %>% 
      tibble::rownames_to_column("Tol")
  
  .rmse.df$Tol <- as.numeric(str_remove(.rmse.df$Tol, "tol"))
  .rmse.df
}

my.read.table <- function(x){
  #' drop-in replacement for base read.table
  read.table(x, header = T, sep = ",")
}

my.write.csv <- function(df, filename){
  write.csv(
    df, 
    file=filename, 
    quote = F, 
    row.names = F
  )
}

plot_cv_rej <- function(rej.df){
  #' Scatterplot for cv4postpr
  #' Percentage misclassified ~ tolerance
  cv.plot <- rej.df %>% 
    ggplot(aes(x=tolerance, y=perc.misclassified)) + 
    geom_point()
  cv.plot
}

gof_for_population <- function (population, tol = 0.025){
  #' Calculate Goodness-of-fit for a given population
  #' This function is partially hardcoded and therefore unreliable
  #' Proceed with caution
  
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

abc_with_distributions <- 
  function (population, col.name, tol = 0.05, method = "rejection"){
  #' hardcoded
  
  estim_x <- stat.1000g[population, ] %>% as.matrix 
  
  abc_x <- abc(target = estim_x,
                 param = params.bott,   #### A confirmer
                 sumstat = sumstats.bott,
                 tol=tol,
                 method = method)
  
  plot_x <- 
    plot_abc_distributions(res.abc = abc_x, params.df = params.bott, col.name)
  
  list(
    abc.obj = abc_x,
    plot = plot_x
  )
}



join_abc_distributions <- function(res.abc, params.df, col.name){
  #' Generate a data frame with three columns :
  #' * One for the prior 
  #' * One for the naive posterior
  #' * One fot the corrected posterior 
  data.frame(
    Prior = params.df[, col.name],
    Posterior = res.abc$unadj.values[, col.name],
    Posterior.corr = res.abc$adj.values[, col.name]
  )
}

plot_abc_distributions <- function(res.abc, params.df, col.name){
  #' Plot prior, posterior and corrected posterior
  .dist.plot <- join_abc_distributions(res.abc, params.df, col.name) %>% 
    reshape2::melt() %>% 
    tibble::as_tibble() %>% 
    dplyr::select(Loi = variable, everything()) %>% 
    ggplot2::ggplot(aes(x=value, colour=Loi)) + 
    geom_density() + labs(title = glue::glue("Lois pour '{col.name}'"), x=col.name)
  .dist.plot
}

normalise_rmse <- function(rmse.df, .melt=T){
  #' Normalise the RMSE so that all parameters
  #' are on the same order of magnitude.
  #' Basically a fancy wrapper for scale(center=F)
  #' IMPORTANT :
  #' Tolerance is assumed to be the first column
  normalised.tb <- rmse.df[, 2:ncol(rmse.df)] %>% 
    scale(center = F) %>% 
    as.data.frame %>% 
    tibble::as_tibble()
  
  normalised.tb %<>% 
    dplyr::mutate(Tolerance = rmse.df[, 1]) %>% 
    dplyr::select(Tolerance, dplyr::everything())
  
  if (.melt) {
    normalised.tb %<>% 
      reshape2::melt(
        id.vars="Tolerance", 
        variable.name = "Parameter",
        value.name = "Scaled Error"
      ) %>% 
      tibble::as_tibble()
  }
  
  normalised.tb
}

plot_abc_cv_rmse <- function(rmse.df){
  #' Chain a call to normalise_rmse and 
  #' some ggplot ;)
  
  .rmse <- normalise_rmse(rmse.df = rmse.df, .melt = T)
  .cv.rmse.plot <- .rmse %>% 
    ggplot(aes(x=Tolerance, y=`Scaled Error`)) + 
    geom_point(aes(colour=Parameter))
  .cv.rmse.plot
}

rmse_weighted_average <- function(rmse.df, .weights=NULL){
  #' Compute a weighted average of parameters' RMSE
  #' in order to obtain a sensible estimate of the global error
  #' IMPORTANT :
  #' Tolerance is assumed to be the first column
  #' 
  #' Weights can be optionally passed as a vector of the
  #' same lenght as the number of sumstats
  #' Otherwise their standard deviation is used 
  #' w ~ sd(stat)
  
  if (is.null(.weights)) .w <- apply(rmse.df[, 2:ncol(rmse.df)], 2, sd)
  else .w <- .weights
    
  data.frame(
    Tolerance=rmse.df[, 1],
    w.RMSE=apply(
      rmse.df[, 2:ncol(rmse.df)], 
      1, 
      weighted.mean, 
      w=.w
    )
  )
}
