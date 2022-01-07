
shhh <- suppressPackageStartupMessages # It's a library, so shhh!

shhh(library(abc))
shhh(library(foreach))
shhh(library(doParallel))
shhh(library(doRNG))

parallel_cv4postpr <- 
  function(index, sumstat, nval, tols, seed = 1234, subset = NULL, nthreads = NULL){
  #' Like cv4postpr, but parallelised
  #' Why would you want to do this sequentially?
  #' see ?cv4postpr
  
  mymodels <- levels(factor(index))
  if (length(colnames(sumstat))) {
    statnames <- colnames(sumstat)
  }
  else {
    warning("No statistics names are given, using S1, S2, ...", 
            call. = F, immediate = T)
    statnames <- paste("S", 1:numstat, sep = "")
  }  
    
  gwt <- rep(TRUE, length(sumstat[, 1]))
  gwt[attributes(na.omit(sumstat))$na.action] <- FALSE
  if (is.null(subset)) 
    subset <- rep(TRUE, length(sumstat[, 1]))
  gwt <- as.logical(gwt * subset)
  cvsamp <- unlist(tapply(c(1:length(index))[gwt], index[gwt], 
                          sample, nval))
  tols <- sort(tols)
  
  if (is.null(nthreads)){
    .ncores <- min(parallel::detectCores(), length(tols))
  } else {
    .ncores <- min(nthreads, parallel::detectCores())
  }

  cl <- makeCluster(.ncores)
  registerDoParallel(cl)
  set.seed(seed)
  tol.iter <- iterators::iter(tols)
tryCatch(
{
  allprobs <- foreach(
    mytol=tols, .inorder = T, .packages = c("abc"), .combine = c, .init = list()
  ) %dorng% {
    res <- matrix(ncol = length(unique(index)), nrow = length(cvsamp))
    for (i in seq_along(cvsamp)){
      mysamp <- cvsamp[i]
      mytrue <- index[mysamp]
      mytarget <- sumstat[mysamp, ]
      myindex <- index[-mysamp]
      mysumstat <- sumstat[-mysamp, ]
      mysubset <- subset[-mysamp]
      subres <- postpr(
        target = mytarget, 
        index = myindex, 
        sumstat = mysumstat, 
        method = "rejection",
        tol = mytol, 
        subset = mysubset
      )
      res[i, ] <- summary(subres, print = F)$Prob
    }
    colnames(res) <- mymodels
    rownames(res) <- index[cvsamp]
    someprobs <- list()
    someprobs[[paste("tol", mytol, sep = "")]] <- res
    
    someprobs
  }
}, 
  finally = stopCluster(cl)
)
  allnames <- 
    lapply(allprobs, apply, 1, function(xx) mymodels[which(xx == max(xx))])
  
  seeds <- attr(allprobs, "rng")
  doRNG_version <- attr(allprobs, "doRNG_version")
  attr(allprobs, "rng") <- NULL
  attr(allprobs, "doRNG_version") <- NULL
  
  cv4postpr.out <- list(
    cvsamples = cvsamp, 
    tols = tols, 
    true = index[cvsamp], 
    estim = allnames, 
    model.probs = allprobs, 
    method = "rejection", 
    names = list(
      models = mymodels, 
      statistics.names = statnames
    ), 
    seeds = seeds,
    doRNG_version = doRNG_version
  )
  
  class(cv4postpr.out) <- "cv4postpr.parallel"
  return(invisible(cv4postpr.out))
}

#  ind.fold <- sample(1:k_folds, nrow(df), replace=T)
#  
#  .splits <- list()
#  for (i in 1:k_folds){
#    .splits[[i]] <- list(
#      train = df[ind.fold != i, ],
#      test = df[ind.fold == i, ]
#    )
#  }