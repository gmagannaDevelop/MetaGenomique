
library(abc)
library(foreach)
library(doParallel)
library(doRNG)

parallel_cv4postpr <- 
  function(index, sumstat, nval, tols, seed = 1234, subset = NULL){
  #' Like cv4postpr, but parallelised
  #' Why would you want to do this sequentially?
  #' see ?cv4postpr
  
  gwt <- rep(TRUE, length(sumstat[, 1]))
  gwt[attributes(na.omit(sumstat))$na.action] <- FALSE
  if (is.null(subset)) 
    subset <- rep(TRUE, length(sumstat[, 1]))
  gwt <- as.logical(gwt * subset)
  cvsamp <- unlist(tapply(c(1:length(index))[gwt], index[gwt], 
                          sample, nval))
  tols <- sort(tols)
  allprobs <- list()
  
  #.ncores <- min(parallel::detectCores(), length(tols))
  #cl <- makeCluster(.ncores)
  #registerDoParallel(cl)
  set.seed(seed)
  tol.iter <- iterators::iter(tols)
tryCatch(
{
  parallel.cv <- foreach(
    tol=tols, .inorder = T, .packages = c("abc"),
    .combine = c
  ) %do% {
    res <- matrix(ncol = length(unique(index)), nrow = length(cvsamp))
    for (i in seq_along(cvsamp)){
      mysamp <- cvsamp[i]
      mytrue <- index[mysamp]
      mytarget <- sumstat[mysamp, ]
      myindex <- index[-mysamp]
      mysumstat <- sumstat[-mysamp, ]
      mysubset <- subset[-mysamp]
      subres <- postpr(target = mytarget, index = myindex, 
                       sumstat = mysumstat, tol = mytol, subset = mysubset)
      res[i, ] <- summary.postpr(subres, print = F, ...)$Prob
    }
    colnames(res) <- mymodels
    rownames(res) <- index[cvsamp]
    allprobs[[paste("tol", mytol, sep = "")]] <- res
    allnames <- 
      lapply(allprobs, apply, 1, function(xx) mymodels[which(xx == max(xx))])
    
    res
}
  }, 
  finally = print("ok")
  #finally = stopCluster(cl)
)
  
  return(parallel.cv)
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