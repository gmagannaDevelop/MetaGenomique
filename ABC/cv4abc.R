shhh <- suppressPackageStartupMessages # It's a library, so shhh!

shhh(library(abc))
shhh(library(foreach))
shhh(library(doParallel))
shhh(library(doRNG))

parallel_cv4abc <- function (param, sumstat, abc.out = NULL, nval, tols, statistic = "median", 
          prior.range = NULL, method, hcorr = TRUE, transf = "none", 
          logit.bounds = c(0, 0), subset = NULL, kernel = "epanechnikov", 
          numnet = 10, sizenet = 5, lambda = c(1e-04, 0.001, 0.01), 
          trace = FALSE, maxit = 500, seed = 1234, nthreads = NULL, ...) {
  mywarn <- options()$warn
  options(warn = -1)
  linout <- TRUE
  if (!any(statistic == c("median", "mean", "mode"))) {
    stop("Statistic has to be mean, median or mode.", call. = F)
  }
  if (is.null(abc.out) && missing(method)) 
    stop("Method must be supplied when 'abc.out' is NULL.", 
         call. = F)
  if (missing(nval)) 
    stop("'nval' must be supplied.", call. = F)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (!is.null(abc.out)) {
    subset <- abc.out$na.action
    method <- abc.out$method
    transf <- abc.out$transf
    logit.bounds <- abc.out$logit.bounds
    kernel <- "epanechnikov"
  }
  if (is.null(dim(param))) {
    np <- 1
    param <- as.data.frame(param)
    names(param) <- "P1"
  }
  else np <- dim(param)[2]
  if (is.null(dim(sumstat))) {
    numstat <- 1
    sumstat <- as.data.frame(sumstat)
    names(sumstat) <- "S1"
  }
  else {
    numstat <- dim(sumstat)[2]
  }
  numsim <- dim(sumstat)[1]
  if (!is.null(abc.out)) {
    if (np != abc.out$numparam || numstat != abc.out$numstat || 
        numsim != length(abc.out$na.action)) {
      stop("The number of parameters, summary statistics, or simulations provided in 'param' or 'sumstat' are not the same as in 'abc.out'.", 
           call. = F)
    }
    else if (!prod(colnames(param) %in% abc.out$names$parameter.names)) {
      stop("Parameters in 'param' are not the same as in 'abc.out', or different names are used.", 
           call. = F)
    }
    else if (!prod(colnames(sumstat) %in% abc.out$names$statistics.names)) {
      stop("Summary statistics in 'sumstat' are not the same as in 'abc.out', or different names are used.", 
           call. = F)
    }
    else {
      paramnames <- abc.out$names$parameter.names
      statnames <- abc.out$names$statistics.names
    }
  }
  else {
    if (length(colnames(param))) {
      paramnames <- colnames(param)
    }
    else {
      paramnames <- paste("P", 1:np, sep = "")
    }
  }
  if (length(colnames(sumstat))) {
    statnames <- colnames(sumstat)
  }
  else {
    statnames <- paste("S", 1:numstat, sep = "")
  }
  gwt <- rep(TRUE, length(sumstat[, 1]))
  gwt[attributes(na.omit(sumstat))$na.action] <- FALSE
  if (is.null(subset)) 
    subset <- rep(TRUE, length(sumstat[, 1]))
  gwt <- as.logical(gwt * subset)
  cvsamp <- sample(1:numsim, nval, prob = gwt/sum(gwt))
  tols <- sort(tols)
  num.panel <- length(tols)
  #alltol <- list()
  mycall <- list()
  
  if (is.null(nthreads)){
    .ncores <- min(parallel::detectCores(), length(tols))
  } else {
    .ncores <- min(nthreads, parallel::detectCores())
  }
  
  cl <- makeCluster(.ncores)
  registerDoParallel(cl)
  set.seed(seed)
tryCatch(
{  
  alltol <- 
  foreach (
    mytol=tols, .inorder = T, .packages = c("abc"), .combine = c, .init = list()
  ) %dorng% {
    res <- matrix(ncol = np, nrow = nval)
    for (i in seq_along(nval)) {
      mysamp <- cvsamp[i]
      mytrue <- param[mysamp, ]
      mytarget <- sumstat[mysamp, ]
      myparam <- param[-mysamp, ]
      mysumstat <- sumstat[-mysamp, ]
      mysubset <- subset[-mysamp]
      subres <- withCallingHandlers(abc(target = mytarget, 
                                        param = myparam, sumstat = mysumstat, tol = mytol, 
                                        subset = mysubset, method = method, transf = transf, 
                                        logit.bounds = logit.bounds, kernel = kernel, 
                                        hcorr = hcorr))
      if (statistic == "median") 
        estim <- invisible(summary(subres, print = F)[3, ])
      
      if (statistic == "mean") 
        estim <- invisible(summary(subres, print = F)[4, ])
      if (statistic == "mode") 
        estim <- invisible(summary(subres, print = F)[5, ])
      res[i, ] <- estim
    }
    if (np == 1) 
      res <- c(res)
    else colnames(res) <- paramnames
    sometol <- list()
    sometol[[paste("tol", mytol, sep = "")]] <- res

    sometol
  }
}, finally = stopCluster(cl)
)
  
  if (np == 1) {
    true <- as.data.frame(param[cvsamp, ])
    names(true) <- paramnames
  }
  else true <- param[cvsamp, ]
  
  seeds <- attr(alltol, "rng")
  doRNG_version <- attr(alltol, "doRNG_version")
  attr(alltol, "rng") <- NULL
  attr(alltol, "doRNG_version") <- NULL
  
  cv4abc.out <- list(
      cvsamples = cvsamp, tols = tols, 
      true = true, estim = alltol, 
      names = list(
        parameter.names = paramnames, 
        statistics.names = statnames
      ), 
      seeds = seeds
  )
  options(warn = mywarn)
  class(cv4abc.out) <- "cv4abc.parallel"
  invisible(cv4abc.out)
}
