
library(stringr)

compute_misclassif_from_cv_model <- function(cv.modsel.reject){
  misclassif.reject <- list()
  for (tol in names(cv.modsel.reject$estim)) {
    misclassif.reject[[tol]] = mean(cv.modsel.reject$estim[[tol]] != cv.modsel.reject$true)
  }
  .tols <- as.numeric(str_remove(names(misclassif.reject), "tol"))
  data.frame(tolerance=.tols, perc.misclassified=as.numeric(misclassif.reject))
}