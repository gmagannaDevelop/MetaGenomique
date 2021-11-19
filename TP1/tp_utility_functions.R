
library(docstring)

clean_data_frame <- function(my.matrix){
  #' Convert a matrix to a dataframe, 
  #' eliminating null columns
  .df <- as.data.frame(my.matrix)
  .null.col.sums <- colSums(.df) == 0
  .to.drop <- .null.col.sums[.null.col.sums]
  .cols.to.drop <- names(.to.drop)
  .clean.df <- .df[, !colnames(.df) %in% .cols.to.drop]
  .retval <- list(
    clean.df=.clean.df,
    null.cols=.cols.to.drop
  )
  return(.retval)
}

H <- function(x){ 
  #' Shannon's entropy of an entry vector
  #' @param x, the vector whose entropy is to be calculated
  p <- x[x > 0] 
  - sum(p * log(p)) 
}

