
require(docstring)

clean_data_frame <- function(my.matrix){
  #' Clean data.frames 
  #' Cast a matrix into a data.frame, eliminating
  #' all columns having a sum equal to zero.
  #' @param my.matrix A matrix to be pruned and converted to data.frame
  #' @return A list with two attributes, clean.df (the pruned data.frame)
  #'         and null.cols (the colnames whose sum was equal to zero)
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
  #' Shannon's entropy
  #' @param x the vector whose entropy is to be calculated
  p <- x[x > 0] 
  - sum(p * log(p)) 
}

