
require(docstring)
library(magrittr)
library(rlang)
library(MASS)

tidy_lda <- function(phylum.matrix, metadata.tb, grouping.var){
  # Defuse grouping var
  grouping.var <- rlang::enquo(grouping.var)
  
  # Select the categorical variable to perform LDA
  .grouper <- metadata.tb %>% 
    dplyr::select(!!grouping.var) %>%
    unlist %>% 
    as.character()
  
  # Find the set of unique labels of grouper
  .levels <- unique(.grouper)
  
  # Perform lda
  LDA <- 
    suppressWarnings(MASS::lda(t(phylum.matrix), grouping=.grouper))
  
  .scale.lda <- function(x){ LDA$scaling * x }
  # Create the list of LDA components
  lda.ls <- sapply(
    .levels,
    function(x){
      phylum.matrix[, .grouper == x] %>% 
        as.matrix %>% 
        apply(2, .scale.lda) %>% 
        colSums
    },
    simplify = T
  )
  
  lda.tb <- lda.ls %>% 
    purrr::map2(names(lda.ls), ~ tibble(LDA_means=.x, Group=.y)) %>% 
    purrr::reduce(~ rbind(.x, .y))
  
  lda.tb
}



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



