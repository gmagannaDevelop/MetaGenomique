cross_validation_split <- function(df, k_folds){
  #' Créer 'k' ensembles disjoints pour faire une 
  #' validation croisée à 'k' plis.
  ind.fold <- sample(1:k_folds, nrow(df), replace=T)
  
  .splits <- list()
  for (i in 1:k_folds){
    .splits[[i]] <- list(
      train = df[ind.fold != i, ],
      test = df[ind.fold == i, ]
    )
  }
  return(.splits)
}
