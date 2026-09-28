predict_ensemble <- function(newdata, ensemble_path = 'final_model_averaging_ensemble.rds') {
  ensemble <- readRDS(ensemble_path)
  models <- ensemble$models
  global_means <- ensemble$global_means
  global_sds <- ensemble$global_sds
  std_vars <- c('age', 'tumor_size', 'lymph_nodes', 'grade')
  for (v in std_vars) {
    newdata[[paste0(v, '_std')]] <- (newdata[[v]] - global_means[v]) / global_sds[v]
  }
  lp <- sapply(models, function(m) predict(m, newdata = newdata, type = 'lp'))
  if (is.null(dim(lp))) {
    mean(lp)
  } else {
    rowMeans(lp)
  }
}
