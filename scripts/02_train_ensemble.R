# scripts/02_train_ensemble.R
# Train four cohort-specific Cox models and save the equal-weight ensemble

suppressPackageStartupMessages({
  library(survival); library(dplyr)
})
setwd("/content/drive/MyDrive/metabric_project")

pooled <- read.csv("pooled_all_cohorts.csv", stringsAsFactors = FALSE)

# Standardise continuous predictors using global means/SDs
std_vars <- c("age", "tumor_size", "lymph_nodes", "grade")
for (v in std_vars) {
  pooled[[paste0(v, "_std")]] <- as.numeric(scale(pooled[[v]]))
}
global_means <- sapply(std_vars, function(v) mean(pooled[[v]], na.rm = TRUE))
global_sds   <- sapply(std_vars, function(v) sd(pooled[[v]], na.rm = TRUE))

cohorts <- split(pooled, pooled$cohort)
features <- c("age_std", "tumor_size_std", "lymph_nodes_std", "grade_std",
              "er_pos", "pr_pos", "her2_pos")

# Train one Cox model per cohort
models <- list()
for (name in names(cohorts)) {
  m <- coxph(as.formula(paste("Surv(time, event) ~", paste(features, collapse = " + "))),
             data = cohorts[[name]])
  models[[name]] <- m
}

# Save ensemble
ensemble <- list(models = models,
                 global_means = global_means,
                 global_sds = global_sds,
                 features = features)
saveRDS(ensemble, "final_model_averaging_ensemble.rds")
cat("Ensemble saved to final_model_averaging_ensemble.rds\n")
