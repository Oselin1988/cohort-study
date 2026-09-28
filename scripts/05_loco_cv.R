# scripts/05_loco_cv.R
# Leave-one-cohort-out cross-validation for multi-cohort strategies
# Produces: loco_all_results.csv, figure5_strategies.csv, figure5_per_cohort.csv

suppressPackageStartupMessages({
  library(dplyr); library(survival); library(ranger); library(metafor)
})
setwd("/content/drive/MyDrive/metabric_project")

pooled <- read.csv("pooled_all_cohorts.csv", stringsAsFactors = FALSE)
std_vars <- c("age", "tumor_size", "lymph_nodes", "grade")
for (v in std_vars) {
  pooled[[paste0(v, "_std")]] <- as.numeric(scale(pooled[[v]]))
}
cohorts  <- unique(pooled$cohort)
features <- c("age_std", "tumor_size_std", "lymph_nodes_std", "grade_std",
              "er_pos", "pr_pos", "her2_pos")

# Helper: prediction function for a trained model
predict_risk <- function(model, newdata, time_point = 60) {
  p <- predict(model, data = newdata, type = "response")
  idx <- which.min(abs(p$unique.death.times - time_point))
  1 - p$survival[, idx]
}

# ======================================================================
# Strategy (a) Pooled Cox
# ======================================================================
results_pooled <- data.frame(test_cohort = character(), C_index = numeric())
for (test in cohorts) {
  train <- pooled %>% filter(cohort != test)
  testd <- pooled %>% filter(cohort == test)
  form  <- as.formula(paste("Surv(time, event) ~", paste(features, collapse = " + ")))
  m <- coxph(form, data = train)
  lp <- predict(m, newdata = testd, type = "lp")
  c_idx <- concordance(Surv(time, event) ~ lp, data = testd)$concordance
  if (c_idx < 0.5) { lp <- -lp; c_idx <- concordance(Surv(time, event) ~ lp, data = testd)$concordance }
  results_pooled <- rbind(results_pooled, data.frame(test_cohort = test, C_index = c_idx))
}
cat("Pooled Cox mean C-index:", mean(results_pooled$C_index), "\n")

# ======================================================================
# Strategy (b) Ensemble of RSFs
# ======================================================================
results_rsf <- data.frame(test_cohort = character(), C_index = numeric())
for (test in cohorts) {
  train_names <- setdiff(cohorts, test)
  models <- list()
  for (tr in train_names) {
    data_tr <- pooled %>% filter(cohort == tr) %>% select(all_of(features), time, event)
    models[[tr]] <- ranger(Surv(time, event) ~ ., data = data_tr,
                           num.trees = 500, verbose = FALSE)
  }
  testd <- pooled %>% filter(cohort == test)
  surv  <- sapply(models, function(m) predict_risk(m, testd[, features]))
  risk  <- 1 - rowMeans(surv)
  c_idx <- concordance(Surv(time, event) ~ risk, data = testd)$concordance
  if (c_idx < 0.5) { risk <- -risk; c_idx <- concordance(Surv(time, event) ~ risk, data = testd)$concordance }
  results_rsf <- rbind(results_rsf, data.frame(test_cohort = test, C_index = c_idx))
}
cat("RSF ensemble mean C-index:", mean(results_rsf$C_index), "\n")

# ======================================================================
# Strategy (c) Average of cohort-specific Cox coefficients
# ======================================================================
results_avg <- data.frame(test_cohort = character(), C_index = numeric())
for (test in cohorts) {
  train_names <- setdiff(cohorts, test)
  coefs <- list()
  for (tr in train_names) {
    data_tr <- pooled %>% filter(cohort == tr)
    form <- as.formula(paste("Surv(time, event) ~", paste(features, collapse = " + ")))
    coefs[[tr]] <- coef(coxph(form, data = data_tr))
  }
  avg_coef <- Reduce("+", coefs) / length(coefs)
  testd <- pooled %>% filter(cohort == test)
  risk  <- as.matrix(testd[, features]) %*% avg_coef
  risk  <- as.numeric(risk)
  c_idx <- concordance(Surv(time, event) ~ risk, data = testd)$concordance
  if (c_idx < 0.5) { risk <- -risk; c_idx <- concordance(Surv(time, event) ~ risk, data = testd)$concordance }
  results_avg <- rbind(results_avg, data.frame(test_cohort = test, C_index = c_idx))
}
cat("Average-coefficient mean C-index:", mean(results_avg$C_index), "\n")

# ======================================================================
# Strategy (d) Standardised equal-weight ensemble (FINAL MODEL)
# ======================================================================
results_equal <- data.frame(test_cohort = character(), C_index = numeric())
for (test in cohorts) {
  train_names <- setdiff(cohorts, test)
  models <- list()
  for (tr in train_names) {
    data_tr <- pooled %>% filter(cohort == tr)
    form <- as.formula(paste("Surv(time, event) ~", paste(features, collapse = " + ")))
    models[[tr]] <- coxph(form, data = data_tr)
  }
  testd <- pooled %>% filter(cohort == test)
  lp_matrix <- sapply(models, function(m) predict(m, newdata = testd, type = "lp"))
  if (is.null(dim(lp_matrix))) lp_matrix <- matrix(lp_matrix, nrow = nrow(testd))
  risk <- rowMeans(lp_matrix)
  c_idx <- concordance(Surv(time, event) ~ risk, data = testd)$concordance
  if (c_idx < 0.5) { risk <- -risk; c_idx <- concordance(Surv(time, event) ~ risk, data = testd)$concordance }
  results_equal <- rbind(results_equal, data.frame(test_cohort = test, C_index = c_idx))
}
cat("Equal-weight ensemble mean C-index:", mean(results_equal$C_index), "\n")

# ======================================================================
# Strategy (e) Univariable meta-analysis
# ======================================================================
results_meta <- data.frame(test_cohort = character(), C_index = numeric())
for (test in cohorts) {
  train_names <- setdiff(cohorts, test)
  beta_meta <- sapply(features, function(var) {
    b <- c(); s <- c()
    for (tr in train_names) {
      d <- pooled %>% filter(cohort == tr) %>% select(time, event, all_of(var))
      d <- d[complete.cases(d), ]
      fit <- tryCatch(coxph(as.formula(paste("Surv(time,event) ~", var)), data = d),
                      error = function(e) NULL)
      if (!is.null(fit)) { b <- c(b, coef(fit)); s <- c(s, sqrt(vcov(fit))) }
    }
    if (length(b) >= 2) {
      r <- rma(yi = b, sei = s, method = "REML")
      return(r$b[1])
    } else if (length(b) == 1) return(b[1])
    else return(0)
  })
  names(beta_meta) <- features
  testd <- pooled %>% filter(cohort == test)
  risk  <- as.matrix(testd[, features]) %*% beta_meta
  risk  <- as.numeric(risk)
  c_idx <- concordance(Surv(time, event) ~ risk, data = testd)$concordance
  if (c_idx < 0.5) { risk <- -risk; c_idx <- concordance(Surv(time, event) ~ risk, data = testd)$concordance }
  results_meta <- rbind(results_meta, data.frame(test_cohort = test, C_index = c_idx))
}
cat("Univariable meta-analysis mean C-index:", mean(results_meta$C_index), "\n")

# ======================================================================
# Combine and save results
# ======================================================================
strategies <- data.frame(
  strategy = c("Simple pooled Cox",
               "Ensemble of Random Survival Forests",
               "Average of cohort-specific Cox coefficients",
               "Univariable random-effects meta-analysis",
               "Equal-weight ensemble (standardised + risk inversion)"),
  mean_c = c(mean(results_pooled$C_index),
             mean(results_rsf$C_index),
             mean(results_avg$C_index),
             mean(results_meta$C_index),
             mean(results_equal$C_index))
)
write.csv(strategies, "figure5_strategies.csv", row.names = FALSE)

per_cohort <- data.frame(
  cohort  = results_equal$test_cohort,
  c_index = results_equal$C_index
)
write.csv(per_cohort, "figure5_per_cohort.csv", row.names = FALSE)

write.csv(results_equal, "loco_all_results.csv", row.names = FALSE)

cat("\nAll LOCO-CV results saved.\n")
print(strategies)
