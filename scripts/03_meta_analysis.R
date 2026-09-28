# scripts/03_meta_analysis.R
# Univariable random-effects meta-analysis + forest plot

suppressPackageStartupMessages({
  library(dplyr); library(survival); library(metafor); library(ggplot2)
})
setwd("/content/drive/MyDrive/metabric_project")

pooled <- read.csv("pooled_all_cohorts.csv", stringsAsFactors = FALSE)
for (v in c("age", "tumor_size", "lymph_nodes", "grade"))
  pooled[[paste0(v, "_std")]] <- as.numeric(scale(pooled[[v]]))

cohorts  <- unique(pooled$cohort)
features <- c("age_std", "tumor_size_std", "lymph_nodes_std", "grade_std",
              "er_pos", "pr_pos", "her2_pos")

# Meta-analysis function for one predictor
get_meta <- function(var) {
  b <- c(); s <- c()
  for (coh in cohorts) {
    d <- pooled %>% filter(cohort == coh) %>% select(time, event, all_of(var))
    d <- d[complete.cases(d), ]
    if (nrow(d) < 3) next
    fit <- tryCatch(coxph(as.formula(paste("Surv(time,event) ~", var)), data = d),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      b <- c(b, coef(fit))
      s <- c(s, sqrt(vcov(fit)))
    }
  }
  if (length(b) >= 2) {
    r <- rma(yi = b, sei = s, method = "REML")
    return(data.frame(predictor = var,
                      HR = exp(r$b[1]),
                      HR_lower = exp(r$ci.lb),
                      HR_upper = exp(r$ci.ub),
                      tau2 = r$tau2, I2 = r$I2, pval = r$pval))
  }
  NULL
}

meta_df <- do.call(rbind, lapply(features, get_meta))
meta_df$label <- c("Age", "Tumor size", "Lymph nodes", "Grade",
                   "ER status", "PR status", "HER2 status")
write.csv(meta_df, "meta_analysis_results.csv", row.names = FALSE)

# Forest plot
p <- ggplot(meta_df, aes(x = HR, y = reorder(label, HR))) +
  geom_point(size = 3, color = "navy") +
  geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2, color = "navy") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_x_log10() +
  labs(x = "Hazard Ratio (95% CI)", y = "",
       title = "Pooled hazard ratios from univariable meta-analysis") +
  theme_minimal(base_size = 14)

ggsave("Figure5_forest_plot.png", p, width = 8, height = 5, dpi = 600, bg = "white")
ggsave("Figure5_forest_plot.tiff", p, width = 8, height = 5, dpi = 600,
       compression = "lzw", bg = "white")
cat("Meta-analysis and forest plot saved.\n")
