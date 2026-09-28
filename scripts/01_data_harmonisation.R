# scripts/01_data_harmonisation.R
# Harmonise four breast cancer cohorts into pooled_all_cohorts.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(GEOquery)
  library(TH.data)
  library(breastCancerNKI)
  library(Biobase)
})
setwd("/content/drive/MyDrive/metabric_project")

# ---- METABRIC ----
df <- read.csv("METABRIC_RNA_Mutation.csv", stringsAsFactors = FALSE)
metabric <- data.frame(
  time = df$overall_survival_months,
  event = ifelse(df$overall_survival == 0, 1, 0),
  tumor_size = df$tumor_size,
  age = df$age_at_diagnosis,
  lymph_nodes = df$lymph_nodes_examined_positive,
  grade = df$neoplasm_histologic_grade,
  er_pos = as.numeric(df$er_status == "Positive"),
  pr_pos = as.numeric(df$pr_status == "Positive"),
  her2_pos = as.numeric(df$her2_status == "Positive")
)
metabric <- metabric[complete.cases(metabric), ]
metabric$cohort <- "METABRIC"

# ---- GBSG2 ----
data("GBSG2", package = "TH.data")
gbsg2 <- data.frame(
  time = GBSG2$time / 30.44,
  event = 1 - GBSG2$cens,
  tumor_size = GBSG2$tsize,
  age = GBSG2$age,
  lymph_nodes = GBSG2$pnodes,
  grade = as.numeric(GBSG2$tgrade),
  er_pos = as.numeric(GBSG2$estrec > 10),
  pr_pos = as.numeric(GBSG2$progrec > 10),
  her2_pos = 0
)
gbsg2 <- gbsg2[complete.cases(gbsg2), ]
gbsg2$cohort <- "GBSG2"

# ---- GSE96058 ----
gse <- getGEO("GSE96058", GSEMatrix = TRUE)
gse_raw <- pData(phenoData(gse[[1]]))
get_col <- function(patterns) {
  for (pat in patterns) {
    col <- grep(pat, names(gse_raw), value = TRUE, ignore.case = TRUE)[1]
    if (!is.na(col)) return(gse_raw[[col]])
  }
  return(rep(NA, nrow(gse_raw)))
}
gse96058 <- data.frame(
  time = as.numeric(get_col(c("overall survival days", "os days"))) / 30.44,
  event = as.numeric(get_col(c("overall survival event", "os event"))),
  tumor_size = as.numeric(get_col(c("tumor size", "size"))) * 10,
  age = as.numeric(get_col(c("age at diagnosis", "age"))),
  lymph_nodes = ifelse(grepl("Positive", get_col(c("lymph node status", "node status")), ignore.case = TRUE), 1, 0),
  grade = as.numeric(factor(get_col(c("nhg", "grade")), levels = c("G1","G2","G3"), labels = 1:3)),
  er_pos = as.numeric(get_col(c("er status", "er")) == 1),
  pr_pos = as.numeric(get_col(c("pgr status", "pr")) == 1),
  her2_pos = as.numeric(get_col(c("her2 status", "her2")) == 1)
)
for (col in c("age","tumor_size","lymph_nodes","grade")) {
  med <- median(gse96058[[col]], na.rm = TRUE); if (is.na(med)) med <- 0
  gse96058[[col]][is.na(gse96058[[col]])] <- med
}
for (col in c("er_pos","pr_pos","her2_pos")) gse96058[[col]][is.na(gse96058[[col]])] <- 0
gse96058 <- gse96058[complete.cases(gse96058), ]
gse96058$cohort <- "GSE96058"

# ---- NKI ----
data("nki", package = "breastCancerNKI")
nki_raw <- pData(nki)
nki <- data.frame(
  time = nki_raw$t.dmfs / 30.44,
  event = nki_raw$e.dmfs,
  tumor_size = nki_raw$size,
  age = nki_raw$age,
  lymph_nodes = as.numeric(nki_raw$node == "1"),
  grade = nki_raw$grade,
  er_pos = as.numeric(nki_raw$er == 1),
  pr_pos = as.numeric(nki_raw$pgr == 1),
  her2_pos = as.numeric(nki_raw$her2 == 1)
)
for (col in c("age","tumor_size","lymph_nodes","grade")) {
  med <- median(nki[[col]], na.rm = TRUE); if (is.na(med)) med <- 0
  nki[[col]][is.na(nki[[col]])] <- med
}
for (col in c("er_pos","pr_pos","her2_pos")) nki[[col]][is.na(nki[[col]])] <- 0
nki <- nki[complete.cases(nki), ]
nki$cohort <- "NKI"

# ---- Combine ----
pooled <- bind_rows(metabric, gbsg2, gse96058, nki)
write.csv(pooled, "pooled_all_cohorts.csv", row.names = FALSE)
cat("Saved pooled_all_cohorts.csv -", nrow(pooled), "rows\n")
