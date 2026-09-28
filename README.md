# Breast cancer survival model transportability

## Overview
Multi-cohort study comparing ensemble strategies for breast cancer survival prediction.

## Cohorts
- METABRIC (n=1,815) - overall survival
- GBSG2 (n=686) - recurrence-free survival
- GSE96058 (n=3,069) - overall survival
- NKI (n=319) - distant metastasis-free survival

Total: 5,889 patients.

## Repository Structure
- data/ - pooled_all_cohorts.csv
- models/ - final_model_averaging_ensemble.rds, predict_ensemble.R
- scripts/ - 01_data_harmonisation.R, 02_train_ensemble.R, 03_meta_analysis.R, 04_figures.py
- results/ - CSV result files
- figures/ - Publication-ready figures (600 DPI)

## Key Results
| Strategy | Mean LOCO-CV C-index |
|----------|----------------------|
| Simple pooled Cox | 0.372 |
| Ensemble of RSFs | 0.492 |
| Average of cohort-specific Cox coefficients | 0.543 |
| Univariable meta-analysis | 0.576 |
| Equal-weight ensemble (final model) | 0.595 |

## Contact
Linda Osaghale - lindaosaghale@gmail.com