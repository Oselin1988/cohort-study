# scripts/04_figures.py
# Generate Figures 1-5 for the breast cancer cohort study

!pip install lifelines -q

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from google.colab import drive, files
from lifelines import KaplanMeierFitter

drive.mount('/content/drive')
drive_path = '/content/drive/MyDrive/metabric_project/'


# ---------- Figure 1: Cohort Overview ----------
df = pd.read_csv(drive_path + 'pooled_all_cohorts.csv')
df['cohort'] = df['cohort'].astype('category')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
sns.boxplot(x='cohort', y='time', data=df, ax=ax1, color='lightgray', showfliers=False)
sns.swarmplot(x='cohort', y='time', data=df, ax=ax1, size=2, alpha=0.6, palette='Set2')
ax1.set_xlabel('Cohort'); ax1.set_ylabel('Survival time (months)')
ax1.text(-0.15, 1.05, 'A', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

kmf = KaplanMeierFitter()
for cohort in df['cohort'].unique():
    mask = df['cohort'] == cohort
    kmf.fit(df.loc[mask, 'time'], event_observed=df.loc[mask, 'event'], label=cohort)
    kmf.plot(ax=ax2, ci_show=False)
ax2.set_xlabel('Time (months)'); ax2.set_ylabel('Survival probability')
ax2.legend(title='Cohort')
ax2.text(-0.15, 1.05, 'B', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

plt.tight_layout()
plt.savefig(drive_path + 'Figure1_cohort_overview.png', dpi=600, bbox_inches='tight')
plt.savefig(drive_path + 'Figure1_cohort_overview.tiff', dpi=600, bbox_inches='tight', format='tiff')
plt.show()


# ---------- Figure 2: Internal Validation ----------
cv_df = pd.read_csv(drive_path + 'figure2_cv_results.csv')

fig, ax = plt.subplots(figsize=(10, 6))
ax.errorbar(x=cv_df['cohort'], y=cv_df['cv_mean'], yerr=cv_df['cv_sd'],
            fmt='o', color='navy', ecolor='lightgray', elinewidth=3,
            capsize=5, markersize=8, label='5-fold CV (mean +/- SD)')
ax.scatter(cv_df['cohort'], cv_df['oob'], color='darkgreen', s=80, label='OOB', zorder=5)
ax.axhline(y=0.5, color='red', linestyle='--', linewidth=1, label='Random (0.5)')
ax.set_ylabel('C-index'); ax.set_xlabel('Cohort'); ax.set_ylim(0, 1)
ax.legend()
plt.tight_layout()
plt.savefig(drive_path + 'Figure2_internal_validation.png', dpi=600, bbox_inches='tight')
plt.savefig(drive_path + 'Figure2_internal_validation.tiff', dpi=600, bbox_inches='tight', format='tiff')
plt.show()


# ---------- Figure 3: METABRIC Model ----------
imp_df = pd.read_csv(drive_path + 'figure3_variable_importance.csv')
risk_df = pd.read_csv(drive_path + 'metabric_risk_groups.csv')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
imp_df = imp_df.sort_values('importance', ascending=True)
ax1.stem(imp_df['importance'], imp_df['variable'], basefmt=" ", linefmt='gray', markerfmt='ro')
ax1.set_xlabel('Permutation importance'); ax1.set_ylabel('Variable')
ax1.text(-0.15, 1.05, 'A', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

kmf = KaplanMeierFitter()
for group in risk_df['risk_group'].unique():
    mask = risk_df['risk_group'] == group
    kmf.fit(risk_df.loc[mask, 'time'], event_observed=risk_df.loc[mask, 'event'], label=group)
    kmf.plot(ax=ax2, ci_show=False)
ax2.set_xlabel('Time (months)'); ax2.set_ylabel('Survival probability')
ax2.legend(title='Risk group')
ax2.text(-0.15, 1.05, 'B', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

plt.tight_layout()
plt.savefig(drive_path + 'Figure3_metabric_model.png', dpi=600, bbox_inches='tight')
plt.savefig(drive_path + 'Figure3_metabric_model.tiff', dpi=600, bbox_inches='tight', format='tiff')
plt.show()


# ---------- Figure 4: External Validation ----------
ext_df = pd.read_csv(drive_path + 'figure4_external_val.csv')
calib_df = pd.read_csv(drive_path + 'calibration_data.csv')
internal_c = 0.701

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
ax1.scatter(ext_df['cohort'], ext_df['c_index'], color='purple', s=80, zorder=5)
ax1.axhline(y=0.5, color='gray', linestyle='--', linewidth=1, label='Random (0.5)')
ax1.axhline(y=internal_c, color='red', linestyle='-', linewidth=1, label=f'METABRIC internal ({internal_c:.3f})')
ax1.set_ylabel('C-index'); ax1.set_xlabel('External cohort'); ax1.set_ylim(0, 1)
ax1.legend()
ax1.text(-0.15, 1.05, 'A', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

for cohort in calib_df['cohort'].unique():
    subset = calib_df[calib_df['cohort'] == cohort]
    ax2.plot(subset['predicted'], subset['observed'], marker='o', label=cohort, linewidth=2, markersize=8)
ax2.plot([0, 1], [0, 1], linestyle='--', color='gray', linewidth=1, label='Perfect calibration')
ax2.set_xlabel('Predicted 5-year survival'); ax2.set_ylabel('Observed 5-year survival')
ax2.set_xlim(0.55, 0.8); ax2.set_ylim(0.85, 1.02)
ax2.legend()
ax2.text(-0.15, 1.05, 'B', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

plt.tight_layout()
plt.savefig(drive_path + 'Figure4_external_validation.png', dpi=600, bbox_inches='tight')
plt.savefig(drive_path + 'Figure4_external_validation.tiff', dpi=600, bbox_inches='tight', format='tiff')
plt.show()


# ---------- Figure 5: Multi-Cohort Strategies ----------
strat_df = pd.read_csv(drive_path + 'figure5_strategies.csv')
per_df = pd.read_csv(drive_path + 'figure5_per_cohort.csv')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
strat_df = strat_df.sort_values('mean_c', ascending=False)
ax1.scatter(strat_df['mean_c'], strat_df['strategy'], color='darkblue', s=80, zorder=5)
ax1.axvline(x=0.5, color='gray', linestyle='--', linewidth=1, label='Random (0.5)')
ax1.set_xlabel('Mean LOCO-CV C-index'); ax1.set_ylabel('Strategy')
ax1.set_xlim(0, 0.7); ax1.legend()
ax1.text(-0.15, 1.05, 'A', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

per_df = per_df.sort_values('cohort')
ax2.scatter(per_df['c_index'], per_df['cohort'], color='darkgreen', s=80, zorder=5)
mean_c = per_df['c_index'].mean()
ax2.axvline(x=0.5, color='gray', linestyle='--', linewidth=1, label='Random (0.5)')
ax2.axvline(x=mean_c, color='red', linestyle='-', linewidth=1, label=f'Mean ({mean_c:.3f})')
ax2.set_xlabel('LOCO-CV C-index (final ensemble)'); ax2.set_ylabel('Cohort')
ax2.set_xlim(0, 0.8); ax2.legend()
ax2.text(-0.15, 1.05, 'B', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

plt.tight_layout()
plt.savefig(drive_path + 'Figure5_strategies.png', dpi=600, bbox_inches='tight')
plt.savefig(drive_path + 'Figure5_strategies.tiff', dpi=600, bbox_inches='tight', format='tiff')
plt.show()

print("All figures saved to Drive.")
