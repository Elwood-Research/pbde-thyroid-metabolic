
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import MultiComparison

# Paths
DATA_DIR = "/home/joshbot/NHANES_BOT/Processed Data/Data/"
OUTPUT_DIR = "/home/joshbot/NHANES_BOT/studies/pbde-thyroid-metabolic-2026-02-02/04-analysis/outputs/"
TABLE_DIR = os.path.join(OUTPUT_DIR, "tables/")
FIGURE_DIR = os.path.join(OUTPUT_DIR, "figures/")

# Cycles to analyze
CYCLES = ['E', 'F', 'G']

def load_data(cycle):
    print(f"Loading data for cycle {cycle}...")
    
    # Core datasets
    demo = pd.read_csv(f"{DATA_DIR}DEMO_{cycle}.csv")
    bmx = pd.read_csv(f"{DATA_DIR}BMX_{cycle}.csv")
    
    # Cotinine (different names)
    if cycle == 'D':
        cot = pd.read_csv(f"{DATA_DIR}COT_D.csv")
        cot = cot.rename(columns={'LBXCOT': 'COT'})
    else:
        cot = pd.read_csv(f"{DATA_DIR}COTNAL_{cycle}.csv")
        cot = cot.rename(columns={'LBXCOT': 'COT'})
        
    thyroid = pd.read_csv(f"{DATA_DIR}THYROD_{cycle}.csv")
    tchol = pd.read_csv(f"{DATA_DIR}TCHOL_{cycle}.csv")
    hdl = pd.read_csv(f"{DATA_DIR}HDL_{cycle}.csv")
    trig = pd.read_csv(f"{DATA_DIR}TRIGLY_{cycle}.csv")
    glu = pd.read_csv(f"{DATA_DIR}GLU_{cycle}.csv")
    
    # Pooled sample linkage
    pool_link = pd.read_csv(f"{DATA_DIR}POOLTF_{cycle}.csv")
    pool_data = pd.read_csv(f"{DATA_DIR}BFRPOL_{cycle}.csv")
    
    # Merge individual data
    df = demo[['SEQN', 'RIAGENDR', 'RIDAGEYR', 'RIDRETH1' if cycle != 'G' else 'RIDRETH3']]
    df = df.rename(columns={'RIDRETH1': 'RACE', 'RIDRETH3': 'RACE'})
    
    df = df.merge(bmx[['SEQN', 'BMXBMI']], on='SEQN', how='left')
    df = df.merge(cot[['SEQN', 'COT']], on='SEQN', how='left')
    df = df.merge(thyroid[['SEQN', 'LBXTSH1', 'LBXT4F']], on='SEQN', how='left')
    df = df.merge(tchol[['SEQN', 'LBDTCSI']], on='SEQN', how='left')
    df = df.merge(hdl[['SEQN', 'LBDHDD']], on='SEQN', how='left')
    df = df.merge(trig[['SEQN', 'LBDTRSI', 'LBDLDL']], on='SEQN', how='left')
    df = df.merge(glu[['SEQN', 'LBXGLU', 'LBXIN']], on='SEQN', how='left')
    
    # Merge with pool linkage
    df = df.merge(pool_link[['SEQN', 'SAMPLEID']], on='SEQN', how='inner')
    
    # Merge with pool exposure data
    pbde_cols = ['SAMPLEID', 'WTSMSMPA', 'LBCBR3', 'LBCBR5', 'LBCBR4', 'LBCBR6', 'LBCBR7']
    df = df.merge(pool_data[pbde_cols], on='SAMPLEID', how='inner')
    
    df['Cycle'] = cycle
    return df

# Combined data
full_df = pd.concat([load_data(c) for c in CYCLES], ignore_index=True)

# Define Variables
# PBDEs (log-transformed)
pbde_names = {
    'LBCBR3': 'BDE-47',
    'LBCBR5': 'BDE-99',
    'LBCBR4': 'BDE-100',
    'LBCBR6': 'BDE-153',
    'LBCBR7': 'BDE-154'
}

# Calculated outcomes
full_df['HOMA_IR'] = (full_df['LBXGLU'] * full_df['LBXIN']) / 405.0

# Initial Count
initial_count = len(full_df)
print(f"Initial linked sample size: {initial_count}")

# Clean missing
# Exclude if no weight or no PBDEs (though merge was inner)
full_df = full_df.dropna(subset=['WTSMSMPA', 'LBCBR3'])
after_exposure_missing = len(full_df)

# STROBE tracking
exclusions = []
exclusions.append(('Linked individual participants', initial_count))

# Screen outcomes
outcomes = ['LBXTSH1', 'LBXT4F', 'LBDTCSI', 'LBDHDD', 'LBDTRSI', 'LBDLDL', 'HOMA_IR']
# For simplicity, we keep rows that have AT LEAST ONE thyroid outcome AND AT LEAST ONE metabolic outcome
# Actually, better to analyze each per outcome, but for Table 1 we need a clean analytic sample.
# Let's drop if missing Thyroid OR Metabolic
full_df = full_df.dropna(subset=['LBXTSH1', 'LBXT4F'])
after_thyroid_missing = len(full_df)
exclusions.append(('Missing thyroid hormones', initial_count - after_thyroid_missing))

full_df = full_df.dropna(subset=['LBDTCSI', 'LBDHDD', 'LBDTRSI', 'LBDLDL', 'HOMA_IR'])
after_metabolic_missing = len(full_df)
exclusions.append(('Missing metabolic outcomes', after_thyroid_missing - after_metabolic_missing))

# Screen Covariates
covars = ['RIAGENDR', 'RIDAGEYR', 'RACE', 'BMXBMI', 'COT']
full_df = full_df.dropna(subset=covars)
after_covar_missing = len(full_df)
exclusions.append(('Missing covariates', after_metabolic_missing - after_covar_missing))

# Outlier removal (|z| > 4)
continuous_vars = outcomes + ['RIDAGEYR', 'BMXBMI', 'COT'] + list(pbde_names.keys())
outlier_mask = pd.Series(False, index=full_df.index)
for var in continuous_vars:
    z = np.abs(stats.zscore(full_df[var], nan_policy='omit'))
    outlier_mask |= (z > 4)

before_outliers = len(full_df)
full_df = full_df[~outlier_mask]
after_outliers = len(full_df)
exclusions.append(('Extreme outliers (|z|>4)', before_outliers - after_outliers))

# Categorical levels (< 5%)
# Check Race and Gender (Gender should be fine)
for col in ['RACE', 'RIAGENDR']:
    counts = full_df[col].value_counts(normalize=True)
    low_membership = counts[counts < 0.05].index
    if len(low_membership) > 0:
        full_df = full_df[~full_df[col].isin(low_membership)]
        
after_categorical = len(full_df)
exclusions.append(('Small categorical levels (<5%)', after_outliers - after_categorical))

final_sample_size = len(full_df)
print(f"Final analytic sample size: {final_sample_size}")

# Save final sample info for STROBE
strobe_df = pd.DataFrame(exclusions, columns=['Step', 'Count'])
strobe_df.to_csv(os.path.join(OUTPUT_DIR, "strobe_counts.csv"), index=False)

# Log-transform PBDEs
for bde in pbde_names.keys():
    full_df[f'log_{bde}'] = np.log10(full_df[bde])

# Table 1: Descriptives
def get_weighted_stats(df, vars, weight_col):
    stats_list = []
    for var in vars:
        if df[var].dtype == 'object' or len(df[var].unique()) < 10:
            # Categorical
            counts = df.groupby(var)[weight_col].sum()
            percent = (counts / counts.sum()) * 100
            for val, pct in percent.items():
                stats_list.append({'Variable': f"{var}_{val}", 'Stat': f"{pct:.1f}%"})
        else:
            # Continuous
            mean = np.average(df[var], weights=df[weight_col])
            std = np.sqrt(np.average((df[var]-mean)**2, weights=df[weight_col]))
            stats_list.append({'Variable': var, 'Stat': f"{mean:.2f} ({std:.2f})"})
    return pd.DataFrame(stats_list)

table1_vars = ['RIDAGEYR', 'RIAGENDR', 'RACE', 'BMXBMI', 'COT', 'LBXTSH1', 'LBXT4F', 'LBDTCSI', 'LBDHDD', 'LBDTRSI', 'LBDLDL', 'HOMA_IR']
table1 = get_weighted_stats(full_df, table1_vars, 'WTSMSMPA')
table1.to_latex(os.path.join(TABLE_DIR, "table1.tex"), index=False)

# Regression Models
results_summary = []

def run_weighted_reg(df, exposure, outcome, covars, weight_col):
    X = df[[exposure] + covars]
    X = sm.add_constant(X)
    y = df[outcome]
    weights = df[weight_col]
    
    model = sm.WLS(y, X, weights=weights).fit()
    coef = model.params[exposure]
    se = model.bse[exposure]
    p = model.pvalues[exposure]
    conf_int = model.conf_int().loc[exposure]
    
    return {
        'Exposure': pbde_names.get(exposure, exposure),
        'Outcome': outcome,
        'Beta': coef,
        'Lower CI': conf_int[0],
        'Upper CI': conf_int[1],
        'P-value': p
    }

for bde in [f'log_{k}' for k in pbde_names.keys()]:
    for outcome in outcomes:
        res = run_weighted_reg(full_df, bde, outcome, ['RIDAGEYR', 'RIAGENDR', 'BMXBMI', 'COT'], 'WTSMSMPA')
        results_summary.append(res)

results_df = pd.DataFrame(results_summary)
results_df.to_csv(os.path.join(OUTPUT_DIR, "regression_results.csv"), index=False)

# Format for LaTeX
results_df['Result'] = results_df.apply(lambda x: f"{x['Beta']:.3f} ({x['Lower CI']:.3f}, {x['Upper CI']:.3f})", axis=1)
pivot_results = results_df.pivot(index='Exposure', columns='Outcome', values='Result')
pivot_results.to_latex(os.path.join(TABLE_DIR, "regression_table.tex"))

# Figures
# Forest Plot for TSH and T4
for outcome in ['LBXTSH1', 'LBXT4F']:
    plt.figure(figsize=(10, 6))
    subset = results_df[results_df['Outcome'] == outcome]
    plt.errorbar(subset['Beta'], subset['Exposure'], xerr=[subset['Beta']-subset['Lower CI'], subset['Upper CI']-subset['Beta']], fmt='o', capsize=5)
    plt.axvline(0, color='red', linestyle='--')
    plt.title(f'Association between PBDEs and {outcome}')
    plt.xlabel('Beta Coefficient (95% CI)')
    plt.savefig(os.path.join(FIGURE_DIR, f"forest_{outcome}.png"))
    plt.close()

# STROBE Diagram
def draw_strobe(exclusions, save_path):
    fig, ax = plt.subplots(figsize=(8, 10))
    y = 0.9
    for step, count in exclusions:
        if step == 'Linked individual participants':
            text = f"{step}\n(n = {count})"
        else:
            text = f"Excluded: {step}\n(n = {count})"
        
        ax.text(0.5, y, text, ha='center', va='center', bbox=dict(boxstyle='round', facecolor='white'))
        if y < 0.9:
            ax.annotate('', xy=(0.5, y+0.05), xytext=(0.5, y+0.1), arrowprops=dict(arrowstyle='->'))
        y -= 0.15
        
    ax.text(0.5, y, f"Final Analytic Sample\n(n = {final_sample_size})", ha='center', va='center', bbox=dict(boxstyle='round', facecolor='lightblue'))
    ax.annotate('', xy=(0.5, y+0.05), xytext=(0.5, y+0.1), arrowprops=dict(arrowstyle='->'))
    
    ax.set_axis_off()
    plt.savefig(save_path)
    plt.close()

draw_strobe(strobe_df.values, os.path.join(FIGURE_DIR, "strobe_diagram.png"))

# Summary Markdown
with open(os.path.join(OUTPUT_DIR, "../results_summary.md"), 'w') as f:
    f.write(f"# Analysis Results Summary\n\n")
    f.write(f"## Sample Size\n")
    f.write(f"- Initial linked sample: {initial_count}\n")
    f.write(f"- Final analytic sample: {final_sample_size}\n\n")
    f.write(f"## Key Findings\n")
    f.write(f"Significant associations (p < 0.05):\n\n")
    sig = results_df[results_df['P-value'] < 0.05]
    if sig.empty:
        f.write("No significant associations found at p < 0.05 level.\n")
    else:
        for _, row in sig.iterrows():
            f.write(f"- {row['Exposure']} vs {row['Outcome']}: Beta={row['Beta']:.3f}, p={row['P-value']:.3f}\n")
    
    f.write(f"\n## Outputs Saved\n")
    f.write(f"- Tables: `04-analysis/outputs/tables/`\n")
    f.write(f"- Figures: `04-analysis/outputs/figures/` (including STROBE diagram)\n")

print("Analysis complete.")
