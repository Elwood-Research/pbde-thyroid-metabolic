# Research Plan: PBDEs and Metabolic Health

## 1. Study Design
This is a cross-sectional study using pooled serum samples from NHANES 2007-2012.

## 2. Study Population
- **Inclusion Criteria**:
  - Participants in NHANES 2007-2012 cycles.
  - Aged 20 years and older.
  - Included in the environmental toxicant pooled samples for PBDEs.
  - Available data for thyroid hormones and metabolic markers.
- **Exclusion Criteria**:
  - Pregnancy.
  - History of thyroid disease or use of thyroid medications (to be refined in Phase 3).
  - Missing key covariates (Age, Sex, BMI).

## 3. Data Integration
- **Pooled Samples**: PBDE levels are measured in pools of individuals. Each individual in a pool will be assigned the concentration measured for that pool.
- **Weights**: Use the specific subsample weights provided for environmental toxicants (`WTSA2YR` or `WTSAF2YR` as appropriate, with adjustments for pooled samples as per NHANES documentation).

## 4. Statistical Analysis
- **Descriptive Statistics**: Summarize participant characteristics and distribution of PBDEs and metabolic markers.
- **Correlation Analysis**: Spearman correlations between PBDE congeners and outcomes.
- **Regression Modeling**:
  - Multiple linear regression to assess associations between PBDEs (log-transformed) and thyroid/metabolic markers.
  - Models will be adjusted for Age, Sex, Race/Ethnicity, BMI, and Smoking status.
- **Mediation Analysis**: Use the `mediation` package in R (or equivalent in Python) to estimate the proportion of the effect of PBDEs on metabolic outcomes mediated by thyroid hormones.
- **Outlier Handling**: Remove observations with |z| > 4 for continuous variables.
- **Categorical Handling**: Exclude levels with <5% membership.

## 5. Software
- Python (Pandas, Statsmodels, Scipy, Seaborn).
