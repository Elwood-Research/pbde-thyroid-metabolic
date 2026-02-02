# Statistical Methods

## Study Population
This study uses data from the National Health and Nutrition Examination Survey (NHANES) cycles from 2005-2006 through 2011-2012, where both pooled serum polybrominated diphenyl ethers (PBDEs) and individual thyroid hormone data were available. Inclusion criteria were defined as:
- Participants aged 20 years and older.
- Non-pregnant individuals (identified by self-report or laboratory pregnancy tests).
- Valid laboratory measurements for serum PBDEs, thyroid hormones (TSH and Free T4), and metabolic markers (lipids and insulin resistance).

A STROBE flow diagram will be constructed to detail the exclusions and the final analytic sample size.

## Laboratory Measurements and Data Handling
### Pooled Serum PBDEs
Serum PBDE levels were measured in pooled samples as part of the Brominated Flame Retardants (BFRs) component. Individual participant data will be linked to pooled results using the pool identifiers provided in the `POOLTF` datasets. Each individual within a pool is assigned the concentration measured for that pool.

### Thyroid Hormones
Thyroid-stimulating hormone (TSH) and Free Thyroxine (Free T4) were measured using standard NHANES laboratory protocols. 

### Metabolic Outcomes
- **Lipid Profile**: Low-density lipoprotein (LDL) cholesterol, high-density lipoprotein (HDL) cholesterol, and triglycerides were measured in fasting subsamples.
- **Insulin Resistance**: HOMA-IR was calculated using the formula: [Fasting Insulin (µU/mL) × Fasting Glucose (mg/dL)] / 405.

### Outlier Screening and Categorical Variables
- **Outliers**: Continuous variables will be screened for extreme outliers. Observations with an absolute z-score > 4 (|z| > 4) will be removed prior to modeling.
- **Categorical Exclusions**: Categorical variable levels with less than 5% membership will be excluded or collapsed into other categories to ensure statistical stability.

## Statistical Analysis
All analyses will account for the complex, multistage, probability sampling design of NHANES. 

### Survey Weights
Analyses will use the specific pooled sample weights (`WTSMSMPA`) provided in the `BFRPOL` datasets to account for the pooled sampling design and provide nationally representative estimates.

### Descriptive Statistics
Descriptive statistics will be calculated for the total population and stratified by demographic characteristics. Continuous variables will be presented as means (standard errors) or medians (interquartile ranges), and categorical variables as frequencies (percentages).

### Regression Models
Multiple linear regression models will be used to evaluate the associations:
1. **PBDEs vs. Thyroid Hormones**: Linear regression models will assess the relationship between log-transformed PBDE concentrations and thyroid hormone levels (log-TSH and Free T4).
2. **PBDEs vs. Metabolic Outcomes**: Linear regression models will assess the relationship between PBDEs and metabolic markers (LDL, HDL, Triglycerides, and HOMA-IR).
3. **Thyroid Hormones vs. Metabolic Outcomes**: The association between thyroid hormones and metabolic markers will also be examined.

All models will be adjusted for potential confounders, including:
- Age (years, continuous)
- Sex (Male, Female)
- Race/Ethnicity (Non-Hispanic White, Non-Hispanic Black, Mexican American, Other)
- Body Mass Index (BMI, kg/m²)
- Serum Cotinine (ng/mL, as a proxy for smoking status)

### Mediation Analysis
If significant associations are found between PBDEs and both thyroid hormones and metabolic outcomes, a mediation analysis will be conducted to explore whether thyroid hormones act as intermediaries in the pathway between PBDE exposure and metabolic dysfunction.

### Software
All statistical analyses will be performed using Python (v3.10+) with `pandas`, `statsmodels`, and `scipy` libraries.
