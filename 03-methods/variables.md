# Study Variables

## Exposure Variables (PBDEs)
Source: Brominated Flame Retardants (BFRs) - Pooled Samples (`BFRPOL_D`, `BFRPOL_E`, `BFRPOL_F`, `BFRPOL_G`).

| Variable Name | Description | Transformation | Source Dataset |
|---------------|-------------|----------------|----------------|
| LBCBR3LA | PBDE 47 (ng/g lipid) | Log10 | BFRPOL |
| LBCBR5LA | PBDE 99 (ng/g lipid) | Log10 | BFRPOL |
| LBCBR6LA | PBDE 100 (ng/g lipid) | Log10 | BFRPOL |
| LBCBR7LA | PBDE 153 (ng/g lipid) | Log10 | BFRPOL |
| LBCBR8LA | PBDE 154 (ng/g lipid) | Log10 | BFRPOL |

## Thyroid Hormone Variables
Source: Thyroid Profile (`THYROD_E`, `THYROD_F`, `THYROD_G`).

| Variable Name | Description | Transformation | Source Dataset |
|---------------|-------------|----------------|----------------|
| LBXTSH1 | Thyroid Stimulating Hormone (mIU/L) | Log10 | THYROD |
| LBXT4F | Free Thyroxine (ng/dL) | None | THYROD |

## Metabolic Outcome Variables
Source: Fasting Glucose & Insulin (`GLU`), Triglycerides & LDL (`TRIGLY`), HDL Cholesterol (`HDL`).

| Variable Name | Description | Transformation | Source Dataset |
|---------------|-------------|----------------|----------------|
| LBDLDL | LDL-Cholesterol (mg/dL) | None | TRIGLY |
| LBDHDD | HDL-Cholesterol (mg/dL) | None | HDL |
| LBXTR | Triglycerides (mg/dL) | Log10 | TRIGLY |
| HOMA-IR | Insulin Resistance (Calculated) | Log10 | GLU |

*Note: HOMA-IR = [LBXIN (µU/mL) × LBXGLU (mg/dL)] / 405*

## Covariates
Source: Demographics (`DEMO`), Body Measures (`BMX`), Cotinine (`COT` / `COTNAL`).

| Variable Name | Description | Type | Source Dataset |
|---------------|-------------|------|----------------|
| RIDAGEYR | Age at screening (years) | Continuous | DEMO |
| RIAGENDR | Gender (1=Male, 2=Female) | Categorical | DEMO |
| RIDRETH1 | Race/Ethnicity (Recode) | Categorical | DEMO |
| BMXBMI | Body Mass Index (kg/m²) | Continuous | BMX |
| LBXCOT | Serum Cotinine (ng/mL) | Continuous | COT / COTNAL |

## Linking and Weights
- **Sample Linkage**: `POOLTF` datasets will be used to link individual `SEQN` to `SAMPLEID` in `BFRPOL`.
- **Survey Weight**: `WTSMSMPA` (Pooled Sample Weight) will be used for all analyses involving PBDEs.
