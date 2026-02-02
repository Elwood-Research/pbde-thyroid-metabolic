# Requirements: Data and Variables

## NHANES Cycles
- 2007-2008 (Cycle E)
- 2009-2010 (Cycle F)
- 2011-2012 (Cycle G)

## Datasets and Variables

### 1. Exposure: PBDEs (Pooled Samples)
- **Dataset**: `BFRPOL_E`, `BFRPOL_F`, `BFRPOL_G`
- **Linkage**: `POOLTF_E`, `POOLTF_F`, `POOLTF_G` (links SEQN to SAMPLEID)
- **Variables**:
  - `LBCBR3LA`: 2,2',4,4'-Tetrabromodiphenyl ether Lipid Adjusted (PBDE 47)
  - `LBCBR5LA`: 2,2',4,4',5-Pentabromodiphenyl ether Lipid Adjusted (PBDE 99)
  - `LBCBR7LA`: 2,2',4,4',5,5'-Hexabromodiphenyl ether Lipid Adjusted (PBDE 153)
  - `LBCBR11L`: Decabromodiphenyl ether Lipid Adjusted (PBDE 209)

### 2. Outcomes: Thyroid Hormones
- **Dataset**: `THYROD_E`, `THYROD_F`, `THYROD_G`
- **Variables**:
  - `LBXTSH1`: Thyroid Stimulating Hormone (TSH)
  - `LBXT4F`: Free Thyroxine (Free T4)
  - `LBXT3F`: Free Triiodothyronine (Free T3)

### 3. Outcomes: Lipids
- **Datasets**: `TCHOL`, `HDL`, `TRIGLY` (Cycles E, F, G)
- **Variables**:
  - `LBXTC`: Total Cholesterol
  - `LBDHDD`: HDL-Cholesterol
  - `LBDLDL`: LDL-Cholesterol
  - `LBXTR`: Triglycerides

### 4. Outcomes: Insulin & Glucose
- **Dataset**: `GLU_E`, `GLU_F`, `GLU_G`
- **Variables**:
  - `LBXGLU`: Fasting Glucose
  - `LBXIN`: Serum Insulin
  - *Calculated*: HOMA-IR = (Glucose * Insulin) / 405

### 5. Covariates
- **Demographics** (`DEMO`): `RIDAGEYR` (Age), `RIAGENDR` (Sex), `RIDRETH1`/`RIDRETH3` (Race/Ethnicity)
- **Body Measures** (`BMX`): `BMXBMI` (BMI)
- **Smoking** (`COT`): `LBXCOT` (Cotinine)
