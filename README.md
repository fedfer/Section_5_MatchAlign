# MatchAlign (Section 5)

This repository demonstrates **MatchAlign**, a post-processing algorithm for Bayesian factor model output. The core steps include:

1. **Data**  
   - Uses NHANES data (2015-2016) with \(4468 \times 107\) measurements.  
   - About 36% of the data are missing, and factor models are fit to impute missing values.

2. **MSGP Prior**  
   - A Bayesian factor model with a multilevel shrinkage prior (MSGP) is fitted.  
   - Missing data are imputed within the MCMC steps.

3. **Post-Processing**  
   - Compares alignment strategies: **WOP**, **RSP-full-SA**, **RSP-partial-SA**, and **MatchAlign**.  
   - Evaluates accuracy (covariance reconstruction), effective sample size (ESS), and runtimes.

4. **Results**  
   - **MatchAlign** shows improved accuracy with higher ESS, offering a more efficient factor alignment.

5. **Usage**  
   - Main R code is provided in `main.md` along with the helper functions (`helper_fcts.R`).  
   - Simply run the code blocks in `main.md` to reproduce the analysis and plots.

For details, see the commented code in `main.md`.
