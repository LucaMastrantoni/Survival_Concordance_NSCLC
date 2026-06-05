# Survival Tools for Randomized Controlled Trials Concordance in NSCLC

Repository for the work: **"Design-Outcome Concordance in Randomized Trials of Targeted Therapies in Non-Small Cell Lung Cancer: A Systematic Review and Meta-Research Study"** *(currently under review)*.

The web app is publicly available at: [https://lucamastrantoni.shinyapps.io/RCT_Concordance/](https://lucamastrantoni.shinyapps.io/RCT_Concordance/)

---

## Overview

This app is intended for researchers working with randomized controlled trials (RCTs) and individual patient data (IPD) or reconstructed KM curves. It allows you to:

- Evaluate **power**, **Type S error**, and **Type M exaggeration** under alternative design assumptions.
- Visualize **Kaplan-Meier curves** and fit **Cox proportional hazards models**.
- Compare **Bayesian posterior hazard ratio distributions** against the original design assumptions, and calculate **Wasserstein-1 distance**.
- Quantify **discrepancies between hazard ratios and restricted mean survival time (RMST)**-based treatment effects.

---

## Tab Description

Tabs are intended to be run in order, as information is shared across them.

### 1. Retrodesign (Single Trial)

Implementation of the retrodesign analysis as proposed by Gelman et al.

- Statistical power as a function of the assumed hazard ratio.
- Type S error (probability of sign error).
- Type M error (exaggeration factor among statistically significant results).

### 2. Survival Analysis (IPD)

Standard frequentist survival analysis using individual patient data from reconstructed KM curves. Uploaded data must include time, event status, and treatment indicator.

- Kaplan-Meier curves with confidence intervals and risk tables.
- Median survival times with 95% confidence intervals.
- Cox proportional hazards regression.
- Assessment of the proportional hazards assumption using Schoenfeld residuals.

### 3. Bayesian Cox (IPD)

Comparison of Bayesian posterior distributions of the treatment log-hazard ratio against the original design assumptions.

- Upload posterior draws from a Bayesian Cox model (log-HR scale). The app does not directly calculate posterior draws; we provide an R script (`bayesian.R`) to do so locally, requiring the same IPD input data, which produces the correct `posterior_b_treat.rds` file.
- Wasserstein distance between posterior and design distributions.
- Density plots for visual comparison.

### 4. RMST (IPD)

Restricted mean survival time (RMST) analysis.

- RMST estimated directly from Kaplan-Meier curves.
- RMST difference and RMST ratio evaluated over increasing τ.
- Quantification of mismatch between RMST under an exponential proportional hazards assumption and the observed RMST ratio from KM curves.

---

## Input Data Requirements

- **IPD files** must contain columns named: `time`, `status` (0 = censored, 1 = event), and `treat` (0 = control, 1 = experimental).
- **Bayesian posterior files** must be `.rds` objects containing a column named `b_treat`. Please use the provided R script `bayesian.R`.
- **Example files** are provided in the `Shiny App files/` folder on GitHub:
  - `ALEX_KM.xlsx`: Excel file for survival analysis (Tab 2) and RMST analysis (Tab 4).
  - `posterior_b_treat.rds`: R file for Bayesian analysis (Tab 3), obtainable using `bayesian.R`.
- This app **does not** include a method to reconstruct IPD from published KM curves. We refer to published methods.

---

## Contact

- **Correspondence**: [luca.mastrantoni01@icatt.it](mailto:luca.mastrantoni01@icatt.it)

---

## Affiliation

- Università Cattolica del Sacro Cuore, Rome, Italy
- Comprehensive Cancer Center, Fondazione Policlinico Universitario Agostino Gemelli IRCCS, Rome, Italy
