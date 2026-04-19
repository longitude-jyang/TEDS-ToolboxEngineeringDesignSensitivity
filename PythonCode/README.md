# Robust FIM & KPI-Based Sensitivity Analysis Framework

This repository contains a Python-based computational framework for performing Monte Carlo simulations, estimating Joint Probability Density Functions (jPDF), calculating the Fisher Information Matrix (FIM), and conducting KPI-based sensitivity analysis. 

Originally translated from MATLAB, this tool has been heavily optimized for Python. It features a fully vectorized mathematical backend specifically designed to prevent memory overloads and hardware-level segmentation faults on modern Apple Silicon (M1/M2/M3) chips.

## 🚀 Key Features
* **Monte Carlo Simulation:** Dynamically samples random variables based on user-defined nominal values and Coefficients of Variation (CoV).
* **N-Dimensional Histogram & jPDF:** Custom implementation of MATLAB's `histcn` and `accumarray` to calculate joint probability density functions and sensitivities.
* **Fisher Information Matrix (FIM):** Evaluates parameter sensitivity and information gain. Includes an Eigen Analysis pipeline to identify the most critical parameter combinations.
* **Failure Probability (KPI) Sensitivity:** Calculates unconditional failure probabilities ($P_f$) against user-defined thresholds (absolute or percentile) and computes the sensitivity of these failures to input parameters.
* **Apple Silicon (M-Series) Safe:** The `cal_jFisher` function has been completely vectorized. It replaces hardware-crashing `for`-loops and $0/0$ division errors with memory-safe broadcasting and contiguous array masking, allowing it to run flawlessly on Apple's Accelerate framework.

## 📂 File Structure

* `main.ipynb`
  The primary Jupyter Notebook. This is the entry point of the application. It initializes the random variables, configures the simulation options, and sequentially executes the Monte Carlo simulation, FIM estimation, Eigen Analysis, and KPI evaluations.
* `utils_TEDS.py`
  The core utility module containing all the heavy lifting. It includes:
  * Variable sampling and transformation (`parList`, `parSampling`, `parTran`)
  * The Blackbox System Model (`design_B3`, `solve4frf_B3`)
  * Statistical tools (`cal_jpdf_hist`, `histcn`)
  * Matrix calculations (`cal_jFisher`, `calSen_KPI`)
  * Data Visualization (`display_jpdf_design_B3`, `display_kpi_sensitivity`)

## 🛠️ Dependencies
To run this project, you will need Python 3.x and the following standard scientific libraries:
```bash
pip install numpy scipy matplotlib jupyter
