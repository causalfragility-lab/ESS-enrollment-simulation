# ESS: Enrollment Strategy Simulation

R code implementing the Enrollment Strategy Simulation (ESS) framework used in:

Hait, S. (2026). A Simulation-Based Decision Framework for University Enrollment Planning Under Uncertainty.

## Description
This repository contains the Monte Carlo simulation used to analyze university marketing allocation strategies under uncertainty in cost-per-lead and conversion rates.

The simulation includes:
- 4-year discounted tuition revenue
- ROI and NPV estimation
- OLS and quantile regression analysis
- marketing allocation optimization
- policy scenario stress testing
- tuition sensitivity analysis

## How to run

Open `ESS_simulation_main.R` in R and run the script from top to bottom.

Required packages:
dplyr  
ggplot2  
quantreg  
purrr  
tidyr  
broom  
konfound

Figures and tables will be saved in the `outputs/` folder.
