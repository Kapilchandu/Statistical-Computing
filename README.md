# Statistical Computing Project: Parameter Estimation using Newton–Raphson and Bootstrap

## Overview
This project implements numerical parameter estimation using the Newton–Raphson (NR) algorithm and evaluates estimator variability through Parametric and Non-Parametric Bootstrap methods.

## Objectives
- Estimate unknown distribution parameters using Maximum Likelihood Estimation (MLE).
- Apply the Newton–Raphson iterative method to solve likelihood equations.
- Construct confidence intervals using bootstrap techniques.
- Compare Parametric and Non-Parametric Bootstrap results.

## Dataset
- Dataset File: `fort.91`
- Sample Size: 36 observations

## Methodology

### 1. Initial Parameter Estimation
Moment-based approximations were used to obtain starting values for the Newton–Raphson algorithm.

### 2. Newton–Raphson Algorithm
- Numerical optimization of the log-likelihood function.
- Convergence criterion:

```text
max|θ(k+1) − θ(k)| < 10⁻⁶
