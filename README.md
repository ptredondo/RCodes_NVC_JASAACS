# RCodes_NVC

This repository contains the R codes needed to reproduce key results of the paper entitled "Nonlinear Coherence for Vector Time Series: Defining Region-to-Region Functional Brain Connectivity" by Paolo Victor Redondo, Raphaël Huser and Hernando Ombao.

The implementation of the proposed nonlinear vector coherence (NVC) measure and other relevant functions is available through the R package [nvc](https://github.com/ptredondo/nvc).

## 01 EEG Analysis.R

The contents of `01 EEG Analysis.R` are as follows:

1. Loading the required R packages.
2. Loading and preprocessing EEG data and defining the regions of interest (ROIs).
3. EEG analysis using the proposed NVC framework, including:
   - Estimation of region-to-region functional connectivity.
   - Testing for differences in connectivity across patient groups (AD, CN and FTD).
   - Sensitivity analysis on the choice of sampling rate, standardization, and block length.
4. EEG analysis using pairwise band coherence (PBC), Kendall coherence (KenCoh), and the nonspectral T measure.
5. Estimation and group comparisons of relative band power (RBP).
6. Visualization of functional connectivity results, including the results of the sensitivity analyses.
7. Classification of CN vs. AD and CN vs. FTD using functional connectivity and relative band power as features, including:
   - Support vector machines (SVM) with different kernels.
   - k-nearest neighbors (k-NN).
   - Random forests (RF).
   - Evaluation of classification performance using leave-one-subject-out cross-validation.

The EEG dataset is publicly available through [OpenNeuro (ds004504, version 1.0.9)](https://openneuro.org/datasets/ds004504/versions/1.0.9).

## 02 Numerical Experiments.R

The contents of `02 Numerical Experiments.R` are as follows:

1. Loading the required R packages.
2. Defining user-defined functions for simulating multivariate data and nonlinear spectral dependence.
3. Numerical experiments under linear spectral dependence, including:
   - Simulation of multivariate time series under five dependence scenarios.
   - Estimation of NVC and competing coherence measures.
   - Evaluation of the empirical size and power of the proposed permutation test.
   - Visualization and comparison of estimation results.
4. Numerical experiments for general random vectors, including:
   - Simulation under independence, linear, quadratic, and sinusoidal dependence.
   - Consideration of different copula families, marginal distributions, sample sizes, and vector dimensions.
   - Evaluation of the empirical size and power of the permutation-invariant T measure.
5. Numerical illustrations under nonlinear spectral vector dependence, including:
   - Simulation of multivariate time series with nonlinear frequency-specific dependence.
   - Comparison of NVC with competing coherence measures.
6. Sensitivity analysis on the choice of the number of response-component permutations used to estimate the permutation-invariant NVC measure.
