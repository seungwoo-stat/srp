# Spherical Random Projection

## Overview

1. Two files contain the functions for spherical random projection and cluster validation.
   - `functions.R` : Basic functions for operations on spheres (`DIST`, `LOG`, `EXP`, `PARALLEL`) and stability validation with spherical k-means clustering using a single core (`cluster_stability`, `perturb_accept_rate`) and multiple cores (`mccluster_stability`, `mcperturb_accept_rate`). Also, a function that validates the result of spectral clustering with multiple cores (`mccluster_stability_spec`). 
   - `spectral.R` :  spectral clustering for spherical data.
2. Three files contain the experiments on synthetic data sets.

   - `experiment-vmf3.R`: three mixtures of von Mises-Fisher distributions on 1000-dimensional sphere.

   - `experiment-vmf5.R`: five mixtures of von Mises-Fisher distributions on 2000-dimensional sphere.

   - `experiment-bulls.R`: two non spherical distributions on 1000-dimensional sphere.
3. Two files contain the experiment on real data sets.
   - `real-face`: analysis of JAFFE ([Japanese Female Facial Expression](https://doi.org/10.5281/zenodo.4029679)) data set.
   - `real-oz`: analysis of the famous Oz books corpus.



## Reference

-   Seungwoo Kang and Hee-Seok Oh. (2024) [Spherical Random Projection](https://doi.org/10.1093/jrsssb/qkae035). *Journal of the Royal Statistical Society Series B: Statistical Methodology* **86**(5), 1364-1382.
