# FGM_Neighborhood
Code downloaded from https://github.com/PercyZhai/FGM_Neighborhood/tree/main.
REF: Zhao, B., Zhai, P. S., Wang, Y. S., & Kolar, M. (2024). High-dimensional functional graphical model structure learning via neighborhood selection approach. Electronic Journal of Statistics, 18(1), 1042-1129.

Note from the same author there is also a paper about differential network: 
- Zhao, B., Wang, Y. S., & Kolar, M. (2022). FuDGE: A method to estimate a functional differential graph in a high-dimensional setting. Journal of Machine Learning Research, 23(82), 1-82.

I also found another reference that do something similar to what we want to do: 
- Tugnait, J. K. (2023). Learning High-Dimensional Differential Graphs From Multiattribute Data. IEEE Transactions on Signal Processing, 72, 415-431.
  
## Folder Function
- A.prelim.func: Code to generate the precision matrix of simulation model A.
- ADMM.new: Code that define the optimization problem to solve and implement the resolution algorithm.
- B.prelim.func: Code to generate the precision matrix of simulation model B.
- C.prelim.func: Code to generate the precision matrix of simulation model C.
- D.prelim.func: Code to generate the precision matrix of simulation model D.
- FPCA.score: Function to calculate FPC score matrix from observation.
- ProxAlg_FGM: Implement Qiao's method from Qiao, X., Guo, S., & James, G. M. (2019). Functional graphical models. Journal of the American Statistical Association, 114(525), 211-222.
- auc: Given an array of TPR and FPR, calculate its AUC ROC.
- bases.func: Implement Fourier basis function definition; Evaualte Fourier, bspline, exponential and power bases in specific obs.time
- prec.rec: Given the true adj matrix and the estimated adjacency matrix calculate TP, FP, TN, FN and return precision, recall, TPR and FPR

## Folder Simulations_with_other_method_comparison
This folder contains codes running multiple methods, namely FPCA-gX, FPCA-gY, FGLasso, PSKL, and FPCA-PSKL (for Model D). Threshold is set to $\epsilon = 0$. Penalty parameter $\lambda_n$ varies. Output results are TPR and FPR for each methods. Runtime of each method under each run is also recorded. A script is given for plotting ROC for all methods under each setting.

## Folder Simulations_CV
This folder contains codes running FPCA-gX method (FGM_Neighborhood with dimentionality reduction computed on the PC of each function). Threshold varies under each setting. Penalty parameter $\lambda_n$ varies under each threshold parameter. Output results are adjacency matrices under AND and OR schemes under the optimal choice of $\lambda_n$, $t_\epsilon$ pairs that maximizes SCV-BIC. Runtime of each run is also recorded. A script is given for computing precision and recall for the optimal adjacency matrices.

## Folder Reak Data Anlysis
This folder contains codes running FPCA-gX method. Two datasets: ABIDE for autism and ADHD. Both raw data and time series extracted from raw data can be found in Data folder. Time series are extracted by "ABIDE.read.data.R" and "ADHD.read.data.R" respectively.

#### ABIDE dataset:
"SCV.autism.R" selects and outputs the optimal adjacency matrices using SCV under FPCA-gX method. "spa.ctrl.autism.R" outputs an adjacency matrix with 2% connection density under FPCA-gX method.

#### ADHD dataset:
"SCV.ADHD.R" selects and outputs the optimal adjacency matrices using SCV under FPCA-gX method. "spa.ctrl.ADHD.R" outputs an adjacency matrix with 2% connection density under FPCA-gX method. "spa.ctrl.ADHD.FGLasso.R" outputs an adjacency matrix estimated by FGLasso method.

