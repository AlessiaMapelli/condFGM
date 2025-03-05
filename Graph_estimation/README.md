# File name: cFGGM_functions
## Function name: FGGReg_cov_estimation
#### Purpose
This function estimates a covariance network from the functional scores on a defined basis using sparse group lasso regression. 
#### Input
1. scores (Dataframe) – Functional score on a defined basis.
    - Rows represent subjects (nrow: n_samples).
    - Columns represent the scores for each of the functions considered (ncol: functions*n_basis)
2. n_basis (Numeric, Default: 1) -  Number of bases considered
3. covariates (Dataframe, Default: NULL) – Additional covariates to regress on nrow: n_samples, ncol: n_covariates)
    - If provided, numeric covariates are standardized.
    - In case of only a grouping factor pass a dataframe with one factor column. 
4. scr (Boolean, Default: TRUE) – Whether to perform correlation-based screening.
    - Speeds up computations by filtering weak interactions.
5. gamma (Numeric, Default: NULL) – Threshold for correlation screening.
    - If NULL, defaults to the 10th percentile of correlation values.
6. lambda (Numeric, Default: NULL) – Penalization term for lasso regression.
    - If NULL, cross-validation determines the optimal value.
7. lambda_type (String, Default: "1se") – Selection criterion for penalization.
    - "1se": More regularized model.
    - "min": Minimizes cross-validation error.
8. asparse (Numeric, Default: 0.75) – Relative weight of L1-norm in sparse group lasso.
9. verbose (Boolean, Default: FALSE) – Whether to print progress updates.
10. eps (Numeric, Default: 1e-08) – Small tolerance for numerical stability





