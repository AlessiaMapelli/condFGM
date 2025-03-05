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

#### Output
The function returns a list with the following elements:
1. Delta_scores (Matrix) – Estimated coefficients of the multiple penlaized regression.
2. Dic_delta_scores (List of Matrices) – Adjacency matrices of the graph between the scores for each covariate.
3. Dic_delta_function (List of Matrices) – Adjacency matrices of the graph between the original function for each covariate.

#### Function Breakdown:
1. Prepares Data:
    - Check if covariates are present and standardize them.
    - Warns if the sample size is too small for reliable estimation.
2.	Optional Screening (scr = TRUE):
    - Computes the correlation matrix of the provided data.
    - If all the scores of one function have a correlation lover then the threshold (gamma) with another function, one is not going to be considered was regressing the other.
3.	Definition of the design matrix:
    - Constructs design matrix and groups for sparse group lasso.
    - Rows represent subjects (samples).
    - Columns represent features (function scores) and  interactions between the features and the covariates (function scores *pat_group)
    - The group is defined based on the covariate the interaction is linked to and the functions the scores are from.
4.	Sparse Group Multivariate Lasso Estimation (ONGOING):
    - Uses parallel computing (foreach, doParallel) to speed up matrix calculations.
    - For each function (group of scores representing the function):
        - Define the specific design matrix excluding any term including the regressed term and those excluded from the scr procedure (if scr=TRUE)
        - Fits a multivariate penalized regression model grouping by scores to determine the coefficient strengths.
        - Store the coefficient of the model in Delta_scores
5.	Constructs Output Matrices (ONGOING):
    - Assembles estimated adjacency nonsimmetric matrix for the score checking if the value of each regression coefficent is > eps.
    - Define the covariates link symmetric adjacency on the scores using either the AND or the OR condition (Dic_delta_scores).
    - Define the covariates link symmetric adjacency on the functions connecting two function is at least two of their scores are connected in  Dic_delta_scores (Dic_delta_function).

## Function name: CV_hyperparam (TO DO)
#### Purpose
This function this function define the best parameters in the sparse group lasso regression exploiting cross validation.










