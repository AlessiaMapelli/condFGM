# Functional Association Network Computational Pipeline

A robust and scalable pipeline for estimating sparse association networks between functional variables, supporting both sequential and parallel processing modes.

## Overview

This pipeline performs neighbourhood selection for functional graphical model estimation. In a graphical model, nodes correspond to elements of a multivariate functional process and edges encode pairwise conditional associations. The pipeline applies a neighbourhood selection approach with efficient ADMM optimisation to identify significant functional connections and their variation with external covariates.

## Key Features

- **Sparse Association Estimation**: Applies group lasso regularisation to identify meaningful function-to-function associations
- **Covariate Support**: Incorporates subject-level covariates into the model (e.g., age or clinical status in brain connectivity analyses); if no covariates are provided, estimates a single population-level network
- **Scalable Processing**: Supports both sequential execution (single machine) and parallel execution (HPC with SLURM)
- **Selective Cross-Validation**: Automated node-wise tuning of the penalty and threshold via K-fold selective cross-validation (SCV)
- **Visualization**: Generates publication-ready adjacency matrix plots
- **Flexible Configuration**: YAML-based configuration for reproducible parameter management

---

## Repository Structure

```
Computation_templates/
├── config_template.yaml          # Annotated configuration template with all available parameters
├── config_example.yaml           # Ready-to-run configuration for the provided example data
├── algorithm_functions.R         # ADMM solver and group lasso core functions (shared library)
├── Data_preprocessing.R          # Extracts functional principal component scores from raw functional data
├── Script_sequential.R           # Single-threaded neighbourhood selection for all nodes
├── Script_sbatch_parallel.R      # Single-node neighbourhood selection for parallel HPC execution
├── Results_evaluation.R          # Assembles adjacency matrices and generates visualisations
├── Sbatch_parallel.sbatch        # SLURM job script for individual node computations
├── Sbatch_parallel_luncher.sh    # Batch job submission manager
├── Example_data.RData            # Raw example data: 400 samples × 10 variables × 100 time points
└── Example_data_score.RData      # Pre-computed fPCA scores from Example_data.RData
```

---

## System Requirements

- R (≥ 4.0)
- R packages (must be installed manually):
  - **Preprocessing**: `yaml`, `fda`, `matrixcalc`, `MASS`, `mvtnorm`, `abind`, `ggplot2`, `reshape2`, `pheatmap`
  - **Estimation**: `yaml`
  - **Results**: `yaml`, `reshape2`, `ggplot2`, `dplyr`, `stringr`, `ggnewscale` (optional: second colour scale for not-in-baseline edges; without it they are drawn in a single colour)
- SLURM workload manager (for parallel processing only)

> **Important**: all scripts must be run from the **repository root** (the folder containing `Computation_templates/`), not from inside `Computation_templates/` itself.

---

## Quick Start

### Step 1 — Configure Parameters

Copy `config_template.yaml` and adapt it to your analysis. A ready-to-run example is provided in `config_example.yaml`.

### Step 2 — Prepare Your Data

Your input should be an `.RData` file containing:

- **`discrete_fun_obs`**: a three-dimensional array of dimensions `n_samples × n_nodes × time_points`, where each entry along the first axis is a matrix `(n_nodes × time_points)` representing the functional observations for one subject.
- **`covariates_df`** *(optional)*: a data frame of subject-level covariates, with factor columns properly defined
  - Rows: subjects/samples
  - Columns: covariate values

```r
# Example data structure
load("Computation_templates/Example_data.RData")
dim(discrete_fun_obs)   # 400 samples × 10 functional variables × 100 time points
head(covariates_df)     # 400 samples × 1 covariate (group)
```

Run the preprocessing script to extract fPCA scores and save them to the path specified by `input_path` in your config:

```bash
Rscript Computation_templates/Data_preprocessing.R Computation_templates/config_example.yaml
```

The output workspace will contain:

- **`scores_df`**: fPCA scores of dimensions `n_samples × (n_nodes × n_basis_for_dim_reduction)`
  - Rows: subjects/samples
  - Columns: basis coefficients for each functional variable (named `f{node}.{basis_index}`)
- **`covariates_df`** *(optional)*: the same covariate data frame, passed through unchanged

> If you have already computed scores with a different dimensionality reduction method, you can provide them directly in this format and skip the preprocessing step.

```r
# Example preprocessed data structure
load("Computation_templates/Example_data_score.RData")
head(scores_df)     # 400 samples × 50 features (10 nodes × 5 basis functions)
head(covariates_df) # 400 samples × 1 covariate (group)
```

### Step 3 — Run the Estimation

Choose the execution mode appropriate for your setup:

**Option A — Sequential (recommended for testing or small analyses)**
```bash
Rscript Computation_templates/Script_sequential.R Computation_templates/config_example.yaml
```
Runs neighbourhood selection for all nodes one at a time on a single machine.

**Option B — Parallel on HPC with SLURM (recommended for large analyses)**
```bash
bash Computation_templates/Sbatch_parallel_luncher.sh Computation_templates/config_example.yaml
```
Submits one SLURM job per node. Before running, edit `Sbatch_parallel.sbatch` to set your email address, partition name, memory, and wall time. The launcher respects a job limit of 70 concurrent jobs by default.

### Step 4 — Assemble and Visualise Results

Once all node computations from step 3 are complete:

```bash
Rscript Computation_templates/Results_evaluation.R Computation_templates/config_example.yaml
```

---

## Configuration Parameters

All parameters are set in the YAML configuration file. See `config_template.yaml` for an annotated reference.

### Data Preprocessing

| Parameter | Description |
|---|---|
| `observed_functional_data_path` | Path to the `.RData` file containing `discrete_fun_obs` |
| `rec_basis_type` | Basis type for fitting the functional data: `"fourier"` or `"bsplines"` |
| `rec_basis_number` | Number of basis functions used to *fit* the raw functional data (e.g. 15) |
| `time_range` | Domain of the functional data, e.g. `[0, 1]` |
| `time_points` | Number of discrete time points per observation |

### Input / Output

| Parameter | Description |
|---|---|
| `input_path` | Path to the `.RData` file with pre-computed fPCA scores (output of preprocessing) |
| `output_path` | Directory where intermediate and final results will be saved |
| `name_output` | Prefix for all output files |
| `n_nodes` | Number of functional variables (nodes) in the multivariate process |
| `n_basis_for_dim_reduction` | Number of fPCA components to *retain* per node (e.g. 5); must be ≤ `rec_basis_number` |

> **`rec_basis_number` vs `n_basis_for_dim_reduction`**: the first controls how many basis functions are used to *smooth* the raw data; the second controls how many principal components to *keep* for modelling. Our choice was `rec_basis_number = 15` and `n_basis_for_dim_reduction = 5`.

### Algorithm Parameters

| Parameter | Description |
|---|---|
| `L` | Number of lambda (regularisation) values to evaluate |
| `K` | Number of folds for cross-validation |
| `thres_ctrl` | Grid of threshold values for edge inclusion (e.g. `[0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0]`) |
| `p_rand_lam` | Proportion of the lambda grid sampled at random instead of exhaustively (e.g. `0.5` = test 50% of lambda values); set to `1` for a full grid search |
| `p_rand_thr` | Same as above but for the threshold grid; set to `1` to test all thresholds |
| `scv_penalty` | SCV criterion on each test fold: `bic` (profiled-variance Gaussian BIC, $n_{test}\sum_m \log(RSS_m/n_{test}) + \log(n_{test})\,\vert\hat N\vert$, invariant to the scale of the scores) or `none` (residual sum of squares only). Default if omitted: `none` |
| `cv_folds` | Fold construction: `random` (rows randomly permuted, seeded per node, then split into K equal blocks) or `contiguous` (consecutive blocks of rows; unsuitable when rows are sorted, e.g. by group). Default if omitted: `contiguous` |
| `refit_ridge` | Ridge term $\kappa I$ of the selective refit: `relative` ($\kappa$ = `refit_ridge_const` × mean diagonal of $X^\top X$ on the training fold) or `fixed` ($\kappa$ = `refit_ridge_const`). Default if omitted: `fixed` |
| `refit_ridge_const` | Optional constant for `refit_ridge`; defaults to `0.001` (`relative`) or `0.1` (`fixed`) |
| `type` | Edge symmetrisation method: `OR` (edge present if detected in *either* direction) or `AND` (edge present only if detected in *both* directions) |
| `verbose` | `TRUE` to print progress messages; `FALSE` for silent execution |

> **SCV defaults**: when `scv_penalty`, `cv_folds` and `refit_ridge` are omitted, the pipeline uses `none` / `contiguous` / `fixed`, which reproduces the results published in the paper. The provided `config_example.yaml` and `config_template.yaml` use the recommended `bic` / `random` / `relative`.

---

## Output Files

### Per-Node Intermediate Results (Step 3)

Written to `output_path` for each node `j`:

- `{name_output}_{j}.rda` — optimal neighbourhood estimates for node `j`
- `{name_output}_{j}coeff.rda` — full coefficient matrices and Frobenius norms for node `j`

### Final Results (Step 4)

- `{name_output}_Adj_estimation.rda` — complete estimated adjacency matrices (see below)
- `{name_output}_Adjacency_matrix_node_{covariate}.png` — binary adjacency matrix plot per covariate level
- `{name_output}_Adjacency_matrix_node_diff_weights_{covariate}.png` — weighted differential matrix plot per covariate level (if covariates are present)
- `{name_output}_Adjacency_matrix_node_group_{covariate}.png` — group-specific binary adjacency matrix (factor covariates only)

The saved `.rda` file contains the following objects:

| Object | Description |
|---|---|
| `G.our.symm` | Named list of symmetrised binary adjacency matrices, one per covariate level (plus `"population"`) |
| `G.our.symm.weighted` | Named list of symmetrised weighted differential matrices for each covariate level |
| `G.our.weighted` | Named list of asymmetric weighted differential matrices for each covariate level |
| `G.our.symm.warning` | Named list of matrices flagging edges (value = 1) where the two directions gave inconsistent differential signals |
| `G.our.symm.not.in.baseline` / `G.our.not.in.baseline` | Named lists of matrices (symmetrised / asymmetric) for differential edges that are *not in the baseline* (the population block is not selected, so the variation rate is not defined): value = $\Vert\hat B^c_k\Vert_F$ (geometric mean over directions in the symmetrised matrix), `NA` elsewhere |
| `G.our.group.symm` / `G.our.group` | Named lists, one per factor dummy (e.g. `group2`), of the symmetrised / asymmetric group-specific networks: the population edges, updated where a differential edge is selected to $\Vert\hat B^0_k + \hat B^c_k\Vert_F > \varepsilon$ |

---

## Interpreting the Results

### Binary Adjacency Matrix

- **1**: a significant conditional association was detected between the two functional variables
- **0**: no detected association

### Differential Analysis (when covariates are included)

The pipeline estimates one adjacency matrix per covariate level:

- **`population`**: baseline connectivity pattern shared across all subjects
- **Covariate matrices**: incremental connectivity contribution associated with each covariate level, relative to the population baseline

The weighted differential matrices quantify the *direction and magnitude* of each connection change. The colour scale is anchored at 1 (no change from baseline):

- **Values > 1** (red in the plot): the connection is *strengthened* in this covariate group relative to the baseline
- **Values < 1** (blue in the plot): the connection is *weakened* relative to the baseline
- **Not in baseline** (white-to-dark-green scale in the plot; `NA` in the weighted matrix, value in `G.our.symm.not.in.baseline`): a differential edge is detected but the population edge is not selected, so the connection is absent in the reference group (or, for a numeric covariate, at its mean) and induced by the covariate; the colour shows its magnitude $\Vert\hat B^c_k\Vert_F$
- **`NA`** elsewhere: no differential edge detected for this pair

The symmetric weight of an edge is the geometric mean of the directional variation rates over the directions in which both the differential and the population block are selected.

---

## Citation

If you use this pipeline, please cite:

> Mapelli, A., Carini, L., Ieva, F., & Sommariva, S. (2026).
> A neighbour selection approach for identifying differential networks in conditional functional graphical models.
> arXiv:2601.02292. https://arxiv.org/abs/2601.02292

---

**Note**: The parallel execution mode is designed for HPC environments with SLURM and Singularity-based R execution. The `Sbatch_parallel.sbatch` script contains placeholder values (`YOUR_EMAIL_ADDRESS`, `cpuq`) that must be updated to match your cluster configuration before use.

For the full R and Python session information used to produce the results in the paper, see [`session_info.txt`](../session_info.txt) in the repository root.
