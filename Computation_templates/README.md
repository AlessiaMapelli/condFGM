# README

## Brain Network Analysis Computational Pipeline
A robust and scalable computational pipeline estimates sparse brain networks and supports both sequential and parallel processing modes.

#### Overview
This pipeline performs neighborhood selection for brain network analysis, estimating sparse connectivity patterns between brain regions. It uses a neighboorhood approach  with efficient optimization algorithms (ADMM) to identify significant connections while controlling for confounding variables.

#### Key Features

- **Sparse Network Estimation**: Uses Sparse Group LASSO regularization to identify meaningful brain connections
- **Covariate Support**: Incorporates demographic, clinical, or other covariates into the analysis
- **Scalable Processing**: Supports both sequential and parallel execution modes
- **Cross-Validation**: Automated parameter tuning using K-fold cross-validation
- **Visualization**: Generates publication-ready adjacency matrix plots
- **Flexible Configuration**: YAML-based configuration for easy parameter management

#### Pipeline Components
###### Core Scripts
| Script | Purpose | Execution Mode |
|--------|---------|----------------|
| `Script_sequential.R` | Single-threaded analysis of all brain nodes | Sequential |
| `Script_sbatch_parallel.R` | Single-node analysis for parallel processing within an HPC environment | Parallel |
| `Results_evaluation.R` | Post-processing and visualization | Post-analysis |

###### Configuration Files

- `config_file.yaml` - Main configuration template

###### Job Management (HPC)

- `Sbatch_parallel.sbatch` - SLURM job script for individual nodes
- `Sbatch_parallel_luncher.sh` - Batch job submission manager

#### System Requirements - Software Dependencies
- R (≥ 4.0)
- R packages: `yaml`
- SLURM (for parallel processing)

#### Quick Start
###### 1. Prepare Your Data
Your input data should be an R workspace (.RData) file containing:

- **`scores_df`**: Functional scores matrix coming from some dimentionality reduction of the function onto a functional basis (e.g., scores from fPCA) (n_samples × n_features)
  - Rows: subjects/samples
  - Columns: brain region features (n_nodes × n_basis_functions)
- **`covariates_df`** (optional): Covariate matrix with proper factor definitions
    - Rows: subjects/samples
    - Columns: covariates values

```r
# Example data structure
load(example1.RData)
scores_df     # 400 samples × 50 features (10 nodes × 5 basis functions)
covariates_df # 400 samples × 1 covariate (group)
```

###### 2. Configure Parameters
Copy and modify `config_file.yaml`:

```yaml
# Essential Parameters
input_path: "your_data.RData"
output_path: "results/"
name_output: "my_analysis"
n_nodes: 10                    # Number of brain nodes
n_basis_for_dim_reduction: 5   # Basis functions per node

# Model Parameters
L: 100                         # Number of lambda values to be tested
K: 5                          # K-fold cross-validation
thres_ctrl: [0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0]  # Threshold grid
p_rand_lam: 0.5               # Proportion of lambda values to test
p_rand_thr: 1                 # Proportion of thresholds to test
type: OR                      # Edge symmetrization method. Options: OR (Edge exists if detected in either direction), AND (Edge exists only if detected in both directions)

verbose: TRUE                 #Enable progress messages. Options: TRUE (Recommended for monitoring), FALSE
```

###### 3. Choose Execution Mode
**Option A: Sequential Processing (Recommended for Testing)**
```bash
# Run complete analysis on single machine
Rscript Script_sequential.R config_file.yaml
```

**Option B: Parallel Processing (Recommended for Large Analyses)**
```bash
# Submit parallel jobs (requires SLURM)
bash Sbatch_parallel_luncher.sh config_file.yaml
```

###### 4. Process Results

```bash
# Generate final adjacency matrices and visualizations
Rscript Results_evaluation.R config_file.yaml
```

## Output Files

### Individual Node Results
- `{name_output}_{node}.rda` - Optimal neighborhoods for each node
- `{name_output}_{node}coeff.rda` - Detailed coefficients and statistics

### Final Results
- `{name_output}_Adj_estimation.rda` - Complete adjacency matrices
- `{name_output}_Adjacency_matrix_node_*.png` - Visualization plots

### Visualization Types
- **Binary adjacency matrices**: Show presence/absence of connections
- **Weighted differential matrices**: Show strengthened/weakened connections by covariates


## Interpreting Results

### Adjacency Matrices
- **Binary matrices**: 1 indicates significant connection, 0 indicates no connection

### Differential Analysis
When covariates are included:
- **Population matrix**: Base connectivity pattern
- **Covariate matrices**: Additional connections for each covariate level
- **Weighted matrices**: Show strengthened (red) vs. weakened (blue) connections

## Citation

If you use this pipeline in your research, please cite:

```
[Your publication details here]
```

---
