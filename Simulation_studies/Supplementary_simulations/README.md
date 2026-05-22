# Supplementary Simulations for condFGM

This folder contains the simulation pipeline for the supplementary analyses presented in our paper. These simulations evaluate the performance of the conditional Functional Graphical Model (condFGM) across a broader and more systematic range of settings, varying both the number of nodes (`p`) and sample sizes (`n`) across multiple iterations.

## Repository Structure

```
Supplementary_simulations/
├── Functions/                              # Helper functions for data generation and algorithm
├── config_template.yaml                    # Template configuration file for batch simulations
├── modify_yaml_config.py                   # Python utility to generate per-iteration config files
├── Data_simulator_full_data.R              # Generates functional data for all simulation settings
├── Step3_node_parallel.R                   # Node-level parallel algorithm execution
├── Sbatch_simulations_luncher_tests.sh     # Main launcher: config generation, data simulation, and job submission
├── Sbatch_parallel_luncher_tests.sh        # Submits node-level parallel SLURM jobs for a single config
├── Sbatch_parallel_tests.sbatch            # SLURM job definition for individual node jobs
├── Check_results_screening_procedure_tests.R  # Validates and collects results per simulation setting
└── Plot_simulations_results.R              # Results visualization
```

## Prerequisites

- R (version 4.0 or higher)
- Python 3 with `pyyaml` package
- SLURM job scheduler (for parallel execution)
- Required R packages (will be installed automatically by the scripts)

## Quick Start

The full pipeline can be launched with a single command from the `Supplementary_simulations/` directory:

```bash
bash Sbatch_simulations_luncher_tests.sh
```

This orchestrates all steps: config file generation, data simulation, parallel job submission, and results checking.

## Detailed Pipeline Instructions

### Configuration Setup

Each simulation iteration has its own configuration file generated from `config_template.yaml`. Key parameters include:

#### Data Configuration
- `save_path`: Directory where simulation results will be saved
- `simulation_name`: Automatically set to `p{p}_n{n_g1}_n{n_g2}`
- `p`: Number of nodes (graph dimension)
- `n_g1`, `n_g2`: Sample sizes for the two groups
- `red_number`: Number of edges affected by the group-specific model (set to `floor(p/3)`)
- `model_g1`, `model_g2`: Model specifications for the two groups
- `model_continuous`: Model specification for the continuous covariate component

#### Algorithm Parameters
- `rec_basis_type`: Reconstruction basis type, either `"fourier"` or `"bsplines"`
- `rec_basis_number`: Number of basis functions used for data reconstruction
- `M`: Number of functional principal components to retain
- `L`: Number of lambda values to test
- `K`: Number of folds in cross-validation optimisation
- `thres_ctrl`: Threshold values to test

### Step 1: Generate Configuration Files and Data

The main launcher `Sbatch_simulations_luncher_tests.sh` iterates over all combinations of `p` and `n` defined at **lines 13–15**, and for each combination generates per-iteration config files using `modify_yaml_config.py`.

The parameter grid is defined as:
- `p_seq`: e.g. `(10 50)`
- `n_g1_seq` / `n_g2_seq`: e.g. `(50 75 100 150 200 400 750 1000)` (equal group sizes)
- `tot_iteration`: number of independent repetitions (default: 10)

After config generation, data is simulated for all settings at once:

```bash
Rscript Data_simulator_full_data.R "$SAVE_PATH"
```

This script loops over all simulation directories under `$SAVE_PATH` and generates the functional data for each iteration using the per-iteration config files.

### Step 2: Run the Algorithm in Parallel

For each per-iteration config file, the launcher submits node-level parallel SLURM jobs:

```bash
bash Sbatch_parallel_luncher_tests.sh "$yaml_iter_file_path"
```

This dispatches one SLURM job per node (up to `p` jobs) via `Sbatch_parallel_tests.sbatch`, each running `Step3_node_parallel.R` for a specific node index. A job limit (`JOBS_LIMIT=70`) is respected to avoid overloading the scheduler.

### Step 3: Verify Job Completion and Check Results

The launcher waits for all submitted jobs to complete before evaluating results for each simulation setting:

```bash
Rscript Check_results_screening_procedure_tests.R "$yaml_iter_file_path"
```

This validates and collects the results across all iterations for the given setting.

### Step 4: Generate Figures

After all simulations and results checks are complete:

```bash
Rscript Plot_simulations_results.R
```

This script generates all figures for the supplementary material.

## Scenario Configuration

### Model Types
The simulation scenarios are defined by specifying models for the two groups (`model_g1`, `model_g2`) and for the continuous covariate component (`model_continuous`) in `config_template.yaml`. Available options are detailed in the README file located in the `Functions/` folder.

### Varying Parameters
The launcher systematically varies `p` and `n` while keeping `n_g1 = n_g2` and deriving `red_number = floor(p/3)`. To change the parameter grid, modify **lines 13–15** of `Sbatch_simulations_luncher_tests.sh`.

## Output Structure

Results are saved under the directory specified by `save_path`, organised by simulation setting and iteration:

```
results/
└── p{p}_n{n_g1}_n{n_g2}/
    └── seed_{iteration}/
        ├── config_single_iter.yaml
        └── results/
            ├── logs/
            │   └── Node_{i}.log
            └── test_results_metrics_{name_output}.csv
```

## Functions Folder

The `Functions/` folder contains utility functions for:
- Adjacency matrix simulation (`Sim.prec.matrix.func.R`)
- Basis function computation (`bases.func.R`)
- FPCA score computation (`FPCA.score.R`)
- Precision matrix reconstruction (`prec.rec.R`)

## Citation

If you use this simulation pipeline in your research, please cite our paper:

---

**Note**: This pipeline is designed for high-performance computing environments with SLURM job scheduling. The `Sbatch_parallel_tests.sbatch` file uses the `cpuq` partition and Singularity-based R execution. Modifications may be needed for other computing environments.
