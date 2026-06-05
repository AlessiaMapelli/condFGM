# Simulation Studies for condFGM

This folder contains the simulation pipeline for reproducing the results presented in our paper. The simulation studies evaluate the performance of the conditional Functional Graphical Model (condFGM) across different settings and scenarios, and compare it with competing algorithms proposed in "Zhao, B., Wang, Y. S., & Kolar, M. (2022). FuDGE: A method to estimate a functional differential graph in a high-dimensional setting. Journal of Machine Learning Research, 23(82), 1-82."

---

## Repository Structure

```
Simulation_studies/
├── config_template_simulation.yaml         # Template configuration file for batch simulations
├── Full_data_simulator.R                   # Generates functional data for all simulation settings
├── Script_sbatch_parallel_simulations.R    # Node-level estimation script (called by SLURM)
├── Sbatch_parallel_luncher_simulations.sh  # Submits node-level SLURM jobs for one config file
├── Results_evaluation_simulations.R        # Assembles results and computes F1/TPR/FPR per setting
├── Plot_simulations_results.R              # Generates all manuscript figures (Figures 3–6, 11)
├── Helper_functions/                       # Shared utility functions and YAML config generator
│   ├── modify_yaml_config.py               # Python utility to generate per-iteration config files
│   ├── Sim.prec.matrix.func.R              # Adjacency matrix simulation functions
│   ├── FPCA.score.R                        # FPCA score computation
│   ├── bases.func.R                        # Basis function computation
│   ├── prec.rec.R                          # Precision/recall evaluation
│   └── plot_binary_adj_matrix.R            # Binary adjacency matrix plotting
├── Step1/                                  # Single-setting analysis with literature comparison
│   ├── Sbatch_simulations_luncher_Step1.sh # Main launcher for Step 1
│   ├── Process_logs_comp_time.R            # Extracts computational times from SLURM logs
│   └── FuDGE/                              # Literature comparison scripts
│       ├── Literature_comparison.R        # Runs FuDGE and other competing methods
│       ├── Results_evaluation_sim_litt_comp.R  # Evaluates competing methods results
│       ├── Sbatch_litt_comp.sbatch         # SLURM job script for literature comparison
│       └── Function_Definitions/           # FuDGE helper functions downloaded directly from https://github.com/boxinz17/FuDGE/tree/master/Function_Definitions
├── Step2/                                  # Performance across different scenarios
│   ├── Sbatch_simulations_luncher_Step2.sh # Main launcher for Step 2
└── Supplementary_simulations/             # See Supplementary_simulations/README.md
```

---

## System Requirements

- R (≥ 4.0)
- Python 3
- SLURM workload manager
- R packages (must be installed manually):

| Script | Packages |
|---|---|
| `Full_data_simulator.R` | `yaml`, `fda`, `matrixcalc`, `MASS`, `mvtnorm`, `abind`, `ggplot2`, `reshape2`, `pheatmap` |
| `Script_sbatch_parallel_simulations.R` | `yaml` (plus packages from `Computation_templates/algorithm_functions.R`) |
| `Results_evaluation_simulations.R` | `yaml`, `reshape2`, `ggplot2`, `dplyr`, `stringr` |
| `Plot_simulations_results.R` | `tidyverse`, `readr`, `latex2exp` |
| `Step1/FuDGE/Literature_comparison.R` | `yaml`, `reshape2`, `ggplot2`, `dplyr`, `Matrix`, `MASS`, `matrixcalc`, `fda`, `quadprog`, `foreach`, `doParallel`, `JGL` |

- Python packages (must be installed manually):

| Script | Packages |
|---|---|
| `Helper_functions/modify_yaml_config.py` | `pyyaml` |

> **Important**: all scripts must be run from the **repository root** (the folder containing `Simulation_studies/`), not from inside any subfolder.

---

## Quick Start

Run all Step 1 simulations (varying `p` and `n`, with literature comparison):

```bash
bash Simulation_studies/Step1/Sbatch_simulations_luncher_Step1.sh
```

Run all Step 2 simulations (varying scenarios and `p`, with fixed  `n`):

```bash
bash Simulation_studies/Step2/Sbatch_simulations_luncher_Step2.sh
```

Generate all simulation related manuscript figures once both steps are complete:

```bash
Rscript Simulation_studies/Plot_simulations_results.R
```

The complete pre-computed simulation results (including `computational_times.csv`) are available at https://osf.io/ta3bq/overview for reproducibility without re-running the full pipeline.

---

## Detailed Pipeline Instructions

Both launchers follow the same four-stage structure: config generation → data simulation → parallel job submission → results evaluation. Step 1 additionally runs competing methods and extracts computational times.

### Step 1 — Single-Setting Analysis and Literature Comparison

Produces **Figures 3, 4, and 5** of the manuscript. Runs condFGM and competing algorithms across a grid of node counts `p` and sample sizes `n`, for 10 independent iterations each.

**Parameter grid** (defined at lines 20–22 of `Sbatch_simulations_luncher_Step1.sh`):
- `p_seq`: `(10 15 25 50)`
- `n_g1_seq` / `n_g2_seq`: `(50 75 100 150 200)` (equal group sizes)
- `tot_iteration`: 10

#### Stage 1.1 — Generate Per-Iteration Configuration Files

For each combination of `p` and `n`, `modify_yaml_config.py` generates a dedicated `config_single_iter.yaml` from `config_template_simulation.yaml`, one per iteration, saved to `Simulation_studies/Step1/simulation_settings/p{p}_n{n}_n{n}/seed_{iteration}/`.

#### Stage 1.2 — Simulate Functional Data

```bash
Rscript Simulation_studies/Full_data_simulator.R "$SAVE_PATH"
```

Generates the functional datasets for all settings and iterations at once. The original simulated data are saved to `Simulation_studies/Step1/simulation_settings/p{p}_n{n}_n{n}/seed_{iteration}/simulated_data` with plots of the simulated adjacency matrix and original data. The preprocessed simulated data are saved to `Simulation_studies/Step1/simulation_settings/p{p}_n{n}_n{n}/seed_{iteration}/`. 

#### Stage 1.3 — Run condFGM and Competing Methods in Parallel

For each per-iteration config, submits:
- One SLURM job per node for **condFGM** via `Simulation_studies/Sbatch_parallel_luncher_simulations.sh`; results are stored in `Simulation_studies/Step1/simulation_settings/p{p}_n{n}_n{n}/seed_{iteration}/results`.
- One SLURM job for the **literature comparison** (FuDGE and other competing methods) via `Simulation_studies/Step1/FuDGE/Sbatch_litt_comp.sbatch`; results are stored in `Simulation_studies/Step1/simulation_settings/p{p}_n{n}_n{n}/seed_{iteration}/results_lit_comparison`.

The launcher waits for all jobs to finish before moving to the next setting.

#### Stage 1.4 — Evaluate Results

Once jobs complete for a given setting:

```bash
Rscript Simulation_studies/Results_evaluation_simulations.R "$yaml_iter_file_path"
Rscript Simulation_studies/Step1/FuDGE/Results_evaluation_sim_litt_comp.R "$yaml_iter_file_path"
```

Per-setting summary CSVs are written to `Simulation_studies/Step1/simulation_settings/p{p}_n{n}_n{n}/`:
- `Results_performance_metrics_results_simulation_step1.csv` — F1, TPR, FPR for condFGM
- `Litt_comp_results__performance_metrics_results_simulation_step1.csv` — same metrics for competing methods

Additional inspectable plots of estimated results are available in `Simulation_studies/Step1/simulation_settings/p{p}_n{n}_n{n}/seed_{iteration}/results` and `Simulation_studies/Step1/simulation_settings/p{p}_n{n}_n{n}/seed_{iteration}/results_lit_comparison`.

#### Stage 1.5 — Extract Computational Times

After all settings are evaluated:

```bash
Rscript Simulation_studies/Step1/Process_logs_comp_time.R "$SAVE_PATH"
```

Parses SLURM log files to extract and compare wall-clock times for condFGM and FuDGE, and stores the results in `Simulation_studies/Step1/simulation_settings/computational_times.csv`.

---

### Step 2 — Performance Across Different Scenarios

Produces **Figure 6** of the manuscript. Fixes `p ∈ {10, 50}` and `n_g1 = n_g2 = 200`, and runs 6 scenarios varying the covariance model for the two groups, for 10 independent iterations each.

**Scenario definitions** (defined at lines 18–20 of `Sbatch_simulations_luncher_Step2.sh`):

| Scenario | `model_g1` | `model_g2` | `red_number` |
|---|---|---|---|
| S1 | `cov.mat.model.red.top` | `cov.mat.model` | `floor(p/3)` |
| S2 | `cov.mat.model` | `cov.mat.model.red.top` | `floor(p/3)` |
| S3 | `cov.mat.model.red.top` | `cov.mat.model.red.bottom` | `floor(p/3)` |
| S4 | `cov.mat.model.red.top` | `cov.mat.model.red.bottom` | `floor(p/3)` |
| S5 | `cov.mat.model.diff.weights.top` | `cov.mat.model` | `floor(p/2)+1` |
| S6 | `cov.mat.model` | `cov.mat.model.diff.weights.top` | `floor(p/3)` |

Available model names are the functions defined in `Helper_functions/Sim.prec.matrix.func.R`.

The pipeline stages are the same as Step 1 (config generation → data simulation → parallel jobs → results evaluation), without the literature comparison. Per-setting summary CSVs are written to `Simulation_studies/Step2/simulation_settings/p{p}_n{n}_n{n}_s{scenario}/`:
- `Results_performance_metrics_results_simulation_step2.csv` — F1, TPR, FPR for condFGM

---

### Step 3 — Generate Figures

After both steps are complete:

```bash
Rscript Simulation_studies/Plot_simulations_results.R
```

Generates Figures 3, 4, 5, 6, and 11 of the manuscript, saved to `Simulation_studies/figures/`.

---

## Configuration Parameters

Per-iteration config files are generated automatically from `config_template_simulation.yaml` by `modify_yaml_config.py`. The table below lists parameters that are **externally set** by the launcher (do not edit these in the template) and those that can be **freely changed**.

| Parameter | Set by | Description |
|---|---|---|
| `save_path` | Launcher | Base path for results; actual output goes to `save_path/simulation_name/` |
| `simulation_name` | Launcher | Auto-set to `p{p}_n{n_g1}_n{n_g2}` (Step 1) or `p{p}_n{n}_n{n}_s{id}` (Step 2) |
| `iteration` | Launcher | Iteration index; sets random seeds for reproducibility |
| `tot_iteration` | Launcher | Total number of iterations; controls result collection |
| `p` | Launcher | Number of nodes in the multivariate functional process |
| `n_g1`, `n_g2` | Launcher | Sample sizes for groups 1 and 2 |
| `model_g1`, `model_g2` | Launcher (Step 2) / Template (Step 1) | Covariance model for each group |
| `red_number` | Launcher | Number of nodes with group-specific connectivity (`floor(p/3)` by default) |
| `rec_basis_type` | **User** | Basis type for data reconstruction: `"fourier"` or `"bsplines"` |
| `rec_basis_number` | **User** | Number of basis functions for smoothing (default: 15) |
| `n_basis_for_dim_reduction` | **User** | fPCA components retained per node (default: 5) |
| `L` | **User** | Number of lambda values to evaluate (default: 100) |
| `K` | **User** | Number of cross-validation folds (default: 5) |
| `thres_ctrl` | **User** | Threshold grid for edge inclusion |
| `seed_g1`, `seed_g2` | **User** | Random seeds for data generation |

---

## Citation

If you use this simulation pipeline in your research, please cite our paper:

> Mapelli, A., Carini, L., Ieva, F., & Sommariva, S. (2026).
> A neighbour selection approach for identifying differential networks in conditional functional graphical models.
> arXiv:2601.02292. https://arxiv.org/abs/2601.02292

---

**Note**: This pipeline is designed for high-performance computing environments with SLURM job scheduling. Modifications may be needed for other computing environments.

For the full R and Python session information used to produce the results in the paper, see [`session_info.txt`](../session_info.txt) in the repository root.
