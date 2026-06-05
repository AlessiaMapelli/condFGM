# Supplementary Simulations for condFGM

This folder contains the simulation pipeline for the supplementary analyses presented in our paper. These simulations extend the main studies by evaluating condFGM in a setting with **both a categorical group covariate and a continuous covariate**, varying the number of nodes (`p`) and sample sizes (`n`) across a wider range than the main Step 1 and Step 2 analyses.

---

## Repository Structure

```
Supplementary_simulations/
├── config_template_supp_simulation.yaml    # Template configuration file for batch simulations
├── Full_data_simulator_supp_simulation.R   # Generates functional data for all simulation settings
├── Results_evaluation_supp_simulation.R    # Assembles results and computes F1/TPR/FPR per setting
├── Plot_simulations_results_supp.R         # Generates supplementary manuscript figures (Figures 9, 10, 12)
└──  Sbatch_simulations_luncher_Supp.sim.sh  # Main launcher: config generation, data simulation, job submission, and results evaluation
```

> **Shared scripts**: the node-level parallel estimation and SLURM job submission reuse the shared scripts from `Simulation_studies/`:
> - `Simulation_studies/Sbatch_parallel_luncher_simulations.sh` — submits one SLURM job per node
> - `Simulation_studies/Script_sbatch_parallel_simulations.R` — node-level estimation script
> - `Simulation_studies/Helper_functions/modify_yaml_config.py` — per-iteration config generator

---

## System Requirements

- R (≥ 4.0)
- Python 3
- SLURM workload manager
- R packages (must be installed manually):

| Script | Packages |
|---|---|
| `Full_data_simulator_supp_simulation.R` | `yaml`, `fda`, `matrixcalc`, `MASS`, `mvtnorm`, `abind`, `ggplot2`, `reshape2`, `pheatmap` |
| `Script_sbatch_parallel_simulations.R` (shared) | `yaml` (plus packages from `Computation_templates/algorithm_functions.R`) |
| `Results_evaluation_supp_simulation.R` | `yaml`, `reshape2`, `ggplot2`, `dplyr`, `stringr` |
| `Plot_simulations_results_supp.R` | `tidyverse`, `readr`, `png`, `grid`, `patchwork` |

- Python packages (must be installed manually):

| Script | Packages |
|---|---|
| `Helper_functions/modify_yaml_config.py` (shared) | `pyyaml` |

> **Important**: all scripts must be run from the **repository root** (the folder containing `Simulation_studies/`), not from inside any subfolder.

---

## Quick Start

```bash
bash Simulation_studies/Supplementary_simulations/Sbatch_simulations_luncher_Supp.sim.sh
```

This orchestrates all stages: config file generation, data simulation, parallel job submission, and results evaluation.

The complete pre-computed simulation results are available at https://osf.io/ta3bq/overview for reproducibility without re-running the full pipeline.

---

## Detailed Pipeline Instructions

The launcher follows the same four-stage structure as the main simulation studies: config generation → data simulation → parallel job submission → results evaluation.

**Parameter grid** (defined at lines 17–19 of `Sbatch_simulations_luncher_Supp.sim.sh`):
- `p_seq`: `(10 50)`
- `n_g1_seq` / `n_g2_seq`: `(50 75 100 150 200 400 750 1000)` (equal group sizes)
- `tot_iteration`: 10

### Stage 1 — Generate Per-Iteration Configuration Files

For each combination of `p` and `n`, `modify_yaml_config.py` generates a dedicated `config_single_iter.yaml` from `config_template_supp_simulation.yaml`, one per iteration, saved to `Simulation_studies/Supplementary_simulations/simulation_settings/p{p}_n{n}_n{n}/seed_{iteration}/`.

### Stage 2 — Simulate Functional Data

```bash
Rscript Simulation_studies/Supplementary_simulations/Full_data_simulator_supp_simulation.R "$SAVE_PATH"
```

Generates functional datasets for all settings and iterations at once. Unlike the main simulations, data are generated under a model with **two covariates**: a categorical group factor and a continuous covariate (controlled by `model_continuous` and `var_cont_type` in the config). The original simulated data are saved to `Simulation_studies/Supplementary_simulations/simulation_settings/p{p}_n{n}_n{n}/seed_{iteration}/simulated_data` with plots of the simulated adjacency matrix and original data. The preprocessed simulated data are saved to `Simulation_studies/Supplementary_simulations/simulation_settings/p{p}_n{n}_n{n}/seed_{iteration}/`. 

### Stage 3 — Run condFGM in Parallel

For each per-iteration config, submits one SLURM job per node using the shared launcher:

```bash
bash Simulation_studies/Sbatch_parallel_luncher_simulations.sh "$yaml_iter_file_path"
```

Each job runs `Simulation_studies/Script_sbatch_parallel_simulations.R` for a specific node index; results are stored in `Simulation_studies/Supplementary_simulations/simulation_settings/p{p}_n{n}_n{n}/seed_{iteration}/results`. The launcher waits for all jobs to complete before moving to the next setting.

### Stage 4 — Evaluate Results

Once jobs complete for a given setting:

```bash
Rscript Simulation_studies/Supplementary_simulations/Results_evaluation_supp_simulation.R "$yaml_iter_file_path"
```

Per-setting summary CSVs are written to `Simulation_studies/Supplementary_simulations/simulation_settings/p{p}_n{n}_n{n}/`:
- `Results_performance_metrics_results_supp_simulations.csv` — F1, TPR, FPR for condFGM

Additional inspectable plots of estimated results are available in `Simulation_studies/Supplementary_simulations/simulation_settings/p{p}_n{n}_n{n}/seed_{iteration}/results`.

### Step 5 — Generate Figures

After all settings are evaluated:

```bash
Rscript Simulation_studies/Supplementary_simulations/Plot_simulations_results_supp.R
```

Generates Figures 9, 10, and 12 of the manuscript (Appendix B.1.2 and B.1.3), saved to `Simulation_studies/Supplementary_simulations/figures/`. Figure 10 requires the seed_3 results for setting `p10_n1000_n1000` to be present (produced by Stage 4 above, or downloaded from https://osf.io/ta3bq/overview).

---

## Configuration Parameters

Per-iteration config files are generated automatically from `config_template_supp_simulation.yaml` by `modify_yaml_config.py`. The table below lists parameters that are **externally set** by the launcher (do not edit these in the template) and those that can be **freely changed**.

| Parameter | Set by | Description |
|---|---|---|
| `save_path` | Launcher | Base path for results; actual output goes to `save_path/simulation_name/` |
| `simulation_name` | Launcher | Auto-set to `p{p}_n{n_g1}_n{n_g2}` |
| `iteration` | Launcher | Iteration index; sets random seeds for reproducibility |
| `tot_iteration` | Launcher | Total number of iterations; controls result collection |
| `p` | Launcher | Number of nodes in the multivariate functional process |
| `n_g1`, `n_g2` | Launcher | Sample sizes for groups 1 and 2 |
| `red_number` | Launcher | Number of nodes with group-specific connectivity (`floor(p/3)`) |
| `model_g1`, `model_g2` | **User** | Covariance model for the categorical group factor (see `Helper_functions/Sim.prec.matrix.func.R`) |
| `model_continuous` | **User** | Covariance model for the continuous covariate component |
| `var_cont_type` | **User** | Variance type for the continuous component: `"pos_cont"` or `"cont"` |
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

For the full R and Python session information used to produce the results in the paper, see [`session_info.txt`](../../session_info.txt) in the repository root.
