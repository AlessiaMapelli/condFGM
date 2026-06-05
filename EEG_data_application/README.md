# EEG Data Application for condFGM

This folder contains the real-world application of the conditional Functional Graphical Model (condFGM) method to electroencephalography (EEG) data, studying brain connectivity changes associated with alcohol use disorder. It follows the computational pipeline described in `Computation_templates/README.md`.

---

## Dataset Description

The EEG data used in this application were obtained from the repository associated with:

> Zhao, B., Wang, Y. S. & Kolar, M. (2022). FuDGE: A method to estimate a functional differential graph in a high-dimensional setting. *Journal of Machine Learning Research* 23(82), 1–82.

**FuDGE data repository**: https://github.com/boxinz17/FuDGE/tree/master/EEG/

The raw EEG recordings were originally downloaded from the UCI KDD EEG database:
https://kdd.ics.uci.edu/databases/eeg/eeg.data.html

### Input Data Files

The pre-processed source files must be downloaded directly from the FuDGE repository (`generated_data/` folder) before running this pipeline:

| File | Source | Description |
|---|---|---|
| `alco_filtered_array.Rdata` | FuDGE repo | 3D array `77 subjects × 64 electrodes × 256 time_points` — filtered EEG for the alcoholic group |
| `contrl_filtered_array.Rdata` | FuDGE repo | 3D array `45 subjects × 64 electrodes × 256 time_points` — filtered EEG for the control group |
| `Position_list.Rdata` | This repo | List of 64 electrode labels (extended 10–20 system, e.g. FP1, CZ, OZ) used for topographic plotting |

> **Note**: `alco_filtered_array.Rdata`, `contrl_filtered_array.Rdata`, and the combined `EEG_application_input_data.RData` are excluded from this repository due to their size (~30 MB total). Download the source files from the FuDGE repository before running Step 1.

---

## Repository Structure

```
EEG_data_application/
├── config_file.yaml                   # Configuration file for this application
├── Input_preprocessing.R              # Merges the two group arrays into a single input file
├── EEG_data_application_plots.R       # Generates topographic and adjacency matrix figures
├── Position_list.Rdata                # Electrode labels and scalp positions for topographic plots
└── EEG_application_score_data.RData   # Pre-computed fPCA scores (ready to use, skips preprocessing)
```

---

## System Requirements

- R (≥ 4.0)
- R packages (must be installed manually):
  - **Input preprocessing**: `yaml`, `abind`
  - **fPCA scores & estimation**: see `Computation_templates/README.md`
  - **Topographic plots**: `igraph`, `reshape2`, `ggplot2`, `dplyr`, `stringr`, `grid`

> **Important**: all scripts must be run from the **repository root** (the folder containing `EEG_data_application/` and `Computation_templates/`), not from inside `EEG_data_application/` itself.

---

## Pipeline

### Step 1 — Merge Input Data

> Skip this step if you use the pre-computed `EEG_application_score_data.RData` provided in this folder and go directly to Step 3.

After downloading `alco_filtered_array.Rdata` and `contrl_filtered_array.Rdata` from the FuDGE repository into `EEG_data_application/`, run:

```bash
Rscript EEG_data_application/Input_preprocessing.R
```

This merges the two arrays into `EEG_application_input_data.RData`, containing:
- **`discrete_fun_obs`**: 3D array of dimensions `122 subjects × 64 electrodes × 256 time_points` (45 controls first, then 77 alcoholics)
- **`covariates_df`**: data frame with column `group` (factor: `"control"` / `"alcoholic"`), dimensions `122 × 1`

### Step 2 — Extract fPCA Scores

> Skip this step if you use the pre-computed `EEG_application_score_data.RData` and go directly to Step 3.

```bash
Rscript Computation_templates/Data_preprocessing.R EEG_data_application/config_file.yaml
```

Produces `EEG_application_score_data.RData`, containing:
- **`scores_df`**: fPCA score matrix of dimensions `122 × 320` (subjects × p × M, where p = 64 electrodes, M = 5 components). Rows are subjects (controls first, then alcoholics); columns are named `f{node}.{index}` (e.g. `f1.0` … `f64.4`)
- **`covariates_df`**: same covariate data frame as above

### Step 3 — Configure Parameters

Parameters are set in `EEG_data_application/config_file.yaml`. Key settings for this application:

| Parameter | Value | Notes |
|---|---|---|
| `rec_basis_type` | `fourier` | Fourier basis to reconstruct EEG signals |
| `rec_basis_number` | `15` | Basis functions used to smooth the discrete observations |
| `n_basis_for_dim_reduction` | `5` | fPCA components retained per electrode |
| `n_nodes` | `64` | Number of EEG electrodes |
| `time_points` | `256` | Time points per trial |
| `type` | `OR` | Edge symmetrisation: present if detected in either direction |

### Step 4 — Run Neighbourhood Selection

```bash
Rscript Computation_templates/Script_sequential.R EEG_data_application/config_file.yaml
```

Writes one result file per electrode to `EEG_data_application/results/`. For faster execution on an HPC cluster, use the parallel launcher instead (see `Computation_templates/README.md`).

### Step 5 — Assemble Adjacency Matrices

```bash
Rscript Computation_templates/Results_evaluation.R EEG_data_application/config_file.yaml
```

Produces the final adjacency matrices and standard visualisation plots in `EEG_data_application/results/`.

### Step 6 — Generate Manuscript Figures

```bash
Rscript EEG_data_application/EEG_data_application_plots.R
```

Produces three figures saved to `EEG_data_application/figures/`:

| Output file | Description |
|---|---|
| `Figure_7_diff_network_AUD.pdf` | Figure 7 of the manuscript — topographic network plot with electrodes at their scalp positions; edge colours (blue → white → red) represent the log₁₀ differential weight |
| `Figure_13_Adjacency_matrix_node_diff_weighted.png` | Supplementary Figure 13 — full weighted differential adjacency matrix with electrode labels; colour scale is log₁₀ with midpoint = 0 |
| `Figure_13_Adjacency_matrix_node_population.png` | Supplementary Figure 13 — binary population adjacency matrix with electrode labels |

> **Note on colour scale**: the differential weight in these plots is transformed to log₁₀ scale (midpoint = 0), so values > 0 (red) indicate strengthened connections and values < 0 (blue) indicate weakened connections in the alcoholic group relative to controls. This differs from the standard `Results_evaluation.R` output, which uses the raw weight scale with midpoint = 1.

---

## Output Files Summary

| File | Written by | Location |
|---|---|---|
| `EEG_application_input_data.RData` | Step 1 | `EEG_data_application/` |
| `EEG_application_score_data.RData` | Step 2 | `EEG_data_application/` |
| `EEG_application_{j}.rda` | Step 4 | `EEG_data_application/results/` |
| `EEG_application_{j}coeff.rda` | Step 4 | `EEG_data_application/results/` |
| `EEG_application_Adj_estimation.rda` | Step 5 | `EEG_data_application/results/` |
| `EEG_application_Adjacency_matrix_node_*.png` | Step 5 | `EEG_data_application/results/` |
| `Figure_7_diff_network_AUD.pdf` | Step 6 | `EEG_data_application/figures/` |
| `Figure_13_Adjacency_matrix_node_diff_weighted.png` | Step 6 | `EEG_data_application/figures/` |
| `Figure_13_Adjacency_matrix_node_population.png` | Step 6 | `EEG_data_application/figures/` |

---

## Citation

If you use this application or pipeline, please cite our paper:

> Mapelli, A., Carini, L., Ieva, F., & Sommariva, S. (2026).
> A neighbour selection approach for identifying differential networks in conditional functional graphical models.
> arXiv:2601.02292. https://arxiv.org/abs/2601.02292

And the original EEG data source:

> Zhao, B., Wang, Y. S. & Kolar, M. (2022). FuDGE: A method to estimate a functional differential graph in a high-dimensional setting. *Journal of Machine Learning Research* 23(82), 1–82.

For the full R and Python session information used to produce the results in the paper, see [`session_info.txt`](../session_info.txt) in the repository root.
