# condFGM — Conditional Functional Graphical Model

This repository provides the full codebase accompanying the paper:

> Mapelli, A., Carini, L., Ieva, F., & Sommariva, S. (2026).
> A neighbour selection approach for identifying differential networks in conditional functional graphical models.
> arXiv:2601.02292. https://arxiv.org/abs/2601.02292

The repository contains reusable computational templates for applying condFGM to any dataset, simulation studies for reproducing all paper results, and a real-data application to EEG brain connectivity data.

---

## Overview

condFGM is a method for estimating sparse conditional functional graphical models in the presence of subject-level covariates (categorical or continuous). It applies a neighbourhood selection approach with Sparse Group LASSO regularisation and ADMM optimisation to identify pairwise functional associations and their differential structure across covariate groups.

The repository is organised into three main components, each with its own detailed README:

| Folder | Contents |
|---|---|
| [`Computation_templates/`](Computation_templates/README.md) | Reusable pipeline scripts for applying condFGM to any dataset, with a toy example |
| [`EEG_data_application/`](EEG_data_application/README.md) | Application to EEG data studying brain connectivity changes in alcohol use disorder |
| [`Simulation_studies/`](Simulation_studies/README.md) | Full simulation pipeline reproducing all paper figures |

---

## Repository Structure

```
cond_fgm/
├── README.md                                   # This file
├── full_pipeline.sh                            # Master script: runs all five pipeline stages in order
│
├── Computation_templates/                      # Reusable templates to apply condFGM to any dataset
│   ├── README.md                               # Step-by-step usage guide
│   └── .....                                   
│
├── EEG_data_application/                       # Real-data application: EEG brain connectivity
│   ├── README.md                               # Dataset provenance and full pipeline guide
│   └── ..... 
│
└── Simulation_studies/                         # Full simulation pipeline
    ├── README.md                               # Pipeline guide for Step 1 and Step 2
    ├── Supplementary_simulations/              # Supplementary analyses (dual-covariate model); see README.md
    ├── Step1/                                  # Single-setting analysis with literature comparison
    ├── Step2/                                  # Performance across scenarios S1–S6
    └── .....
```

---

## Full Reproduction Pipeline

All paper results can be reproduced in sequence by running:

```bash
bash full_pipeline.sh
```

This executes five stages in order. **Steps 1–3 submit SLURM jobs and may take hours on an HPC cluster.** Each sub-launcher waits for its jobs to complete before proceeding.

| Stage | Script | Figures produced |
|---|---|---|
| 1 — Single-setting simulation + literature comparison | `Simulation_studies/Step1/Sbatch_simulations_luncher_Step1.sh` | 3, 4, 5, 11 (intermediate data) |
| 2 — Scenario simulations S1–S6 | `Simulation_studies/Step2/Sbatch_simulations_luncher_Step2.sh` | 6 (intermediate data) |
| 3 — Supplementary simulations (dual-covariate model) | `Simulation_studies/Supplementary_simulations/Sbatch_simulations_luncher_Supp.sim.sh` | 9, 10, 12 (intermediate data) |
| 4 — Generate all simulation figures | `Simulation_studies/Plot_simulations_results.R`, `Simulation_studies/Supplementary_simulations/Plot_simulations_results_supp.R` | 3, 4, 5, 6, 9, 10, 11, 12 |
| 5 — EEG data application | `Computation_templates/Script_sequential.R`, `Computation_templates/Results_evaluation.R`, `EEG_data_application/EEG_data_application_plots.R` | 7, 13 |

> **Figure 10** (estimated adjacency structure for one simulation replicate) is assembled automatically by `Simulation_studies/Supplementary_simulations/Plot_simulations_results_supp.R` from three intermediate images in `Simulation_studies/Supplementary_simulations/simulation_settings/p10_n1000_n1000/seed_3/results/`. Those images are produced by `Results_evaluation_supp_simulation.R` or can be downloaded from the OSF archive.
>
> **Figures 2 and 8** (graphical model diagrams) are hand-drawn; no generating script exists for these.

Pre-computed simulation results (including `computational_times.csv`) are available at https://osf.io/ta3bq/overview if you wish to regenerate figures without re-running the full HPC pipeline. Download all CSVs into the corresponding `simulation_settings/` directories, then run Stage 4 directly.

> **EEG shortcut**: `EEG_application_score_data.RData` (pre-computed fPCA scores) is included in the repository. `full_pipeline.sh` starts Stage 5 from these scores, skipping the raw data merging and fPCA extraction steps (Steps 1–2 of the EEG pipeline in `EEG_data_application/README.md`). To run Stage 5 from raw data instead, first complete Steps 1–2 manually.

---

## System Requirements

- R (≥ 4.0)
- Python 3 (≥ 3.8; tested with Python 3.11.5 and PyYAML 6.0.1)
- SLURM workload manager (for parallel execution of Stages 1–3 only)

See each folder's README for the full list of required R and Python packages per script.

The full R session information (R 4.4.0 on Red Hat Enterprise Linux 8.8, with all package versions) is recorded in [`session_info.txt`](session_info.txt).

> **SLURM note**: Stages 1–3 submit SLURM batch jobs and are designed for HPC clusters. They use placeholder values (`YOUR_EMAIL_ADDRESS`, `cpuq`) in the `.sbatch` files that must be updated to match your cluster before running. Stages 4–5 run locally and do not require SLURM.

---

## EEG Dataset

The EEG data in `EEG_data_application/` were obtained from the repository accompanying:

> Zhao, B., Wang, Y. S. & Kolar, M. (2022). FuDGE: A method to estimate a functional differential graph in a high-dimensional setting. *Journal of Machine Learning Research*, 23(82), 1–82.
> Repository: https://github.com/boxinz17/FuDGE/tree/master/EEG/

The raw recordings were originally published in the UCI KDD EEG database: https://kdd.ics.uci.edu/databases/eeg/eeg.data.html

The source data files (`alco_filtered_array.Rdata`, `contrl_filtered_array.Rdata`) are excluded from this repository due to size and must be downloaded directly from the FuDGE repository before running the EEG application. See [`EEG_data_application/README.md`](EEG_data_application/README.md) for instructions.

---

## Citation

If you use this repository in your research, please cite:

> Mapelli, A., Carini, L., Ieva, F., & Sommariva, S. (2026).
> A neighbour selection approach for identifying differential networks in conditional functional graphical models.
> arXiv:2601.02292. https://arxiv.org/abs/2601.02292
