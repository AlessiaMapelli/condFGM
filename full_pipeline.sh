#!/bin/bash
# Full reproduction pipeline for all paper figures and results.
#
# Paper: "A neighbour selection approach for identifying differential
#         networks in conditional functional graphical models"
#
# Produced outputs
#   Figure 3  – F1 performance (single-covariate setting)
#   Figure 4  – Literature comparison (FuDGE, FGL, FFGL, FFGL2)
#   Figure 5  – Computational time comparison
#   Figure 6  – Performance across scenarios S1–S6
#   Figure 7  – EEG differential network (topographic plot)
#   Figure 9  – F1 performance (supplementary, multiple covariates)
#   Figure 10 – Estimated adjacency structure for one simulation replicate (extracted from the results Simulation_studies/Supplementary_simulations/simulation_settings/p10_n1000_n1000/seed_3/results/ combing Estimated_adjacency_matrix_node_pop.png, Estimated_adjacency_matrix_node_diff_signed_weights.png and Estimated_adjacency_matrix_node_cont_signed_weights.png)
#   Figure 11 - FPR performance (single-covariate setting)
#   Figure 12 – FPR performance (supplementary, multiple covariates)
#   Figure 13 – EEG adjacency matrices
#
# Prerequisites:
#   R >= 4.0 with all packages listed in README.md
#   Python 3 with all packages listed in README.md
#   SLURM job scheduler
#
# Usage (from the repository root):
#   bash full_pipeline.sh
#
# NOTE: Steps 1–3 submit SLURM jobs and may take hours on HPC.
# Each sub-launcher waits for its jobs to complete before returning.

# ──────────────────────────────────────────────────────────────────
# STEP 1  Single-setting simulation + literature comparison
# Intermediate results for Figures 3, 4, 5
# ──────────────────────────────────────────────────────────────────
echo "========== STEP 1: single-setting simulation =========="
bash Simulation_studies/Step1/Sbatch_simulations_luncher_Step1.sh


# ──────────────────────────────────────────────────────────────────
# STEP 2  Scenario-varying simulation (S1–S6, p=10 and p=50)
# Intermediate results for Figure 6
# ──────────────────────────────────────────────────────────────────
echo ""
echo "========== STEP 2: scenario simulations S1–S6 =========="
bash Simulation_studies/Step2/Sbatch_simulations_luncher_Step2.sh

# ──────────────────────────────────────────────────────────────────
# STEP 3  Supplementary simulations (multiple covariates)
# Intermediate results for Figures 9, 12 (Appendix B)
# ──────────────────────────────────────────────────────────────────
echo ""
echo "========== STEP 3: supplementary simulations =========="
bash Simulation_studies/Supplementary_simulations/Sbatch_simulations_luncher_Supp.sim.sh

# ──────────────────────────────────────────────────────────────────
# STEP 4  Generate all simulation figures
# Saves PNG files to Simulation_studies/ and
# Simulation_studies/Supplementary_simulations/
# Figure 10 is assembled by Plot_simulations_results_supp.R from
# intermediate images in simulation_settings/p10_n1000_n1000/seed_3/results/.
# ──────────────────────────────────────────────────────────────────
echo ""
echo "========== STEP 4: generating simulation figures =========="
Rscript Simulation_studies/Plot_simulations_results.R
Rscript Simulation_studies/Supplementary_simulations/Plot_simulations_results_supp.R

# ──────────────────────────────────────────────────────────────────
# STEP 5  EEG data application
# Saves PNG files to EEG_data_application/
# Starts from pre-computed fPCA scores (EEG_application_score_data.RData);
# Steps 1–2 of the EEG pipeline (raw data merging and fPCA extraction) are skipped.
# ──────────────────────────────────────────────────────────────────
echo ""
echo "========== STEP 5: EEG data application =========="
echo "--- Fitting the model ---"
Rscript Computation_templates/Script_sequential.R EEG_data_application/config_file.yaml

echo "--- Evaluating results ---"
Rscript Computation_templates/Results_evaluation.R EEG_data_application/config_file.yaml

echo "--- Plotting topographic network (Figure 7 and Figure 13) ---"
Rscript EEG_data_application/EEG_data_application_plots.R

# ──────────────────────────────────────────────────────────────────
echo ""
echo "Pipeline complete."
echo "Figures saved to:"
echo "  Simulation_studies/         (Figures 3, 4, 5, 6, 11)"
echo "  Simulation_studies/Supplementary_simulations/  (Figures 9, 10, 12)"
echo "  EEG_data_application/       (Figures 7, 13)"
