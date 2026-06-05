#!/bin/bash

# Configuration
TEMPLATE_CONFIG="Simulation_studies/config_template_simulation.yaml"
YAML_MODIFIER="Simulation_studies/Helper_functions/modify_yaml_config.py"
SAVE_PATH="Simulation_studies/Step2/simulation_settings/"
JOBS_LIMIT=100
tot_iteration=10

# Simulation parameter grid
p_seq=(10 50)
n_g1=200
n_g2=200
scenario_ids=(1 2 3 4 5 6)
model_g1_seq=("cov.mat.model.red.top" "cov.mat.model" "cov.mat.model.red.top" "cov.mat.model.red.top" "cov.mat.model.diff.weights.top" "cov.mat.model")
model_g2_seq=("cov.mat.model" "cov.mat.model.red.top" "cov.mat.model.red.bottom" "cov.mat.model.red.bottom" "cov.mat.model" "cov.mat.model.diff.weights.top")
name_output="results_simulation_step2"

# Loop through the parameter grid and create configuration files
for i in "${!scenario_ids[@]}"; do
  scenario_id="${scenario_ids[i]}"
  model_g1="${model_g1_seq[i]}"
  model_g2="${model_g2_seq[i]}"

  for p in "${p_seq[@]}"; do
    if [ "$scenario_id" = "5" ]; then
      red_number=$((p / 2+1))
    else
      red_number=$((p / 3))
    fi
    simulation_name="p${p}_n${n_g1}_n${n_g2}_s${scenario_id}"

    echo "=========================================="
    echo "Starting simulation: $simulation_name"
    echo "p=$p, n_g1=$n_g1, n_g2=$n_g2, scenario=S${scenario_id}, model_g1=$model_g1, model_g2=$model_g2, red_number=$red_number"
    echo "=========================================="

    for iteration in $(seq 1 "$tot_iteration"); do
      echo "Processing iteration $iteration/$tot_iteration for $simulation_name"
      config_dir="${SAVE_PATH}${simulation_name}/seed_${iteration}"
      mkdir -p "$config_dir"
      current_config="${config_dir}/config_single_iter.yaml"

      echo "Step 1: Creating config file: $current_config"
      python3 "$YAML_MODIFIER" "$TEMPLATE_CONFIG" "$current_config" \
              "p=$p" \
              "save_path=$SAVE_PATH" \
              "n_g1=$n_g1" \
              "n_g2=$n_g2" \
              "model_g1=$model_g1" \
              "model_g2=$model_g2" \
              "red_number=$red_number" \
              "simulation_name=$simulation_name" \
              "iteration=$iteration" \
              "tot_iteration=$tot_iteration" \
              "name_output=$name_output"
    done
  done
done

echo "  Step 2: Running data simulator"

# Run the data simulator for all iterations
Rscript Simulation_studies/Full_data_simulator.R $SAVE_PATH

# Loop through the generated simulation directories and submit parallel jobs for each iteration
for sim_dir in "$SAVE_PATH"/*/; do
  for iter_sim_dir in "$sim_dir"*/; do
    echo "Step 3: Launching parallel jobs for $iter_sim_dir"
    yaml_iter_file_path="${iter_sim_dir}config_single_iter.yaml"

    # Submit the jobs to perform our algorithm for this iteration
    bash Simulation_studies/Sbatch_parallel_luncher_simulations.sh "$yaml_iter_file_path"
  done

  echo "All iterations submitted for $sim_dir"
  echo "Waiting for all jobs to complete before running results checker..."
      
  while [ "$(squeue -u $USER | wc -l)" -ge 2 ]; do
    echo "Jobs limit reached, sleeping for a while..."
    sleep 240
  done
  echo "Step 4: Evaluating results for $sim_dir"
  # Evaluate results
  Rscript Simulation_studies/Results_evaluation_simulations.R "$yaml_iter_file_path"
done