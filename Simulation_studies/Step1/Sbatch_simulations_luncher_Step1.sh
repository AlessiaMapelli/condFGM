#!/bin/bash

# Configuration
TEMPLATE_CONFIG="Simulation_studies/config_template_simulation.yaml"
YAML_MODIFIER="Simulation_studies/Helper_functions/modify_yaml_config.py"
SAVE_PATH="Simulation_studies/Step1/simulation_settings/"
tot_iteration=10

LIT_COMP_SCRIPT="Simulation_studies/Step1/FuDGE/Literature_comparison.R"
CURRENT_FOLDER=$(pwd)

# Simulation parameter grid
p_seq=(10 15 25 50)
n_g1_seq=(50 75 100 150 200)
n_g2_seq=(50 75 100 150 200)  # Assuming same values as n_g1

# Loop through the parameter grid and create configuration files
for p in "${p_seq[@]}"; do
  for n_g1 in "${n_g1_seq[@]}"; do
    #for n_g2 in "${n_g2_seq[@]}"; do
      n_g2=$n_g1
      red_number=$((p / 3))
      simulation_name="p${p}_n${n_g1}_n${n_g2}"

      echo "=========================================="
      echo "Starting simulation: $simulation_name"
      echo "p=$p, n_g1=$n_g1, n_g2=$n_g2, red_number=$red_number"
      echo "=========================================="

      for iteration in $(seq 1 "$tot_iteration"); do
        echo "Processing iteration $iteration/$tot_iteration for $simulation_name"
        config_dir=$SAVE_PATH$simulation_name/seed_$iteration
        mkdir -p "$config_dir"
        current_config="$config_dir/config_single_iter.yaml"
        echo "Step 1: Creating config file: $current_config"
        python3 "$YAML_MODIFIER" "$TEMPLATE_CONFIG" "$current_config" \
                    "p=$p" \
                    "save_path=$SAVE_PATH" \
                    "n_g1=$n_g1" \
                    "n_g2=$n_g2" \
                    "red_number=$red_number" \
                    "simulation_name=$simulation_name" \
                    "iteration=$iteration" \
                    "tot_iteration=$tot_iteration"
      done
    #done
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

    # Submit the jobs to perform the comparison with the literature methods for this iteration
    LOG_FILE="${iter_sim_dir}results_lit_comparison/logs/results_litt.log"
    RES=$(sbatch --parsable --output="$LOG_FILE" "Simulation_studies/Step1/FuDGE/Sbatch_litt_comp.sbatch" "$LIT_COMP_SCRIPT" "$yaml_iter_file_path" "$CURRENT_FOLDER")
    echo "For literature comparison running job ID: $RES"

  done

  echo "All iterations submitted for $sim_dir"
  echo "Waiting for all jobs to complete before running results checker..."
      
  while [ "$(squeue -u $USER | wc -l)" -ge 2 ]; do
    echo "Jobs limit reached, sleeping for a while..."
    sleep 240
  done
  echo "Step 4: Evaluating results for $sim_dir"
  # Evaluate results for our method and FuDGE
  Rscript Simulation_studies/Results_evaluation_simulations.R "$yaml_iter_file_path"
  Rscript Simulation_studies/Step1/FuDGE/Results_evaluation_sim_litt_comp.R "$yaml_iter_file_path"
done

# Process log files to extract computational times for both FuDGE and our method
Rscript Simulation_studies/Step1/Process_logs_comp_time.R $SAVE_PATH

