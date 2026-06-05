#!/bin/bash

# Configuration
TEMPLATE_CONFIG="Simulation_studies/Supplementary_simulations/config_template_supp_simulation.yaml"
YAML_MODIFIER="Simulation_studies/Helper_functions/modify_yaml_config.py"
SAVE_PATH="Simulation_studies/Supplementary_simulations/simulation_settings/"
tot_iteration=10

# Simulation parameter grid
p_seq=(10 50)
n_g1_seq=(50 75 100 150 200 400 750 1000)
n_g2_seq=(50 75 100 150 200 400 750 1000)  # Assuming same values as n_g1

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

        
echo "  Step 2: Running data simulator (iteration $iteration)"

# Run the data simulator for all iterations
Rscript Simulation_studies/Supplementary_simulations/Full_data_simulator_supp_simulation.R $SAVE_PATH

# Loop through the generated simulation directories and submit parallel jobs for each iteration
for sim_dir in "$SAVE_PATH"/*/; do
  for iter_sim_dir in "$sim_dir"*/; do
    echo "Step 3: Launching parallel jobs for $iter_sim_dir"
    yaml_iter_file_path="${iter_sim_dir}config_single_iter.yaml"

    # Submit the jobs to perform our algorithm for this iteration
    bash Simulation_studies/Sbatch_parallel_luncher_simulations.sh "$yaml_iter_file_path"
  done

  echo "All iterations submitted for $sim_dir"
  echo "Waiting for all screening jobs to complete before running results checker..."
  
  while [ "$(squeue -u $USER | wc -l)" -ge 2 ]; do
    echo "Jobs limit reached, sleeping for a while..."
    sleep 60
  done
  
  echo "Step 4: Evaluating results for $sim_dir"
  # Evaluate results
  Rscript Simulation_studies/Supplementary_simulations/Results_evaluation_supp_simulation.R "$yaml_iter_file_path"
done


