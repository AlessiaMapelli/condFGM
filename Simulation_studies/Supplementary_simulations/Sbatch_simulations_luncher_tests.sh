#!/bin/bash

cd /Simulation_studies/Step3

# Configuration
TEMPLATE_CONFIG="config_template.yaml"
YAML_MODIFIER="modify_yaml_config.py"
SAVE_PATH="/Simulation_studies/Step3/Simulations_results/"
JOBS_LIMIT=200
tot_iteration=10

# Simulation parameter arrays (proper bash syntax)
p_seq=(10 50)
n_g1_seq=(50 75 100 150 200 400 750 1000)
n_g2_seq=(50 75 100 150 200 400 750 1000)  # Assuming same values as n_g1

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

Rscript Data_simulator_full_data.R $SAVE_PATH

for sim_dir in "$SAVE_PATH"/*/; do
  for iter_sim_dir in "$sim_dir"*/; do
    echo "Step 3: Launching parallel screening jobs for $iter_sim_dir"
    yaml_iter_file_path="${iter_sim_dir}config_single_iter.yaml"

    while [ "$(squeue -u $USER | wc -l)" -ge "$JOBS_LIMIT" ]; do
      echo "Jobs limit reached, sleeping for a while..."
      sleep 60
    done
    bash Sbatch_parallel_luncher_tests.sh "$yaml_iter_file_path"
  done

  echo "All iterations submitted for $simulation_name"
  echo "Waiting for all screening jobs to complete before running results checker..."
  
  while [ "$(squeue -u $USER | wc -l)" -ge 2 ]; do
    echo "Jobs limit reached, sleeping for a while..."
    sleep 60
  done
  
  echo "Step 4: Evaluating results for $sim_dir"
  Rscript Check_results_screening_procedure_tests.R "$yaml_iter_file_path"

done


