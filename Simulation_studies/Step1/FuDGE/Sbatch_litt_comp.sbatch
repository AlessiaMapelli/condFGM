#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUR_EMAIL_ADDRESS
#SBATCH --partition=cpuq
#SBATCH --mem=3GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=Brain_sim_litt_comp
#SBATCH --time=12:00:00

r_script=$1
conf_file=$2

cd $3

Rscript "${r_script}" "${conf_file}" 
