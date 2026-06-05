#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUR_EMAIL_ADDRESS
#SBATCH --partition=cpuq
#SBATCH --mem=1GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=Node_computation
#SBATCH --time=00:20:00

input=$1
r_script=$2
conf_file=$3

cd $4

Rscript "${r_script}" "${conf_file}" "${input}"
