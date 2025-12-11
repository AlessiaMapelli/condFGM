#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alessia-mapelli@fht.org
#SBATCH --output=/logs/log_name.log
#SBATCH --partition=cpuq
#SBATCH --mem=16GB
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --job-name=Brain_sim_1
#SBATCH --time=07:00:00

source /center/healthds/singularity_functions
cd $HOME
Rscript Scrtipt_sbatch_sequential.R



