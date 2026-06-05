# Produces individual node results for simulations. Run one node per job in a SLURM array.
# Usage: Rscript Computation_templates/Script_sbatch_parallel_simulations.R config_file.yaml

# Required packages
suppressPackageStartupMessages(library(yaml))

# Load parameters from YAML configuration file
args <- commandArgs(trailingOnly = TRUE)
yaml_file_path <- if (length(args) > 0) args[[1]] else stop("Please provide the path to the YAML configuration file as an argument.")
config <- yaml.load_file(yaml_file_path)
time.start <- Sys.time()

save_path       <- config$save_path
simulation_name <- config$simulation_name
iteration       <- config$iteration
name_output     <- config$name_output
input_path      <- paste0(save_path, simulation_name, "/", "seed_", iteration, "/",
                          "net_est_input_", name_output,".rda")
output_path     <- paste0(save_path, simulation_name, "/", "seed_", iteration, "/",
                          "results/")

n_basis              <- config$n_basis_for_dim_reduction
L                    <- config$L
K                    <- config$K
thres.ctrl           <- unlist(config$thres_ctrl)
tol.abs              <- 0.0001
tol.rel              <- 0.0001
verbose              <- config$verbose
p.rand.lam           <- config$p_rand_lam
p.rand.thr           <- config$p_rand_thr
pre_screen           <- FALSE
pre_screen_threshold <- NULL
distance_path        <- NULL

cat("Parameter source:", yaml_file_path, "\n")

#######################################
#  NODE-WISE NEIGHBORHOOD ESTIMATION  #
#######################################

j           <- if (length(args) > 1) as.numeric(args[[2]]) else stop("Please provide the node index.")
load(input_path)
if (!(exists("covariates_df", envir = .GlobalEnv) && is.data.frame(covariates_df))) {
  covariates_df <- NULL
  cat("The covariate file was not included or invalid. Carrying on analysis to estimate the unique population network.")
}
if (pre_screen) {
  screening_matrix <- get(load(distance_path))
  diag(screening_matrix) <- 0
} else {
  screening_matrix <- NULL
}
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

source(file.path("Computation_templates", "algorithm_functions.R"))

run_node_computation(
  j = j, scores_df = scores_df, covariates_df = covariates_df,
  n_basis = n_basis, L = L, K = K, thres.ctrl = thres.ctrl,
  p.rand.lam = p.rand.lam, p.rand.thr = p.rand.thr,
  tol.abs = tol.abs, tol.rel = tol.rel,
  iteration = iteration,
  output_path = output_path, name_output = name_output,
  pre_screen = pre_screen, pre_screen_threshold = pre_screen_threshold,
  screening_matrix = screening_matrix, verbose = verbose
)

cat("\nComputational time:", as.numeric(difftime(Sys.time(), time.start, units = "secs")), "seconds\n")
