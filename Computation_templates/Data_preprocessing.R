# Required packages
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(matrixcalc))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(abind))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pheatmap))

# Load parameters from YAML configuration file
args <- commandArgs(trailingOnly = TRUE)
yaml_file_path <- if (length(args) > 0) args[[1]] else
  stop("Please provide the path to the YAML configuration file as an argument.")
config <- yaml.load_file(yaml_file_path)

time.start <- Sys.time()
print("STARTING: Preprocessing of the data")
input_path <- config$input_path
observed_functional_data_path <- config$observed_functional_data_path
n_nodes <- config$n_nodes
tau <- config$time_points 
time_range <- config$time_range
rec_basis_type <- config$rec_basis_type
rec_basis_number <- config$rec_basis_number
M <- config$n_basis_for_dim_reduction

if (!dir.exists(dirname(input_path))) {
  dir.create(dirname(input_path), recursive = TRUE)
}

####################################
#     DIM REDUCTION: FPC SCORES    #
####################################

load(observed_functional_data_path)
if (rec_basis_type == "bsplines"){
  basis <- create.bspline.basis(rangeval=c(time_range[1], time_range[2]), nbasis=rec_basis_number)
} else{
  basis <- create.fourier.basis(rangeval=c(time_range[1], time_range[2]),nbasis=rec_basis_number)
}

step <- (time_range[2] - time_range[1]) / tau
obs.time <- seq(time_range[1] + step, time_range[2], step)

fpc.score <- numeric(0)
for(j in 1:n_nodes){
  obs.val.matrix <- matrix(0, nrow=tau, ncol=dim(discrete_fun_obs)[1])
  for (i in c(1:dim(discrete_fun_obs)[1])){
    obs.val.vec <- as.vector(discrete_fun_obs[i, j, ])
    obs.val.matrix[, i] <- obs.val.vec
  }
  fd.object.array <- Data2fd(argvals=obs.time, y=obs.val.matrix, basisobj=basis)
  fpc.score <- cbind(fpc.score, pca.fd(fd.object.array, nharm=M)$scores)
}

n_nodes <- ncol(fpc.score)/M
n_samples <- nrow(fpc.score)
names <- rep(NA, ncol(fpc.score))
for(l in 1:ncol(fpc.score)){
  names[l] <- paste("f",ceiling(l/M) ,".",l%%M, sep ="")
}
colnames(fpc.score) <- names

scores_df <- as.data.frame(fpc.score)
save(scores_df,covariates_df, file=input_path)
print(paste0("END: Preprocessing of the data completed in: ", Sys.time()-time.start))
