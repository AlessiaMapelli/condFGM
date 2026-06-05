# Required packages
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(matrixcalc))
suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(quadprog))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(JGL))

# Load parameters from YAML configuration file
args <- commandArgs(trailingOnly = TRUE)
yaml_file_path <- if (length(args) > 0) args[[1]] else
  stop("Please provide the path to the YAML configuration file as an argument.")
config <- yaml.load_file(yaml_file_path)
time.start <- Sys.time()

save_path = config$save_path
simulation_name = config$simulation_name
iteration = config$iteration
name_output = config$name_output
input_path      <- paste0(save_path, simulation_name, "/", "seed_", iteration, "/",
                          "net_est_input_", name_output,".rda")
output_path     <- paste0(save_path, simulation_name, "/", "seed_", iteration, "/",
                          "results_lit_comparison/")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
p <- config$p
n_basis <- config$n_basis_for_dim_reduction
func.path <- "Simulation_studies/Step1/FuDGE/Function_Definitions"
source(file.path(func.path, "ProxAlg_DFGM.R"))
source(file.path(func.path, "FFGL_ADMM.R"))
source(file.path(func.path, "FFGL2_ADMM.R"))

load(input_path)
n_nodes <- ncol(scores_df)/n_basis
n_samples <- nrow(scores_df)
full_data <- cbind(covariates_df,scores_df)


#################################################
## RUN COMPARISON ALGORITHMS on our data following the code provided in
## https://github.com/boxinz17/FuDGE/tree/master/Simulation/Simulation_Code
#################################################


################### FuDGE #########################
# Data Preparation
principle.score.X <- full_data[full_data$group == "1",][, -1]
principle.score.Y <- full_data[full_data$group == "2",][, -1]

principle.score.X.cen <- scale(principle.score.X, center=T, scale=F)
principle.score.Y.cen <- scale(principle.score.Y, center=T, scale=F)

estimated.cov.pc.X <- (t(principle.score.X.cen) %*% principle.score.X.cen) / (n_samples-1)
estimated.cov.pc.Y <- (t(principle.score.Y.cen) %*% principle.score.Y.cen) / (n_samples-1)

# Setting hyperparams
lambdamax <- 2
lambdamin <- 0
steplength <- 0.05
lambda.v.fudge <- seq(lambdamax, lambdamin, -steplength)
Delta.initial <- matrix(0, nrow=(p*n_basis), ncol=(p*n_basis))

# Estimate the Differential Graph
AdjMats_lambda_FuDGE <- list()
for (lambda.ind in c(1:length(lambda.v.fudge)) ){
  lambda.choose <- lambda.v.fudge[lambda.ind]
  DFGM_g <- ProxAlg_DFGM(estimated.cov.pc.X, estimated.cov.pc.Y, p, n_basis, 
                         lambda.choose, Intialization=Delta.initial)
  
  # Update initialization (warm start)
  Delta.initial <- DFGM_g$DelatMathat
  
  FrobMat <- DFGM_g$blockFrob
  
  Edgeshat <- matrix(0, nrow=p, ncol=p)
  Edgeshat[which(FrobMat > 0)] <- 1
  AdjMats_lambda_FuDGE[[lambda.ind]]<- list(Diff=Edgeshat)
  
}

cat("\n FuDGE computational time:",
    as.numeric(difftime(Sys.time(), time.start, units = "secs")), "seconds\n")

################### FGL #########################

# Data Preparation
prin.score <- list(X=principle.score.X, Y=principle.score.Y)

# Setting hyperparams
lambda1.choose.fgl <- 0.01
lambda2.v.fgl <- rev(seq(0, 1, length.out = length(lambda.v.fudge)))

# Estimate the Differential Graph
AdjMats_lambda_FGL <- list()
for (lambda.ind in c(1:length(lambda2.v.fgl)) ){
  lambda2.choose <- lambda2.v.fgl[lambda.ind]
  JGL.result <- JGL(prin.score, penalty="fused", lambda1=lambda1.choose.fgl, 
                    lambda2=lambda2.choose, penalize.diagonal=FALSE, 
                    return.whole.theta=TRUE)
  
  ThetaX <- JGL.result$theta[[1]]
  ThetaY <- JGL.result$theta[[2]]
  Edgeshat <- matrix(0, nrow=p, ncol=p)
  Edgeshat.pop <- matrix(0, nrow=p, ncol=p)
  Edgeshat.group <- matrix(0, nrow=p, ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      diff <- frobenius.norm(ThetaX[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)] - 
                               ThetaY[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      pop <- frobenius.norm(ThetaX[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      group <- frobenius.norm(ThetaY[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      if (diff>0){
        Edgeshat[i,j] <- 1
      }
      if (pop>0){
        Edgeshat.pop[i,j] <- 1
      }
      if (group>0){
        Edgeshat.group[i,j] <- 1
      }
      
      
    }
  }
  AdjMats_lambda_FGL[[lambda.ind]]<- list(Pop=Edgeshat.pop,Group=Edgeshat.group, Diff=Edgeshat)
}

################### FFGL/FFGL2 #########################

# Data Preparation
cov.list <- list(X=estimated.cov.pc.X, Y=estimated.cov.pc.Y)
n.sam <- c(n_samples, n_samples)

# Setting hyperparams
lambda2.v.ffgl <- rev(seq(0, 100, length.out = length(lambda.v.fudge)))
lambda1.choose.ffgl <- 0.1

# Estimate the Differential Graph 
AdjMats_lambda_FFGL <- list()
AdjMats_lambda_FFGL2 <- list()
for (lambda.ind in c(1:length(lambda2.v.ffgl)) ){

  ################### FFGL #########################
  lambda2.choose <- lambda2.v.ffgl[lambda.ind]
  FFGL.result <- FFGL_ADMM(cov.list, n.sam, p, n_basis, lambda1.choose.ffgl, lambda2.choose)
  ThetaX <- FFGL.result$est[[1]]
  ThetaY <- FFGL.result$est[[2]]
  Edgeshat <- matrix(0, nrow=p, ncol=p)
  Edgeshat.pop <- matrix(0, nrow=p, ncol=p)
  Edgeshat.group <- matrix(0, nrow=p, ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      diff <- frobenius.norm(ThetaX[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)] - 
                               ThetaY[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      pop <- frobenius.norm(ThetaX[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      group <- frobenius.norm(ThetaY[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      if (diff>0.02){ 
        Edgeshat[i,j] <- 1
      }
      if (pop>0.02){ 
         Edgeshat.pop[i,j] <- 1
      }
      if (group>0.02){
         Edgeshat.group[i,j] <- 1
      }

      
    }
  }
  AdjMats_lambda_FFGL[[lambda.ind]]<- list(Pop=Edgeshat.pop,Group=Edgeshat.group, Diff=Edgeshat)

  ################### FFGL2 #########################
  
  FFGL2.result <- FFGL2_ADMM(cov.list, n.sam, p, n_basis, lambda1.choose.ffgl, lambda2.choose)
  ThetaX <- FFGL2.result$est[[1]]
  ThetaY <- FFGL2.result$est[[2]]
  Edgeshat <- matrix(0, nrow=p, ncol=p)
  Edgeshat.pop <- matrix(0, nrow=p, ncol=p)
  Edgeshat.group <- matrix(0, nrow=p, ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      diff <- frobenius.norm(ThetaX[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)] - 
                               ThetaY[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      pop <- frobenius.norm(ThetaX[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      group <- frobenius.norm(ThetaY[((i-1)*n_basis+1):(i*n_basis), ((j-1)*n_basis+1):(j*n_basis)])
      if (diff>0.02){ 
        Edgeshat[i,j] <- 1
      }
      if (pop>0.02){
         Edgeshat.pop[i,j] <- 1
      }
      if (group>0.02){
         Edgeshat.group[i,j] <- 1
      }

      
    }
  }
  AdjMats_lambda_FFGL2[[lambda.ind]]<- list(Pop=Edgeshat.pop,Group=Edgeshat.group, Diff=Edgeshat)
  
}

full_result_path = paste(output_path, name_output, "_litt_comparison", ".rda", sep ="")
save(AdjMats_lambda_FuDGE, AdjMats_lambda_FGL, AdjMats_lambda_FFGL, AdjMats_lambda_FFGL2,
     file = full_result_path)
cat("\n Output saved to: ", full_result_path,"\n" )
