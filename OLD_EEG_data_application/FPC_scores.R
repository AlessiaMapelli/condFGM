# Required packages
rm(list=ls(all=TRUE))
packages <- c('yaml')
install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(abind))

#################################################
## USER DEFINED PARAMETERS (MODIFY THE PATH TO THE CORRECT YAML FILE)
#################################################
save.path <- "EEG_dataset/"
dir.create(save.path)
dir.create(paste0(save.path))
iteration = 1 
  
foldname = paste0(save.path, "/", "seed_", iteration, "/")
dir.create(foldname)
args <- commandArgs(trailingOnly = TRUE)
yaml_file_path = "EEG_dataset/config_template.yaml"
config <- yaml.load_file(yaml_file_path)
load("EEG_dataset/alco_filtered_array.Rdata")
load("EEG_dataset/contrl_filtered_array.Rdata")

name_output = "Fudge_dataset"

p <- config$p
n_g1 <- config$n_g1
n_g2 <- config$n_g2
mu <- config$mu
tau <- config$tau 
group=c(rep(1,n_g1), rep(2,n_g2))

rec_basis_type <- config$rec_basis_type
rec_basis_number <- config$rec_basis_number
M <- config$M

suppressPackageStartupMessages(library(fda))

obs.time <- seq(1/tau, 1, 1/tau)
h_g1 <- contrl.filtered.array
h_g2 <- alco.filtered.array
h <- abind(h_g1, h_g2, along = 1)
save(h,group, file=paste(foldname, "/Original_data_generated_", name_output, ".RData", sep=""))


####################################
#     PART 2: GAIN FPC SCORE       #
####################################

if (rec_basis_type == "bsplines"){
  basis <- create.bspline.basis(rangeval=c(0, 1), nbasis=rec_basis_number)
} else{
  basis <- create.fourier.basis(rangeval=c(0,1),nbasis=rec_basis_number)
}


fpc.score <- numeric(0)
for(j in 1:p){
  obs.val.matrix <- matrix(0, nrow=tau, ncol=dim(h)[1])
  for (i in c(1:dim(h)[1])){
    obs.val.vec <- as.vector(h[i, j, ])
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

write.csv(fpc.score, paste(foldname, "fpc_scores_", name_output, ".csv", sep=""))
write.csv(data.frame(group=group),paste(foldname, "grouping_factor_", name_output, ".csv", sep=""))

save(fpc.score, group, file =paste(foldname, "/fpc_scores_and_group", name_output, ".RData", sep=""))
# save(group, paste(foldname, "/fpc_scores_and_group", name_output, ".RData", sep=""))

fpc.score <- as.data.frame(fpc.score)
scoredata_matrix <- as.matrix(fpc.score)
Scor <- as.data.frame(cor(scoredata_matrix))
colnames(Scor) <- colnames(fpc.score)
rownames(Scor) <- colnames(fpc.score)

save(Scor, file=paste(foldname,"Score_corr_matrix_", name_output, ".rda", sep=""))

node_Scor <- matrix(NA, nrow = p, ncol = p)
for (i in 1:p) {
  for (j in 1:p) {
    row_idx <- ((i - 1) * M + 1):(i * M)
    col_idx <- ((j - 1) * M + 1):(j * M)
    
    block <- abs(Scor[row_idx, col_idx])
    node_Scor[i, j] <- max(block, na.rm = TRUE)
  }
}

names_node_corr <- rep(NA, p)
for(l in 1:p){
  names_node_corr[l] <- paste("f",l, sep ="")
}
colnames(node_Scor) <- names_node_corr
rownames(node_Scor) <- names_node_corr

save(node_Scor, file=paste(foldname,"Node_corr_matrix_", name_output, ".rda", sep=""))
print("END: Simulation of the data completed in: ")

output_path = paste0(foldname, "results/")
dir.create(output_path)


# fpc_comp = pca.fd(fd.object.array, nharm=M)
# fpc.score <- cbind(fpc.score, fpc_comp$scores)
# 
# # var_explained[[j]] <- fpc_comp$varprop
# # matplot(do.call(rbind, var_explained), type = "l", main = "Variance Explained per Node")
# 
# 
# fd_reconstructed <- fd(fpc_comp$harmonics$coefs %*% t(fpc_comp$scores), fpc_comp$harmonics$basis)
# matplot(obs.time, obs.val.matrix[,1], main = "Original curve - node 61", type = "l")
# plot(fd.object.array, main = "bsline basis fitted curve - node 61")
# plot(fd_reconstructed, main = "FPC fitted curve - node 61")
# 
# plotfit.fd(argvals = obs.time, y = obs.val.matrix, fdobj = fd.object.array, index = c(1), rng = c(0,1) , main = "Bspline basis approximated signal - Node 61, subject 1" )
# plotfit.fd(argvals = obs.time, y = obs.val.matrix, fdobj = fd_reconstructed, index = c(1), rng = c(0,1), main = "FPC basis approximated signal - Node 61, subject 1"   )


covariates_df = data.frame(group = group)
scores_df = fpc.score
save(scores_df, covariates_df, file = "EEG_dataset/seed_1/Fudge_dataset.RData")
