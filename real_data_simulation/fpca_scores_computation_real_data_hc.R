# rm(list=ls(all=TRUE))
# setwd("conditional_neurofgm/")
library(fda)
library(RcppCNPy)

n = 20
p  = 61 
M=5

idx_hc = c(0:n)

mat_t = seq(0, length.out = n)
k = 1
for(i in idx_hc){
  print( paste('Processing patient ', i))
  signal <- npyLoad(file = paste('real_data/hc_data/array_fd_data/fd_data_', as.character(i), '.npy', sep=''))  ##da modificare x prendere soggetti giusti
  mat_t[k] = dim(signal)[2]
  k = k + 1
}

#qui uno può capire max e min di mat_t che contiene la lunghezza di ogni segnale e di conseguenza capire quanto tempo si vuole considerare

##
fs = 0.01
tau = min(mat_t)
T_end = fs*tau
time <- seq(from = 0, to = T_end, length.out = tau)  
##


fpc.score <- numeric(0)
for(j in 1:p){
  print( paste('Processing node ', j))
  obs.val.matrix <- matrix(0, nrow= tau, ncol=n)  ##
  kk = 1
  for (i in idx_hc){
    print( paste('Processing patient ', i))
    signal <- npyLoad(file = paste('real_data/hc_data/array_fd_data/fd_data_', as.character(i), '.npy', sep=''))
    obs.val.vec = as.vector(signal[j, 1:tau])
    obs.val.matrix[, kk] = obs.val.vec
    kk = kk+1
  }
  bspline.basis <- create.bspline.basis(rangeval=c(0, T_end), nbasis=M)
  fd.object.array <- Data2fd(argvals=time, y=obs.val.matrix, basisobj=bspline.basis)
  # FPCA process
  fpc.score <- cbind(fpc.score, pca.fd(fd.object.array, nharm=M)$scores)
}

write.csv(fpc.score, "real_data/fpca_scores_hc.csv")