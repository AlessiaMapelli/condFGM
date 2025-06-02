# rm(list=ls(all=TRUE))
# setwd("conditional_neurofgm/")
library(fda)
library(RcppCNPy)

n = 33
p  = 61 ##
M = 5

idx_dlb = c(0, 1, 2, 3, 4, '5_0', 6, 7, 8, 9, 
              11, '12_1', '14_1', 15, '16_0', 17, 18, '19_1',
              20, 21, 23, 24, 25, 26, 27,
              28, 29, '30_0',31 ,'32_2' ,33, 34, 35)

mat_t = c(0,(n-1))
k = 1
for(i in idx_dlb){
  print( paste('Processing patient ', i))
  signal <- npyLoad(file = paste('real_data/dlb_data/array_fd_data/filtered_indiv_data/fd_data_fil_indiv_', as.character(i), '.npy', sep=''))  ##da modificare x prendere soggetti giusti
  mat_t[k] = dim(signal)[2]
  k = k + 1
}

#qui uno può capire max e min di mat_t che contiene la lunghezza di ogni segnale e di conseguenza capire quanto tempo si vuole considerare

##
fs =  512
tau = min(mat_t) 
T_end = tau/fs
time <- seq(from = 0, to = T_end, length.out = tau)
##

fpc.score <- numeric(0)
for(j in 1:p){
  print( paste('Processing node ', j))
  obs.val.matrix <- matrix(0, nrow= tau, ncol=n)  ##
  k = 1
  for (i in idx_dlb){
    #print( paste('Processing patient ', i))
    signal <- npyLoad(file = paste('real_data/dlb_data/array_fd_data/filtered_indiv_data/fd_data_fil_indiv_' ,as.character(i), '.npy', sep=''))
    obs.val.vec = as.vector(signal[j, 1:tau])
    obs.val.matrix[, k] = obs.val.vec
    k = k+1
  }
  bspline.basis <- create.bspline.basis(rangeval=c(0, T_end), nbasis=10)  #nbasis = M
  fd.object.array <- Data2fd(argvals=time, y=obs.val.matrix, basisobj=bspline.basis)
  # FPCA process
  #fpc.score <- cbind(fpc.score, pca.fd(fd.object.array, nharm=M)$scores)
}

#write.csv(fpc.score, "real_data_scores/fpca_scores_filtered_individ_data/fpca_scores_individ_dlb.csv")


