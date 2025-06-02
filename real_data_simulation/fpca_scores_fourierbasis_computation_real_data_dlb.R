library(fda)
library(RcppCNPy)
traceback()

n = 33
p  = 61 ##
M = 5


idx_dlb = c(0, 1, 2, 3, 4, '5_0', 6, 7, 8, 9,
            11, '12_1', '14_1', 15, '16_0', 17, 18, '19_1',
            20, 21, 23, 24, 25, 26, 27,
            28, 29, '30_0',31 ,'32_2' ,33, 34, 35)


fs =  512
tau = 73219
T = round(tau/fs)
time <- seq(from = 0, to = T, length.out = tau)
n_scores = round(T*2*12) + 100 
w = 100

fpc_fourier.score = numeric(0)
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
  
  drop_indx = c(c(1:(T*4*2-1)),c((T*4*2+2):(T*5*2-1)),c((T*5*2+2):(T*6*2-1)), c((T*6*2+2):(T*7*2-1)) ,c((T*7*2+2):(T*8*2-1)),c((T*8*2+2):(T*9*2-1)),c((T*9*2+2):(T*10*2-1)), c((T*10*2+2):(T*11*2-1)), c((T*11*2+2):(T*12*2-1)), c(c(T*12*2+1):(n_scores))) 
  fourier.basis = create.fourier.basis(rangeval = c(0,T), nbasis = n_scores, dropind = drop_indx)
  fd_fourier.object.array = Data2fd(argvals= time, y=obs.val.matrix, basisobj=fourier.basis, lambda = 0)
  # plot(fd_fourier.object.array[1]) 

  
  # # FPCA process
  #MM = n_scores- length(drop_indx)+1
  fpc_comp = pca.fd(fd_fourier.object.array, nharm=M)
  fpc_fourier.score <- cbind(fpc_fourier.score, fpc_comp$scores)
  
}

# matplot(signal, type="l")
# plot(fd.object.array)
plot(fd_fourier.object.array)
# plot(fd_fromfdPar)
# matplot(fpc.score, type="l")
#matplot(fpc_fourier.score, type="l")


#write.csv(fpc_fourier.score, "real_data_scores/fpca_scores_fourierbasis_filtered_individ_data/fpca5_scores_fourierbasis_8_12_individ_dlb.csv")
