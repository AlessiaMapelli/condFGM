####################################
#  PART 0: PRELIMINARY FUNCTIONS   #
####################################

# Generating M*M Tridiagonal
tridiag <- function(M){
  result <- diag(M)
  for(i in 1:M){
    for(j in 1:M){
      if(abs(i-j)==1) result[i,j] <- 0.5
    }
  }
  return(result)
}

toeplitz <- function(M){
  result <- diag(M)
  for(i in 1:M){
    for(j in 1:M){
      if(i != j){
        result[i,j] <- 0.5 / abs(i-j)
      }
    }
  }
  return(result)
}

####################################
# 0. Base model for precision matrix of delta #
####################################

cov.mat.model <- function(p,M,diff_nodes){
  # Input:
  #   p, number of covariates
  #   M, number of basis functions
  # Output:
  #   the precision matrix that generates delta via MVN
  Theta <- matrix(nrow=p*M, ncol=p*M)
  for(i in 1:p){
    for(j in 1:p){
      if(i==j) Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- toeplitz(M)
      else if(abs(i-j)==1) Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- tridiag(M) * 0.4
      else if(abs(i-j)==2) Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- diag(M) * 0.2
      else Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- matrix(0, M, M)
    }
  }
  eig.val <- eigen(Theta)$values
  if(min(eig.val)<0.05){
    Theta = Theta + (abs(min(eig.val)) + 0.05) * diag(p*M)
  }
  
  return(Theta)
}

####################################
# 1. Possible modifications due to a binary covariate#
####################################

cov.mat.model.red.top <- function(p,M,diff_nodes){
  # Input:
  #   p, number of covariates
  #   M, number of basis functions
  # Output:
  #   the precision matrix that generates delta via MVN
  Theta.full <- cov.mat.model(p,M)
  Theta.red <- Theta.full
  for(i in 1:diff_nodes){
    for(j in c(i-2, i-1, i+1, i+2)){
      if (j >= 1 && j <= diff_nodes && j != i) {  # Ensure valid index and exclude self
        Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- 0
      }
    }
  }
  return(Theta.red)
}

cov.mat.model.red.bottom <- function(p,M,diff_nodes) {
  # Input:
  #   p           - number of covariates (nodes)
  #   M           - number of basis functions per node
  #   diff_nodes  - number of nodes to disconnect (at the end)
  # Output:
  #   Theta.red   - precision matrix with last diff_nodes disconnected

  Theta.full <- cov.mat.model(p, M)
  Theta.red <- Theta.full

  # Start from the last 'diff_nodes' indices
  for (i in (p - diff_nodes + 1):p) {
    for (j in c(i - 2, i - 1, i + 1, i + 2)) {
      if (j >= (p - diff_nodes + 1) && j <= p && j != i) {
        Theta.red[((i - 1) * M + 1):(i * M), ((j - 1) * M + 1):(j * M)] <- 0
      }
    }
  }

  return(Theta.red)
}

cov.mat.model.red.top.bottom <- function(p,M,diff_nodes) {
  # Input:
  #   p           - number of covariates (nodes)
  #   M           - number of basis functions per node
  #   diff_nodes  - number of nodes to disconnect (at the end)
  # Output:
  #   Theta.red   - precision matrix with last diff_nodes disconnected

  Theta.full <- cov.mat.model(p, M)
  Theta.red <- Theta.full

  for(i in 1:diff_nodes){
    for(j in c(i-2, i-1, i+1, i+2)){
      if (j >= 1 && j <= diff_nodes && j != i) {  # Ensure valid index and exclude self
        Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- 0
      }
    }
  }

  for (i in (p - diff_nodes + 1):p) {
    for (j in c(i - 2, i - 1, i + 1, i + 2)) {
      if (j >= (p - diff_nodes + 1) && j <= p && j != i) {
        Theta.red[((i - 1) * M + 1):(i * M), ((j - 1) * M + 1):(j * M)] <- 0
      }
    }
  }

  return(Theta.red)
}

####################################
# 2. Possible modifications due to a positive continuous covariate #
####################################

cov.mat.model.bottom <- function(p, M, diff_nodes) {
  # Input:
  #   p           - number of covariates (nodes)
  #   M           - number of basis functions per node
  #   diff_nodes  - number of nodes to disconnect (at the end)
  # Output:
  #   Theta.red   - precision matrix with last diff_nodes disconnected
  Theta.red <- matrix(0, nrow=p*M, ncol=p*M)
  
  # Start from the last 'diff_nodes' indices
  for (i in (p - diff_nodes + 1):p) {
    for (j in c(i - 2, i - 1, i + 1, i + 2)) {
      if (j >= (p - diff_nodes + 1) && j <= p && j != i) {
        if(i==j) Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- toeplitz(M)
        else if(abs(i-j)==1) Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- tridiag(M) * 0.4
        else if(abs(i-j)==2) Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- diag(M) * 0.2
        else Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- matrix(0, M, M)
      }
    }
  }
  
  return(Theta.red)
}


####################################
# 3. Possible modifications due to a continuous covariate #
####################################

cov.mat.model.bottom.continuous <- function(p, M, diff_nodes) {
  # Input:
  #   p           - number of covariates (nodes)
  #   M           - number of basis functions per node
  #   diff_nodes  - number of nodes to disconnect (at the end)
  # Output:
  #   Theta.red   - precision matrix with last diff_nodes disconnected
  Theta.red <- matrix(0, nrow=p*M, ncol=p*M)
  
  # Start from the last 'diff_nodes' indices
  for (i in (p - diff_nodes + 1):p) {
    for (j in c(i - 2, i - 1, i, i + 1, i + 2)) {
      if (j >= (p - diff_nodes + 1) && j <= p) {
        #if(i==j) Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- toeplitz(M)*0.25
        # else 
	if(abs(i-j)==1) Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- tridiag(M) * 0.1
        else if(abs(i-j)==2) Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- diag(M) * 0.05
        else Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- matrix(0, M, M)
      }
    }
  }
  
  return(Theta.red)
}
