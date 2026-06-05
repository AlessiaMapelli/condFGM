# Precision matrix generation functions for all simulation settings.
# All model functions follow the uniform signature (p, M, diff_nodes = NULL)
# so they can be dispatched via match.fun() from Data_simulator scripts.

####################################
# Helper matrix constructors
####################################

# @param M  integer; matrix dimension
# @return   M x M tridiagonal matrix with 0.5 on super- and sub-diagonals
tridiag <- function(M) {
  result <- diag(M)
  for (i in 1:M) {
    for (j in 1:M) {
      if (abs(i - j) == 1) result[i, j] <- 0.5
    }
  }
  return(result)
}

# @param M  integer; matrix dimension
# @return   M x M Toeplitz matrix: off-diagonal (i,j) entry = 0.5/|i-j|
toeplitz <- function(M) {
  result <- diag(M)
  for (i in 1:M) {
    for (j in 1:M) {
      if (i != j) result[i, j] <- 0.5 / abs(i - j)
    }
  }
  return(result)
}

####################################
# Base model
####################################

# @param p          integer; number of nodes
# @param M          integer; number of basis functions per node
# @param diff_nodes unused; present for uniform match.fun() signature
# @return           (p*M) x (p*M) Toeplitz-structured block precision matrix
cov.mat.model <- function(p, M, diff_nodes = NULL) {
  Theta <- matrix(nrow = p * M, ncol = p * M)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j)              Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- toeplitz(M)
      else if (abs(i-j) == 1) Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- tridiag(M) * 0.4
      else if (abs(i-j) == 2) Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- diag(M) * 0.2
      else                     Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- matrix(0, M, M)
    }
  }
  eig.val <- eigen(Theta)$values
  if (min(eig.val) < 0.05) Theta <- Theta + (abs(min(eig.val)) + 0.05) * diag(p * M)
  return(Theta)
}

####################################
# Models with zeroed differential edges (binary covariate / group difference)
####################################

# @param p          integer; number of nodes
# @param M          integer; number of basis functions per node
# @param diff_nodes integer; number of nodes to modify at the top of the network
# @return           (p*M) x (p*M) precision matrix with edges among first diff_nodes nodes zeroed
cov.mat.model.red.top <- function(p, M, diff_nodes) {
  Theta.full <- cov.mat.model(p, M)
  Theta.red <- Theta.full
  for (i in 1:diff_nodes) {
    for (j in c(i-2, i-1, i+1, i+2)) {
      if (j >= 1 && j <= diff_nodes && j != i)
        Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- 0
    }
  }
  return(Theta.red)
}

# @param p          integer; number of nodes
# @param M          integer; number of basis functions per node
# @param diff_nodes integer; number of nodes to modify at the bottom of the network
# @return           (p*M) x (p*M) precision matrix with edges among last diff_nodes nodes zeroed
cov.mat.model.red.bottom <- function(p, M, diff_nodes) {
  Theta.full <- cov.mat.model(p, M)
  Theta.red <- Theta.full
  for (i in (p - diff_nodes + 1):p) {
    for (j in c(i-2, i-1, i+1, i+2)) {
      if (j >= (p - diff_nodes + 1) && j <= p && j != i)
        Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- 0
    }
  }
  return(Theta.red)
}

# @param p          integer; number of nodes
# @param M          integer; number of basis functions per node
# @param diff_nodes integer; number of nodes to modify at both ends of the network
# @return           (p*M) x (p*M) precision matrix with edges at both ends zeroed
cov.mat.model.red.top.bottom <- function(p, M, diff_nodes) {
  Theta.full <- cov.mat.model(p, M)
  Theta.red <- Theta.full
  for (i in 1:diff_nodes) {
    for (j in c(i-2, i-1, i+1, i+2)) {
      if (j >= 1 && j <= diff_nodes && j != i)
        Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- 0
    }
  }
  for (i in (p - diff_nodes + 1):p) {
    for (j in c(i-2, i-1, i+1, i+2)) {
      if (j >= (p - diff_nodes + 1) && j <= p && j != i)
        Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- 0
    }
  }
  return(Theta.red)
}

####################################
# Models with weakened (not zeroed) differential edges
####################################

# @param p          integer; number of nodes
# @param M          integer; number of basis functions per node
# @param diff_nodes integer; number of nodes to modify at the top of the network
# @return           (p*M) x (p*M) precision matrix with edges among first diff_nodes nodes weakened
cov.mat.model.diff.weights.top <- function(p, M, diff_nodes) {
  Theta.full <- cov.mat.model(p, M)
  Theta.red <- Theta.full
  for (i in 1:diff_nodes) {
    for (j in c(i-2, i-1, i+1, i+2)) {
      if (j >= 1 && j <= diff_nodes && j != i) {
        if      (abs(i-j) == 1) Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- tridiag(M) * 0.1
        else if (abs(i-j) == 2) Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- diag(M) * 0.05
        else                     Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- matrix(0, M, M)
      }
    }
  }
  return(Theta.red)
}

# @param p          integer; number of nodes
# @param M          integer; number of basis functions per node
# @param diff_nodes integer; number of nodes to modify at the bottom of the network
# @return           (p*M) x (p*M) precision matrix with edges among last diff_nodes nodes weakened
cov.mat.model.diff.weights.bottom <- function(p, M, diff_nodes) {
  Theta.full <- cov.mat.model(p, M)
  Theta.red <- Theta.full
  for (i in (p - diff_nodes + 1):p) {
    for (j in c(i-2, i-1, i+1, i+2)) {
      if (j >= (p - diff_nodes + 1) && j <= p && j != i) {
        if      (abs(i-j) == 1) Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- tridiag(M) * 0.1
        else if (abs(i-j) == 2) Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- diag(M) * 0.05
        else                     Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- matrix(0, M, M)
      }
    }
  }
  return(Theta.red)
}

####################################
# Empty / diagonal-only model
####################################

# @param p          integer; number of nodes
# @param M          integer; number of basis functions per node
# @param diff_nodes unused; present for uniform match.fun() signature
# @return           (p*M) x (p*M) block-diagonal precision matrix (no cross-node edges)
cov.mat.model.empty <- function(p, M, diff_nodes = NULL) {
  Theta <- matrix(nrow = p * M, ncol = p * M)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- toeplitz(M)
      else        Theta[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- matrix(0, M, M)
    }
  }
  eig.val <- eigen(Theta)$values
  if (min(eig.val) < 0.05) Theta <- Theta + (abs(min(eig.val)) + 0.05) * diag(p * M)
  return(Theta)
}

####################################
# Continuous covariate models (Supplementary simulations)
####################################

# @param p          integer; number of nodes
# @param M          integer; number of basis functions per node
# @param diff_nodes integer; number of bottom nodes contributing to the covariate effect
# @return           (p*M) x (p*M) sparse matrix encoding positive-continuous covariate edges
cov.mat.model.bottom <- function(p, M, diff_nodes) {
  Theta.red <- matrix(0, nrow = p * M, ncol = p * M)
  for (i in (p - diff_nodes + 1):p) {
    for (j in c(i-2, i-1, i+1, i+2)) {
      if (j >= (p - diff_nodes + 1) && j <= p && j != i) {
        if (i == j)              Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- toeplitz(M)
        else if (abs(i-j) == 1) Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- tridiag(M) * 0.4
        else if (abs(i-j) == 2) Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- diag(M) * 0.2
        else                     Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- matrix(0, M, M)
      }
    }
  }
  return(Theta.red)
}

# @param p          integer; number of nodes
# @param M          integer; number of basis functions per node
# @param diff_nodes integer; number of bottom nodes contributing to the covariate effect
# @return           (p*M) x (p*M) sparse matrix encoding signed-continuous covariate edges
cov.mat.model.bottom.continuous <- function(p, M, diff_nodes) {
  Theta.red <- matrix(0, nrow = p * M, ncol = p * M)
  for (i in (p - diff_nodes + 1):p) {
    for (j in c(i-2, i-1, i, i+1, i+2)) {
      if (j >= (p - diff_nodes + 1) && j <= p) {
        if      (abs(i-j) == 1) Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- tridiag(M) * 0.1
        else if (abs(i-j) == 2) Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- diag(M) * 0.05
        else                     Theta.red[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)] <- matrix(0, M, M)
      }
    }
  }
  return(Theta.red)
}
