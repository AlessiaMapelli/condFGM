# Shared ADMM / group-LASSO functions used to solve the optimization problem
# (Section 3.2 of the paper)

# @param A.X  n x (n_blocks*M) design matrix
# @param A.Y  n x M response matrix
# @return scalar: supremum lambda value
lambda.sup <- function(A.X, A.Y){
  n <- nrow(A.X)
  M <- ncol(A.Y)
  p <- ncol(A.X)/M
  candidates <- rep(0,p)
  for(j in 1:p){
    A.X.j <- A.X[,(j-1)*M + (1:M)]
    candidates[j] <- norm(t(A.X.j) %*% A.Y,"F") / n
  }
  return(max(candidates))
}


# @param V.Grp  list of p matrices (Vk), each M x M
# @param d      vector of length p*M
# @param lambda penalty parameter of group lasso
# @param rho    penalty parameter of augmented Lagrangian
# @return list of p matrices (Pk), each M x M
grp.soft.thres.two.groups <- function(V.Grp, d, lambda, rho){
  M <- nrow(V.Grp[[1]])
  n_blocks <- length(V.Grp)
  D.Grp <- list()
  P.Grp <- list()

  for(k in 1:n_blocks){
    D.Grp[[k]] <- diag(d[(k-1)*M+(1:M)])
    left.plus <- matrix(0, M, M)
    norm.F <- norm(D.Grp[[k]] %*% V.Grp[[k]], "F")
    for(l in 1:M){
      left.plus[l,l] <- max(0, 1 - lambda / rho * D.Grp[[k]][l,l]^2 / norm.F )
    }

    P.Grp[[k]] <- left.plus %*% V.Grp[[k]]
  }

  return(P.Grp)
}

# @param A.X    n x (n_blocks*M) design matrix
# @param A.Y    n x M response matrix
# @param d      vector of length n_blocks*M
# @param B      (n_blocks*M) x M coefficient matrix
# @param lambda penalty parameter
# @return scalar: objective function value
objective.function.two.groups <- function(A.X, A.Y, d, B, lambda){
  n <- nrow(A.X)
  M <- ncol(A.Y)
  n_blocks <- ncol(A.X) / M

  D <- diag(d)
  term.1 <- norm(A.X %*% D %*% B - A.Y, "F")^2 / (2*n)
  term.2 <- 0
  for(k in 1:n_blocks){
    B.k <- B[(k-1)*M + (1:M), ]
    Dkk <- diag(d[(k-1)*M + (1:M)])
    term.2 <- term.2 + lambda * norm(Dkk %*% B.k ,"F")
  }
  obj <- term.1 + term.2
  return(obj)
}

# @param A.X      n x (2*(p-1)*M) design matrix
# @param A.Y      n x M response matrix
# @param d        vector of length 2*(p-1)*M
# @param lambda   group-lasso penalty
# @param rho.init initial ADMM penalty (default 1)
# @param maxiter  maximum iterations (default 2000)
# @param mu       rho adaptation trigger ratio (default 10)
# @param tau      rho adaptation factor (default 2)
# @param P.in, Q.in, U.in  warm-start matrices
# @param tol.rel, tol.abs  convergence tolerances
# @return list: P, Q, U, prim.res, dual.res, rho.path, tol.prim.path, tol.dual.path, obj.func.path
ADMM.grplasso.two.groups <- function(A.X, A.Y, d, lambda,
                                     rho.init=1, maxiter=2000, mu=10, tau=2,
                                     P.in, Q.in, U.in, tol.rel=1e-4, tol.abs=1e-4){
  n <- nrow(A.Y)
  M <- ncol(A.Y)
  n_blocks <- ncol(A.X) / M

  D <- diag(d)
  AX.D <- A.X %*% D
  DXXD.n <- (t(AX.D) %*% AX.D) / n
  DXY.n <- (t(AX.D) %*% A.Y)/n

  P.old <- P.in
  Q.old <- Q.in
  U.old <- U.in

  prim.resid.vec <- numeric(0)
  dual.resid.vec <- numeric(0)
  rho.path <- rho.init
  tol.prim.path <- numeric(0)
  tol.dual.path <- numeric(0)
  obj.func.path <- numeric(0)

  V <- Q.old - U.old
  V.Grp <- list()
  for(k in 1:n_blocks){
    V.Grp[[k]] <- V[(k-1)*M + (1:M), ]
  }

  rho <- rho.init
  rho.changed <- TRUE

  for(iter in 1:maxiter){

    P.new <- matrix(NA, nrow=n_blocks*M, ncol=M)
    P.Grp <- grp.soft.thres.two.groups(V.Grp, d, lambda, rho)
    for(k in 1:n_blocks){
      P.new[(k-1)*M + (1:M), ] <- P.Grp[[k]]
    }

    if(rho.changed){
      left.inv <- tryCatch(
        solve(DXXD.n + rho * diag(nrow(DXXD.n))),
        error = function(e) {
          stop("Matrix inversion error: Check dimensions of DXXD.n and rho term")
        }
      )
    }
    Q.new <- left.inv %*% (DXY.n + rho * P.new + rho * U.old)

    U.new <- U.old + P.new - Q.new

    prim.resid <- norm( (P.new - Q.new) , "F")
    dual.resid <- norm( rho * (Q.new - Q.old) , "F")
    prim.resid.vec <- c(prim.resid.vec, prim.resid)
    dual.resid.vec <- c(dual.resid.vec, dual.resid)

    tol.prim <- tol.abs * sqrt(n_blocks*M*M) + tol.rel * max(norm(P.new, "F"),
                                                             norm(Q.new, "F"))
    tol.dual <- tol.abs * sqrt(n_blocks*M*M) + tol.rel * norm(U.new, "F")

    if(prim.resid < tol.prim && dual.resid < tol.dual){
      break
    }

    V <- Q.new - U.new
    V.Grp <- list()
    for(k in 1:n_blocks){
      V.Grp[[k]] <- V[((k-1)*M + 1):(k*M), ]
    }

    if(prim.resid > mu * dual.resid){
      rho <- tau * rho
      rho.changed <- TRUE
    }else if(dual.resid > mu * prim.resid){
      rho <- rho / tau
      rho.changed <- TRUE
    }else{
      rho.changed <- FALSE
    }
    rho.path <- c(rho.path, rho)

    P.old <- P.new
    Q.old <- Q.new
    U.old <- U.new

    obj <- tryCatch(
      objective.function.two.groups(A.X=A.X, A.Y=A.Y, d=d, B=P.new, lambda=lambda),
      error = function(e) {
        stop("Error in objective function calculation: Check objective.function.two.groups()")
      }
    )

    obj.func.path <- c(obj.func.path, obj)
  }

  if(iter == maxiter){
    cat("ADMM did not converge.")
  }

  result <- list(P=P.new, Q=Q.new, U=U.new,
                 prim.res=prim.resid.vec, dual.res=dual.resid.vec,
                 rho.path=rho.path,
                 tol.prim.path=tol.prim.path, tol.dual.path=tol.dual.path,
                 obj.func.path=obj.func.path)
  return(result)
}


# @param j Node index (integer, 1-based)
# @param scores_df FPC score matrix (n_samples x p*n_basis)
# @param covariates_df Covariate data frame or NULL for population-only model
# @param n_basis Number of basis functions per node
# @param L Number of lambda values in grid
# @param K Number of cross-validation folds
# @param thres.ctrl Vector of threshold control values
# @param p.rand.lam Proportion of lambda values to randomly select
# @param p.rand.thr Proportion of thresholds to randomly select
# @param tol.abs ADMM absolute tolerance
# @param tol.rel ADMM relative tolerance
# @param iteration Simulation iteration index; NULL for non-simulation analyses
# @param output_path Directory path for saving results
# @param name_output Prefix for saved result files
# @param pre_screen Logical; apply pre-screening using screening_matrix
# @param pre_screen_threshold Threshold for pre-screening; auto-computed if NULL
# @param screening_matrix Pre-computed node correlation matrix; required if pre_screen=TRUE
# @param verbose Logical; print progress messages
# @return Invisibly returns list(N.hat.optimal, P.frob, P.values, threshold)
run_node_computation <- function(
  j,
  scores_df,
  covariates_df        = NULL,
  n_basis,
  L,
  K,
  thres.ctrl,
  p.rand.lam,
  p.rand.thr,
  tol.abs              = 0.0001,
  tol.rel              = 0.0001,
  iteration            = NULL,
  output_path,
  name_output,
  pre_screen           = FALSE,
  pre_screen_threshold = NULL,
  screening_matrix     = NULL,
  verbose              = TRUE
) {
  len.t <- length(thres.ctrl)
  n     <- nrow(scores_df)
  M     <- n_basis
  p     <- ceiling(ncol(scores_df) / M)

  if (is.null(covariates_df)) {
    iU <- data.frame(rep(1, n))
    colnames(iU) <- "(Intercept)"
    q <- 0
  } else {
    numeric_columns <- sapply(covariates_df, is.numeric)
    covariates_df[, numeric_columns] <- scale(covariates_df[, numeric_columns])
    C  <- model.matrix(~ ., covariates_df)
    q  <- ncol(C) - 1
    iU <- C
  }

  interM      <- data.frame(Intercept = rep(1, n))
  temp_groups <- c(0)

  if (verbose) cat("Computation of the design matrix\n")
  for (c in 1:(q + 1)) {
    product <- scores_df * iU[, c]
    if (c != 1) {
      colnames(product) <- paste(colnames(iU)[c], ":", colnames(product), sep = "")
    }
    for (i in 1:p) temp_groups <- c(temp_groups, rep(i + (p * (c - 1)), M))
    interM <- cbind(interM, product)
  }
  interM      <- interM[, -1]
  temp_groups <- temp_groups[-1]

  cat(paste("Processing node", j, "\n"))
  jth.range.y <- (j - 1) * M + (1:M)

  scr_index <- rep(FALSE, p)
  if (pre_screen) {
    if (is.null(pre_screen_threshold)) {
      pre_screen_threshold <- quantile(abs(screening_matrix), probs = 0.2)
    }
    for (it in 1:p) {
      if (screening_matrix[j, it] > pre_screen_threshold) scr_index[it] <- TRUE
    }
  } else {
    scr_index    <- rep(TRUE, p)
    scr_index[j] <- FALSE
  }

  num_nodes_included_j <- sum(scr_index)
  num_cov_j            <- num_nodes_included_j * M * (q + 1)

  d.array <- matrix(1, nrow = p, ncol = num_cov_j)
  d_out   <- list()
  for (k in 1:p) d_out[[k]] <- d.array[k, ]

  P.def <- matrix(0,    num_cov_j, M)
  Q.def <- matrix(0.1,  num_cov_j, M)
  U.def <- matrix(0.01, num_cov_j, M)

  A.Y         <- as.matrix(interM[, jth.range.y])
  jth.range.x <- rep(scr_index, each = M)
  jth.range.x <- rep(jth.range.x, times = (q + 1))
  A.X         <- as.matrix(interM[, jth.range.x])
  groups      <- temp_groups[jth.range.x]

  P <- P.def; Q <- Q.def; U <- U.def

  SCV.mat    <- matrix(NA, ceiling(L * p.rand.lam), ceiling(len.t * p.rand.thr))
  lambda.max <- lambda.sup(A.X, A.Y)
  lambdas    <- exp(seq(log(lambda.max), log(1e-4), length.out = L))

  seed_base <- if (is.null(iteration)) 1L else as.integer(iteration)

  if (p.rand.lam == 1) {
    random.sel.lambdas <- lambdas
  } else {
    set.seed(123 * seed_base + 10 * j)
    random.sel.indexes <- sample(seq_along(lambdas), ceiling(L * p.rand.lam))
    random.sel.lambdas <- lambdas[random.sel.indexes[order(random.sel.indexes)]]
  }
  L_random <- length(random.sel.lambdas)

  for (l in 1:L_random) {
    lambda <- random.sel.lambdas[l]
    if (l %% 10 == 0) cat(paste0("Lambda ", l, " of node ", j, "\n"))

    grp.lasso.result <- tryCatch({
      ADMM.grplasso.two.groups(A.X = A.X, A.Y = A.Y, d = d_out[[j]],
                               lambda = lambda, rho.init = 1,
                               P.in = P, Q.in = Q, U.in = U,
                               tol.abs = tol.abs, tol.rel = tol.rel,
                               maxiter = 400)
    }, error = function(e) {
      message("Error in ADMM.grplasso.two.groups: ", e$message)
      return(NULL)
    })
    if (is.null(grp.lasso.result)) next
    P <- grp.lasso.result$P
    Q <- grp.lasso.result$Q
    U <- grp.lasso.result$U

    n_blocks <- nrow(P) / M
    P.frob   <- list()
    for (k in 1:n_blocks) {
      key <- unique(groups[(k - 1) * M + (1:M)])
      if (length(key) > 1) {
        message("Error in group definition for block ", k)
        next
      } else {
        P.frob[[key]] <- norm(P[(k - 1) * M + (1:M), ], "F")
      }
    }

    if (p.rand.thr == 1) {
      random.sel.thresholds <- thres.ctrl
    } else {
      set.seed(234 * seed_base + 10 * j)
      random.sel.thresholds <- thres.ctrl[sample(seq_along(thres.ctrl),
                                                  ceiling(len.t * p.rand.thr))]
    }
    len.t.random <- length(random.sel.thresholds)

    for (ind.t in 1:len.t.random) {
      threshold <- lambda * random.sel.thresholds[ind.t]
      N.hat.jlt <- rep(FALSE, length(P.frob))
      for (n_block in 1:length(P.frob)) {
        if (!is.null(P.frob[[n_block]])) {
          if (P.frob[[n_block]] > threshold) N.hat.jlt[n_block] <- TRUE
        }
      }
      if (length(N.hat.jlt) < max(temp_groups)) {
        N.hat.jlt <- c(N.hat.jlt, rep(FALSE, max(temp_groups) - length(N.hat.jlt)))
      }
      card.N.hat.jlt <- sum(N.hat.jlt)
      slice.pM       <- rep(N.hat.jlt, each = M)
      A.X.eff        <- as.matrix(interM[, slice.pM])

      SCV.single <- rep(NA, K)
      for (k in 1:K) {
        fold.k.ind  <- seq(floor((k - 1) * n / K) + 1, floor(k * n / K))
        fold.k.size <- length(fold.k.ind)
        A.Y.cv    <- A.Y[fold.k.ind, ]
        A.Y.train <- A.Y[-fold.k.ind, ]

        if (card.N.hat.jlt > 0) {
          A.X.cv    <- A.X.eff[fold.k.ind, ]
          A.X.train <- A.X.eff[-fold.k.ind, ]
          B.tilde   <- solve(t(A.X.train) %*% A.X.train +
                             0.1 * diag(card.N.hat.jlt * M)) %*%
                       t(A.X.train) %*% A.Y.train
          residual  <- A.Y.cv - A.X.cv %*% B.tilde
        } else {
          residual  <- A.Y.cv
        }
        SCV.single[k] <- norm(residual, "F")^2
      }
      SCV.mat[l, ind.t] <- mean(SCV.single)
    }
  }

  scv.min        <- which(SCV.mat == min(SCV.mat), arr.ind = TRUE)
  index.optimal  <- scv.min[nrow(scv.min), ]
  l.optimal      <- index.optimal[1]
  ind.t.optimal  <- index.optimal[2]
  lambda.optimal <- random.sel.lambdas[l.optimal]
  t.optimal      <- random.sel.thresholds[ind.t.optimal]
  cat("Optimal lambda:", lambda.optimal, "\n")
  cat("Optimal threshold:", t.optimal, "\n")

  cat("Computing optimal neighbours\n")
  grp.lasso.result <- tryCatch({
    ADMM.grplasso.two.groups(A.X = A.X, A.Y = A.Y, d = d_out[[j]],
                             lambda = lambda.optimal, rho.init = 1,
                             P.in = P, Q.in = Q, U.in = U,
                             tol.abs = tol.abs, tol.rel = tol.rel,
                             maxiter = 400)
  }, error = function(e) {
    message("Error in ADMM.grplasso.two.groups: ", e$message)
    return(NULL)
  })
  P <- grp.lasso.result$P
  Q <- grp.lasso.result$Q
  U <- grp.lasso.result$U

  n_blocks <- nrow(P) / M
  P.frob   <- list()
  for (k in 1:n_blocks) {
    key <- unique(groups[(k - 1) * M + (1:M)])
    if (length(key) > 1) {
      message("Error in group definition for block ", k)
      next
    } else {
      P.frob[[key]] <- norm(P[(k - 1) * M + (1:M), ], "F")
    }
  }

  P.values <- list()
  for (k in 1:n_blocks) {
    key <- unique(groups[(k - 1) * M + (1:M)])
    if (length(key) > 1) {
      message("Error in group definition for block ", k)
      next
    } else {
      P.values[[key]] <- P[(k - 1) * M + (1:M), ]
    }
  }

  N.hat.optimal <- rep(FALSE, length(P.frob))
  for (n_block in 1:length(P.frob)) {
    if (!is.null(P.frob[[n_block]])) {
      threshold <- t.optimal * lambda.optimal
      if (P.frob[[n_block]] > threshold) N.hat.optimal[n_block] <- TRUE
    }
  }
  if (length(N.hat.optimal) < max(temp_groups)) {
    N.hat.optimal <- c(N.hat.optimal,
                       rep(FALSE, max(temp_groups) - length(N.hat.optimal)))
  }
  cat("Optimal neighbours of node", j, ":\n")
  cat(N.hat.optimal)
  cat("\n")

  full_result_path      <- paste(output_path, name_output, "_", j, ".rda", sep = "")
  full_result_path_coef <- paste(output_path, name_output, "_", j, "coeff.rda", sep = "")
  save(N.hat.optimal, file = full_result_path)
  save(N.hat.optimal, P.frob, P.values, threshold, file = full_result_path_coef)
  cat("Output saved to:", full_result_path, "\n")

  invisible(list(N.hat.optimal = N.hat.optimal, P.frob = P.frob,
                 P.values = P.values, threshold = threshold))
}
