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


args <- commandArgs(trailingOnly = TRUE)
save.path <- if (length(args) > 0) args[[1]] else
  stop("Please provide the path where results are saved.")

time.start <- Sys.time()
print("STARTING: Simulation of the data")

#############################################
###### PART 1: Create a dataframe for each sim setting
#############################################

sim_dirs <- list.dirs(save.path, recursive = FALSE)
for(sim_dir in sim_dirs){
  iter_sim_dirs <- list.dirs(sim_dir, recursive = FALSE)
  for(iter_sim_dir in iter_sim_dirs){
    cat(paste0("Processing directory: ", iter_sim_dir, "\n"))
    yaml_iter_file_path = file.path(iter_sim_dir, "config_single_iter.yaml")
    config <- yaml.load_file(yaml_iter_file_path)
    
    func.path <- config$func_path
    source(file.path(func.path, "bases.func.R"))
    source(file.path(func.path, "Sim.prec.matrix.func.R"))
    
    iteration <- config$iteration
    seed_g1 <- config$seed_g1
    seed_g2 <- config$seed_g2
    name_output <- config$name_output
    
    p <- config$p
    n_g1 <- config$n_g1
    n_g2 <- config$n_g2
    mu <- config$mu
    tau <- config$tau 
    model_g1 <- config$model_g1
    model_g2 <- config$model_g2
    model_continuous <- config$model_continuous
    var_cont_type <- config$var_cont_type
    red.number <- config$red_number
    
    rec_basis_type <- config$rec_basis_type
    rec_basis_number <- config$rec_basis_number
    M <- config$n_basis_for_dim_reduction
    
    ####################################
    #     PART 1: DATA GENERATION      
    ####################################
    foldname <- file.path(iter_sim_dir, "simulated_data")
    if (!dir.exists(foldname)) {
      dir.create(foldname, recursive = TRUE)
    }
    
    # Generate Random Functions and Observation h_ijk       
    # h_ijk = Fourier basis func(t_k) %*% delta_ij + iid error
    
    ####################################
    #     PART 1.1: PRECISION MATRIX GENERATION      
    ####################################
    
    # Precision matrix lined to the continuous variable
    set.seed((seed_g1+seed_g2)*iteration)
    Theta.c <- match.fun(model_continuous)(p, mu, red.number)
    adj_df <- melt(Theta.c)
    colnames(adj_df) <- c("Row", "Col", "Value")
    
    plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "black", limits = c(0, 1)) +
      theme_minimal() +
      coord_fixed() +
      scale_y_reverse() + 
      labs(title = "Score adjacency Matrix", x = "", y = "")
    ggsave(plot,
           filename = paste0(foldname, "/Adjacency_matrix_scores_cont_",
                             name_output, ".png"),
           width = 8, height = 8, dpi = 300)
    
    G.true.cont <- matrix(0, p, p)
    for(i in 1:p){
      for(j in 1:p){
        if(sum(abs(Theta.c[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])) > 0)
          G.true.cont[i,j] <- frobenius.norm(Theta.c[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])
      }
    }
    adj_df <- melt(G.true.cont)
    colnames(adj_df) <- c("Row", "Col", "Value")
    plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "black") +
      theme_minimal() +
      coord_fixed() +
      scale_y_reverse() + 
      labs(title = "Node adjacency Matrix", x = "", y = "")
    ggsave(plot,
           filename = paste0(foldname, "/Adjacency_matrix_node_cont_",
                             name_output, ".png"),
           width = 8, height = 8, dpi = 300)
    
    # Precision matrix lined to the group G1
    Theta.g1 <- match.fun(model_g1)(p, mu, red.number)
    adj_df <- melt(Theta.g1)
    colnames(adj_df) <- c("Row", "Col", "Value")
    
    plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "black") +
      theme_minimal() +
      coord_fixed() +
      scale_y_reverse() + 
      labs(title = "Score adjacency Matrix", x = "", y = "")
    ggsave(plot,
           filename = paste0(foldname, "/Adjacency_matrix_scores_g1_",
                             name_output, ".png"),
           width = 8, height = 8, dpi = 300)
    
    G.true.g1 <- matrix(0, p, p) 
    for(i in 1:p){
      for(j in 1:p){
        if(sum(abs(Theta.g1[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])) > 0)
          G.true.g1[i,j] <- frobenius.norm(Theta.g1[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])
      }
    }
    adj_df <- melt(G.true.g1)
    colnames(adj_df) <- c("Row", "Col", "Value")
    plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "black") +
      theme_minimal() +
      coord_fixed() +
      scale_y_reverse() + 
      labs(title = "Node adjacency Matrix", x = "", y = "")
    ggsave(plot,
           filename = paste0(foldname, "/Adjacency_matrix_node_g1_",
                             name_output, ".png"),
           width = 8, height = 8, dpi = 300)
    
    # Precision matrix lined to the group G2
    Theta.g2 <- match.fun(model_g2)(p, mu, red.number)
    adj_df <- melt(Theta.g2)
    colnames(adj_df) <- c("Row", "Col", "Value")

    plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "black") +
      theme_minimal() +
      coord_fixed() +
      scale_y_reverse() +  
      labs(title = "Score adjacency Matrix", x = "", y = "")
    ggsave(plot,
           filename = paste0(foldname, "/Adjacency_matrix_scores_g2_",
                             name_output, ".png"),
           width = 8, height = 8, dpi = 300)
    
    G.true.g2 <- matrix(0, p, p) 
    for(i in 1:p){
      for(j in 1:p){
        if(sum(abs(Theta.g2[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])) > 0)
          G.true.g2[i,j] <- frobenius.norm(Theta.g2[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])
      }
    }
    adj_df <- melt(G.true.g2)
    colnames(adj_df) <- c("Row", "Col", "Value")
    plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "black") +
      theme_minimal() +
      coord_fixed() +
      scale_y_reverse() + 
      labs(title = "Node adjacency Matrix", x = "", y = "")
    ggsave(plot,
           filename = paste0(foldname, "/Adjacency_matrix_node_g2_",
                             name_output, ".png"),
           width = 8, height = 8, dpi = 300)
    
    ####################################
    #     PART 1.2: DATA GENERATION      
    ####################################  
    set.seed(seed_g1*iteration)
    if(var_cont_type=="pos_cont"){
      var_cont_g1 <- runif(n_g1, min = 0, max=1)
    }else{
      var_cont_g1 <- runif(n_g1, min = -1, max=1)
    }
    delta <- matrix(nrow=n_g1, ncol = ncol(Theta.g1))
    for(sample in 1:n_g1){
      delta[sample,] <- rmvnorm(1, sigma = solve(Theta.g1 + var_cont_g1[sample]*Theta.c))
    }
    
    obs.time <- seq(1/tau, 1, 1/tau) # vector of observation time points of delta
    b.mat.list <- list()
    for(j in 1:p){
      b.mat.list[[j]] <- fda.fourier.mat(obs.time, mu)
    }
    
    h <- array(0, c(n_g1, p, tau))
    for(i in 1:n_g1){
      for(j in 1:p){
        set.seed(seed_g1*iteration + j)
        h[i,j,] <- b.mat.list[[j]] %*% 
          matrix(delta[i, ((j-1)*mu+1) : (j*mu)], ncol=1) + rnorm(tau, 0, 0.5)
      }
    }
    
    h_g1 <- h
    png(paste0(foldname, "/Data_simulated_g1_", name_output, ".png"),
        width = 1200, height = 800, res = 150)
    matplot(t(h_g1[,1,]), type = "l", lty = 1, col = 1:6)
    dev.off()
    
    set.seed(seed_g2*iteration)
    if(var_cont_type=="pos_cont"){
      var_cont_g2 <- runif(n_g2, min = 0, max=1)
    }else{
      var_cont_g2 <- runif(n_g2, min = -1, max=1)
    }
    delta <- matrix(nrow=n_g2, ncol = ncol(Theta.g2))
    for(sample in 1:n_g2){
      sample_sig <- solve(Theta.g2 + var_cont_g2[sample]*Theta.c)
      delta[sample,] <- rmvnorm(1, sigma = sample_sig)
    }
    
    obs.time <- seq(1/tau, 1, 1/tau) # vector of observation time points of delta
    b.mat.list <- list()
    for(j in 1:p){
      b.mat.list[[j]] <- fda.fourier.mat(obs.time, mu)
    }
    
    h <- array(0, c(n_g2, p, tau))
    for(i in 1:n_g2){
      for(j in 1:p){
        set.seed(seed_g2*iteration + j)
        h[i,j,] <- b.mat.list[[j]] %*% 
          matrix(delta[i, ((j-1)*mu+1) : (j*mu)], ncol=1) + rnorm(tau, 0, 0.5)
      }
    }
    
    h_g2 <- h
    png(paste0(foldname, "/Data_simulated_g2_", name_output, ".png"),
        width = 1200, height = 800, res = 150)
    matplot(t(h_g2[,1,]), type = "l", lty = 1, col = 1:6)
    dev.off()
    
    h <- abind(h_g1, h_g2, along = 1)
    group=c(rep(1,n_g1), rep(2,n_g2))
    var_cont <- c(var_cont_g1, var_cont_g2)
    save(h, group, var_cont,
         file = paste0(foldname, "/Original_data_generated_", name_output, ".RData"))
    
    save(G.true.g1, G.true.g2, G.true.cont,
         file = paste0(iter_sim_dir, "/Ground_truth_", name_output, ".rda"))
    
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
    
    covariates_df <- data.frame(group = group, var_cont = var_cont)
    covariates_df$group <- factor(covariates_df$group, levels = c(1, 2))
    scores_df <- as.data.frame(fpc.score)
    save(scores_df, covariates_df,
         file = paste0(iter_sim_dir, "/net_est_input_", name_output, ".rda"))
    output_path = paste0(iter_sim_dir, "/results/")
    if (!dir.exists(output_path)) {
      dir.create(output_path)
    }
  }
}
print("END: Full simulation of the data completed in: ")
Sys.time()-time.start
