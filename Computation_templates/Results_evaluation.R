# Required packages
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))

# Load parameters from YAML configuration file
args <- commandArgs(trailingOnly = TRUE)
yaml_file_path <- if (length(args) > 0) args[[1]] else
  stop("Please provide the path to the YAML configuration file as an argument.")
config <- yaml.load_file(yaml_file_path)

output_path = config$output_path
name_output = config$name_output
input_path = config$input_path
type = config$type

num_nodes = config$n_nodes
n_comp <- config$n_basis_for_dim_reduction

####################################
#   RESULTS EVALUATION and PLOTS   #
####################################

load(input_path)

if ((exists("covariates_df", envir = .GlobalEnv) && is.data.frame(covariates_df))){
  C <- model.matrix(~ ., covariates_df)
  n_groups <- ncol(C)
  cov_names <- colnames(C)
  cov_names[1] <- "population"
} else{ n_groups <- 1
cov_names <- "population" }


G.our <- matrix(NA, nrow = num_nodes, ncol = num_nodes*n_groups)
for(j in 1: num_nodes){
    load(paste(output_path, name_output,"_" ,j, ".rda", sep=""))
    G.our[j, ] <- N.hat.optimal
}

G.our.symm = list()
G.our.symm.weighted = list()
G.our.symm.warning= list()
G.our.weighted= list()


for(i in 1:n_groups){
  G.cov <- G.our[,(num_nodes*(i-1) +1):(num_nodes*i)]
  G.cov.simm <- matrix(0, nrow = num_nodes, ncol = num_nodes)
  if(type=="AND"){
    for(j in 1:(num_nodes-1)){
      for(k in (i+1):num_nodes){
        if(G.cov[j,k] == 1 & G.cov[k,j] == 1){
          G.cov.simm[j,k] = 1
          G.cov.simm[k,j] = 1
        }
      }
    }
  }else if(type=="OR"){
    for(j in 1:(num_nodes-1)){
      for(k in (i+1):num_nodes){
        if(G.cov[j,k] == 1 | G.cov[k,j] == 1){
          G.cov.simm[j,k] = 1
          G.cov.simm[k,j] = 1
        }
      }
    }
  }else{
    G.cov.simm = G.cov
  }
  
  if(i==1){
    plot_title <- paste0("Node adjacency Matrix - ", cov_names[i])
  } else {
    plot_title <- paste0("Node adjacency Matrix - Differential contribution ", cov_names[i])
  }
  

  G.our.symm[[cov_names[i]]] <- G.cov.simm

  adj_df <- melt(G.cov.simm)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "grey92", linewidth = 0.15) +
    coord_fixed() +
    scale_y_reverse() +
    labs(
      title = str_wrap(plot_title, width = 42),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(
        size = 20,
        hjust = 0.5,
        vjust =-1.5,
        lineheight = 0.95 
      ),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.35, "cm"),
      plot.margin = margin(4, 4, 4, 4)
    ) +
    aes(fill = factor(Value, levels = c(0, 1))) + scale_fill_manual(
      values = c( "0" = "white","1" = "black"),
      labels = c("0" = "Absent", "1" = "Present"),
      name = "Edge"
    )
  fname <- paste0(output_path, name_output, "_Adjacency_matrix_node_", cov_names[i], ".png")
  ggsave(plot, filename = fname, width = 8, height = 8, dpi = 300)


  if(i !=1){
    G.cov.weighted <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
    for(j in 1: num_nodes){
      for(k in 1: num_nodes){
        if(G.cov[j,k]){
          load(paste(output_path, name_output,"_" ,j, "coeff.rda", sep=""))
          beta.g2.forbs <- norm(P.values[[k]] + P.values[[k+(num_nodes*(i-1))]], "F")
          if(P.frob[[k]] == 0){
            cat("Warning: Zero norm at denominator. Setting to the numerator value")
            G.cov.weighted[j,k] <- beta.g2.forbs
          }else{
            G.cov.weighted[j,k] <- beta.g2.forbs/(P.frob[[k]])
          }
        }
      }
    }
    G.cov.simm.weighted <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
    G.cov.simm.warning <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
    for (j in 1:(num_nodes - 1)) {
      for (k in (j + 1):num_nodes) {
        if(G.cov.simm[j,k]){
          load(paste(output_path, name_output,"_" ,j, "coeff.rda", sep=""))
          beta.g2.forbs.jk <- norm(P.values[[k]] + P.values[[k+(num_nodes*(i-1))]], "F")
          if(P.frob[[k]] == 0){
            cat("Warning: Zero norm at denominator. Setting to the numerator value")
            ratio.beta.forbs.jk <- beta.g2.forbs.jk
          }else{
            ratio.beta.forbs.jk <- beta.g2.forbs.jk/(P.frob[[k]])
          }
          load(paste(output_path, name_output,"_" ,k, "coeff.rda", sep=""))
          beta.g2.forbs.kj <- norm(P.values[[j]] + P.values[[j+(num_nodes*(i-1))]], "F")
          if(P.frob[[j]] == 0){
            cat("Warning: Zero norm at denominator. Setting to the numerator value")
            ratio.beta.forbs.kj <- beta.g2.forbs.kj
          }else{
            ratio.beta.forbs.kj <- beta.g2.forbs.kj/(P.frob[[j]])
          }
          if(G.cov[j,k] & G.cov[k,j]){
            
            if(ratio.beta.forbs.jk > 1 & ratio.beta.forbs.kj > 1){
              G.cov.simm.warning[j,k] <- 0
              G.cov.simm.warning[k,j] <- 0
            } else if (ratio.beta.forbs.jk < 1 & ratio.beta.forbs.kj < 1){
              G.cov.simm.warning[j,k] <- 0
              G.cov.simm.warning[k,j] <- 0
            } else{
              G.cov.simm.warning[j,k] <- 1
              G.cov.simm.warning[k,j] <- 1
            }
            x <- c(ratio.beta.forbs.jk, ratio.beta.forbs.kj)
            G.cov.simm.weighted[j,k] <- exp(mean(log(x)))
            G.cov.simm.weighted[k,j] <- exp(mean(log(x)))
            
          }else if (G.cov[j,k]){
            G.cov.simm.warning[j,k] <- 0
            G.cov.simm.warning[k,j] <- 0
            G.cov.simm.weighted[j,k] <- ratio.beta.forbs.jk
            G.cov.simm.weighted[k,j] <- ratio.beta.forbs.jk
          }else{
            G.cov.simm.warning[j,k] <- 0
            G.cov.simm.warning[k,j] <- 0
            G.cov.simm.weighted[j,k] <- ratio.beta.forbs.kj
            G.cov.simm.weighted[k,j] <- ratio.beta.forbs.kj
          }
          
        }
      }
    }
  
  G.our.symm.weighted[[cov_names[i]]] <- G.cov.simm.weighted
  G.our.symm.warning[[cov_names[i]]] <- G.cov.simm.warning
  G.our.weighted[[cov_names[i]]] <- G.cov.weighted

  adj_df <- melt(G.cov.simm.weighted, varnames = c("Row","Col"), value.name = "Value")
  plot<- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "grey92", linewidth = 0.15) +
    coord_fixed() +
    scale_y_reverse() +
    labs(
      title = str_wrap(
        paste0("Direction of changes in the estimated adjacency matrix - ", cov_names[i]),
        width = 42
      ),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(
        size = 20,
        hjust = 0.5,
        vjust = -1.5,
        lineheight = 0.95 
      ),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.35, "cm"),
      plot.margin = margin(4, 4, 4, 4)
    ) +
    scale_fill_gradient2(
      low = "darkblue",
      mid = "white",
      high = "darkred",
      midpoint = 1,
      na.value = "white",
      name = "Value"
    )
  fname_w <- paste0(output_path, name_output,
                    "_Adjacency_matrix_node_diff_weights_", cov_names[i], ".png")
  ggsave(plot, filename = fname_w, width = 8, height = 8, dpi = 300)
}
  
}

save(G.our.symm, G.our.symm.weighted, G.our.symm.warning, G.our.weighted,
     file = paste0(output_path, name_output, "_Adj_estimation.rda"))
cat("Results saved to ", paste(output_path, name_output, "_Adj_estimation.rda", sep=""), "\n")