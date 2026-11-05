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

source(file.path("Computation_templates", "network_functions.R"))

####################################
#   RESULTS EVALUATION and PLOTS   #
####################################

load(input_path)

# Columns of the design matrix: 1 = population (intercept), then one column per
# numeric covariate / factor dummy. Group-specific networks are only computed
# for factor dummies (for a numeric covariate B0 + Bc is the network at +1 SD).
if ((exists("covariates_df", envir = .GlobalEnv) && is.data.frame(covariates_df))){
  C <- model.matrix(~ ., covariates_df)
  n_groups <- ncol(C)
  cov_names <- colnames(C)
  cov_names[1] <- "population"
  cov_assign <- attr(C, "assign")
  is_factor_col <- sapply(seq_len(n_groups), function(i)
    cov_assign[i] > 0 && !is.numeric(covariates_df[[cov_assign[i]]]))
} else{ n_groups <- 1
cov_names <- "population"
is_factor_col <- FALSE }


G.our <- matrix(NA, nrow = num_nodes, ncol = num_nodes*n_groups)
for(j in 1: num_nodes){
    load(paste(output_path, name_output,"_" ,j, ".rda", sep=""))
    G.our[j, ] <- N.hat.optimal
}
G.pop <- G.our[, 1:num_nodes]

G.our.symm = list()
G.our.symm.weighted = list()
G.our.symm.warning= list()
G.our.symm.not.in.baseline = list()
G.our.weighted= list()
G.our.not.in.baseline = list()
G.our.group = list()
G.our.group.symm = list()

plot_binary <- function(G, plot_title, fname) {
  adj_df <- melt(G)
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
      name = "Edge",
      drop = FALSE
    )
  ggsave(plot, filename = fname, width = 8, height = 8, dpi = 300)
}

for(i in 1:n_groups){
  G.cov <- G.our[,(num_nodes*(i-1) +1):(num_nodes*i)]
  G.cov.simm <- symmetrise_adjacency(G.cov, type)

  if(i==1){
    plot_title <- paste0("Node adjacency Matrix - ", cov_names[i])
  } else {
    plot_title <- paste0("Node adjacency Matrix - Differential contribution ", cov_names[i])
  }

  G.our.symm[[cov_names[i]]] <- G.cov.simm
  plot_binary(G.cov.simm, plot_title,
              paste0(output_path, name_output, "_Adjacency_matrix_node_", cov_names[i], ".png"))

  if(i !=1){
    # Directional VR (Eq. 5.16): defined only where the baseline block is
    # selected (||B0_k||_F > threshold); otherwise the edge is "not in baseline".
    G.cov.weighted <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
    G.cov.nib      <- matrix(FALSE, nrow = num_nodes, ncol = num_nodes)
    G.cov.group    <- G.pop
    for(j in 1: num_nodes){
      load(paste(output_path, name_output,"_" ,j, "coeff.rda", sep=""))
      for(k in 1: num_nodes){
        if(G.cov[j,k]){
          d <- directional_vr(P.values, P.frob, k, num_nodes, i, threshold)
          G.cov.weighted[j,k] <- d$vr
          G.cov.nib[j,k]      <- !d$baseline_selected
          # Group network, as in the simulation code: where the differential
          # block is selected, the edge is present in the group if ||B0+Bc|| > threshold
          G.cov.group[j,k]    <- ifelse(d$group_norm > threshold, 1, 0)
        }
      }
    }

    # Symmetric VR (Eq. 5.17): geometric mean of the directional VRs over the
    # directions where both the differential and the baseline block are
    # selected; "not in baseline" if there is no such direction.
    G.cov.simm.weighted <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
    G.cov.simm.warning  <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
    G.cov.simm.nib      <- matrix(FALSE, nrow = num_nodes, ncol = num_nodes)
    for (j in 1:(num_nodes - 1)) {
      for (k in (j + 1):num_nodes) {
        if(G.cov.simm[j,k]){
          x <- c(G.cov.weighted[j,k], G.cov.weighted[k,j])
          x <- x[!is.na(x)]
          if (length(x) == 0) {
            G.cov.simm.nib[j,k] <- G.cov.simm.nib[k,j] <- TRUE
            G.cov.simm.warning[j,k] <- G.cov.simm.warning[k,j] <- 0
          } else {
            w <- exp(mean(log(x)))
            G.cov.simm.weighted[j,k] <- G.cov.simm.weighted[k,j] <- w
            # Warning: the two directions disagree on the sign of the change
            disagree <- length(x) == 2 && !((x[1] > 1 & x[2] > 1) | (x[1] < 1 & x[2] < 1))
            G.cov.simm.warning[j,k] <- G.cov.simm.warning[k,j] <- as.numeric(disagree)
          }
        }
      }
    }

    G.our.symm.weighted[[cov_names[i]]] <- G.cov.simm.weighted
    G.our.symm.warning[[cov_names[i]]] <- G.cov.simm.warning
    G.our.symm.not.in.baseline[[cov_names[i]]] <- G.cov.simm.nib
    G.our.weighted[[cov_names[i]]] <- G.cov.weighted
    G.our.not.in.baseline[[cov_names[i]]] <- G.cov.nib

    plot_vr_matrix(
      G.cov.simm.weighted, G.cov.simm.nib,
      title = paste0("Direction of changes in the estimated adjacency matrix - ", cov_names[i]),
      filename = paste0(output_path, name_output,
                        "_Adjacency_matrix_node_diff_weights_", cov_names[i], ".png")
    )

    if (is_factor_col[i]) {
      G.our.group[[cov_names[i]]] <- G.cov.group
      G.cov.group.simm <- symmetrise_adjacency(G.cov.group, type)
      G.our.group.symm[[cov_names[i]]] <- G.cov.group.simm
      plot_binary(G.cov.group.simm,
                  paste0("Node adjacency Matrix - Group ", cov_names[i]),
                  paste0(output_path, name_output, "_Adjacency_matrix_node_group_", cov_names[i], ".png"))
    }
  }
}

save(G.our.symm, G.our.symm.weighted, G.our.symm.warning, G.our.symm.not.in.baseline,
     G.our.weighted, G.our.not.in.baseline, G.our.group, G.our.group.symm,
     file = paste0(output_path, name_output, "_Adj_estimation.rda"))
cat("Results saved to ", paste(output_path, name_output, "_Adj_estimation.rda", sep=""), "\n")
