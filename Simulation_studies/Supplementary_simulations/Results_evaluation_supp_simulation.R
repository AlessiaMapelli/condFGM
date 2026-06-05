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

save_path = config$save_path
simulation_name = config$simulation_name
tot_iteration = config$tot_iteration
name_output = config$name_output
num_nodes <- config$p
n_comp <- config$n_basis_for_dim_reduction
n_groups <- 2
n_cont <- 1
n_g1 <- config$n_g1
n_g2 <- config$n_g2
func.path <- config$func_path
source(file.path(func.path, "prec.rec.R"))
source(file.path(func.path, "plot_binary_adj_matrix.R"))  

####################################
#   RESULTS EVALUATION and PLOTS   #
####################################
results_metrics <- data.frame(network=NA,	symm=NA,	prec=NA,	TPR=NA,	FPR=NA,	F1=NA,
                              Baseline=NA,	Differential=NA)
for(iteration in 1:tot_iteration){
  cat("Processing iteration ", iteration, "\n")
  output_path = paste0(save_path, simulation_name,"/", "seed_", iteration, "/", "results/")
  G.our <- matrix(NA, nrow = num_nodes, ncol = num_nodes*(n_groups+n_cont))
  for(j in 1: num_nodes){
    load(paste(output_path, name_output,"_" ,j, ".rda", sep=""))
    G.our[j, ] <- N.hat.optimal
  }
  
  foldname = paste0(save_path, simulation_name,"/", "seed_", iteration)
  load(paste(foldname,"/Ground_truth_",name_output,".rda", sep=""))

  ####### Population network #######
  G.pop <- (G.true.g1+G.true.cont>0)
  diag(G.pop) <- FALSE
  G.pop <-ifelse(G.pop, 1, 0)
  G.our.pop <- G.our[,1:num_nodes]
  G.our.pop <- ifelse(G.our.pop, 1, 0)

  plot_binary_adj_matrix(G.our.pop,
                      title= "Estimated adjacency Matrix - Population",
                      output_path=output_path,
                      plot_name="Estimated_adjacency_matrix_node_pop.png")

  
  res_and <- prec.rec(G.pop, G.our.pop, type="AND")
  results_metrics <- rbind(results_metrics,c("POP","AND",res_and$prec,
                                             res_and$TPR, res_and$FPR, res_and$F1,
                                             n_g1, n_g2))
  res_or <- prec.rec(G.pop, G.our.pop, type="OR")
  results_metrics <- rbind(results_metrics,c("POP","OR",res_or$prec,
                                             res_or$TPR, res_or$FPR, res_or$F1,
                                             n_g1, n_g2))

  #########  Differential network - group contribution ######### 
  G.diff <- (abs(G.true.g2 - G.true.g1))>0
  diag(G.diff) <- FALSE
  G.diff <-ifelse(G.diff, 1, 0)
  
  G.our.diff <- G.our[,(num_nodes+1):(num_nodes*2)]

  plot_binary_adj_matrix(G.our.diff,
                    title= "Estimated adjacency Matrix - Differential - categorical",
                    output_path=output_path,
                    plot_name="Estimated_adjacency_matrix_node_diff_categorical.png")  

  res_and <- prec.rec(G.diff, G.our.diff, type="AND")
  results_metrics <- rbind(results_metrics,c("DIFF","AND",res_and$prec,
                                             res_and$TPR, res_and$FPR, res_and$F1,
                                             n_g1, n_g2))
  res_or <- prec.rec(G.diff, G.our.diff, type="OR")
  results_metrics <- rbind(results_metrics,c("DIFF","OR",res_or$prec,
                                             res_or$TPR, res_or$FPR, res_or$F1,
                                             n_g1, n_g2))
  
  ######### Differential network - group contribution - weights ######### 
  G.our.diff.ratio.forb <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  G.our.diff.sum.forb <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  G.our.diff.sum.forb.logic <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  G.our.diff.ratio.forb.sum <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  G.our.group <- G.our.pop
  for(j in 1: num_nodes){
    for(k in 1: num_nodes){
      if(G.our.diff[j,k]){
        load(paste(output_path, name_output,"_" ,j, "coeff.rda", sep=""))
        beta.g2.forbs <- norm(P.values[[k]] + P.values[[k+num_nodes]], "F")
        G.our.group[j,k] <- ifelse(beta.g2.forbs > threshold,1,0)
        if(P.frob[[k]] == 0){
          cat("Warning: Zero norm at denominator. Setting to the numerator value")
          G.our.diff.ratio.forb[j,k] <- (P.frob[[k+num_nodes]])
          G.our.diff.ratio.forb.sum[j,k] <- beta.g2.forbs
        }else{
          G.our.diff.ratio.forb[j,k] <- (P.frob[[k+num_nodes]])/(P.frob[[k]])
          G.our.diff.ratio.forb.sum[j,k] <- beta.g2.forbs/(P.frob[[k]])
        }
        G.our.diff.sum.forb[j,k] <- beta.g2.forbs
        G.our.diff.sum.forb.logic[j,k] <- beta.g2.forbs>threshold
      }
    }
  }

  adj_df <- melt(G.our.diff.ratio.forb.sum, varnames = c("Row","Col"), value.name = "Value")
  
  p <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "grey92", linewidth = 0.15) +
    coord_fixed() +
    scale_y_reverse() +
    labs(
      title = str_wrap(
        "Direction of changes in the estimated differential adjacency matrix - categorical",
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
  
  ggsave(p,
         filename = paste0(output_path, "Estimated_adjacency_matrix_node_diff_categorical_signed_weights.png"),
         width = 8, height = 8, dpi = 300)

  ######### Group network ######### 
  
  G.group <- (G.true.g2+G.true.cont>0)
  diag(G.group) <- FALSE
  G.group <-ifelse(G.group, 1, 0)

  plot_binary_adj_matrix(G.our.group,
                  title= "Estimated adjacency Matrix - Group",
                  output_path=output_path,
                  plot_name="Estimated_adjacency_matrix_node_group_g2.png")
  
  
  res_and <- prec.rec(G.group, G.our.group, type="AND")
  results_metrics <- rbind(results_metrics,c("GROUP","AND",res_and$prec,
                                             res_and$TPR, res_and$FPR, res_and$F1,
                                             n_g1, n_g2))
  res_or <- prec.rec(G.group, G.our.group, type="OR")
  results_metrics <- rbind(results_metrics,c("GROUP","OR",res_or$prec,
                                             res_or$TPR, res_or$FPR, res_or$F1,
                                             n_g1, n_g2))
  
  #########  Differential network - continuous contribution ######### 
  G.mod <- (G.true.cont>0)
  diag(G.mod) <- FALSE
  G.mod <-ifelse(G.mod, 1, 0)
  
  G.our.mod <- G.our[,(num_nodes*2+1):(num_nodes*3)]

  plot_binary_adj_matrix(G.our.mod,
                  title= "Estimated adjacency Matrix - Differential - continuous",
                  output_path=output_path,
                  plot_name="Estimated_adjacency_matrix_node_diff_continuous.png")  
  
  res_and <- prec.rec(G.mod, G.our.mod, type="AND")
  results_metrics <- rbind(results_metrics,c("CONT","AND",res_and$prec,
                                             res_and$TPR, res_and$FPR, res_and$F1,
                                             n_g1, n_g2))
  res_or <- prec.rec(G.mod, G.our.mod, type="OR")
  results_metrics <- rbind(results_metrics,c("CONT","OR",res_or$prec,
                                             res_or$TPR, res_or$FPR, res_or$F1,
                                             n_g1, n_g2))
  
  ######### Differential network - continuous contribution - weights ######### 
  G.our.mod.ratio.forb <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  G.our.mod.sum.forb <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  G.our.mod.sum.forb.logic <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  G.our.mod.ratio.forb.sum <- matrix(NA, nrow = num_nodes, ncol = num_nodes)
  for(j in 1: num_nodes){
    for(k in 1: num_nodes){
      if(G.our.mod[j,k]){
        load(paste(output_path, name_output,"_" ,j, "coeff.rda", sep=""))
        beta.g2.forbs <- norm(P.values[[k]] + P.values[[k+2*num_nodes]], "F")
        if(P.frob[[k]] == 0){
          cat("Warning: Zero norm at denominator. Setting to the numerator value")
          G.our.mod.ratio.forb[j,k] <- (P.frob[[k+2*num_nodes]])
          G.our.mod.ratio.forb.sum[j,k] <- beta.g2.forbs
        }else{
          G.our.mod.ratio.forb[j,k] <- (P.frob[[k+2*num_nodes]])/(P.frob[[k]])
          G.our.mod.ratio.forb.sum[j,k] <- beta.g2.forbs/(P.frob[[k]])
        }
        G.our.mod.sum.forb[j,k] <- beta.g2.forbs
        G.our.mod.sum.forb.logic[j,k] <- beta.g2.forbs>threshold
      }
    }
  }
  
  adj_df <- melt(G.our.mod.ratio.forb.sum, varnames = c("Row","Col"), value.name = "Value")
  
  p <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "grey92", linewidth = 0.15) +
    coord_fixed() +
    scale_y_reverse() +
    labs(
      title = str_wrap(
        "Direction of changes in the estimated differential adjacency matrix - continuous",
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
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 1,
      na.value = "white",
      name = "Value"
    )
  ggsave(p,
         filename = paste0(output_path, "Estimated_adjacency_matrix_node_diff_continuous_signed_weights.png"),
         width = 8, height = 8, dpi = 300)
  
  save(G.our.pop, G.our.diff, G.our.mod, G.our.mod.ratio.forb.sum,
       G.our.group, G.our.diff.ratio.forb.sum,
       file = paste0(output_path, "Estimated_adjacency_matrices.RData"))
  
}
results_metrics <- results_metrics[-1,]
out_csv <- paste0(save_path, simulation_name,
                  "/Results_performance_metrics_", name_output, ".csv")
write.csv(results_metrics, out_csv)
cat("Results saved to ", out_csv, "\n")
