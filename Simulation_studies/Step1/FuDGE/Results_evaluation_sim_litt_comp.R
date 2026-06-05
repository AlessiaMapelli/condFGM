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
n_g1 <- config$n_g1
n_g2 <- config$n_g2
func.path <- config$func_path
source(file.path(func.path, "prec.rec.R"))
source(file.path(func.path, "plot_binary_adj_matrix.R"))

####################################
#   RESULTS EVALUATION and PLOTS   #
####################################


results_metrics <- data.frame(network=NA,	method=NA, hyper=NA,	prec=NA,	TPR=NA,	FPR=NA,	F1=NA,
                              Baseline=NA,	Differential=NA)
for(iteration in 1:tot_iteration){
  cat("Processing iteration ", iteration, "\n")
  output_path <- paste0(save_path, simulation_name,
                        "/seed_", iteration, "/results_lit_comparison/")
  load(paste(output_path, name_output,"_litt_comparison", ".rda", sep=""))
  
  foldname = paste0(save_path, simulation_name,"/", "seed_", iteration)
  load(paste(foldname,"/Ground_truth_",name_output,".rda", sep=""))
  
  ###### Evaluate over all tested hyperparameters + select best and worst case ######
  G.pop <- (G.true.g1>0)
  diag(G.pop) <- FALSE
  G.pop <-ifelse(G.pop, 1, 0)
  
  G.diff <- (abs(G.true.g2 - G.true.g1))>0
  diag(G.diff) <- FALSE
  G.diff <-ifelse(G.diff, 1, 0)
  
  G.group <- (G.true.g2>0)
  diag(G.group) <- FALSE
  G.group <-ifelse(G.group, 1, 0)
  
  # Store FuDGE results
  G.FuDGE.diff.opt <- NULL
  G.FuDGE.diff.worse <- NULL
  F1.FuDGE.opt <- 0
  F1.FuDGE.worse <- 1
  
  # Store FGL results
  G.FGL.diff.opt <- NULL
  G.FGL.diff.worse <- NULL
  F1.FGL.diff.opt <- 0
  F1.FGL.diff.worse <- 1
  
  G.FGL.pop.opt <- NULL
  G.FGL.pop.worse <- NULL
  F1.FGL.pop.opt <- 0
  F1.FGL.pop.worse <- 1

  G.FGL.group.opt <- NULL
  G.FGL.group.worse <- NULL
  F1.FGL.group.opt <- 0
  F1.FGL.group.worse <- 1

  # Store FFGL results
  G.FFGL.diff.opt <- NULL
  G.FFGL.diff.worse <- NULL
  F1.FFGL.diff.opt <- 0
  F1.FFGL.diff.worse <- 1
  
  G.FFGL.pop.opt <- NULL
  G.FFGL.pop.worse <- NULL
  F1.FFGL.pop.opt <- 0
  F1.FFGL.pop.worse <- 1

  G.FFGL.group.opt <- NULL
  G.FFGL.group.worse <- NULL
  F1.FFGL.group.opt <- 0
  F1.FFGL.group.worse <- 1

  # Store FFGL2 results
  G.FFGL2.diff.opt <- NULL
  G.FFGL2.diff.worse <- NULL
  F1.FFGL2.diff.opt <- 0
  F1.FFGL2.diff.worse <- 1
  
  G.FFGL2.pop.opt <- NULL
  G.FFGL2.pop.worse <- NULL
  F1.FFGL2.pop.opt <- 0
  F1.FFGL2.pop.worse <- 1

  G.FFGL2.group.opt <- NULL
  G.FFGL2.group.worse <- NULL
  F1.FFGL2.group.opt <- 0
  F1.FFGL2.group.worse <- 1


  for(lam in 1: length(AdjMats_lambda_FuDGE)){

    ################### FuDGE #########################
    G.FuDGE <- AdjMats_lambda_FuDGE[[lam]]
    G.FuDGE.diff <- G.FuDGE$Diff
    res_FuDGE <- prec.rec(G.diff, G.FuDGE.diff, type="AND")
    results_metrics <- rbind(results_metrics,c("DIFF","FuDGE", lam, res_FuDGE$prec,
                                               res_FuDGE$TPR, res_FuDGE$FPR, res_FuDGE$F1,
                                               n_g1, n_g2))
    if(res_FuDGE$F1 > F1.FuDGE.opt){
      F1.FuDGE.opt = res_FuDGE$F1
      G.FuDGE.diff.opt <- G.FuDGE.diff
    }
    
    if(res_FuDGE$F1 < F1.FuDGE.worse){
      F1.FuDGE.worse = res_FuDGE$F1
      G.FuDGE.diff.worse <- G.FuDGE.diff
    }
    
    ################### FGL #########################
    G.FGL <- AdjMats_lambda_FGL[[lam]]
    
    ## DIFERENTIAL
    G.FGL.diff <- G.FGL$Diff
    res_FGL <- prec.rec(G.diff, G.FGL.diff, type="AND")
    results_metrics <- rbind(results_metrics,c("DIFF","FGL", lam, res_FGL$prec,
                                               res_FGL$TPR, res_FGL$FPR, res_FGL$F1,
                                               n_g1, n_g2))
    
    if(res_FGL$F1 > F1.FGL.diff.opt){
      F1.FGL.diff.opt = res_FGL$F1
      G.FGL.diff.opt <- G.FGL.diff
    }
    
    if(res_FGL$F1 < F1.FGL.diff.worse){
      F1.FGL.diff.worse = res_FGL$F1
      G.FGL.diff.worse <- G.FGL.diff
    }
    
    ## POPULATION
    G.FGL.pop <- G.FGL$Pop
    res_FGL <- prec.rec(G.pop, G.FGL.pop, type="AND")
    results_metrics <- rbind(results_metrics,c("POP","FGL", lam, res_FGL$prec,
                                               res_FGL$TPR, res_FGL$FPR, res_FGL$F1,
                                               n_g1, n_g2))
    
    if(res_FGL$F1 > F1.FGL.pop.opt){
      F1.FGL.pop.opt = res_FGL$F1
      G.FGL.pop.opt <- G.FGL.pop
    }
    
    if(res_FGL$F1 < F1.FGL.pop.worse){
      F1.FGL.pop.worse = res_FGL$F1
      G.FGL.pop.worse <- G.FGL.pop
    }
    
    ## GROUP
    G.FGL.group <- G.FGL$Group
    res_FGL <- prec.rec(G.group, G.FGL.group, type="AND")
    results_metrics <- rbind(results_metrics,c("Group","FGL", lam, res_FGL$prec,
                                               res_FGL$TPR, res_FGL$FPR, res_FGL$F1,
                                               n_g1, n_g2))
    
    if(res_FGL$F1 > F1.FGL.group.opt){
      F1.FGL.group.opt = res_FGL$F1
      G.FGL.group.opt <- G.FGL.group
    }
    
    if(res_FGL$F1 < F1.FGL.group.worse){
      F1.FGL.group.worse = res_FGL$F1
      G.FGL.group.worse <- G.FGL.group
    }

    ################### FFGL #########################
    G.FFGL <- AdjMats_lambda_FFGL[[lam]]
    
    ## DIFERENTIAL
    G.FFGL.diff <- G.FFGL$Diff
    res_FFGL <- prec.rec(G.diff, G.FFGL.diff, type="AND")
    results_metrics <- rbind(results_metrics,c("DIFF","FFGL", lam, res_FFGL$prec,
                                               res_FFGL$TPR, res_FFGL$FPR, res_FFGL$F1,
                                               n_g1, n_g2))
    
    if(res_FFGL$F1 > F1.FFGL.diff.opt){
      F1.FFGL.diff.opt = res_FFGL$F1
      G.FFGL.diff.opt <- G.FFGL.diff
    }
    
    if(res_FFGL$F1 < F1.FFGL.diff.worse){
      F1.FFGL.diff.worse = res_FFGL$F1
      G.FFGL.diff.worse <- G.FFGL.diff
    }
    
    ## POPULATION
    G.FFGL.pop <- G.FFGL$Pop
    res_FFGL <- prec.rec(G.pop, G.FFGL.pop, type="AND")
    results_metrics <- rbind(results_metrics,c("POP","FFGL", lam, res_FFGL$prec,
                                               res_FFGL$TPR, res_FFGL$FPR, res_FFGL$F1,
                                               n_g1, n_g2))
    
    if(res_FFGL$F1 > F1.FFGL.pop.opt){
      F1.FFGL.pop.opt = res_FFGL$F1
      G.FFGL.pop.opt <- G.FFGL.pop
    }
    
    if(res_FFGL$F1 < F1.FFGL.pop.worse){
      F1.FFGL.pop.worse = res_FFGL$F1
      G.FFGL.pop.worse <- G.FFGL.pop
    }
    
    ## GROUP
    G.FFGL.group <- G.FFGL$Group
    res_FFGL <- prec.rec(G.group, G.FFGL.group, type="AND")
    results_metrics <- rbind(results_metrics,c("Group","FFGL", lam, res_FFGL$prec,
                                               res_FFGL$TPR, res_FFGL$FPR, res_FFGL$F1,
                                               n_g1, n_g2))
    
    if(res_FFGL$F1 > F1.FFGL.group.opt){
      F1.FFGL.group.opt = res_FFGL$F1
      G.FFGL.group.opt <- G.FFGL.group
    }
    
    if(res_FFGL$F1 < F1.FFGL.group.worse){
      F1.FFGL.group.worse = res_FFGL$F1
      G.FFGL.group.worse <- G.FFGL.group
    }
    
    ################### FFGL2 #########################
    G.FFGL2 <- AdjMats_lambda_FFGL2[[lam]]
    
    ## DIFERENTIAL
    G.FFGL2.diff <- G.FFGL2$Diff
    res_FFGL2 <- prec.rec(G.diff, G.FFGL2.diff, type="AND")
    results_metrics <- rbind(results_metrics,c("DIFF","FFGL2", lam, res_FFGL2$prec,
                                               res_FFGL2$TPR, res_FFGL2$FPR, res_FFGL2$F1,
                                               n_g1, n_g2))
    
    if(res_FFGL2$F1 > F1.FFGL2.diff.opt){
      F1.FFGL2.diff.opt = res_FFGL2$F1
      G.FFGL2.diff.opt <- G.FFGL2.diff
    }
    
    if(res_FFGL2$F1 < F1.FFGL2.diff.worse){
      F1.FFGL2.diff.worse = res_FFGL2$F1
      G.FFGL2.diff.worse <- G.FFGL2.diff
    }
    
    ## POPULATION
    G.FFGL2.pop <- G.FFGL2$Pop
    res_FFGL2 <- prec.rec(G.pop, G.FFGL2.pop, type="AND")
    results_metrics <- rbind(results_metrics,c("POP","FFGL2", lam, res_FFGL2$prec,
                                               res_FFGL2$TPR, res_FFGL2$FPR, res_FFGL2$F1,
                                               n_g1, n_g2))
    
    if(res_FFGL2$F1 > F1.FFGL2.pop.opt){
      F1.FFGL2.pop.opt = res_FFGL2$F1
      G.FFGL2.pop.opt <- G.FFGL2.pop
    }
    
    if(res_FFGL2$F1 < F1.FFGL2.pop.worse){
      F1.FFGL2.pop.worse = res_FFGL2$F1
      G.FFGL2.pop.worse <- G.FFGL2.pop
    }
    
    ## GROUP
    G.FFGL2.group <- G.FFGL2$Group
    res_FFGL2 <- prec.rec(G.group, G.FFGL2.group, type="AND")
    results_metrics <- rbind(results_metrics,c("Group","FFGL2", lam, res_FFGL2$prec,
                                               res_FFGL2$TPR, res_FFGL2$FPR, res_FFGL2$F1,
                                               n_g1, n_g2))
    
    if(res_FFGL2$F1 > F1.FFGL2.group.opt){
      F1.FFGL2.group.opt = res_FFGL2$F1
      G.FFGL2.group.opt <- G.FFGL2.group
    }
    
    if(res_FFGL2$F1 < F1.FFGL2.group.worse){
      F1.FFGL2.group.worse = res_FFGL2$F1
      G.FFGL2.group.worse <- G.FFGL2.group
    }

    
  }
  
  ########## Plot the best and worst case scenario ##########

  
  # FuDGE
  plot_binary_adj_matrix(G.FuDGE.diff.opt,
                         title= "Estimated adjacency Matrix - Differential - FuDGE (best)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_diff_FuDGE_best.png")
  
  plot_binary_adj_matrix(G.FuDGE.diff.worse,
                         title= "Estimated adjacency Matrix - Differential - FuDGE (worse)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_diff_FuDGE_worse.png")
  
  
  # FGL
  plot_binary_adj_matrix(G.FGL.diff.opt,
                         title= "Estimated adjacency Matrix - Differential - FGL (best)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_diff_FGL_best.png")
  
  plot_binary_adj_matrix(G.FGL.diff.worse,
                         title= "Estimated adjacency Matrix - Differential - FGL (worse)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_diff_FGL_worse.png")
  
  plot_binary_adj_matrix(G.FGL.pop.opt,
                         title= "Estimated adjacency Matrix - Population - FGL (best)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_pop_FGL_best.png")
  
  plot_binary_adj_matrix(G.FGL.pop.worse,
                         title= "Estimated adjacency Matrix - Population - FGL (worse)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_pop_FGL_worse.png")
  
  plot_binary_adj_matrix(G.FGL.group.opt,
                         title= "Estimated adjacency Matrix - Group - FGL (best)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_group_FGL_best.png")
  
  plot_binary_adj_matrix(G.FGL.group.worse,
                         title= "Estimated adjacency Matrix - Group - FGL (worse)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_group_FGL_worse.png")
  
  #  FFGL
  plot_binary_adj_matrix(G.FFGL.diff.opt,
                         title= "Estimated adjacency Matrix - Differential - FFGL (best)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_diff_FFGL_best.png")
  
  plot_binary_adj_matrix(G.FFGL.diff.worse,
                         title= "Estimated adjacency Matrix - Differential - FFGL (worse)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_diff_FFGL_worse.png")
  
  plot_binary_adj_matrix(G.FFGL.pop.opt,
                         title= "Estimated adjacency Matrix - Population - FFGL (best)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_pop_FFGL_best.png")
  
  plot_binary_adj_matrix(G.FFGL.pop.worse,
                         title= "Estimated adjacency Matrix - Population - FFGL (worse)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_pop_FFGL_worse.png")
  
  plot_binary_adj_matrix(G.FFGL.group.opt,
                         title= "Estimated adjacency Matrix - Group - FFGL (best)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_group_FFGL_best.png")
  
  plot_binary_adj_matrix(G.FFGL.group.worse,
                         title= "Estimated adjacency Matrix - Group - FFGL (worse)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_group_FFGL_worse.png")
  
  # FFGL2
  plot_binary_adj_matrix(G.FFGL2.diff.opt,
                         title= "Estimated adjacency Matrix - Differential - FFGL2 (best)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_diff_FFGL2_best.png")
  
  plot_binary_adj_matrix(G.FFGL2.diff.worse,
                         title= "Estimated adjacency Matrix - Differential - FFGL2 (worse)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_diff_FFGL2_worse.png")
  
  plot_binary_adj_matrix(G.FFGL2.pop.opt,
                         title= "Estimated adjacency Matrix - Population - FFGL2 (best)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_pop_FFGL2_best.png")
  
  plot_binary_adj_matrix(G.FFGL2.pop.worse,
                         title= "Estimated adjacency Matrix - Population - FFGL2 (worse)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_pop_FFGL2_worse.png")
  
  plot_binary_adj_matrix(G.FFGL2.group.opt,
                         title= "Estimated adjacency Matrix - Group - FFGL2 (best)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_group_FFGL2_best.png")
  
  plot_binary_adj_matrix(G.FFGL2.group.worse,
                         title= "Estimated adjacency Matrix - Group - FFGL2 (worse)",
                         output_path=output_path,
                         plot_name="Estimated_adjacency_matrix_node_group_FFGL2_worse.png")
  
}
results_metrics <- results_metrics[-1,]
out_csv <- paste0(save_path, simulation_name,
                  "/Litt_comp_results__performance_metrics_", name_output, ".csv")
write.csv(results_metrics, out_csv)
cat("Results saved to ", out_csv, "\n")

