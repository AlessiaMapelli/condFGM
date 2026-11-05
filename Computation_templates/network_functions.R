# Helper functions for assembling networks from the node-wise estimates.
# Shared by Computation_templates/Results_evaluation.R and the simulation
# evaluation scripts, so that the template and the simulations compute the
# symmetrised graphs, the group-specific networks and the variation rate (VR)
# in exactly the same way.

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

# Symmetrise a directional adjacency matrix (Eqs. 5.14-5.15).
# @param G    p x p 0/1 (or logical) matrix; G[j,k] = 1 if k is in the estimated
#             neighbourhood of j
# @param type "OR", "AND" or anything else (returned unchanged)
# @return p x p symmetric 0/1 matrix
symmetrise_adjacency <- function(G, type = "OR") {
  G <- (G == 1) & !is.na(G)
  if (type == "AND") {
    S <- G & t(G)
  } else if (type == "OR") {
    S <- G | t(G)
  } else {
    return(G * 1)
  }
  diag(S) <- FALSE
  S * 1
}

# Directional quantities for the covariate block in position `cov_index`
# (1 = population, c + 1 = c-th covariate column of the design matrix) of the
# regression of node j on node k. P.values / P.frob / threshold are the objects
# saved in "<name>_<j>coeff.rda" by run_node_computation().
# @return list(baseline_selected, group_norm, vr)
#   baseline_selected: ||B0_k||_F > threshold, i.e. k is in the population
#                      neighbourhood of j at the selected (lambda, threshold)
#   group_norm:        ||B0_k + Bc_k||_F
#   vr:                group_norm / ||B0_k||_F when the baseline is selected,
#                      NA otherwise ("not in baseline", Eq. 5.16)
directional_vr <- function(P.values, P.frob, k, num_nodes, cov_index, threshold) {
  b0_norm    <- P.frob[[k]]
  group_norm <- norm(P.values[[k]] + P.values[[k + num_nodes * (cov_index - 1)]], "F")
  baseline_selected <- !is.null(b0_norm) && b0_norm > threshold
  vr <- if (baseline_selected) group_norm / b0_norm else NA_real_
  list(baseline_selected = baseline_selected, group_norm = group_norm, vr = vr)
}

# Variation-rate heat map. Edges whose baseline coefficient is not selected
# are drawn in their own colour and labelled "Not in baseline".
# @param vr              p x p matrix of VR values (NA where not defined)
# @param not_in_baseline p x p logical matrix
plot_vr_matrix <- function(vr, not_in_baseline, title, filename,
                           low = "darkblue", high = "darkred",
                           nib_colour = "#E69F00") {
  adj_df <- melt(vr, varnames = c("Row", "Col"), value.name = "Value")
  nib_df <- melt(not_in_baseline, varnames = c("Row", "Col"), value.name = "NIB")
  nib_df <- nib_df[!is.na(nib_df$NIB) & nib_df$NIB, , drop = FALSE]

  plot <- ggplot(adj_df, aes(x = Col, y = Row)) +
    geom_tile(aes(fill = Value), color = "grey92", linewidth = 0.15) +
    coord_fixed() +
    scale_y_reverse() +
    labs(title = str_wrap(title, width = 42), x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 20, hjust = 0.5, vjust = -1.5,
                                lineheight = 0.95),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.35, "cm"),
      plot.margin = margin(4, 4, 4, 4)
    ) +
    scale_fill_gradient2(low = low, mid = "white", high = high, midpoint = 1,
                         na.value = "white", name = "VR")
  if (nrow(nib_df) > 0) {
    plot <- plot +
      geom_tile(data = nib_df, aes(x = Col, y = Row, colour = "Not in baseline"),
                fill = nib_colour, linewidth = 0.15) +
      scale_colour_manual(
        values = c("Not in baseline" = "grey92"), name = NULL,
        guide = guide_legend(override.aes = list(fill = nib_colour),
                             order = 2)
      ) +
      guides(fill = guide_colourbar(order = 1))
  }
  ggsave(plot, filename = filename, width = 8, height = 8, dpi = 300)
  invisible(plot)
}
