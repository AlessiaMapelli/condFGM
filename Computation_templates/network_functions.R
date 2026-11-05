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
# @return list(baseline_selected, group_norm, diff_norm, vr, nib)
#   baseline_selected: ||B0_k||_F > threshold, i.e. k is in the population
#                      neighbourhood of j at the selected (lambda, threshold)
#   group_norm:        ||B0_k + Bc_k||_F
#   diff_norm:         ||Bc_k||_F
#   vr:                group_norm / ||B0_k||_F if the baseline is selected,
#                      NA otherwise (Eq. 5.16, first branch)
#   nib:               diff_norm if the baseline is NOT selected ("not in
#                      baseline", Eq. 5.16, second branch), NA otherwise
directional_vr <- function(P.values, P.frob, k, num_nodes, cov_index, threshold) {
  b0_norm    <- P.frob[[k]]
  bc         <- P.values[[k + num_nodes * (cov_index - 1)]]
  group_norm <- norm(P.values[[k]] + bc, "F")
  diff_norm  <- norm(bc, "F")
  baseline_selected <- !is.null(b0_norm) && b0_norm > threshold
  list(baseline_selected = baseline_selected,
       group_norm = group_norm,
       diff_norm  = diff_norm,
       vr  = if (baseline_selected) group_norm / b0_norm else NA_real_,
       nib = if (baseline_selected) NA_real_ else diff_norm)
}

# Symmetric edge weights (Eq. 5.17) from the directional matrices.
# For each edge of the symmetrised differential graph:
#   - VR weight: geometric mean of the directional VRs over the directions
#     where both the differential and the baseline block are selected;
#   - otherwise the edge is "not in baseline" and its magnitude is the
#     geometric mean of ||Bc||_F over the directions where the differential
#     block is selected.
# @param vr, nib  p x p directional matrices from directional_vr()
# @param G_symm   p x p symmetrised differential adjacency matrix
# @return list(vr, nib, warning): symmetric matrices (NA where not defined);
#   warning = 1 when two qualifying directions disagree on the sign of change
symmetric_weights <- function(vr, nib, G_symm) {
  p <- nrow(G_symm)
  W   <- matrix(NA_real_, p, p)
  NIB <- matrix(NA_real_, p, p)
  WRN <- matrix(NA_real_, p, p)
  for (j in 1:(p - 1)) {
    for (k in (j + 1):p) {
      if (G_symm[j, k]) {
        x <- c(vr[j, k], vr[k, j]); x <- x[!is.na(x)]
        if (length(x) > 0) {
          W[j, k] <- W[k, j] <- exp(mean(log(x)))
          disagree <- length(x) == 2 && !((x[1] > 1 & x[2] > 1) | (x[1] < 1 & x[2] < 1))
          WRN[j, k] <- WRN[k, j] <- as.numeric(disagree)
        } else {
          y <- c(nib[j, k], nib[k, j]); y <- y[!is.na(y)]
          NIB[j, k] <- NIB[k, j] <- exp(mean(log(y)))
          WRN[j, k] <- WRN[k, j] <- 0
        }
      }
    }
  }
  list(vr = W, nib = NIB, warning = WRN)
}

# Variation-rate heat map with two colour scales:
#   - edges in the baseline: VR, diverging blue - white (VR = 1, no change) - red;
#   - edges not in the baseline: ||Bc||_F, sequential white (0) - dark green.
# The second scale requires the 'ggnewscale' package; without it the
# not-in-baseline edges are drawn in a single green colour.
# @param vr  p x p matrix of VR values (NA where not defined)
# @param nib p x p matrix of ||Bc||_F for not-in-baseline edges (NA elsewhere)
plot_vr_matrix <- function(vr, nib, title, filename,
                           low = "darkblue", high = "darkred",
                           nib_high = "darkgreen") {
  adj_df <- melt(vr, varnames = c("Row", "Col"), value.name = "Value")
  nib_df <- melt(nib, varnames = c("Row", "Col"), value.name = "NIB")
  nib_df <- nib_df[!is.na(nib_df$NIB), , drop = FALSE]

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
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.35, "cm"),
      plot.margin = margin(4, 4, 4, 4)
    ) +
    scale_fill_gradient2(low = low, mid = "white", high = high, midpoint = 1,
                         limits = range(c(vr, 1), na.rm = TRUE),
                         na.value = "white", name = "VR w.r.t.\nbaseline",
                         guide = if (all(is.na(vr))) "none" else
                           guide_colourbar(order = 1))
  if (nrow(nib_df) > 0) {
    if (requireNamespace("ggnewscale", quietly = TRUE)) {
      plot <- plot +
        ggnewscale::new_scale_fill() +
        geom_tile(data = nib_df, aes(x = Col, y = Row, fill = NIB),
                  color = "grey92", linewidth = 0.15) +
        scale_fill_gradient(low = "white", high = nib_high,
                            limits = c(0, max(nib_df$NIB)),
                            name = "VR not in\nbaseline",
                            guide = guide_colourbar(order = 2))
    } else {
      warning("Package 'ggnewscale' not available: not-in-baseline edges ",
              "are drawn in a single colour without magnitude.")
      plot <- plot +
        geom_tile(data = nib_df, aes(x = Col, y = Row, colour = "VR not in baseline"),
                  fill = nib_high, linewidth = 0.15) +
        scale_colour_manual(values = c("VR not in baseline" = "grey92"), name = NULL,
                            guide = guide_legend(override.aes = list(fill = nib_high),
                                                 order = 2))
    }
  }
  ggsave(plot, filename = filename, width = 8, height = 8, dpi = 300)
  invisible(plot)
}
