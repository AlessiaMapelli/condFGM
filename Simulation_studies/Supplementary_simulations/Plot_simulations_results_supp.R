# Produces Figures 9, 10, 12 of the manuscript (Appendix B.1.2, B.1.3 supplementary results).
# Run from the repository root:
#   Rscript Simulation_studies/Supplementary_simulations/Plot_simulations_results_supp.R
# Figure 10 requires seed_3 results for p10_n1000_n1000
# (from Results_evaluation_supp_simulation.R or downloaded from https://osf.io/ta3bq/overview).

# Required packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(png))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(patchwork))

base.dir <- "Simulation_studies/Supplementary_simulations/"
figures.dir <- file.path(base.dir, "figures")
if (!dir.exists(figures.dir)) {
  dir.create(figures.dir, recursive = TRUE)
}

#################### 
# SUPPLEMENTARY SETTING 
####################

sim_supp_path <- paste0(base.dir, "simulation_settings")
results_dirs_supp <- list.dirs(sim_supp_path, recursive = FALSE)
results_list_supp <- list()

for (dir_path in results_dirs_supp) {
  file_list <- list.files(dir_path,
                          pattern = "^Results_performance_metrics.*\\.csv$",
                          full.names = TRUE)
  file_path <- file_list[1]
  if (file.exists(file_path) ) {
      df <- read_csv(file_path)[,-1]
      df <- df %>%
        group_by(network, symm) %>%
        summarise(
          med_prec = mean(prec, na.rm = TRUE),
          sd_prec = sd(prec, na.rm = TRUE),
          med_TPR = mean(TPR, na.rm = TRUE),
          sd_TPR = sd(TPR, na.rm = TRUE),
          med_FPR = mean(FPR, na.rm = TRUE),
          sd_FPR = sd(FPR, na.rm = TRUE),
          max_FPR = max(FPR, na.rm = TRUE),
          min_FPR = min(FPR, na.rm = TRUE),
          med_F1 = mean(F1, na.rm = TRUE),
          sd_F1 = sd(F1, na.rm = TRUE),
          max_F1 = max(F1, na.rm = TRUE),
          min_F1 = min(F1, na.rm = TRUE),
          Baseline = first(Baseline),
          Differential = first(Differential),
          .groups = "drop"
        )
      df$simulation <- basename(dir_path)
      df$p <- as.numeric(str_extract(df$simulation, "(?<=p)\\d+(?=_)"))
      results_list_supp[[length(results_list_supp) + 1]] <- df
  }else{
    cat("Results not found for simulation setting: ", dir_path, " \n")
  }
}

combined_data <- bind_rows(results_list_supp)

combined_data <- combined_data %>%
  mutate(p = as.numeric(p))

# Figure 9
temp_data <- combined_data %>%
  filter(network != "GROUP")

network_levels <- c("POP", "DIFF", "CONT")
temp_data$network <- factor(temp_data$network, levels = network_levels)

symm_levels <- c("OR", "AND")
temp_data$symm <- factor(temp_data$symm, levels = symm_levels)

network_cols <- c(
  "POP" = "#38A3A5",
  "DIFF" = "#0072B2",
  "CONT" = "#5C5CFF"
)

method_lt <- c(
  "OR" = "solid",
  "AND" = "dashed"
)

plot_obj <- ggplot(
  temp_data,
  aes(
    x = Baseline, y = med_F1,
    color = network,
    linetype = symm,
    group = interaction(network, symm)
  )
) +
  geom_line(linewidth = 1.10) +
  geom_point(size = 2.3) +
  geom_errorbar(
    aes(ymin = min_F1, ymax = max_F1),
    width = 2,
    linewidth = 0.3
  ) +
  facet_wrap(~ p, scales = "free_x", nrow = 1, labeller = label_both) +
  scale_color_manual(values = network_cols, drop = FALSE, labels = c(
    "POP"  = "Population",
    "DIFF" =  "Differential - categorical",
    "CONT" = "Differential - continuous")) +
  scale_linetype_manual(values = method_lt, drop = FALSE) +
  labs(
    x = "Sample size for group",
    y = "F1 score",
    color = "Network",
    linetype = "Symmetrization"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(
    breaks = c(50, 100, 200, 400, 750, 1000)
  ) + 
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )

ggsave(file.path(figures.dir, "Figure_9_F1_performance_supplementary.png"),
       plot_obj, width=15, height=6, dpi=300)

# Figure 12
plot_obj <- ggplot(
  temp_data,
  aes(
    x = Baseline, y = med_FPR,
    color = network,
    linetype = symm,
    group = interaction(network, symm)
  )
) +
  geom_line(linewidth = 1.10) +
  geom_point(size = 2.3) +
  geom_errorbar(
    aes(ymin = min_FPR, ymax = max_FPR),
    width = 2,
    linewidth = 0.3
  ) +
  facet_wrap(~ p, scales = "free_x", nrow = 1, labeller = label_both) +
  scale_color_manual(values = network_cols, drop = FALSE, labels = c(
    "POP"  = "Population",
    "DIFF" =  "Differential - categorical",
    "CONT" = "Differential - continuous")) +
  scale_linetype_manual(values = method_lt, drop = FALSE) +
  labs(
    x = "Sample size for group",
    y = "F1 score",
    color = "Network",
    linetype = "Symmetrization"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(
    breaks = c(50, 100, 200, 400, 750, 1000)
  ) + 
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )

ggsave(file.path(figures.dir, "Figure_12_FPR_performance_supplementary.png"),
       plot_obj, width=15, height=6, dpi=300)

# Figure 10 — Estimated adjacency structure for one simulation replicate (p=10, n=1000, seed 3)
seed3_dir <- file.path(
  base.dir,
  "simulation_settings", "p10_n1000_n1000", "seed_3", "results"
)
fig10_paths <- c(
  file.path(seed3_dir, "Estimated_adjacency_matrix_node_pop.png"),
  file.path(seed3_dir,
            "Estimated_adjacency_matrix_node_diff_categorical_signed_weights.png"),
  file.path(seed3_dir,
            "Estimated_adjacency_matrix_node_diff_continuous_signed_weights.png")
)

img_to_gg <- function(path) {
  img <- rasterGrob(readPNG(path), interpolate = TRUE)
  ggplot() + annotation_custom(img) + theme_void()
}
p_pop  <- img_to_gg(fig10_paths[1])
p_diff <- img_to_gg(fig10_paths[2])
p_cont <- img_to_gg(fig10_paths[3])
fig10 <- p_pop + p_diff + p_cont
ggsave(
  file.path(figures.dir, "Figure_10_adjacency_replicate.png"),
  fig10, width = 21, height = 7, dpi = 300
)

