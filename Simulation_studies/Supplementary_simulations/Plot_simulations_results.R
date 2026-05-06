# Load required packages
install.packages(setdiff(c("tidyverse", "readr"), rownames(installed.packages())))
library(tidyverse)
library(readr)

# Results directory
base_path <- "/Simulation_studies/Step3/Simulations_results/"
setwd(base_path)

# Get list of folders in the results directory
results_dirs <- list.dirs(base_path, recursive = FALSE)

# Collect data from each folder
results_list <- list()

#################### 
# SETTING 3 
####################

for (dir_path in results_dirs) {
  file_list <- list.files(dir_path, pattern = "^test_results_metices.*\\.csv$", full.names = TRUE)
  file_path <- file_list[1]
  if (file.exists(file_path) ) {
    tryCatch({
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
      df$screening <- "Yes"
      results_list[[length(results_list) + 1]] <- df
    }, error = function(e) {
      message("Error reading: ", file_path)
    })
  }
}

# Combine all results
combined_data <- bind_rows(results_list)

combined_data <- combined_data %>%
  mutate(p = as.numeric(p))

# Figure 2
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
    "DIFF" =  "Differential - group modulated",
    "CONT" = "Differential - continuously modulated")) +
  scale_linetype_manual(values = method_lt, drop = FALSE) +
  labs(
    title = "Simulation Performance Metrics - F1 score",
    x = "Sample size for group",
    y = "F1 score",
    color = "Network",
    linetype = "Symmetrization"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
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

plot_obj
ggsave(plot_obj, filename=paste(base_path, "/F1_behaviour_varying_p_and_n.png", sep=""), width=8, height=8, dpi=300)

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
    aes(ymin = med_FPR - sd_FPR, ymax = med_FPR + sd_FPR),
    width = 2,
    linewidth = 0.3
  ) +
  facet_wrap(~ p, scales = "free_x", nrow = 1, labeller = label_both) +
  scale_color_manual(values = network_cols, drop = FALSE, labels = c(
    "POP"  = "Population",
    "DIFF" = "Differential - group modulated",
    "CONT" = "Differential - continuously modulated")) +
  scale_linetype_manual(values = method_lt, drop = FALSE) +
  labs(
    title = "Simulation Performance Metrics - FPR",
    x = "Sample size for group",
    y = "FPR",
    color = "Network",
    linetype = "Symmetrization"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
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

plot_obj
ggsave(plot_obj, filename=paste(base_path, "/FPR_behaviour_varying_p_and_n.png", sep=""), width=8, height=8, dpi=300)

plot_obj <- ggplot(
  temp_data,
  aes(
    x = Baseline, y = med_TPR,
    color = network,
    linetype = symm,
    group = interaction(network, symm)
  )
) +
  geom_line(linewidth = 1.10) +
  geom_point(size = 2.3) +
  geom_errorbar(
    aes(ymin = med_TPR - sd_TPR, ymax = med_TPR + sd_TPR),
    width = 2,
    linewidth = 0.3
  ) +
  facet_wrap(~ p, scales = "free_x", nrow = 1, labeller = label_both) +
  scale_color_manual(values = network_cols, drop = FALSE, labels = c(
    "POP"  = "Population",
    "DIFF" = "Differential - group modulated",
    "CONT" = "Differential - continuously modulated")) +
  scale_linetype_manual(values = method_lt, drop = FALSE) +
  labs(
    title = "Simulation Performance Metrics - TPR",
    x = "Sample size for group",
    y = "TPR",
    color = "Network",
    linetype = "Symmetrization"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
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

plot_obj
ggsave(plot_obj, filename=paste(base_path, "/TPR_behaviour_varying_p_and_n.png", sep=""), width=8, height=8, dpi=300)
