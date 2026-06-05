# Produces Figures 3, 4, 5, 6, 11 of the manuscript
# (simulation performance and literature comparison).
# Run from the repository root: Rscript Simulation_studies/Plot_simulations_results.R

# Required packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(latex2exp))

base.dir <- "Simulation_studies/"
figures.dir <- file.path(base.dir, "figures")
if (!dir.exists(figures.dir)) {
  dir.create(figures.dir, recursive = TRUE)
}

#################### 
# SETTING 1
####################

sim_step1_path <- paste0(base.dir, "Step1/simulation_settings")
results_dirs_step1 <- list.dirs(sim_step1_path, recursive = FALSE)
results_list_step1 <- list()

for (dir_path in results_dirs_step1) {
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
          max_prec = max(prec, na.rm = TRUE),
          min_prec = min(prec, na.rm = TRUE),
          med_TPR = mean(TPR, na.rm = TRUE),
          sd_TPR = sd(TPR, na.rm = TRUE),
          max_TPR = max(TPR, na.rm = TRUE),
          min_TPR = min(TPR, na.rm = TRUE),
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
      results_list_step1[[length(results_list_step1) + 1]] <- df
  }else{
    cat("Results not found for simulation setting: ", dir_path, " \n")
  }
}

combined_data <- bind_rows(results_list_step1)
combined_data <- combined_data %>%
  mutate(p = as.numeric(p))

# Figure 3
temp_data <- combined_data %>%
  filter(network != "GROUP")

network_levels <- c("POP", "DIFF" )
temp_data$network <- factor(temp_data$network, levels = network_levels)

symm_levels <- c("OR", "AND")
temp_data$symm <- factor(temp_data$symm, levels = symm_levels)

network_cols <- c(
  "POP" = "#38A3A5",
  "DIFF" = "#0072B2"
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
    "DIFF" = "Differential")) +
  scale_linetype_manual(values = method_lt, drop = FALSE) +
  labs(
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

ggsave(file.path(figures.dir, "Figure_3_F1_performance.png"),
       plot_obj, width = 10, height = 5, dpi = 300)

#Figure 11
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
    "DIFF" = "Differential")) +
  scale_linetype_manual(values = method_lt, drop = FALSE) +
  labs(
    x = "Sample size for group",
    y = "False Positive Rate",
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
ggsave(file.path(figures.dir, "Figure_11_FPR_performance.png"),
       plot_obj, width = 10, height = 5, dpi = 300)


####################
# SETTING 1 - Litt Comparison
####################
litt_comp_results_list <- list()

for (dir_path in results_dirs_step1) {
  file_list <- list.files(dir_path,
                          pattern = "^Litt_comp_results__performance_metrics.*\\.csv$",
                          full.names = TRUE)
  file_path <- file_list[1]
  if (file.exists(file_path) ) {
    df <- read_csv(file_path, show_col_types = FALSE)
    df <- df %>%
      group_by(network, method, hyper) %>%
      summarise(
        med_prec = mean(prec, na.rm = TRUE),
        sd_prec = sd(prec, na.rm = TRUE),
        max_prec = max(prec, na.rm = TRUE),
        min_prec = min(prec, na.rm = TRUE),
        med_TPR = mean(TPR, na.rm = TRUE),
        sd_TPR = sd(TPR, na.rm = TRUE),
        max_TPR = max(TPR, na.rm = TRUE),
        min_TPR = min(TPR, na.rm = TRUE),
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
    df <- df%>%
      group_by(network, method) %>%
      slice_max(med_F1, n = 1, with_ties = FALSE) %>%
      dplyr::select(-hyper)
    litt_comp_results_list[[length(litt_comp_results_list) + 1]] <- df
  }else{
    cat("Results of litt comparison not found for simulation setting: ", dir_path, " \n")
  }
}

litt_comp_combined_data <- bind_rows(litt_comp_results_list)

combined_data <- combined_data %>%
  rename(method = symm)

litt_comp_combined_data <- rbind(combined_data,litt_comp_combined_data)

# Figure 4
litt_comp_combined_data_diff <- litt_comp_combined_data %>%
  filter(network=="DIFF" & method != "AND")

litt_comp_combined_data_diff[litt_comp_combined_data_diff$method == "OR", ]$method <- "Our Method"

method_levels <- c("Our Method", "FuDGE", "FGL", "FFGL", "FFGL2")
litt_comp_combined_data_diff$method <- factor(litt_comp_combined_data_diff$method,
                                              levels = method_levels)

method_cols <- c(
  "FFGL"            = "#E69F00",
  "FFGL2"           = "#F4D06F",
  "FGL"             = "#D55E00",
  "FuDGE"           = "#C10A0A",
  "Our Method" = "#0072B2"
)

method_lw <- c(
  "FFGL"            = 0.45,
  "FFGL2"           = 0.45,
  "FGL"             = 0.70,
  "FuDGE"           = 0.70,
  "Our Method" = 1.10
)

plot_obj <- ggplot(
  litt_comp_combined_data_diff,
  aes(x = Baseline, y = med_F1, group = method, color = method, linewidth = method)
) +
  geom_line() +
  geom_point(size = 2.3) +
  facet_wrap(~ p, scales = "free_x", nrow = 1, labeller = label_both) +
  scale_color_manual(values = method_cols, drop = FALSE) +
  scale_linewidth_manual(values = method_lw, guide = "none") +  
  labs(
    title = NULL,  
    x = "Sample size for group",
    y = "F1 score",
    color = "Method"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 14)  +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold")
  )

ggsave(file.path(figures.dir, "Figure_4_literature_comparison.png"),
       plot_obj, width = 10, height = 5, dpi = 300)


#################### 
# SETTING 1 - Comp time
####################

computational_times <- read_csv(file.path(sim_step1_path, "computational_times.csv"))
computational_times <- computational_times %>%
  mutate(comp_time_sec = comp_time) %>%
  dplyr::select(p, n_g1, method, iteration, comp_time_sec) %>%
  rename(n = n_g1) %>%
  filter(n == 200) %>%
  group_by(iteration, method) %>%
  mutate(
    baseline_time = comp_time_sec[p == 10][1],
    scaled_comp_time = comp_time_sec / baseline_time
  ) %>%
  ungroup() %>%
  group_by(p, n, method) %>%
  summarise(
    med_scaled_comp_time = mean(scaled_comp_time, na.rm = TRUE),
    sd_scaled_comp_time  = sd(scaled_comp_time, na.rm = TRUE),
    max_scaled_comp_time = max(scaled_comp_time, na.rm = TRUE),
    min_scaled_comp_time = min(scaled_comp_time, na.rm = TRUE),
    .groups = "drop"
  )

method_levels <- c("Our Method", "FuDGE")
computational_times$method <- factor(computational_times$method,
                                              levels = method_levels)

method_cols <- c(
  "FuDGE"           = "#C10A0A",
  "Our Method" = "#0072B2"
)

method_lw <- c(
  "FuDGE"           = 0.70,
  "Our Method" = 1.10
)


# Figure 5
plot_obj <- ggplot(
  computational_times,
  aes(x = p, y = med_scaled_comp_time, group = method, color = method, linewidth = method)
) +
  geom_line() +
  geom_point(size = 2.3) +
  geom_errorbar(
    aes(ymin = min_scaled_comp_time, ymax = max_scaled_comp_time),
    width = 2, linewidth = 0.3
  ) +
  scale_color_manual(values = method_cols, drop = FALSE) +
  scale_linewidth_manual(values = method_lw, guide = "none") +  #
  labs(
    title = NULL,  
    x = "Network size",
    y = "Computational time proportional to p=10",
    color = "Method"
  ) +
  theme_minimal(base_size = 14)  +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold")
  )

ggsave(file.path(figures.dir, "Figure_5_computational_time.png"),
       plot_obj, width = 10, height = 6, dpi = 300)



####################
# SETTING 2
####################

sim_step2_path <- paste0(base.dir, "Step2/simulation_settings")
sel_dirs <- list.dirs(sim_step2_path, recursive = FALSE)
results_list_step2 <- list()

for (dir_path in sel_dirs) {
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
        max_prec = max(prec, na.rm = TRUE),
        min_prec = min(prec, na.rm = TRUE),
        med_TPR = mean(TPR, na.rm = TRUE),
        sd_TPR = sd(TPR, na.rm = TRUE),
        max_TPR = max(TPR, na.rm = TRUE),
        min_TPR = min(TPR, na.rm = TRUE),
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
    df$simulation_type <- str_sub(df$simulation, -2)
    results_list_step2[[length(results_list_step2) + 1]] <- df
  }else{
    cat("Results not found for simulation setting: ", dir_path, " \n")
  }
}

combined_data <- bind_rows(results_list_step2)
combined_data <- combined_data %>%
  filter(symm == "OR")

p_levels <- c("10", "50")
combined_data$p <- factor(combined_data$p, levels = p_levels)

combined_data$network <- factor(combined_data$network)

p_cols <- c(
  "10" = "#89CFF1",
  "50" = "#0072B2"
)

network_labels <- c(
  DIFF  = "Differential",
  GROUP = "Group $\\textit{G}_1$",
  POP   = "Population"
)

pd <- position_dodge(width = 0.55)

# Figure 6
plot_obj <- ggplot(
  combined_data,
  aes(x = network, y = med_F1, color = p, group = p)
) +
  geom_pointrange(
    aes(ymin = min_F1, ymax = max_F1),
    position = pd,
    linewidth = 0.35
  ) +
  geom_point(position = pd, size = 2.6) +
  facet_wrap(
    ~ simulation_type,
    scales = "fixed"
  ) +
  scale_x_discrete(
    labels = function(x) latex2exp::TeX(network_labels[x])
  ) +
  scale_color_manual(
    values = p_cols,
    drop   = FALSE,
    labels = c(
      "10" = "p = 10",
      "50" = "p = 50"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title  = NULL,
    x      = "Network",
    y      = "F1 score",
    color  = "Nodes number"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position   = "right",
    legend.title      = element_text(face = "bold"),
    strip.background  = element_rect(fill = "grey95", color = NA),
    strip.text        = element_text(face = "bold"),
    axis.title        = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(figures.dir, "Figure_6_scenario_comparison.png"),
       plot_obj, width = 10, height = 5, dpi = 300)
