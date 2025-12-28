# Load required packages
install.packages(setdiff(c("tidyverse", "readr"), rownames(installed.packages())))
library(tidyverse)
library(readr)

base_path <- "Simulation_studies/Step1/"
setwd(base_path)

# Get list of folders (hotspots)
results_dirs <- list.dirs(base_path, recursive = FALSE)

# Collect data from each folder
results_list <- list()

#################### 
# SETTING 1
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


#################### 
# SETTING 1 - Litt Comparison
####################

litt_comp_results_list <- list()
for (dir_path in results_dirs) {
  file_list <- list.files(dir_path, pattern = "^test_results_metices.*\\.csv$", full.names = TRUE)
  file_path <- file_list[1]
  df <- read_csv(file_path, show_col_types = FALSE)
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
  df <- df %>% rename(method=symm)
  litt_comp_results_list[[length(litt_comp_results_list) + 1]] <- df
  
  file_list_litt <- list.files(dir_path, pattern = "^litt_comp_test_results_metices.*\\.csv$", full.names = TRUE)
  file_path_litt <- file_list_litt[1]
  df_litt <- read_csv(file_path_litt, show_col_types = FALSE)
  
  df_litt <- df_litt%>%
    group_by(network, method, hyper) %>%
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
  df_litt$simulation <- basename(dir_path)
  df_litt$p <- as.numeric(str_extract(df_litt$simulation, "(?<=p)\\d+"))
  
  df_litt <- df_litt%>%
    group_by(network, method) %>%
    slice_max(med_F1, n = 1, with_ties = FALSE) %>%
    dplyr::select(-hyper)
  litt_comp_results_list[[length(litt_comp_results_list) + 1]] <- df_litt
  
}


litt_comp_combined_data <- bind_rows(litt_comp_results_list)

litt_comp_combined_data <- litt_comp_combined_data %>%
  mutate(p = as.numeric(p))

# Figure 3

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
  # If you want error bars back:
  # geom_errorbar(aes(ymin = min_F1, ymax = max_F1), width = 2, linewidth = 0.3) +
  facet_wrap(~ p, scales = "free_x", nrow = 1, labeller = label_both) +
  scale_color_manual(values = method_cols, drop = FALSE) +
  scale_linewidth_manual(values = method_lw, guide = "none") +  # hide linewidth legend (cleaner)
  labs(
    title = NULL,  # journals often prefer title in caption, not inside the plot
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

plot_obj

# Supp Figure 1
sel_dirs <- c("Simulation_studies/Step1/p25_n100_n100") 
dir_path <- sel_dirs[1]
results_list <- list()
file_list <- list.files(dir_path, pattern = "^test_results_metices.*\\.csv$", full.names = TRUE)
file_path <- file_list[1]
df <- read_csv(file_path, show_col_types = FALSE)
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

file_list_litt <- list.files(dir_path, pattern = "^litt_comp_test_results_metices.*\\.csv$", full.names = TRUE)
file_path_litt <- file_list_litt[1]
df_litt <- read_csv(file_path_litt, show_col_types = FALSE)
View(df_litt)
try <- df_litt%>%
  group_by(network, method, hyper) %>%
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
try$simulation <- basename(dir_path)
try$p <- as.numeric(str_extract(try$simulation, "(?<=p)\\d+"))

df_expanded <- df %>%
  rename(method = symm) %>%        
  mutate(row_id = row_number()) %>%  
  uncount(weights = 41) %>%          
  group_by(row_id) %>%               
  mutate(hyper = row_number()) %>%   
  ungroup() %>%
  dplyr::select(-row_id) %>%
  relocate(hyper, .after = 2)

df_expanded

try3 <- rbind(try, df_expanded)

try3 <- try3[try3$method != "AND", ]
try3[try3$method == "OR", ]$method <- "Proposed Method"


method_levels <- c("Proposed Method", "FuDGE", "FGL", "FFGL", "FFGL2")
try3$method <- factor(try3$method, levels = method_levels)

method_cols <- c(
  "FFGL"            = "#E69F00",
  "FFGL2"           = "#F4D06F",
  "FGL"             = "#D55E00",
  "FuDGE"           = "#C10A0A",
  "Proposed Method" = "#0072B2"
)

method_lw <- c(
  "FFGL"            = 0.45,
  "FFGL2"           = 0.45,
  "FGL"             = 0.70,
  "FuDGE"           = 0.70,
  "Proposed Method" = 1.10
)

plot_obj <- ggplot(
  try3[try3$network=="DIFF", ],
  aes(x = hyper, y = med_F1, group = method, color = method, linewidth = method)
) +
  geom_line() +
  geom_point(size = 2.3) +
  # If you want error bars back:
  # geom_errorbar(aes(ymin = min_F1, ymax = max_F1), width = 2, linewidth = 0.3) +
  facet_wrap(~ p, scales = "free_x", nrow = 1, labeller = label_both) +
  scale_color_manual(values = method_cols, drop = FALSE) +
  scale_linewidth_manual(values = method_lw, guide = "none") +  # hide linewidth legend (cleaner)
  labs(
    title = "Simulation F1 Performance Metrics different hyper metrics - p=25, n=10",  # journals often prefer title in caption, not inside the plot
    x = "Hyperparam",
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

plot_obj


#################### 
# SETTING 1 - Comp time
####################


computational_times <- read_csv("computational_times.csv")
computational_times$comp_time_sec <- computational_times$comp_time
computational_times <- computational_times %>%
  dplyr::select(p,n,comp_time_sec) %>%
  group_by(p, n) %>%
  summarise(
    med_comp_time_sec= mean(comp_time_sec, na.rm = TRUE),
    sd_comp_time_sec = sd(comp_time_sec, na.rm = TRUE),
    max_comp_time_sec = max(comp_time_sec, na.rm = TRUE),
    min_comp_time_sec = min(comp_time_sec, na.rm = TRUE),
    .groups = "drop"
  )

computational_times$method <- "Our Method"

computational_times_litt_comp <- read_csv("computational_times_litt_comp.csv")
computational_times_litt_comp$comp_time_sec <- computational_times_litt_comp$comp_time
computational_times_litt_comp[computational_times_litt_comp$comp_time_sec<5, ]$comp_time_sec <- computational_times_litt_comp[computational_times_litt_comp$comp_time_sec<5, ]$comp_time_sec*60
computational_times_litt_comp <- computational_times_litt_comp %>%
  dplyr::select(p,n,comp_time_sec) %>%
  group_by(p, n) %>%
  summarise(
    med_comp_time_sec= mean(comp_time_sec, na.rm = TRUE),
    sd_comp_time_sec = sd(comp_time_sec, na.rm = TRUE),
    max_comp_time_sec = max(comp_time_sec, na.rm = TRUE),
    min_comp_time_sec = min(comp_time_sec, na.rm = TRUE),
    .groups = "drop"
  )

computational_times_litt_comp$method  <- "FuDGE"

computational_times_full <- rbind(computational_times, computational_times_litt_comp)

str(computational_times_full)

method_levels <- c("Our Method", "FuDGE")
computational_times_full$method <- factor(computational_times_full$method,
                                              levels = method_levels)

# 2) Colorblind-friendly palette (edit as you like)
method_cols <- c(
  "FuDGE"           = "#C10A0A",
  "Our Method" = "#0072B2"
)

method_lw <- c(
  "FuDGE"           = 0.70,
  "Our Method" = 1.10
)

computational_times_red <- computational_times_full %>%
  filter(n == 200)


# Figure 4
plot_obj <- ggplot(
  computational_times_red,
  aes(x = p, y = med_comp_time_sec, group = method, color = method, linewidth = method)
) +
  geom_line() +
  geom_point(size = 2.3) +
  # If you want error bars back:
  geom_errorbar(aes(ymin = min_comp_time_sec, ymax = max_comp_time_sec), width = 2, linewidth = 0.3) +
  #facet_wrap(~ p, scales = "free_x", nrow = 1, labeller = label_both) +
  scale_color_manual(values = method_cols, drop = FALSE) +
  scale_linewidth_manual(values = method_lw, guide = "none") +  # hide linewidth legend (cleaner)
  labs(
    title = NULL,  # journals often prefer title in caption, not inside the plot
    x = "Network size",
    y = "Computational time",
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

plot_obj



#################### 
# SETTING 2 
####################

base_path <- "Simulation_studies/Step2/"
setwd(base_path)

# Get list of folders (hotspots)
sel_dirs <- list.dirs(base_path, recursive = FALSE)

results_list <- list()

for (dir_path in sel_dirs) {
  file_list <- list.files(dir_path, pattern = "^test_results_metices.*\\.csv$", full.names = TRUE)
  if (length(file_list) == 1) {
    file_path <- file_list[1]
  } else if (length(file_list) > 1) {
    warning("Multiple matching files found in ", dir_path, "; using the first.")
    file_path <- file_list[1]
  } else {
    next  # No matching file, skip this directory
  }
  if (file.exists(file_path) ) {
    tryCatch({
      df <- read_csv(file_path, show_col_types = FALSE)
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
      df$simulation_type <- str_sub(df$simulation, -2)
      df$p <- as.numeric(str_extract(df$simulation, "(?<=p)\\d+(?=_)"))
      results_list[[length(results_list) + 1]] <- df
    }, error = function(e) {
      message("Error reading: ", file_path)
    })
  }
}

# Combine all results
combined_data <- bind_rows(results_list)

combined_data <- combined_data %>%
  filter(symm == "OR")

p_levels <- c("10", "50")
combined_data$p <- factor(combined_data$p, levels = p_levels)

combined_data$network <- factor(combined_data$network)

p_cols <- c(
  "10" = "#89CFF1",
  "50" = "#0072B2"
)

library(latex2exp)

network_labels <- c(
  DIFF  = "Differential",
  GROUP = "Group $\\textit{G}_1$",
  POP   = "Population"
)

pd <- position_dodge(width = 0.55)

# Figure 5

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

plot_obj


