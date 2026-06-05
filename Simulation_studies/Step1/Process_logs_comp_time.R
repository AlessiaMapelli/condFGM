# This script processes the log files generated during the simulation studies to extract computational times for both FuDGE and our method, and saves the results in a CSV file for further analysis.
# Usage: Rscript Simulation_studies/Step1/Process_logs_comp_time.R <save_path>

args <- commandArgs(trailingOnly = TRUE)
save_path <- if (length(args) > 0) args[[1]] else stop("Please provide the path where results are saved.")

logs <- list.files(save_path, "\\.log$", recursive = TRUE, full.names = TRUE)

extract_log_info <- function(log_file_path) {
  sim_name_iter <- regmatches(log_file_path, regexpr("(p\\d+)_n(\\d+)_n\\d+/seed_(\\d+)", log_file_path))
  nums <- regmatches(sim_name_iter, gregexpr("\\d+", sim_name_iter))[[1]]
  p <- as.numeric(nums[1])
  n_g1 <- as.numeric(nums[2])
  n_g2 <- as.numeric(nums[3])
  iteration <- as.numeric(nums[4])

  is_fudge <- grepl("results_lit_comparison", log_file_path)
  pattern <- if (is_fudge) "FuDGE computational time:" else "Computational time:"
  method  <- if (is_fudge) "FuDGE" else "Our Method"

  compt_time_line <- grep(pattern, readLines(log_file_path, warn = FALSE), value = TRUE)
  comp_time <- as.numeric(sub(".* ([0-9.]+) seconds.*", "\\1", compt_time_line[1]))

  if (length(nums) == 4 && !is.na(comp_time)) {
    return(data.frame(p=p, n_g1=n_g1, n_g2=n_g2, iteration=iteration, method=method, comp_time=comp_time))
  } else {
    cat("Warning: Could not extract parameters or computational time from", log_file_path, "\n")
    return(NULL)
  }
}

rows <- lapply(logs, extract_log_info)
result <- aggregate(comp_time ~ p + n_g1 + n_g2 + iteration + method, do.call(rbind, Filter(Negate(is.null), rows)), max)
write.csv(result, file.path(save_path, "computational_times.csv"), row.names = FALSE)

