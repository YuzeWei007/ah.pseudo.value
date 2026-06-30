# ============================================================
# main.R
# Main script for PO-only / AH-Cox-only / Stacked-GMM-AH-Cox simulation
#
# All helper functions are stored in utility.R.
# Put main.R and utility.R in the same folder, then run main.R.
# ============================================================

rm(list = ls())

library(survival)
library(pseudo)
library(rootSolve)
library(eventglm)

# Load utility functions from the same folder as this script when possible.
this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)

if (!is.na(this_file)) {
  script_dir <- dirname(this_file)
  source(file.path(script_dir, "Version_0604_utility.R"))
} else {
  source("Version_0604_utility.R")
}

target_dir <- "/Users/JasonWei/Desktop/Spring lab/529"
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
setwd(target_dir)

# ============================================================
# Main simulation settings
# ============================================================

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_prefix <- "po_ahcox_gmm_ahcox_uno8_standalone"

cat("\nRun ID:\n")
print(run_id)

cat("\nWorking directory:\n")
print(getwd())

cat("\nAH implementation:\n")
print("Functions loaded from utility.R")

set.seed(2026)

beta_true <- c(
  Intercept = -1.2340,
  Age       =  0.0387,
  LogBili   =  0.8371,
  Albumin   = -1.1590
)

tau <- 7
n_sample <- 300

# Debug: use 50 or 100.
# Final: use 1000 or 5000.
nsim <- 1000

conf_level <- 0.95
z_crit <- qnorm(1 - (1 - conf_level) / 2)

# ============================================================
# Method flags
# ============================================================

run_PO_regular <- FALSE
run_PO_group_specific <- FALSE
run_PO_covariate_dependent <- FALSE
run_AHCox <- FALSE
run_GMM <- FALSE
run_GMM_PO_GS <- TRUE

#po_cd_ipcw_method <- "hajek"
po_cd_ipcw_method <- "binder"

methods <- c(
  if (run_PO_regular) "PO-only",
  if (run_PO_group_specific) "PO_group_specific",
  if (run_PO_covariate_dependent) "PO_covariate_dependent",
  if (run_AHCox) "AH-Cox-only",
  if (run_GMM) "Stacked-GMM-AH-Cox",
  if (run_GMM_PO_GS) "Stacked-GMM-POGS-AH-Cox"
)


# ============================================================
# PBC covariate pool
# ============================================================

pbc_pool <- prepare_pbc_pool()

cat("\nPBC covariate pool size:\n")
print(nrow(pbc_pool$data))

cat("\nAge range after possible conversion:\n")
print(range(pbc_pool$data$age))


# ============================================================
# Full setting grid
# ============================================================

setting_grid <- data.frame(
  Model = rep(c("Model 1", "Model 2"), each = 4),
  Weibull_shape = rep(c(1.0, 2.0), each = 4),
  Censoring = rep(
    c("none", "independent", "group_specific", "covariate_dependent"),
    times = 2
  ),
  Censoring_label = rep(
    c(
      "(a) no censoring",
      "(b) independent censoring",
      "(c) group-specific censoring",
      "(d) covariate-dependent censoring"
    ),
    times = 2
  ),
  IPCW_method = rep("AH-Cox via embedded fun-ahreg.R", 8),
  nsim = nsim,
  stringsAsFactors = FALSE
)

cat("\nFull setting grid:\n")
print(setting_grid)

# ============================================================
# Run all settings
# ============================================================

cat("\n############################################################\n")
cat("START PO / AH-COX / GMM-AH-COX 8-SETTING SIMULATION\n")
cat("############################################################\n")

full_start_time <- Sys.time()

all_results <- vector("list", nrow(setting_grid))
setting_times <- vector("list", nrow(setting_grid))
summary_list <- vector("list", nrow(setting_grid))

names(all_results) <- paste0(
  setting_grid$Model,
  "_shape", setting_grid$Weibull_shape,
  "_", setting_grid$Censoring
)

for (ii in seq_len(nrow(setting_grid))) {
  this_setting <- setting_grid[ii, , drop = FALSE]
  
  cat("\n============================================================\n")
  cat("Running setting", ii, "of", nrow(setting_grid), "\n")
  cat(
    this_setting$Model,
    ", shape=", this_setting$Weibull_shape,
    ", ", this_setting$Censoring_label,
    ", ipcw=", this_setting$IPCW_method,
    ", nsim=", this_setting$nsim,
    "\n",
    sep = ""
  )
  cat("============================================================\n")
  
  time_one <- system.time({
    res_one <- run_sim(
      nsim = this_setting$nsim,
      n = n_sample,
      beta = beta_true,
      tau = tau,
      pbc_pool = pbc_pool,
      shape = this_setting$Weibull_shape,
      censoring = this_setting$Censoring,
      run_PO_regular = run_PO_regular,
      run_PO_group_specific = run_PO_group_specific,
      run_PO_covariate_dependent = run_PO_covariate_dependent,
      run_AHCox = run_AHCox,
      run_GMM = run_GMM,
      run_GMM_PO_GS = run_GMM_PO_GS,
      po_cd_ipcw_method = po_cd_ipcw_method
    )
  })
  
  all_results[[ii]] <- res_one
  
  summary_one <- summ_all(res_one, this_setting)
  summary_list[[ii]] <- summary_one
  
  setting_times[[ii]] <- data.frame(
    Model = this_setting$Model,
    Weibull_shape = this_setting$Weibull_shape,
    Censoring = this_setting$Censoring,
    Censoring_label = this_setting$Censoring_label,
    IPCW_method = this_setting$IPCW_method,
    nsim = this_setting$nsim,
    user = unname(time_one["user.self"]),
    system = unname(time_one["sys.self"]),
    elapsed = unname(time_one["elapsed"]),
    Mean_observed_random_censor_rate = mean(
      res_one$res$observed_random_censor_rate,
      na.rm = TRUE
    ),
    Mean_potential_censor_rate_at_tau = mean(
      res_one$res$potential_censor_rate_at_tau,
      na.rm = TRUE
    ),
    Mean_total_censor_rate = mean(
      res_one$res$total_censor_rate,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
  
  partial_summary <- do.call(rbind, summary_list[seq_len(ii)])
  partial_times <- do.call(rbind, setting_times[seq_len(ii)])
  
  write.csv(
    partial_summary,
    paste0(output_prefix, "_summary_partial_", run_id, ".csv"),
    row.names = FALSE
  )
  
  write.csv(
    partial_times,
    paste0(output_prefix, "_setting_times_partial_", run_id, ".csv"),
    row.names = FALSE
  )
  
  saveRDS(
    all_results,
    paste0(output_prefix, "_checkpoint_", run_id, ".rds")
  )
  
  cat("\nFinished setting", ii, "of", nrow(setting_grid), "\n")
  cat("Runtime:\n")
  print(time_one)
  
  cat("\nPartial summary for this setting:\n")
  print(summary_one)
  
  cat("\nCheckpoint saved.\n")
}

full_end_time <- Sys.time()

# ============================================================
# Combine and save final results
# ============================================================

final_summary <- do.call(rbind, summary_list)
final_setting_times <- do.call(rbind, setting_times)

final_summary_rounded <- final_summary
num_cols <- sapply(final_summary_rounded, is.numeric)

final_summary_rounded[num_cols] <- lapply(
  final_summary_rounded[num_cols],
  round,
  4
)

rel_bias_wide <- make_wide_table(final_summary_rounded, "Rel_Bias")
coverage_wide <- make_wide_table(final_summary_rounded, "Coverage")
estimate_wide <- make_wide_table(final_summary_rounded, "Estimate")
sd_wide <- make_wide_table(final_summary_rounded, "Empirical_SD")
mean_se_wide <- make_wide_table(final_summary_rounded, "Mean_SE")

comparison_table <- make_comparison_table(final_summary_rounded)

target_check_table <- make_target_check_table(
  final_summary_rounded,
  conf_level = conf_level,
  bias_tolerance = 0.03,
  coverage_tolerance = 0.03
)

censoring_check <- final_setting_times[, c(
  "Model",
  "Weibull_shape",
  "Censoring",
  "Censoring_label",
  "IPCW_method",
  "nsim",
  "Mean_observed_random_censor_rate",
  "Mean_potential_censor_rate_at_tau",
  "Mean_total_censor_rate"
)]

cat("\n############################################################\n")
cat("FINISHED PO / AH-COX / GMM-AH-COX 8-SETTING SIMULATION\n")
cat("############################################################\n")

cat("\nStart time:\n")
print(full_start_time)

cat("\nEnd time:\n")
print(full_end_time)

cat("\nTotal elapsed time:\n")
print(full_end_time - full_start_time)

cat("\nFinal summary dimension:\n")
print(dim(final_summary_rounded))
cat("Expected output slots: 8 settings x 5 estimator slots x 4 parameters = 160 rows\n")
cat("Only active methods are actually fit; inactive methods appear as NA with Success = 0.\n")

cat("\nWide table dimensions:\n")
cat("Relative bias wide:\n")
print(dim(rel_bias_wide))
cat("Expected: 8 settings x 4 parameters = 32 rows\n")

cat("\nCensoring check:\n")
print(censoring_check)

cat("\nRelative bias wide table:\n")
print(rel_bias_wide)

cat("\nCoverage wide table:\n")
print(coverage_wide)

cat("\nEmpirical SD wide table:\n")
print(sd_wide)

cat("\nTarget check table:\n")
print(target_check_table)

# ============================================================
# Save final files
# ============================================================

write.csv(
  final_summary_rounded,
  paste0(output_prefix, "_summary_long_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  final_setting_times,
  paste0(output_prefix, "_setting_times_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  rel_bias_wide,
  paste0(output_prefix, "_relative_bias_wide_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  coverage_wide,
  paste0(output_prefix, "_coverage_wide_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  estimate_wide,
  paste0(output_prefix, "_estimate_wide_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  sd_wide,
  paste0(output_prefix, "_empirical_sd_wide_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  mean_se_wide,
  paste0(output_prefix, "_mean_se_wide_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  comparison_table,
  paste0(output_prefix, "_comparison_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  target_check_table,
  paste0(output_prefix, "_target_check_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  censoring_check,
  paste0(output_prefix, "_censoring_check_", run_id, ".csv"),
  row.names = FALSE
)

saveRDS(
  all_results,
  paste0(output_prefix, "_raw_results_", run_id, ".rds")
)

cat("\nSaved final files:\n")
cat("  ", paste0(output_prefix, "_summary_long_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_setting_times_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_relative_bias_wide_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_coverage_wide_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_estimate_wide_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_empirical_sd_wide_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_mean_se_wide_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_comparison_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_target_check_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_censoring_check_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_raw_results_", run_id, ".rds\n"), sep = "")
cat("  ", paste0(output_prefix, "_checkpoint_", run_id, ".rds\n"), sep = "")

cat("\nFiles saved in:\n")
print(getwd())

cat("\nFiles generated in this run:\n")
print(list.files(pattern = paste0(output_prefix, ".*", run_id)))