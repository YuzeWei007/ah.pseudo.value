# ============================================================
# POPLAR data analysis for our method
# PO-only / AH-only / Stacked-GMM-AH
#
# Purpose:
#   Apply our proposed real-data analysis pipeline to the POPLAR OS data.
#
# Methods:
#   1. PO-only
#   2. AH-only, independent-censoring AH regression
#   3. Stacked-GMM-PO-AH
#
# Important:
#   - This script DOES NOT modify utility.R.
#   - It sources your utility file and uses the core functions defined there:
#       ahreg()
#       get_PO()
#       fit_PO_with_se()
#       moment_PO()
#       moment_AHCox()
#       standardize_moment()
#       numerical_jacobian()
#       safe_solve()
#       is_bad_beta()
#   - The real-data wrapper functions are defined inside this script.
#
# Required files in target_dir:
#   1. Version_0604_utility.R
#   2. 41591_2018_134_MOESM3_ESM.xlsx
#
# Main outputs:
#   1. KM curve PDF
#   2. risk-set table
#   3. censoring association check table
#   4. Cox reference table
#   5. full method comparison table
#   6. treatment-effect table
#   7. average-subject AH table
#   8. saved fit object
# ============================================================

rm(list = ls())

library(openxlsx)
library(survival)
library(pseudo)
library(rootSolve)

# ============================================================
# 0. Paths and settings
# ============================================================

target_dir <- "/Users/JasonWei/Desktop/文件整理"
setwd(target_dir)

utility_file <- "Version_0604_utility.R"

if (!file.exists(utility_file)) {
  stop(
    paste0(
      "\nCannot find utility file: ", utility_file,
      "\nCurrent working directory: ", getwd(),
      "\nFiles containing 'utility' in this folder:\n",
      paste(list.files(pattern = "utility", ignore.case = TRUE), collapse = "\n"),
      "\n\nPlease edit utility_file or move the utility file into target_dir.\n"
    )
  )
}

source(utility_file)

excel_file <- "41591_2018_134_MOESM3_ESM.xlsx"
excel_sheet <- "POPLAR_Clinical_Data"

tau <- 20
conf_level <- 0.95
z_crit <- qnorm(1 - (1 - conf_level) / 2)

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_prefix <- "poplar_our_method_DA"

# ============================================================
# 1. Read POPLAR data
# ============================================================

D <- as.data.frame(read.xlsx(excel_file, sheet = excel_sheet))

cat("\n============================================================\n")
cat("Raw POPLAR data\n")
cat("============================================================\n")
cat("Dimension:\n")
print(dim(D))

cat("\nVariable names:\n")
print(names(D))

# Do not print head(D), because it may display patient-level rows.

# ============================================================
# 2. Helper functions
# ============================================================

upper_trim <- function(x) {
  toupper(trimws(as.character(x)))
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

extract_first_number <- function(x) {
  x <- upper_trim(x)
  out <- suppressWarnings(as.numeric(sub(".*?([0-9]+\\.?[0-9]*).*", "\\1", x)))
  out[!grepl("[0-9]", x)] <- NA_real_
  out
}

format_p <- function(p) {
  ifelse(p < 0.001, "<.001", sub("^0", "", sprintf("%.3f", p)))
}

round_numeric_cols <- function(dat, digits = 4) {
  num_cols <- sapply(dat, is.numeric)
  dat[num_cols] <- lapply(dat[num_cols], round, digits)
  dat
}

print_raw_coding_tables <- function(D) {
  cat("\n============================================================\n")
  cat("Raw variable tables for coding check\n")
  cat("============================================================\n")
  
  cat("\nTRT01P:\n")
  print(table(D$TRT01P, useNA = "always"))
  
  cat("\nSEX:\n")
  print(table(D$SEX, useNA = "always"))
  
  cat("\nHIST:\n")
  print(table(D$HIST, useNA = "always"))
  
  cat("\nECOGGR:\n")
  print(table(D$ECOGGR, useNA = "always"))
  
  cat("\nrace2:\n")
  print(table(D$race2, useNA = "always"))
  
  cat("\nOS.CNSR:\n")
  print(table(D$OS.CNSR, useNA = "always"))
  
  cat("\nOS summary:\n")
  print(summary(D$OS))
}

select_estimable_covariates <- function(dat, candidate_covariates, force_keep = "trt") {
  selected <- candidate_covariates
  
  unique_counts <- sapply(dat[, selected, drop = FALSE], function(x) length(unique(x)))
  constant_vars <- names(unique_counts)[unique_counts < 2]
  
  if (length(constant_vars) > 0) {
    cat("\nWARNING: Dropping constant covariates:\n")
    print(constant_vars)
    selected <- setdiff(selected, constant_vars)
  }
  
  if (!(force_keep %in% selected)) {
    stop(
      paste0(
        "\nTreatment variable ", force_keep, " is not estimable.",
        "\nCheck treatment coding."
      )
    )
  }
  
  repeat {
    Z_tmp <- cbind(Intercept = 1, as.matrix(dat[, selected, drop = FALSE]))
    storage.mode(Z_tmp) <- "numeric"
    
    qr_rank <- qr(Z_tmp)$rank
    
    if (qr_rank == ncol(Z_tmp)) {
      break
    }
    
    drop_candidates <- setdiff(selected, force_keep)
    
    if (length(drop_candidates) == 0) {
      stop("\nDesign matrix remains rank deficient even after keeping only treatment.")
    }
    
    var_to_drop <- tail(drop_candidates, 1)
    
    cat("\nWARNING: Design matrix is rank deficient.\n")
    cat("Dropping covariate to restore full rank:\n")
    print(var_to_drop)
    
    selected <- setdiff(selected, var_to_drop)
  }
  
  selected
}

make_Z <- function(dat, selected_covariates) {
  Z <- cbind(Intercept = 1, as.matrix(dat[, selected_covariates, drop = FALSE]))
  storage.mode(Z) <- "numeric"
  Z
}

make_display_labels <- function(selected_covariates) {
  label_map <- c(
    trt = "Atezolizumab vs. docetaxel",
    age = "Age (year)",
    male = "Male vs. female",
    squamous = "Squamous vs. other",
    ecog1 = "ECOG PS (0 vs. >=1)",
    white = "White vs. non-White"
  )
  
  c("Intercept", unname(label_map[selected_covariates]))
}

make_km_plot <- function(dat, output_file) {
  fit <- survfit(Surv(time, status) ~ trt, data = dat)
  
  pdf(output_file, width = 7, height = 5)
  plot(
    fit,
    xlab = "Months",
    ylab = "Overall survival probability",
    lty = 1:2,
    mark.time = TRUE
  )
  legend(
    "topright",
    legend = c("Docetaxel", "Atezolizumab"),
    lty = 1:2,
    bty = "n"
  )
  dev.off()
  
  invisible(fit)
}

make_risk_set_table <- function(dat, tau_values = c(12, 18, 20, 21, 24, 25)) {
  fit <- survfit(Surv(time, status) ~ 1, data = dat)
  ss <- summary(fit, times = tau_values)
  
  data.frame(
    Tau = ss$time,
    N_risk = ss$n.risk,
    N_event = ss$n.event,
    N_censor = ss$n.censor
  )
}

fit_cox_reference <- function(dat, selected_covariates) {
  fmla <- as.formula(
    paste0(
      "Surv(time, status) ~ ",
      paste(selected_covariates, collapse = " + ")
    )
  )
  
  fit <- coxph(fmla, data = dat)
  s <- summary(fit)
  
  out <- data.frame(
    Variable = rownames(s$coef),
    HR = s$coef[, "exp(coef)"],
    Lower_95 = s$conf.int[, "lower .95"],
    Upper_95 = s$conf.int[, "upper .95"],
    P_value = s$coef[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  
  list(fit = fit, table = out)
}

check_censoring_association <- function(dat, selected_covariates) {
  dat$censor_status <- 1 - dat$status
  
  fmla <- as.formula(
    paste0(
      "Surv(time, censor_status) ~ ",
      paste(selected_covariates, collapse = " + ")
    )
  )
  
  fit <- coxph(fmla, data = dat)
  s <- summary(fit)
  
  out <- data.frame(
    Variable = rownames(s$coef),
    HR_for_censoring = s$coef[, "exp(coef)"],
    Lower_95 = s$conf.int[, "lower .95"],
    Upper_95 = s$conf.int[, "upper .95"],
    P_value = s$coef[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  
  list(fit = fit, table = out)
}

check_PH_for_treatment <- function(dat) {
  fit <- coxph(Surv(time, status) ~ trt, data = dat)
  zph <- cox.zph(fit)
  
  out <- data.frame(
    Test = rownames(zph$table),
    Chisq = zph$table[, "chisq"],
    P_value = zph$table[, "p"],
    stringsAsFactors = FALSE
  )
  
  list(fit = fit, zph = zph, table = out)
}

# ============================================================
# 3. Generic GMM functions for real data
# ============================================================

get_center_generic <- function(Z) {
  as.numeric(colMeans(Z[, -1, drop = FALSE], na.rm = TRUE))
}

make_Z_centered_generic <- function(Z_raw, center_vec) {
  Zc <- Z_raw
  Zc[, -1] <- sweep(Z_raw[, -1, drop = FALSE], 2, center_vec, "-")
  colnames(Zc) <- colnames(Z_raw)
  Zc
}

raw_beta_to_centered_generic <- function(beta_raw, center_vec) {
  beta_c <- beta_raw
  beta_c[1] <- beta_raw[1] + sum(center_vec * beta_raw[-1])
  beta_c[-1] <- beta_raw[-1]
  names(beta_c) <- names(beta_raw)
  beta_c
}

centered_beta_to_raw_generic <- function(beta_c, center_vec) {
  beta_raw <- beta_c
  beta_raw[1] <- beta_c[1] - sum(center_vec * beta_c[-1])
  beta_raw[-1] <- beta_c[-1]
  names(beta_raw) <- names(beta_c)
  beta_raw
}

centered_var_to_raw_generic <- function(V_c, center_vec) {
  p <- nrow(V_c)
  L <- diag(p)
  L[1, -1] <- -center_vec
  L %*% V_c %*% t(L)
}

fit_GMM_real_with_se <- function(Z,
                                 po,
                                 ipcw,
                                 start,
                                 center_vec = NULL,
                                 use_standardization = TRUE,
                                 beta_limit = 20) {
  p <- length(start)
  
  if (is.null(po) || is.null(ipcw)) {
    return(list(beta = rep(NA_real_, p), se = rep(NA_real_, p)))
  }
  
  if (is.null(center_vec)) {
    center_vec <- get_center_generic(Z)
  }
  
  Zc <- make_Z_centered_generic(Z, center_vec)
  start_c <- raw_beta_to_centered_generic(start, center_vec)
  
  build_M <- function(b, scale_use = NULL) {
    M_po <- moment_PO(b, Zc, po)
    M_ah <- moment_AHCox(b, Zc, ipcw)
    M <- cbind(M_po, M_ah)
    
    if (!use_standardization) {
      return(list(M = M, scale = rep(1, ncol(M))))
    }
    
    if (is.null(scale_use)) {
      standardize_moment(M)
    } else {
      standardize_moment(M, scale_use)
    }
  }
  
  M0 <- build_M(start_c)
  scale0 <- M0$scale
  
  Q_identity <- function(b) {
    if (length(b) != length(start_c)) return(1e20)
    if (is_bad_beta(b, limit = beta_limit)) return(1e20)
    
    M <- build_M(b, scale0)$M
    g <- colMeans(M)
    
    if (any(!is.finite(g))) return(1e20)
    
    sum(g^2)
  }
  
  fit1 <- try(
    optim(
      par = start_c,
      fn = Q_identity,
      method = "BFGS",
      control = list(maxit = 1000, reltol = 1e-10)
    ),
    silent = TRUE
  )
  
  if (inherits(fit1, "try-error")) {
    return(list(beta = rep(NA_real_, p), se = rep(NA_real_, p)))
  }
  
  b_c <- as.numeric(fit1$par)
  names(b_c) <- names(start)
  
  if (is_bad_beta(b_c, limit = beta_limit)) {
    return(list(beta = rep(NA_real_, p), se = rep(NA_real_, p)))
  }
  
  b_raw <- centered_beta_to_raw_generic(b_c, center_vec)
  
  if (is_bad_beta(b_raw, limit = beta_limit)) {
    return(list(beta = rep(NA_real_, p), se = rep(NA_real_, p)))
  }
  
  g_fun <- function(b) colMeans(build_M(b, scale0)$M)
  D <- try(numerical_jacobian(g_fun, b_c), silent = TRUE)
  
  if (inherits(D, "try-error") || any(!is.finite(D))) {
    return(list(beta = b_raw, se = rep(NA_real_, p)))
  }
  
  M_hat <- build_M(b_c, scale0)$M
  n <- nrow(M_hat)
  S_hat <- crossprod(scale(M_hat, center = TRUE, scale = FALSE)) / n
  
  bread <- safe_solve(t(D) %*% D, ridge = 1e-8)
  
  if (is.null(bread)) {
    return(list(beta = b_raw, se = rep(NA_real_, p)))
  }
  
  V_c <- bread %*% t(D) %*% S_hat %*% D %*% bread / n
  V_raw <- centered_var_to_raw_generic(V_c, center_vec)
  
  se_raw <- sqrt(diag(V_raw))
  se_raw[!is.finite(se_raw)] <- NA_real_
  
  names(b_raw) <- names(start)
  names(se_raw) <- names(start)
  
  list(beta = b_raw, se = se_raw)
}

# ============================================================
# 4. Fit PO, AH, and GMM
# ============================================================

fit_AH_independent_real <- function(dat, covariates, tau) {
  fit <- ahreg(
    time = dat$time,
    status = dat$status,
    covariates = covariates,
    tau = tau,
    covariates4cens = NULL
  )
  
  beta <- as.numeric(fit$beta)
  se <- as.numeric(fit$result$SE)
  
  list(
    fit = fit,
    beta = beta,
    se = se
  )
}

make_independent_ipcw_from_ahfit <- function(ah_fit, dat, tau) {
  wt <- ah_fit$fit$ipcw_information$ipcw
  
  if (is.null(wt) || length(wt) != nrow(dat) || any(!is.finite(wt))) {
    return(NULL)
  }
  
  list(
    weight = as.numeric(wt),
    event_tau = as.numeric(dat$time < tau & dat$status == 1),
    time_tau = pmin(dat$time, tau)
  )
}

fit_all_methods_real <- function(dat, selected_covariates, tau) {
  Z <- make_Z(dat, selected_covariates)
  covariates <- as.matrix(dat[, selected_covariates, drop = FALSE])
  storage.mode(covariates) <- "numeric"
  
  po <- get_PO(
    X = dat$time,
    Delta = dat$status,
    tau = tau
  )
  
  ah <- fit_AH_independent_real(
    dat = dat,
    covariates = covariates,
    tau = tau
  )
  
  names(ah$beta) <- colnames(Z)
  names(ah$se) <- colnames(Z)
  
  ipcw <- make_independent_ipcw_from_ahfit(
    ah_fit = ah,
    dat = dat,
    tau = tau
  )
  
  if (all(is.finite(ah$beta))) {
    start <- ah$beta
  } else {
    start <- rep(0, ncol(Z))
    names(start) <- colnames(Z)
  }
  
  po_fit <- fit_PO_with_se(
    Z = Z,
    po = po,
    start = start
  )
  
  names(po_fit$beta) <- colnames(Z)
  names(po_fit$se) <- colnames(Z)
  
  gmm_fit <- fit_GMM_real_with_se(
    Z = Z,
    po = po,
    ipcw = ipcw,
    start = start,
    center_vec = get_center_generic(Z)
  )
  
  names(gmm_fit$beta) <- colnames(Z)
  names(gmm_fit$se) <- colnames(Z)
  
  list(
    Z = Z,
    covariates = covariates,
    po = po,
    ipcw = ipcw,
    fits = list(
      "PO-only" = po_fit,
      "AH-only" = list(beta = ah$beta, se = ah$se, fit = ah$fit),
      "Stacked-GMM-PO-AH" = gmm_fit
    )
  )
}

make_method_regression_table <- function(fit_obj, display_labels, conf_level = 0.95) {
  z <- qnorm(1 - (1 - conf_level) / 2)
  
  out <- data.frame()
  
  for (method_name in names(fit_obj$fits)) {
    beta <- fit_obj$fits[[method_name]]$beta
    se <- fit_obj$fits[[method_name]]$se
    
    lower <- beta - z * se
    upper <- beta + z * se
    pval <- 2 * (1 - pnorm(abs(beta / se)))
    
    tmp <- data.frame(
      Method = method_name,
      Variable = display_labels,
      Beta = as.numeric(beta),
      SE = as.numeric(se),
      AH_ratio = exp(as.numeric(beta)),
      Lower_95 = exp(as.numeric(lower)),
      Upper_95 = exp(as.numeric(upper)),
      CI_95 = paste0(
        "(",
        sprintf("%.3f", exp(as.numeric(lower))),
        "-",
        sprintf("%.3f", exp(as.numeric(upper))),
        ")"
      ),
      P_value = as.numeric(pval),
      P_value_print = format_p(pval),
      stringsAsFactors = FALSE
    )
    
    out <- rbind(out, tmp)
  }
  
  out
}

make_average_subject_AH_table <- function(fit_obj, dat, selected_covariates) {
  Z <- fit_obj$Z
  
  x_doc <- rep(NA_real_, ncol(Z))
  names(x_doc) <- colnames(Z)
  
  x_doc["Intercept"] <- 1
  
  for (v in selected_covariates) {
    if (v == "trt") {
      x_doc[v] <- 0
    } else {
      x_doc[v] <- mean(dat[[v]], na.rm = TRUE)
    }
  }
  
  x_atezo <- x_doc
  x_atezo["trt"] <- 1
  
  out <- data.frame()
  
  for (method_name in names(fit_obj$fits)) {
    beta <- fit_obj$fits[[method_name]]$beta
    
    ah_doc <- as.numeric(exp(sum(x_doc * beta)))
    ah_atezo <- as.numeric(exp(sum(x_atezo * beta)))
    
    tmp <- data.frame(
      Method = method_name,
      AH_Docetaxel = ah_doc,
      AH_Atezolizumab = ah_atezo,
      AH_difference_Atezolizumab_minus_Docetaxel = ah_atezo - ah_doc,
      AH_ratio_Atezolizumab_vs_Docetaxel = ah_atezo / ah_doc,
      stringsAsFactors = FALSE
    )
    
    out <- rbind(out, tmp)
  }
  
  out
}

# ============================================================
# 5. Recode analysis variables
# ============================================================

print_raw_coding_tables(D)

trt_raw <- upper_trim(D$TRT01P)
sex_raw <- upper_trim(D$SEX)
hist_raw <- upper_trim(D$HIST)
race_raw <- upper_trim(D$race2)

ecog_num <- safe_numeric(D$ECOGGR)
if (all(is.na(ecog_num))) {
  ecog_num <- extract_first_number(D$ECOGGR)
}

dat <- data.frame(
  time = safe_numeric(D$OS),
  status = 1 - safe_numeric(D$OS.CNSR),
  trt = as.numeric(trt_raw == "MPDL3280A"),
  age = safe_numeric(D$BAGE),
  male = as.numeric(sex_raw %in% c("MALE", "M")),
  squamous = as.numeric(grepl("^SQUAM", hist_raw)),
  ecog1 = as.numeric(ecog_num >= 1),
  white = as.numeric(grepl("WHITE", race_raw)),
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("Analysis variables before complete-case filtering\n")
cat("============================================================\n")
cat("\nMissing values:\n")
print(colSums(is.na(dat)))

dat <- dat[complete.cases(dat), ]

cat("\n============================================================\n")
cat("Analysis data after complete-case filtering\n")
cat("============================================================\n")
cat("Dimension:\n")
print(dim(dat))

cat("\nCoding checks:\n")
cat("status: 1 = death/event, 0 = censored\n")
print(table(dat$status, useNA = "ifany"))

cat("\ntrt: 1 = atezolizumab / MPDL3280A, 0 = docetaxel\n")
print(table(dat$trt, useNA = "ifany"))

cat("\nmale: 1 = male, 0 = female\n")
print(table(dat$male, useNA = "ifany"))

cat("\nsquamous: 1 = squamous, 0 = other\n")
print(table(dat$squamous, useNA = "ifany"))

cat("\necog1: 1 = ECOG >= 1, 0 = ECOG 0\n")
print(table(dat$ecog1, useNA = "ifany"))

cat("\nwhite: 1 = White, 0 = non-White\n")
print(table(dat$white, useNA = "ifany"))

stopifnot(all(dat$time >= 0))
stopifnot(all(dat$status %in% c(0, 1)))
stopifnot(all(dat$trt %in% c(0, 1)))
stopifnot(all(dat$male %in% c(0, 1)))
stopifnot(all(dat$squamous %in% c(0, 1)))
stopifnot(all(dat$ecog1 %in% c(0, 1)))
stopifnot(all(dat$white %in% c(0, 1)))

# ============================================================
# 6. Select covariates and run diagnostics
# ============================================================

candidate_covariates <- c("trt", "age", "male", "squamous", "ecog1", "white")

selected_covariates <- select_estimable_covariates(
  dat = dat,
  candidate_covariates = candidate_covariates,
  force_keep = "trt"
)

display_labels <- make_display_labels(selected_covariates)

cat("\n============================================================\n")
cat("Selected covariates used in our method DA\n")
cat("============================================================\n")
print(selected_covariates)

cat("\nDisplay labels:\n")
print(display_labels)

if (!("squamous" %in% selected_covariates)) {
  cat("\nNOTE: squamous was not included because it was not estimable in this imported data.\n")
  cat("This keeps the model estimable but means the adjusted model is not exactly the same as Uno Table 4.\n")
}

Z_check <- make_Z(dat, selected_covariates)
cat("\nDesign matrix rank:\n")
print(qr(Z_check)$rank)
cat("Number of columns:\n")
print(ncol(Z_check))

# ============================================================
# 7. Descriptive / diagnostic outputs
# ============================================================

km_file <- paste0(output_prefix, "_km_curve_", run_id, ".pdf")
km_fit <- make_km_plot(dat, km_file)

risk_set_table <- make_risk_set_table(
  dat,
  tau_values = c(12, 18, 20, 21, 24, 25)
)

cat("\n============================================================\n")
cat("Risk set table\n")
cat("============================================================\n")
print(risk_set_table)

cox_ref <- fit_cox_reference(dat, selected_covariates)
cox_table <- round_numeric_cols(cox_ref$table, 4)

cat("\n============================================================\n")
cat("Adjusted Cox reference table\n")
cat("============================================================\n")
print(cox_table)

ph_check <- check_PH_for_treatment(dat)
ph_check_table <- round_numeric_cols(ph_check$table, 4)

cat("\n============================================================\n")
cat("Treatment-only PH test\n")
cat("============================================================\n")
print(ph_check_table)

censor_check <- check_censoring_association(dat, selected_covariates)
censoring_check_table <- round_numeric_cols(censor_check$table, 4)

cat("\n============================================================\n")
cat("Censoring association check table\n")
cat("============================================================\n")
print(censoring_check_table)

# ============================================================
# 8. Run our methods
# ============================================================

cat("\n============================================================\n")
cat("Fitting PO-only / AH-only / Stacked-GMM-PO-AH\n")
cat("============================================================\n")

fit_methods <- fit_all_methods_real(
  dat = dat,
  selected_covariates = selected_covariates,
  tau = tau
)

method_table <- make_method_regression_table(
  fit_obj = fit_methods,
  display_labels = display_labels,
  conf_level = conf_level
)

method_table_rounded <- round_numeric_cols(method_table, 4)

cat("\n============================================================\n")
cat("Full method comparison table\n")
cat("============================================================\n")
print(method_table_rounded, row.names = FALSE)

treatment_effect_table <- method_table_rounded[
  method_table_rounded$Variable == "Atezolizumab vs. docetaxel",
]

cat("\n============================================================\n")
cat("Treatment effect table\n")
cat("============================================================\n")
print(treatment_effect_table, row.names = FALSE)

average_AH_table <- make_average_subject_AH_table(
  fit_obj = fit_methods,
  dat = dat,
  selected_covariates = selected_covariates
)

average_AH_table_rounded <- round_numeric_cols(average_AH_table, 4)

cat("\n============================================================\n")
cat("Average-subject AH table\n")
cat("============================================================\n")
print(average_AH_table_rounded, row.names = FALSE)

# ============================================================
# 9. Save outputs
# ============================================================

write.csv(
  risk_set_table,
  paste0(output_prefix, "_risk_set_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  cox_table,
  paste0(output_prefix, "_cox_reference_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  ph_check_table,
  paste0(output_prefix, "_ph_check_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  censoring_check_table,
  paste0(output_prefix, "_censoring_check_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  method_table_rounded,
  paste0(output_prefix, "_full_method_table_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  treatment_effect_table,
  paste0(output_prefix, "_treatment_effect_table_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  average_AH_table_rounded,
  paste0(output_prefix, "_average_AH_table_", run_id, ".csv"),
  row.names = FALSE
)

saveRDS(
  list(
    n = nrow(dat),
    tau = tau,
    selected_covariates = selected_covariates,
    display_labels = display_labels,
    km_fit = km_fit,
    risk_set_table = risk_set_table,
    cox_ref = cox_ref,
    ph_check = ph_check,
    censor_check = censor_check,
    fit_methods = fit_methods,
    method_table = method_table,
    treatment_effect_table = treatment_effect_table,
    average_AH_table = average_AH_table
  ),
  paste0(output_prefix, "_fit_object_", run_id, ".rds")
)

cat("\n============================================================\n")
cat("Saved files\n")
cat("============================================================\n")
cat("  ", km_file, "\n", sep = "")
cat("  ", paste0(output_prefix, "_risk_set_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_cox_reference_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_ph_check_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_censoring_check_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_full_method_table_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_treatment_effect_table_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_average_AH_table_", run_id, ".csv\n"), sep = "")
cat("  ", paste0(output_prefix, "_fit_object_", run_id, ".rds\n"), sep = "")

cat("\nFiles saved in:\n")
print(getwd())

cat("\nDone.\n")
