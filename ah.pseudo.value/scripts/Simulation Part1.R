# Load required libraries
library(survival)
library(pseudo)   # Used to compute pseudo-observations [cite: 101]
library(rootSolve) # Used to solve estimating equations [cite: 534]
library(knitr)
library(kableExtra)
library(gt)
library(knitr)

# ==========================================
# 1. Simulation setup (based on Uno 2024 DGP)
# ==========================================

# True regression coefficients: log(AH) = beta0 + beta1*age + beta2*log_bili + beta3*albumin [cite: 728]
true_beta <- c(-1.2340, 0.0387, 0.8371, -1.1590)

tau <- 7           # Truncation time [cite: 727]
n <- 300           # Sample size per replication [cite: 752]
nsim <- 1000       # Number of simulation repetitions

# Load and preprocess PBC data as the covariate pool [cite: 725, 731]
data(pbc)
pbc_clean <- na.omit(pbc[, c("time", "status", "age", "bili", "albumin")])
pbc_clean$log_bili <- log(pbc_clean$bili)
pbc_clean$time <- pbc_clean$time / 365.25 # Convert to years

# Pre-allocate storage
estimates <- matrix(NA, nrow = nsim, ncol = 4)
variances <- matrix(NA, nrow = nsim, ncol = 4)

#set.seed(2026) # Ensure reproducibility

# ==========================================
# 2. Simulation loop
# ==========================================
for (s in 1:nsim) {
  
  # Step 1: Resample covariates Z from original data [cite: 731]
  idx <- sample(1:nrow(pbc_clean), n, replace = TRUE)
  Z_matrix <- as.matrix(cbind(1, pbc_clean[idx, c("age", "log_bili", "albumin")]))
  
  # Step 2: Generate event time T (Model 1: Weibull shape=1.0, i.e., exponential) [cite: 734]
  # AH = 1/scale. Under exponential distribution, AH = lambda.
  eta_true <- exp(Z_matrix %*% true_beta)
  T_event <- rexp(n, rate = eta_true)
  
  # Step 3: Generate independent censoring time C (Pattern B) [cite: 739]
  C_cens <- rexp(n, rate = 0.1)
  C_cens <- pmin(C_cens, tau) # Administrative censoring [cite: 738]
  
  # Observed data
  X <- pmin(T_event, C_cens)
  Delta <- as.numeric(T_event <= C_cens)
  
  # ------------------------------------------
  # Step 4: Compute pseudo-observations
  # ------------------------------------------
  
  # A. Compute F(tau) pseudo-values: F = 1 - S
  # Key fix: use unlist(as.data.frame()) to safely extract numeric values from complex objects
  surv_obj <- pseudosurv(X, Delta, tmax = tau)
  surv_df <- as.data.frame(surv_obj)
  pseudo_F <- 1 - as.numeric(unlist(surv_df[, ncol(surv_df)]))
  
  # B. Compute RMST (R_i) pseudo-values
  mean_obj <- pseudomean(X, Delta, tmax = tau)
  pseudo_R <- as.numeric(unlist(mean_obj))
  
  # ------------------------------------------
  # Step 5: Solve estimating equations [cite: 534]
  # ------------------------------------------
  # U(theta) = sum{ Z * (F - exp(theta'Z) * R) } = 0
  estimating_eq <- function(beta) {
    # Add pmin protection to prevent exponential overflow during root finding
    linear_pred <- exp(pmin(Z_matrix %*% beta, 20))
    
    residual <- pseudo_F - (linear_pred * pseudo_R)
    return(colMeans(Z_matrix * as.vector(residual)))
  }
  
  # Root finding. If starting from zero fails, try starting from true_beta [cite: 534]
  fit <- try(multiroot(estimating_eq, start = true_beta, maxiter = 100), silent = TRUE)
  
  if(inherits(fit, "try-error")) next
  
  beta_hat <- as.numeric(fit$root)
  estimates[s, ] <- beta_hat
  
  # ------------------------------------------
  # Step 6: Sandwich variance estimator
  # ------------------------------------------
  linear_pred_hat <- exp(Z_matrix %*% beta_hat)
  residual_hat <- pseudo_F - (linear_pred_hat * pseudo_R)
  
  # A. Empirical sensitivity matrix (Gamma) [cite: 199]
  Gamma_hat <- matrix(0, 4, 4)
  for(i in 1:n) {
    # Derivative term: -Z * Z' * exp(beta'Z) * R
    Gamma_hat <- Gamma_hat - (Z_matrix[i, ] %*% t(Z_matrix[i, ])) *
      as.numeric(linear_pred_hat[i] * pseudo_R[i])
  }
  Gamma_hat <- Gamma_hat / n
  
  # B. Empirical score variance (B) [cite: 200]
  B_hat <- matrix(0, 4, 4)
  for(i in 1:n) {
    U_i <- Z_matrix[i, ] * as.numeric(residual_hat[i])
    B_hat <- B_hat + (U_i %*% t(U_i))
  }
  B_hat <- B_hat / n
  
  # C. Compute Sigma = Gamma_inv * B * Gamma_inv
  Gamma_inv <- solve(Gamma_hat)
  Sigma_hat <- (Gamma_inv %*% B_hat %*% Gamma_inv) / n
  variances[s, ] <- diag(Sigma_hat)
}

# ==========================================
# 3. Results summary and evaluation
# ==========================================

# Count successful convergences
valid_idx <- which(!is.na(estimates[,1]))
success_count <- length(valid_idx)

if(success_count > 0) {
  
  # Compute mean estimates and relative bias
  mean_est <- colMeans(estimates[valid_idx, , drop=FALSE])
  rel_bias <- (mean_est - true_beta) / true_beta
  
  # Compute 95% confidence interval coverage [cite: 192]
  lower_ci <- estimates[valid_idx, ] - 1.96 * sqrt(variances[valid_idx, ])
  upper_ci <- estimates[valid_idx, ] + 1.96 * sqrt(variances[valid_idx, ])
  
  true_mat <- matrix(rep(true_beta, each = success_count), nrow = success_count)
  coverage <- colMeans(lower_ci <= true_mat & upper_ci >= true_mat)
  
  # Create summary table
  results_table <- data.frame(
    Parameter = c("Intercept", "Age", "Log(bili)", "Albumin"),
    True = true_beta,
    Estimate = round(mean_est, 4),
    Rel_Bias = round(rel_bias, 4),
    Coverage = round(coverage, 4)
  )
  
  cat("\nSimulation completed successfully! Number of converged runs:", success_count, "/", nsim, "\n")
  print(results_table)
  
  latex_table <- kable(results_table, format = "latex", booktabs = TRUE)
  
  tex <- paste0(
    "\\documentclass{article}
\\usepackage{booktabs}
\\begin{document}

", latex_table, "

\\end{document}"
  )
  
  writeLines(tex, "table1.tex")
  
  tinytex::pdflatex("table1.tex")
  
} else {
  stop("Simulation failed: no iterations converged. Please check data scaling or initial values.")
}