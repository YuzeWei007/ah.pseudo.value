# ============================================================
# FULL R CODE
# 3 settings for validating:
#   J_i ≈ theta + phi_i
# under single-event right-censored survival data
#
# For each setting, this code does:
#   1) specify T distribution and C distribution
#   2) generate data
#   3) compute true F(tau) and true RMST(tau)
#   4) compute pseudo-observations for F and RMST
#   5) compute theoretical KM-based IF plug-in values
#   6) compare:
#        J_F,i - F(tau)   vs   phi_F,i
#        J_R,i - R(tau)   vs   phi_R,i
#   7) summarize correlation, MSE, intercept, slope
#
# Install once if needed:
# install.packages(c("survival", "pseudo"))
# ============================================================

library(survival)
library(pseudo)

set.seed(12345)

# ============================================================
# Global settings
# ============================================================
tau <- 5

n_small <- 200
n_mid   <- 500
n_large <- 5000

nrep_small <- 20
nrep_mid   <- 20
nrep_large <- 20

# ============================================================
# 1. SETTINGS
# ============================================================
# Setting 1: Exponential survival + Exponential censoring
# Setting 2: Weibull survival + Uniform censoring
# Setting 3: Log-logistic survival + Exponential censoring

get_setting_info <- function(setting) {
  
  if (setting == 1) {
    return(list(
      name = "Setting 1: Exponential survival + Exponential censoring",
      surv_name = "Exponential(rate=0.10)",
      cens_name = "Exponential(rate=0.05)",
      lambda = 0.10,
      lambda_c = 0.05
    ))
  }
  
  if (setting == 2) {
    return(list(
      name = "Setting 2: Weibull survival + Uniform censoring",
      surv_name = "Weibull(shape=2, scale=8)",
      cens_name = "Uniform(0, 12)",
      shape = 2,
      scale = 8,
      cmax = 12
    ))
  }
  
  if (setting == 3) {
    return(list(
      name = "Setting 3: Log-logistic survival + Exponential censoring",
      surv_name = "Log-logistic(shape=2, scale=4)",
      cens_name = "Exponential(rate=0.12)",
      shape = 2,
      scale = 4,
      lambda_c = 0.12
    ))
  }
  
  stop("Unknown setting.")
}

# ============================================================
# 2. TRUE SURVIVAL FUNCTION S(t), F(tau), RMST(tau)
# ============================================================

# ---------- Setting 1: Exponential ----------
S_exp <- function(t, lambda) {
  exp(-lambda * t)
}

F_exp <- function(tau, lambda) {
  1 - exp(-lambda * tau)
}

RMST_exp <- function(tau, lambda) {
  (1 - exp(-lambda * tau)) / lambda
}

# ---------- Setting 2: Weibull ----------
# S(t) = exp(-(t/scale)^shape)
S_weib <- function(t, shape, scale) {
  exp(- (t / scale)^shape)
}

F_weib <- function(tau, shape, scale) {
  1 - exp(- (tau / scale)^shape)
}

# RMST = ∫_0^tau exp(-(t/scale)^shape) dt
#      = scale/shape * gamma(1/shape) * pgamma((tau/scale)^shape, shape=1/shape)
RMST_weib <- function(tau, shape, scale) {
  (scale / shape) * gamma(1 / shape) *
    pgamma((tau / scale)^shape, shape = 1 / shape, scale = 1)
}

# ---------- Setting 3: Log-logistic ----------
# S(t) = 1 / (1 + (t/scale)^shape)
S_llogis <- function(t, shape, scale) {
  1 / (1 + (t / scale)^shape)
}

F_llogis <- function(tau, shape, scale) {
  1 - 1 / (1 + (tau / scale)^shape)
}

# For shape = 2:
# RMST = ∫_0^tau 1/(1 + (t/scale)^2) dt = scale * arctan(tau/scale)
RMST_llogis <- function(tau, shape, scale) {
  if (shape != 2) stop("Closed form here is coded for shape = 2.")
  scale * atan(tau / scale)
}

true_S_tau <- function(tau, setting) {
  pars <- get_setting_info(setting)
  
  if (setting == 1) return(S_exp(tau, pars$lambda))
  if (setting == 2) return(S_weib(tau, pars$shape, pars$scale))
  if (setting == 3) return(S_llogis(tau, pars$shape, pars$scale))
}

true_F <- function(tau, setting) {
  pars <- get_setting_info(setting)
  
  if (setting == 1) return(F_exp(tau, pars$lambda))
  if (setting == 2) return(F_weib(tau, pars$shape, pars$scale))
  if (setting == 3) return(F_llogis(tau, pars$shape, pars$scale))
}

true_RMST <- function(tau, setting) {
  pars <- get_setting_info(setting)
  
  if (setting == 1) return(RMST_exp(tau, pars$lambda))
  if (setting == 2) return(RMST_weib(tau, pars$shape, pars$scale))
  if (setting == 3) return(RMST_llogis(tau, pars$shape, pars$scale))
}

# ============================================================
# 3. DATA GENERATION
# ============================================================

rloglogis <- function(n, shape, scale) {
  u <- runif(n)
  scale * (u / (1 - u))^(1 / shape)
}

simulate_data <- function(n, setting) {
  pars <- get_setting_info(setting)
  
  if (setting == 1) {
    T <- rexp(n, rate = pars$lambda)
    C <- rexp(n, rate = pars$lambda_c)
  }
  
  if (setting == 2) {
    T <- rweibull(n, shape = pars$shape, scale = pars$scale)
    C <- runif(n, min = 0, max = pars$cmax)
  }
  
  if (setting == 3) {
    T <- rloglogis(n, shape = pars$shape, scale = pars$scale)
    C <- rexp(n, rate = pars$lambda_c)
  }
  
  Y <- pmin(T, C)
  Delta <- as.integer(T <= C)
  
  data.frame(
    id = seq_len(n),
    T = T,
    C = C,
    Y = Y,
    Delta = Delta
  )
}

# ============================================================
# 4. PSEUDO-OBSERVATIONS
# ============================================================

compute_pseudo <- function(data, tau) {
  ps_obj <- pseudosurv(time = data$Y, event = data$Delta, tmax = tau)
  J_S <- as.numeric(ps_obj$pseudo[, 1])
  J_F <- 1 - J_S
  
  J_R <- as.numeric(pseudomean(time = data$Y, event = data$Delta, tmax = tau))
  
  list(
    J_S = J_S,
    J_F = J_F,
    J_R = J_R
  )
}

# ============================================================
# 5. KM STEP INFORMATION
# ============================================================

km_step_info <- function(Y, Delta, tau) {
  n <- length(Y)
  
  event_times <- sort(unique(Y[Delta == 1 & Y <= tau]))
  m <- length(event_times)
  
  if (m == 0) {
    return(list(
      n = n,
      event_times = numeric(0),
      m = 0,
      risk = numeric(0),
      deaths = numeric(0),
      dA = numeric(0),
      S_after = numeric(0),
      theta_hat_S = 1,
      theta_hat_F = 0,
      theta_hat_R = tau
    ))
  }
  
  risk <- sapply(event_times, function(tt) sum(Y >= tt))
  deaths <- sapply(event_times, function(tt) sum(Y == tt & Delta == 1))
  dA <- deaths / risk
  
  S_after <- cumprod(1 - dA)
  
  theta_hat_S <- S_after[m]
  theta_hat_F <- 1 - theta_hat_S
  
  interval_starts <- c(0, event_times)
  interval_ends   <- c(event_times, tau)
  widths <- interval_ends - interval_starts
  surv_levels <- c(1, S_after)
  
  theta_hat_R <- sum(widths * surv_levels)
  
  list(
    n = n,
    event_times = event_times,
    m = m,
    risk = risk,
    deaths = deaths,
    dA = dA,
    S_after = S_after,
    theta_hat_S = theta_hat_S,
    theta_hat_F = theta_hat_F,
    theta_hat_R = theta_hat_R
  )
}

# ============================================================
# 6. THEORETICAL KM-BASED INFLUENCE FUNCTION PLUG-IN
# ============================================================
# This is the corrected part.
# We compute phi_F,i and phi_R,i from the KM influence representation.

compute_km_if <- function(data, tau) {
  Y <- data$Y
  Delta <- data$Delta
  n <- nrow(data)
  
  info <- km_step_info(Y, Delta, tau)
  
  if (info$m == 0) {
    return(list(
      phi_S = rep(0, n),
      phi_F = rep(0, n),
      phi_R = rep(0, n),
      theta_hat_S = info$theta_hat_S,
      theta_hat_F = info$theta_hat_F,
      theta_hat_R = info$theta_hat_R
    ))
  }
  
  ev <- info$event_times
  m <- info$m
  risk <- info$risk
  dA <- info$dA
  S_after <- info$S_after
  
  # Hhat(u) = P_hat(Y >= u)
  Hhat <- risk / n
  
  # Y_ij = I(Y_i >= t_j)
  Y_ind <- outer(Y, ev, FUN = ">=") * 1
  
  # dN_ij = I(Y_i = t_j, Delta_i = 1)
  dN_ind <- outer(Y, ev, FUN = "==") * outer(Delta, rep(1, m), FUN = "*")
  
  # Increment: (dN_i(t_j) - Y_i(t_j)dA_j) / Hhat_j
  inc_mat <- matrix(0, nrow = n, ncol = m)
  for (j in seq_len(m)) {
    inc_mat[, j] <- (dN_ind[, j] - Y_ind[, j] * dA[j]) / Hhat[j]
  }
  
  # cumulative sum up to each event time
  cum_inc_mat <- t(apply(inc_mat, 1, cumsum))
  
  # phi_S(t_j) = -S_hat(t_j) * cumulative increment up to t_j
  phi_S_mat <- matrix(0, nrow = n, ncol = m)
  for (j in seq_len(m)) {
    phi_S_mat[, j] <- - S_after[j] * cum_inc_mat[, j]
  }
  
  phi_S_tau <- phi_S_mat[, m]
  phi_F_tau <- -phi_S_tau
  
  # RMST IF = ∫_0^tau phi_S(u) du
  widths_after_event <- c(diff(ev), tau - ev[m])
  
  phi_R_tau <- rowSums(
    phi_S_mat * matrix(widths_after_event, nrow = n, ncol = m, byrow = TRUE)
  )
  
  list(
    phi_S = phi_S_tau,
    phi_F = phi_F_tau,
    phi_R = phi_R_tau,
    theta_hat_S = info$theta_hat_S,
    theta_hat_F = info$theta_hat_F,
    theta_hat_R = info$theta_hat_R
  )
}

# ============================================================
# 7. COMPARISON METRICS
# ============================================================

compare_vectors <- function(J, theta, phi) {
  centered <- J - theta
  diff <- centered - phi
  
  fit <- lm(centered ~ phi)
  
  c(
    correlation = suppressWarnings(cor(centered, phi)),
    mse = mean(diff^2),
    intercept = unname(coef(fit)[1]),
    slope = unname(coef(fit)[2])
  )
}

# ============================================================
# 8. PLOT FUNCTION
# ============================================================

plot_compare <- function(phi, centered, main_title = "") {
  plot(
    x = phi,
    y = centered,
    pch = 16,
    xlab = "Theoretical IF plug-in value",
    ylab = "Pseudo-observation minus true value",
    main = main_title
  )
  abline(a = 0, b = 1, col = "red", lwd = 2)
}

# ============================================================
# 9. ONE SIMULATION RUN
# ============================================================

run_once <- function(n, setting, tau, make_plot = FALSE) {
  data <- simulate_data(n = n, setting = setting)
  
  pseudo_obj <- compute_pseudo(data = data, tau = tau)
  if_obj <- compute_km_if(data = data, tau = tau)
  
  theta_true_F <- true_F(tau = tau, setting = setting)
  theta_true_R <- true_RMST(tau = tau, setting = setting)
  
  # This is the key theoretical validation
  res_F_true <- compare_vectors(
    J = pseudo_obj$J_F,
    theta = theta_true_F,
    phi = if_obj$phi_F
  )
  
  res_R_true <- compare_vectors(
    J = pseudo_obj$J_R,
    theta = theta_true_R,
    phi = if_obj$phi_R
  )
  
  # Optional numerical consistency check
  res_F_hat <- compare_vectors(
    J = pseudo_obj$J_F,
    theta = if_obj$theta_hat_F,
    phi = if_obj$phi_F
  )
  
  res_R_hat <- compare_vectors(
    J = pseudo_obj$J_R,
    theta = if_obj$theta_hat_R,
    phi = if_obj$phi_R
  )
  
  censoring_rate <- mean(data$Delta == 0)
  
  if (make_plot) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    
    par(mfrow = c(1, 2))
    plot_compare(
      phi = if_obj$phi_F,
      centered = pseudo_obj$J_F - theta_true_F,
      main_title = paste0("F(tau), setting ", setting, ", n=", n)
    )
    plot_compare(
      phi = if_obj$phi_R,
      centered = pseudo_obj$J_R - theta_true_R,
      main_title = paste0("RMST, setting ", setting, ", n=", n)
    )
  }
  
  list(
    data = data,
    censoring_rate = censoring_rate,
    theta_true_F = theta_true_F,
    theta_true_R = theta_true_R,
    theta_hat_F = if_obj$theta_hat_F,
    theta_hat_R = if_obj$theta_hat_R,
    J_F = pseudo_obj$J_F,
    J_R = pseudo_obj$J_R,
    phi_F = if_obj$phi_F,
    phi_R = if_obj$phi_R,
    result_true_F = res_F_true,
    result_true_R = res_R_true,
    result_hat_F = res_F_hat,
    result_hat_R = res_R_hat
  )
}

# ============================================================
# 10. REPEATED SIMULATION
# ============================================================

run_simulation <- function(n, nrep, setting, tau, plot_first_rep = FALSE) {
  out_true_F <- matrix(NA_real_, nrow = nrep, ncol = 4)
  out_true_R <- matrix(NA_real_, nrow = nrep, ncol = 4)
  out_hat_F  <- matrix(NA_real_, nrow = nrep, ncol = 4)
  out_hat_R  <- matrix(NA_real_, nrow = nrep, ncol = 4)
  censoring_vec <- numeric(nrep)
  
  cn <- c("correlation", "mse", "intercept", "slope")
  colnames(out_true_F) <- cn
  colnames(out_true_R) <- cn
  colnames(out_hat_F)  <- cn
  colnames(out_hat_R)  <- cn
  
  first_run <- NULL
  
  for (i in seq_len(nrep)) {
    ans <- run_once(
      n = n,
      setting = setting,
      tau = tau,
      make_plot = (plot_first_rep && i == 1)
    )
    
    if (i == 1) first_run <- ans
    
    out_true_F[i, ] <- ans$result_true_F
    out_true_R[i, ] <- ans$result_true_R
    out_hat_F[i, ]  <- ans$result_hat_F
    out_hat_R[i, ]  <- ans$result_hat_R
    censoring_vec[i] <- ans$censoring_rate
  }
  
  summarize_mat <- function(M) {
    rbind(
      mean = colMeans(M, na.rm = TRUE),
      sd   = apply(M, 2, sd, na.rm = TRUE)
    )
  }
  
  list(
    n = n,
    nrep = nrep,
    setting = setting,
    censoring_mean = mean(censoring_vec),
    censoring_sd   = sd(censoring_vec),
    true_F = summarize_mat(out_true_F),
    true_R = summarize_mat(out_true_R),
    hat_F  = summarize_mat(out_hat_F),
    hat_R  = summarize_mat(out_hat_R),
    true_F_all = out_true_F,
    true_R_all = out_true_R,
    hat_F_all  = out_hat_F,
    hat_R_all  = out_hat_R,
    first_run = first_run
  )
}

# ============================================================
# 11. PRINT HELPERS
# ============================================================

print_setting_header <- function(setting) {
  pars <- get_setting_info(setting)
  cat("\n========================================================\n")
  cat(pars$name, "\n")
  cat("Survival:", pars$surv_name, "\n")
  cat("Censoring:", pars$cens_name, "\n")
  cat("tau =", tau, "\n")
  cat("True F(tau)   =", round(true_F(tau, setting), 6), "\n")
  cat("True RMST(tau)=", round(true_RMST(tau, setting), 6), "\n")
  cat("========================================================\n")
}

print_simulation_summary <- function(obj) {
  cat("\n--------------------------------------------\n")
  cat("Setting =", obj$setting, " | n =", obj$n, " | nrep =", obj$nrep, "\n")
  cat("Mean censoring rate =", round(obj$censoring_mean, 4),
      " (sd =", round(obj$censoring_sd, 4), ")\n")
  cat("--------------------------------------------\n")
  
  cat("\nTHEORETICAL VALIDATION: J_i - theta_true vs phi_hat\n")
  cat("\nF(tau):\n")
  print(round(obj$true_F, 4))
  cat("\nRMST(tau):\n")
  print(round(obj$true_R, 4))
  
  cat("\nNUMERICAL CHECK: J_i - theta_hat vs phi_hat\n")
  cat("\nF(tau):\n")
  print(round(obj$hat_F, 4))
  cat("\nRMST(tau):\n")
  print(round(obj$hat_R, 4))
}

# ============================================================
# 12. RUN ALL SAMPLE SIZES FOR ONE SETTING
# ============================================================

run_one_setting_all_n <- function(setting) {
  print_setting_header(setting)
  
  res_small <- run_simulation(
    n = n_small, nrep = nrep_small,
    setting = setting, tau = tau,
    plot_first_rep = TRUE
  )
  
  res_mid <- run_simulation(
    n = n_mid, nrep = nrep_mid,
    setting = setting, tau = tau,
    plot_first_rep = TRUE
  )
  
  res_large <- run_simulation(
    n = n_large, nrep = nrep_large,
    setting = setting, tau = tau,
    plot_first_rep = TRUE
  )
  
  print_simulation_summary(res_small)
  print_simulation_summary(res_mid)
  print_simulation_summary(res_large)
  
  list(
    small = res_small,
    mid = res_mid,
    large = res_large
  )
}

# ============================================================
# 13. RUN ALL 3 SETTINGS
# ============================================================

all_res_1 <- run_one_setting_all_n(setting = 1)
all_res_2 <- run_one_setting_all_n(setting = 2)
all_res_3 <- run_one_setting_all_n(setting = 3)

# ============================================================
# 14. PAPER-STYLE SUMMARY TABLES
# ============================================================

extract_mean_row <- function(obj, which_block = c("true_F", "true_R", "hat_F", "hat_R")) {
  which_block <- match.arg(which_block)
  obj[[which_block]]["mean", ]
}

make_summary_table <- function(setting_res, which_block) {
  rbind(
    small_mean = extract_mean_row(setting_res$small, which_block),
    mid_mean   = extract_mean_row(setting_res$mid, which_block),
    large_mean = extract_mean_row(setting_res$large, which_block)
  )
}

cat("\n========================================================\n")
cat("SETTING 1: TRUE-CENTERED F\n")
print(round(make_summary_table(all_res_1, "true_F"), 4))

cat("\nSETTING 1: TRUE-CENTERED RMST\n")
print(round(make_summary_table(all_res_1, "true_R"), 4))

cat("\nSETTING 2: TRUE-CENTERED F\n")
print(round(make_summary_table(all_res_2, "true_F"), 4))

cat("\nSETTING 2: TRUE-CENTERED RMST\n")
print(round(make_summary_table(all_res_2, "true_R"), 4))

cat("\nSETTING 3: TRUE-CENTERED F\n")
print(round(make_summary_table(all_res_3, "true_F"), 4))

cat("\nSETTING 3: TRUE-CENTERED RMST\n")
print(round(make_summary_table(all_res_3, "true_R"), 4))
cat("\n========================================================\n")

############ Make plot

library(knitr)

make_pdf <- function(df, filename, title) {
  
  tex_file <- paste0(filename, ".tex")
  
  
  latex_table <- kable(
    df,
    format = "latex",
    booktabs = TRUE,
    digits = 4,
    align = "lcccc",
    caption = title
  )
  
  tex <- paste0(
    "\\documentclass{article}
\\usepackage{booktabs}
\\usepackage[margin=1in]{geometry}
\\begin{document}

\\begin{center}
\\large ", title, "
\\end{center}

", latex_table, "

\\end{document}"
  )
  
  writeLines(tex, tex_file)
  
  tinytex::pdflatex(tex_file)
  
  unlink(c(
    tex_file,
    paste0(filename, ".log"),
    paste0(filename, ".aux"),
    paste0(filename, ".out")
  ))
}


make_pdf(table_s1_F, "simulation2_table1", "Setting 1: F(tau)")
make_pdf(table_s1_R, "simulation2_table2", "Setting 1: RMST")

make_pdf(table_s2_F, "simulation2_table3", "Setting 2: F(tau)")
make_pdf(table_s2_R, "simulation2_table4", "Setting 2: RMST")

make_pdf(table_s3_F, "simulation2_table5", "Setting 3: F(tau)")
make_pdf(table_s3_R, "simulation2_table6", "Setting 3: RMST")




library(knitr)


s1_F <- round(make_summary_table(all_res_1, "true_F"), 4)
s1_R <- round(make_summary_table(all_res_1, "true_R"), 4)

s2_F <- round(make_summary_table(all_res_2, "true_F"), 4)
s2_R <- round(make_summary_table(all_res_2, "true_R"), 4)

s3_F <- round(make_summary_table(all_res_3, "true_F"), 4)
s3_R <- round(make_summary_table(all_res_3, "true_R"), 4)


make_all_pdf <- function() {
  
  tex <- paste0(
    "\\documentclass{article}
\\usepackage{booktabs}
\\usepackage[margin=1in]{geometry}
\\begin{document}

\\section*{Simulation 2 Results}

\\subsection*{Setting 1: TRUE-CENTERED F}
", kable(s1_F, "latex", booktabs=TRUE), "

\\subsection*{Setting 1: TRUE-CENTERED RMST}
", kable(s1_R, "latex", booktabs=TRUE), "

\\subsection*{Setting 2: TRUE-CENTERED F}
", kable(s2_F, "latex", booktabs=TRUE), "

\\subsection*{Setting 2: TRUE-CENTERED RMST}
", kable(s2_R, "latex", booktabs=TRUE), "

\\subsection*{Setting 3: TRUE-CENTERED F}
", kable(s3_F, "latex", booktabs=TRUE), "

\\subsection*{Setting 3: TRUE-CENTERED RMST}
", kable(s3_R, "latex", booktabs=TRUE), "

\\end{document}"
  )
  
  writeLines(tex, "simulation2_tables.tex")
  tinytex::pdflatex("simulation2_tables.tex")
  
  unlink(c("simulation2_tables.tex",
           "simulation2_tables.log",
           "simulation2_tables.aux"))
}


make_all_pdf()


save_plot_pdf <- function(n, setting, tau, filename) {
  
  pdf(paste0(filename, ".pdf"), width = 10, height = 5)
  
  run_once(
    n = n,
    setting = setting,
    tau = tau,
    make_plot = TRUE
  )
  
  dev.off()
}

save_plot_pdf(200, 1, tau, "simulation2_fig1")
save_plot_pdf(200, 2, tau, "simulation2_fig2")
save_plot_pdf(200, 3, tau, "simulation2_fig3")

save_plot_pdf(500, 1, tau, "simulation2_fig4")
save_plot_pdf(500, 2, tau, "simulation2_fig5")
save_plot_pdf(500, 3, tau, "simulation2_fig6")



full_output <- capture.output({
  
  # ===== Setting 1 =====
  print_setting_header(1)
  print_simulation_summary(all_res_1$small)
  print_simulation_summary(all_res_1$mid)
  print_simulation_summary(all_res_1$large)
  
  # ===== Setting 2 =====
  print_setting_header(2)
  print_simulation_summary(all_res_2$small)
  print_simulation_summary(all_res_2$mid)
  print_simulation_summary(all_res_2$large)
  
  # ===== Setting 3 =====
  print_setting_header(3)
  print_simulation_summary(all_res_3$small)
  print_simulation_summary(all_res_3$mid)
  print_simulation_summary(all_res_3$large)
  
})

tex <- paste0(
  "\\documentclass{article}
\\usepackage[margin=1in]{geometry}
\\begin{document}

\\section*{Simulation 2 Full Output}

\\begin{verbatim}
",
  paste(full_output, collapse = "\n"),
  "
\\end{verbatim}

\\end{document}"
)

writeLines(tex, "simulation2_full.tex")
tinytex::pdflatex("simulation2_full.tex")

unlink(c(
  "simulation2_full.tex",
  "simulation2_full.log",
  "simulation2_full.aux"
))