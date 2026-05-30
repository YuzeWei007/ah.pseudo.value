# ============================================================
# PO-only / IPCW-only / Stacked GMM simulation
# Based on Uno 2024 AH simulation settings
#
# This is NOT the original Uno 5-method Table 2 reproduction.
# This is our new 3-method output using Uno-style settings:
#
# Methods:
#   1. PO-only
#   2. IPCW-only
#   3. Stacked GMM
#
# Settings:
#   Model 1: Weibull shape = 1.0
#   Model 2: Weibull shape = 2.0
#
# Under each model:
#   (a) no censoring
#   (b) independent censoring
#   (c) group-specific censoring
#   (d) covariate-dependent censoring
#
# Total:
#   2 x 4 = 8 settings
#
# Each setting:
#   nsim = 1000
#   n = 300
#   tau = 7
#
# Outputs:
#   - Long summary table
#   - Relative bias wide table
#   - Coverage wide table
#   - Method comparison table
#   - Setting/censoring check table
#   - Raw results RDS
#   - Checkpoint RDS
#
# Output folder:
#   /Users/JasonWei/Desktop/Spring lab/529
# ============================================================


# ============================================================
# 0. Packages and output folder
# ============================================================

rm(list = ls())

library(survival)
library(pseudo)
library(rootSolve)

target_dir <- "/Users/JasonWei/Desktop/Spring lab/529"
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
setwd(target_dir)

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

cat("\nRun ID:\n")
print(run_id)

cat("\nWorking directory:\n")
print(getwd())


# ============================================================
# 1. Global settings
# ============================================================

# The paper does not report a Monte Carlo seed.
# This seed is only for reproducibility of this run.
set.seed(2026)

# Uno 2024 Section 3.1 / Equation 7 coefficients
# age must be in years.
beta_true <- c(
  Intercept = -1.2340,
  Age       =  0.0387,
  LogBili   =  0.8371,
  Albumin   = -1.1590
)

tau <- 7
n_sample <- 300

# Current project setting.
# Use smaller number only for debugging.
nsim <- 1000

conf_level <- 0.95
z_crit <- qnorm(1 - (1 - conf_level) / 2)

methods <- c("PO-only", "IPCW-only", "Stacked-GMM")

cat("\nGlobal settings:\n")
print(list(
  n_sample = n_sample,
  nsim = nsim,
  tau = tau,
  conf_level = conf_level
))


# ============================================================
# 2. Prepare covariate pool from PBC
# ============================================================
# Uno paper uses 312 randomized PBC participants.
# In survival::pbc, randomized participants have non-missing trt.
# The model uses age in years, log(bilirubin), and albumin.

prepare_pbc_pool <- function() {
  data(pbc, package = "survival")
  
  pbc_clean <- pbc[!is.na(pbc$trt), ]
  pbc_clean <- pbc_clean[
    complete.cases(pbc_clean[, c("age", "bili", "albumin")]),
  ]
  
  # Some versions of survival::pbc store age in days.
  # Convert to years when needed.
  if (max(pbc_clean$age, na.rm = TRUE) > 120) {
    pbc_clean$age <- pbc_clean$age / 365.25
  }
  
  pbc_clean$log_bili <- log(pbc_clean$bili)
  
  cov_raw <- as.matrix(pbc_clean[, c("age", "log_bili", "albumin")])
  
  list(
    data = pbc_clean,
    cov_raw = cov_raw
  )
}

pbc_pool <- prepare_pbc_pool()

cat("\nPBC covariate pool size:\n")
print(nrow(pbc_pool$data))

cat("\nAge range after possible conversion:\n")
print(range(pbc_pool$data$age))


# ============================================================
# 3. Weibull helper
# ============================================================

weibull_AH <- function(scale, shape, tau) {
  S_tau <- exp(- (tau / scale)^shape)
  F_tau <- 1 - S_tau
  
  R_tau <- integrate(
    f = function(u) exp(- (u / scale)^shape),
    lower = 0,
    upper = tau,
    subdivisions = 100L,
    rel.tol = 1e-8
  )$value
  
  F_tau / R_tau
}

solve_weibull_scale <- function(eta, shape, tau) {
  if (!is.finite(eta) || eta <= 0) return(NA_real_)
  
  # For shape = 1, Weibull is exponential and AH equals rate = 1 / scale.
  if (abs(shape - 1) < 1e-12) {
    return(1 / eta)
  }
  
  f <- function(scale) weibull_AH(scale, shape, tau) - eta
  
  out <- try(
    uniroot(f, lower = 1e-6, upper = 1e6, extendInt = "yes", tol = 1e-8),
    silent = TRUE
  )
  
  if (inherits(out, "try-error")) return(NA_real_)
  
  out$root
}


# ============================================================
# 4. Data generation under Uno setting
# ============================================================

gen_data_uno <- function(n, beta, tau, pbc_pool,
                         shape = 1.0,
                         censoring = c(
                           "none",
                           "independent",
                           "group_specific",
                           "covariate_dependent"
                         )) {
  
  censoring <- match.arg(censoring)
  
  idx <- sample(seq_len(nrow(pbc_pool$cov_raw)), n, replace = TRUE)
  
  covariates <- pbc_pool$data[idx, c("age", "bili", "log_bili", "albumin")]
  
  Z <- as.matrix(cbind(
    Intercept = 1,
    Age = covariates$age,
    LogBili = covariates$log_bili,
    Albumin = covariates$albumin
  ))
  
  eta <- as.numeric(exp(pmin(Z %*% beta, 20)))
  
  scale_T <- vapply(
    eta,
    FUN = solve_weibull_scale,
    FUN.VALUE = numeric(1),
    shape = shape,
    tau = tau
  )
  
  if (any(!is.finite(scale_T))) {
    stop("Failed to solve Weibull scale for at least one subject.")
  }
  
  T_event <- rweibull(n, shape = shape, scale = scale_T)
  
  if (censoring == "none") {
    C_potential <- rep(Inf, n)
  }
  
  if (censoring == "independent") {
    # Paper setting:
    # common exponential censoring so that Pr(C <= tau) = 0.50.
    rate_c <- -log(0.50) / tau
    C_potential <- rexp(n, rate = rate_c)
  }
  
  if (censoring == "group_specific") {
    # Paper setting:
    # bilirubin > 1.35: exponential scale 19.61
    # bilirubin <= 1.35: exponential scale 5.83
    high_bili <- covariates$bili > 1.35
    scale_C <- ifelse(high_bili, 19.61, 5.83)
    C_potential <- rexp(n, rate = 1 / scale_C)
  }
  
  if (censoring == "covariate_dependent") {
    # Paper setting:
    # scale = exp{1.4386 + 0.0151 age + 0.2120 log(bilirubin)}
    scale_C <- exp(
      1.4386 +
        0.0151 * covariates$age +
        0.2120 * covariates$log_bili
    )
    C_potential <- rexp(n, rate = 1 / scale_C)
  }
  
  X <- pmin(T_event, C_potential, tau)
  
  Delta <- as.numeric(T_event <= C_potential & T_event <= tau)
  
  gs_group <- as.integer(covariates$bili > 1.35)
  
  list(
    Z = Z,
    X = X,
    Delta = Delta,
    T = T_event,
    C = C_potential,
    covariates = covariates,
    gs_group = gs_group,
    observed_random_censor_rate = mean(Delta == 0 & X < tau),
    potential_censor_rate_at_tau = mean(C_potential <= tau),
    total_censor_rate = mean(Delta == 0)
  )
}


# ============================================================
# 5. Pseudo-observations for PO-only
# ============================================================

get_PO <- function(X, Delta, tau) {
  out <- try({
    S_obj <- pseudosurv(X, Delta, tmax = tau)
    S_df <- as.data.frame(S_obj)
    S_hat <- as.numeric(unlist(S_df[, ncol(S_df)]))
    
    F_hat <- 1 - S_hat
    R_hat <- as.numeric(unlist(pseudomean(X, Delta, tmax = tau)))
    
    list(F = F_hat, R = R_hat)
  }, silent = TRUE)
  
  if (inherits(out, "try-error")) return(NULL)
  if (any(!is.finite(out$F)) || any(!is.finite(out$R))) return(NULL)
  
  out
}


# ============================================================
# 6. IPCW quantities for Uno AH estimating equation
# ============================================================

get_ipcw_ind <- function(X, Delta, tau, G_floor = 1e-8) {
  # Delta = 0 at X = tau is administrative truncation, not random censoring.
  # Only count censoring events before tau for the censoring KM.
  censor_ind <- as.numeric(Delta == 0 & X < tau)
  
  fit <- survfit(Surv(X, censor_ind) ~ 1)
  times_eval <- pmin(X, tau)
  
  Ghat <- summary(fit, times = times_eval, extend = TRUE)$surv
  Ghat[is.na(Ghat)] <- 1
  Ghat <- pmax(Ghat, G_floor)
  
  V <- as.numeric((X <= tau & Delta == 1) | (X >= tau))
  weight <- V / Ghat
  weight[!is.finite(weight)] <- 0
  
  list(
    weight = weight,
    event_tau = as.numeric(X < tau & Delta == 1),
    time_tau = pmin(X, tau)
  )
}

get_ipcw_gs <- function(X, Delta, tau, strata, G_floor = 1e-8) {
  weight <- rep(NA_real_, length(X))
  
  for (g in sort(unique(strata))) {
    idx <- strata == g
    tmp <- get_ipcw_ind(X[idx], Delta[idx], tau, G_floor = G_floor)
    weight[idx] <- tmp$weight
  }
  
  weight[!is.finite(weight)] <- 0
  
  list(
    weight = weight,
    event_tau = as.numeric(X < tau & Delta == 1),
    time_tau = pmin(X, tau)
  )
}

get_ipcw_cox <- function(X, Delta, tau, covariates4cens, G_floor = 1e-8) {
  censoring <- as.numeric(Delta == 0 & X < tau)
  
  # If no random censoring before tau, Cox censoring model cannot be fitted.
  # Weights reduce to 1 in this degenerate case.
  if (sum(censoring) == 0) {
    return(list(
      weight = rep(1, length(X)),
      event_tau = as.numeric(X < tau & Delta == 1),
      time_tau = pmin(X, tau)
    ))
  }
  
  dat <- data.frame(
    X = X,
    censoring = censoring,
    covariates4cens
  )
  
  fit <- try(coxph(Surv(X, censoring) ~ ., data = dat), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  
  bh <- basehaz(fit, centered = FALSE)
  
  times_eval <- pmin(X, tau)
  
  Lam0 <- approx(
    x = bh$time,
    y = bh$hazard,
    xout = times_eval,
    method = "constant",
    rule = 2,
    f = 0
  )$y
  
  Lam0[is.na(Lam0)] <- 0
  
  lp <- as.numeric(predict(fit, newdata = dat, type = "lp"))
  
  Ghat <- exp(-Lam0 * exp(lp))
  Ghat <- pmax(Ghat, G_floor)
  
  V <- as.numeric((X <= tau & Delta == 1) | (X >= tau))
  weight <- V / Ghat
  weight[!is.finite(weight)] <- 0
  
  list(
    weight = weight,
    event_tau = as.numeric(X < tau & Delta == 1),
    time_tau = pmin(X, tau)
  )
}

get_IPCW_AH <- function(dat, tau, method = c("ind", "gs", "cox")) {
  method <- match.arg(method)
  
  if (method == "ind") {
    return(get_ipcw_ind(dat$X, dat$Delta, tau))
  }
  
  if (method == "gs") {
    return(get_ipcw_gs(dat$X, dat$Delta, tau, strata = dat$gs_group))
  }
  
  if (method == "cox") {
    covariates4cens <- data.frame(
      Age = dat$covariates$age,
      LogBili = dat$covariates$log_bili
    )
    return(get_ipcw_cox(dat$X, dat$Delta, tau, covariates4cens))
  }
}


# ============================================================
# 7. Moment functions
# ============================================================

moment_PO <- function(beta, Z, po) {
  mu <- as.numeric(exp(pmin(Z %*% beta, 20)))
  residual <- as.numeric(po$F - mu * po$R)
  Z * matrix(residual, nrow = nrow(Z), ncol = ncol(Z))
}

moment_IPCW_AH <- function(beta, Z, ipcw) {
  mu <- as.numeric(exp(pmin(Z %*% beta, 20)))
  residual <- ipcw$weight * (ipcw$event_tau - mu * ipcw$time_tau)
  Z * matrix(residual, nrow = nrow(Z), ncol = ncol(Z))
}

standardize_moment <- function(M, s = NULL) {
  if (is.null(s)) {
    s <- apply(M, 2, sd)
    s[!is.finite(s) | s < 1e-6] <- 1
  }
  
  list(
    M = sweep(M, 2, s, "/"),
    scale = s
  )
}


# ============================================================
# 8. Safe solvers and sandwich SE
# ============================================================

is_bad_beta <- function(b, limit = 20) {
  any(!is.finite(b)) || any(abs(b) > limit)
}

solve_root_safe <- function(fun, start, limit = 20) {
  out <- try(multiroot(fun, start = start, maxiter = 100), silent = TRUE)
  
  if (inherits(out, "try-error")) return(rep(NA_real_, length(start)))
  
  b <- as.numeric(out$root)
  
  if (length(b) != length(start)) return(rep(NA_real_, length(start)))
  if (is_bad_beta(b, limit = limit)) return(rep(NA_real_, length(start)))
  
  names(b) <- names(start)
  b
}

sandwich_just_identified <- function(M, A) {
  n <- nrow(M)
  
  B <- crossprod(scale(M, center = TRUE, scale = FALSE)) / n
  
  A_inv <- try(solve(A), silent = TRUE)
  if (inherits(A_inv, "try-error")) {
    return(rep(NA_real_, ncol(M)))
  }
  
  V <- A_inv %*% B %*% t(A_inv) / n
  se <- sqrt(diag(V))
  
  se[!is.finite(se)] <- NA_real_
  se
}

fit_PO_with_se <- function(Z, po, start) {
  if (is.null(po)) {
    return(list(beta = rep(NA_real_, ncol(Z)), se = rep(NA_real_, ncol(Z))))
  }
  
  f <- function(b) colMeans(moment_PO(b, Z, po))
  
  b <- solve_root_safe(f, start)
  
  if (is_bad_beta(b)) {
    return(list(beta = rep(NA_real_, ncol(Z)), se = rep(NA_real_, ncol(Z))))
  }
  
  mu <- as.numeric(exp(pmin(Z %*% b, 20)))
  
  A <- - t(Z) %*% (Z * as.numeric(mu * po$R)) / nrow(Z)
  
  M <- moment_PO(b, Z, po)
  se <- sandwich_just_identified(M, A)
  
  list(beta = b, se = se)
}

fit_IPCW_with_se <- function(Z, ipcw, start) {
  if (is.null(ipcw)) {
    return(list(beta = rep(NA_real_, ncol(Z)), se = rep(NA_real_, ncol(Z))))
  }
  
  f <- function(b) colMeans(moment_IPCW_AH(b, Z, ipcw))
  
  b <- solve_root_safe(f, start)
  
  if (is_bad_beta(b)) {
    return(list(beta = rep(NA_real_, ncol(Z)), se = rep(NA_real_, ncol(Z))))
  }
  
  mu <- as.numeric(exp(pmin(Z %*% b, 20)))
  
  A <- - t(Z) %*% (Z * as.numeric(ipcw$weight * mu * ipcw$time_tau)) / nrow(Z)
  
  M <- moment_IPCW_AH(b, Z, ipcw)
  se <- sandwich_just_identified(M, A)
  
  list(beta = b, se = se)
}


# ============================================================
# 9. Centering helpers for stable GMM
# ============================================================

get_center_from_pbc_pool <- function(pbc_pool_obj) {
  as.numeric(colMeans(
    pbc_pool_obj$data[, c("age", "log_bili", "albumin")],
    na.rm = TRUE
  ))
}

make_Z_centered <- function(Z_raw, center_vec) {
  Zc <- Z_raw
  Zc[, 2:4] <- sweep(Z_raw[, 2:4, drop = FALSE], 2, center_vec, "-")
  colnames(Zc) <- colnames(Z_raw)
  Zc
}

raw_beta_to_centered <- function(beta_raw, center_vec) {
  beta_c <- beta_raw
  beta_c[1] <- beta_raw[1] + sum(center_vec * beta_raw[2:4])
  beta_c[2:4] <- beta_raw[2:4]
  names(beta_c) <- names(beta_raw)
  beta_c
}

centered_beta_to_raw <- function(beta_c, center_vec) {
  beta_raw <- beta_c
  beta_raw[1] <- beta_c[1] - sum(center_vec * beta_c[2:4])
  beta_raw[2:4] <- beta_c[2:4]
  names(beta_raw) <- names(beta_c)
  beta_raw
}

centered_var_to_raw <- function(V_c, center_vec) {
  p <- nrow(V_c)
  
  L <- diag(p)
  L[1, 2:4] <- -center_vec
  
  L %*% V_c %*% t(L)
}

numerical_jacobian <- function(fun, b, eps = 1e-5) {
  g0 <- fun(b)
  q <- length(g0)
  p <- length(b)
  
  J <- matrix(NA_real_, nrow = q, ncol = p)
  
  for (j in seq_len(p)) {
    step <- eps * max(abs(b[j]), 1)
    
    b_plus <- b
    b_minus <- b
    
    b_plus[j] <- b_plus[j] + step
    b_minus[j] <- b_minus[j] - step
    
    g_plus <- fun(b_plus)
    g_minus <- fun(b_minus)
    
    J[, j] <- (g_plus - g_minus) / (2 * step)
  }
  
  J
}


# ============================================================
# 10. Stacked GMM estimator with sandwich SE
# ============================================================

fit_GMM_with_se <- function(Z, po, ipcw, start,
                            center_vec = NULL,
                            use_standardization = TRUE,
                            beta_limit = 20) {
  
  p <- length(start)
  
  if (is.null(po) || is.null(ipcw)) {
    return(list(beta = rep(NA_real_, p), se = rep(NA_real_, p)))
  }
  
  if (is.null(center_vec)) {
    center_vec <- colMeans(Z[, 2:4, drop = FALSE], na.rm = TRUE)
  }
  
  Zc <- make_Z_centered(Z, center_vec)
  start_c <- raw_beta_to_centered(start, center_vec)
  
  build_M <- function(b, scale_use = NULL) {
    M_po <- moment_PO(b, Zc, po)
    M_ip <- moment_IPCW_AH(b, Zc, ipcw)
    M <- cbind(M_po, M_ip)
    
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
  
  Q <- function(b) {
    if (length(b) != length(start_c)) return(1e20)
    if (is_bad_beta(b, limit = beta_limit)) return(1e20)
    
    M <- build_M(b, scale0)$M
    g <- colMeans(M)
    
    if (any(!is.finite(g))) return(1e20)
    
    sum(g^2)
  }
  
  fit <- try(
    optim(
      par = start_c,
      fn = Q,
      method = "BFGS",
      control = list(maxit = 1000, reltol = 1e-10)
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(list(beta = rep(NA_real_, p), se = rep(NA_real_, p)))
  }
  
  b_c <- as.numeric(fit$par)
  names(b_c) <- names(start)
  
  if (is_bad_beta(b_c, limit = beta_limit)) {
    return(list(beta = rep(NA_real_, p), se = rep(NA_real_, p)))
  }
  
  b_raw <- centered_beta_to_raw(b_c, center_vec)
  
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
  
  S <- crossprod(scale(M_hat, center = TRUE, scale = FALSE)) / n
  
  bread <- try(solve(t(D) %*% D), silent = TRUE)
  
  if (inherits(bread, "try-error")) {
    return(list(beta = b_raw, se = rep(NA_real_, p)))
  }
  
  V_c <- bread %*% t(D) %*% S %*% D %*% bread / n
  
  V_raw <- centered_var_to_raw(V_c, center_vec)
  
  se_raw <- sqrt(diag(V_raw))
  se_raw[!is.finite(se_raw)] <- NA_real_
  
  list(beta = b_raw, se = se_raw)
}


# ============================================================
# 11. One simulation replicate
# ============================================================

one_run <- function(n, beta, tau, pbc_pool,
                    shape = 1.0,
                    censoring = "independent",
                    ipcw_method = "ind",
                    run_gmm = TRUE,
                    center_vec = NULL) {
  
  dat <- gen_data_uno(
    n = n,
    beta = beta,
    tau = tau,
    pbc_pool = pbc_pool,
    shape = shape,
    censoring = censoring
  )
  
  po <- get_PO(dat$X, dat$Delta, tau)
  ipcw <- get_IPCW_AH(dat, tau, method = ipcw_method)
  
  b_po <- fit_PO_with_se(dat$Z, po, beta)
  b_ipcw <- fit_IPCW_with_se(dat$Z, ipcw, beta)
  
  if (run_gmm) {
    b_gmm <- fit_GMM_with_se(dat$Z, po, ipcw, beta, center_vec = center_vec)
  } else {
    b_gmm <- list(
      beta = rep(NA_real_, length(beta)),
      se = rep(NA_real_, length(beta))
    )
  }
  
  list(
    beta = list(
      PO = b_po$beta,
      IPCW = b_ipcw$beta,
      GMM = b_gmm$beta
    ),
    se = list(
      PO = b_po$se,
      IPCW = b_ipcw$se,
      GMM = b_gmm$se
    ),
    observed_random_censor_rate = dat$observed_random_censor_rate,
    potential_censor_rate_at_tau = dat$potential_censor_rate_at_tau,
    total_censor_rate = dat$total_censor_rate
  )
}


# ============================================================
# 12. Simulation loop for one setting
# ============================================================

run_sim <- function(nsim, n, beta, tau, pbc_pool,
                    shape = 1.0,
                    censoring = "independent",
                    ipcw_method = "ind",
                    run_gmm = TRUE) {
  
  p <- length(beta)
  center_vec <- get_center_from_pbc_pool(pbc_pool)
  
  res <- list(
    beta = list(
      PO = matrix(NA_real_, nsim, p),
      IPCW = matrix(NA_real_, nsim, p),
      GMM = matrix(NA_real_, nsim, p)
    ),
    se = list(
      PO = matrix(NA_real_, nsim, p),
      IPCW = matrix(NA_real_, nsim, p),
      GMM = matrix(NA_real_, nsim, p)
    ),
    observed_random_censor_rate = rep(NA_real_, nsim),
    potential_censor_rate_at_tau = rep(NA_real_, nsim),
    total_censor_rate = rep(NA_real_, nsim)
  )
  
  for (nm in names(res$beta)) {
    colnames(res$beta[[nm]]) <- names(beta)
    colnames(res$se[[nm]]) <- names(beta)
  }
  
  for (s in seq_len(nsim)) {
    cat("sim", s, "of", nsim, "\n")
    
    out <- try(
      one_run(
        n = n,
        beta = beta,
        tau = tau,
        pbc_pool = pbc_pool,
        shape = shape,
        censoring = censoring,
        ipcw_method = ipcw_method,
        run_gmm = run_gmm,
        center_vec = center_vec
      ),
      silent = TRUE
    )
    
    if (inherits(out, "try-error")) {
      next
    }
    
    res$beta$PO[s, ] <- out$beta$PO
    res$beta$IPCW[s, ] <- out$beta$IPCW
    res$beta$GMM[s, ] <- out$beta$GMM
    
    res$se$PO[s, ] <- out$se$PO
    res$se$IPCW[s, ] <- out$se$IPCW
    res$se$GMM[s, ] <- out$se$GMM
    
    res$observed_random_censor_rate[s] <- out$observed_random_censor_rate
    res$potential_censor_rate_at_tau[s] <- out$potential_censor_rate_at_tau
    res$total_censor_rate[s] <- out$total_censor_rate
  }
  
  list(
    beta_true = beta,
    res = res,
    settings = list(
      n = n,
      tau = tau,
      shape = shape,
      censoring = censoring,
      ipcw_method = ipcw_method,
      nsim = nsim
    )
  )
}


# ============================================================
# 13. Summary functions
# ============================================================

summ_one <- function(beta_mat, se_mat, beta_true, estimator, setting_row) {
  parameter_names <- names(beta_true)
  
  out <- data.frame()
  
  for (j in seq_along(beta_true)) {
    bhat <- beta_mat[, j]
    se <- se_mat[, j]
    
    valid <- is.finite(bhat) & is.finite(se)
    
    if (sum(valid) == 0) {
      row <- data.frame(
        Model = setting_row$Model,
        Weibull_shape = setting_row$Weibull_shape,
        Censoring = setting_row$Censoring,
        Censoring_label = setting_row$Censoring_label,
        IPCW_method = setting_row$IPCW_method,
        nsim = setting_row$nsim,
        Estimator = estimator,
        Parameter = parameter_names[j],
        True = as.numeric(beta_true[j]),
        Estimate = NA_real_,
        Bias = NA_real_,
        Rel_Bias = NA_real_,
        Empirical_SD = NA_real_,
        Mean_SE = NA_real_,
        Coverage = NA_real_,
        Success = 0,
        stringsAsFactors = FALSE
      )
    } else {
      lower <- bhat[valid] - z_crit * se[valid]
      upper <- bhat[valid] + z_crit * se[valid]
      
      row <- data.frame(
        Model = setting_row$Model,
        Weibull_shape = setting_row$Weibull_shape,
        Censoring = setting_row$Censoring,
        Censoring_label = setting_row$Censoring_label,
        IPCW_method = setting_row$IPCW_method,
        nsim = setting_row$nsim,
        Estimator = estimator,
        Parameter = parameter_names[j],
        True = as.numeric(beta_true[j]),
        Estimate = mean(bhat[valid]),
        Bias = mean(bhat[valid]) - beta_true[j],
        Rel_Bias = (mean(bhat[valid]) - beta_true[j]) / beta_true[j],
        Empirical_SD = sd(bhat[valid]),
        Mean_SE = mean(se[valid]),
        Coverage = mean(lower <= beta_true[j] & upper >= beta_true[j]),
        Success = sum(valid),
        stringsAsFactors = FALSE
      )
    }
    
    out <- rbind(out, row)
  }
  
  out
}

summ_all <- function(res, setting_row) {
  rbind(
    summ_one(res$res$beta$PO, res$res$se$PO, res$beta_true, "PO-only", setting_row),
    summ_one(res$res$beta$IPCW, res$res$se$IPCW, res$beta_true, "IPCW-only", setting_row),
    summ_one(res$res$beta$GMM, res$res$se$GMM, res$beta_true, "Stacked-GMM", setting_row)
  )
}

make_wide_table <- function(summary_table, value_col) {
  tmp <- summary_table[, c(
    "Model",
    "Weibull_shape",
    "Censoring_label",
    "Parameter",
    "Estimator",
    value_col
  )]
  
  wide <- reshape(
    tmp,
    idvar = c("Model", "Weibull_shape", "Censoring_label", "Parameter"),
    timevar = "Estimator",
    direction = "wide"
  )
  
  names(wide) <- gsub(paste0(value_col, "\\."), "", names(wide))
  
  desired_cols <- c(
    "Model",
    "Weibull_shape",
    "Censoring_label",
    "Parameter",
    "PO-only",
    "IPCW-only",
    "Stacked-GMM"
  )
  
  wide <- wide[, desired_cols]
  
  model_order <- c("Model 1", "Model 2")
  censor_order <- c(
    "(a) no censoring",
    "(b) independent censoring",
    "(c) group-specific censoring",
    "(d) covariate-dependent censoring"
  )
  param_order <- names(beta_true)
  
  wide <- wide[
    order(
      match(wide$Model, model_order),
      match(wide$Censoring_label, censor_order),
      match(wide$Parameter, param_order)
    ),
  ]
  
  rownames(wide) <- NULL
  wide
}

make_comparison_table <- function(summary_table) {
  tmp <- summary_table
  
  tmp$Abs_Rel_Bias <- abs(tmp$Rel_Bias)
  tmp$Coverage_Distance <- abs(tmp$Coverage - conf_level)
  
  key_cols <- c("Model", "Weibull_shape", "Censoring_label", "Parameter")
  
  wide_bias <- reshape(
    tmp[, c(key_cols, "Estimator", "Abs_Rel_Bias")],
    idvar = key_cols,
    timevar = "Estimator",
    direction = "wide"
  )
  names(wide_bias) <- gsub("Abs_Rel_Bias\\.", "AbsRelBias_", names(wide_bias))
  
  wide_cov <- reshape(
    tmp[, c(key_cols, "Estimator", "Coverage_Distance")],
    idvar = key_cols,
    timevar = "Estimator",
    direction = "wide"
  )
  names(wide_cov) <- gsub("Coverage_Distance\\.", "CovDist_", names(wide_cov))
  
  out <- merge(wide_bias, wide_cov, by = key_cols)
  
  out$GMM_better_than_IPCW_bias <- out$`AbsRelBias_Stacked-GMM` < out$`AbsRelBias_IPCW-only`
  out$GMM_better_than_IPCW_coverage <- out$`CovDist_Stacked-GMM` < out$`CovDist_IPCW-only`
  
  out
}


# ============================================================
# 14. Full setting grid: 8 settings
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
  IPCW_method = rep(c("ind", "ind", "gs", "cox"), times = 2),
  nsim = nsim,
  stringsAsFactors = FALSE
)

cat("\nFull setting grid:\n")
print(setting_grid)


# ============================================================
# 15. Run all settings
# ============================================================

cat("\n############################################################\n")
cat("START PO / IPCW / GMM 8-SETTING SIMULATION\n")
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
      ipcw_method = this_setting$IPCW_method,
      run_gmm = TRUE
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
    paste0("po_ipcw_gmm_uno8_summary_partial_", run_id, ".csv"),
    row.names = FALSE
  )
  
  write.csv(
    partial_times,
    paste0("po_ipcw_gmm_uno8_setting_times_partial_", run_id, ".csv"),
    row.names = FALSE
  )
  
  saveRDS(
    all_results,
    paste0("po_ipcw_gmm_uno8_checkpoint_", run_id, ".rds")
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
# 16. Combine and save final results
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
cat("FINISHED PO / IPCW / GMM 8-SETTING SIMULATION\n")
cat("############################################################\n")

cat("\nStart time:\n")
print(full_start_time)

cat("\nEnd time:\n")
print(full_end_time)

cat("\nTotal elapsed time:\n")
print(full_end_time - full_start_time)

cat("\nFinal summary dimension:\n")
print(dim(final_summary_rounded))
cat("Expected: 8 settings x 3 estimators x 4 parameters = 96 rows\n")

cat("\nWide table dimensions:\n")
cat("Relative bias wide:\n")
print(dim(rel_bias_wide))
cat("Expected: 8 settings x 4 parameters = 32 rows\n")

cat("\nCensoring check:\n")
print(censoring_check)

cat("\nIndependent censoring check:\n")
cat("Potential censoring rate at tau should be close to 0.50 for independent censoring.\n")
print(censoring_check[censoring_check$Censoring == "independent", ])

cat("\nRelative bias wide table:\n")
print(rel_bias_wide)

cat("\nCoverage wide table:\n")
print(coverage_wide)

cat("\nGMM vs IPCW comparison table:\n")
print(comparison_table)


# ============================================================
# 17. Save final files
# ============================================================

write.csv(
  final_summary_rounded,
  paste0("po_ipcw_gmm_uno8_summary_long_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  final_setting_times,
  paste0("po_ipcw_gmm_uno8_setting_times_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  rel_bias_wide,
  paste0("po_ipcw_gmm_uno8_relative_bias_wide_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  coverage_wide,
  paste0("po_ipcw_gmm_uno8_coverage_wide_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  estimate_wide,
  paste0("po_ipcw_gmm_uno8_estimate_wide_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  sd_wide,
  paste0("po_ipcw_gmm_uno8_empirical_sd_wide_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  mean_se_wide,
  paste0("po_ipcw_gmm_uno8_mean_se_wide_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  comparison_table,
  paste0("po_ipcw_gmm_uno8_gmm_vs_ipcw_comparison_", run_id, ".csv"),
  row.names = FALSE
)

write.csv(
  censoring_check,
  paste0("po_ipcw_gmm_uno8_censoring_check_", run_id, ".csv"),
  row.names = FALSE
)

saveRDS(
  all_results,
  paste0("po_ipcw_gmm_uno8_raw_results_", run_id, ".rds")
)

cat("\nSaved final files:\n")
cat("  ", paste0("po_ipcw_gmm_uno8_summary_long_", run_id, ".csv\n"), sep = "")
cat("  ", paste0("po_ipcw_gmm_uno8_setting_times_", run_id, ".csv\n"), sep = "")
cat("  ", paste0("po_ipcw_gmm_uno8_relative_bias_wide_", run_id, ".csv\n"), sep = "")
cat("  ", paste0("po_ipcw_gmm_uno8_coverage_wide_", run_id, ".csv\n"), sep = "")
cat("  ", paste0("po_ipcw_gmm_uno8_estimate_wide_", run_id, ".csv\n"), sep = "")
cat("  ", paste0("po_ipcw_gmm_uno8_empirical_sd_wide_", run_id, ".csv\n"), sep = "")
cat("  ", paste0("po_ipcw_gmm_uno8_mean_se_wide_", run_id, ".csv\n"), sep = "")
cat("  ", paste0("po_ipcw_gmm_uno8_gmm_vs_ipcw_comparison_", run_id, ".csv\n"), sep = "")
cat("  ", paste0("po_ipcw_gmm_uno8_censoring_check_", run_id, ".csv\n"), sep = "")
cat("  ", paste0("po_ipcw_gmm_uno8_raw_results_", run_id, ".rds\n"), sep = "")
cat("  ", paste0("po_ipcw_gmm_uno8_checkpoint_", run_id, ".rds\n"), sep = "")

cat("\nFiles saved in:\n")
print(getwd())

cat("\nFiles generated in this run:\n")
print(list.files(pattern = paste0("po_ipcw_gmm_uno8_.*", run_id)))