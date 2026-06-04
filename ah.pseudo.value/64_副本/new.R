# ============================================================
# Standalone PO-only / AH-Cox-only / Stacked-GMM-AH-Cox simulation
#
# No source("fun-ahreg.R") needed.
# Embedded functions:
#   kmcens2()
#   cox1()
#   vtm()
#   ahreg()
#
# Main methods:
#   1. PO-only
#   2. AH-Cox-only
#   3. Stacked-GMM-AH-Cox
#
# All IPCW components use AH-Cox via embedded fun-ahreg.R logic.
# ============================================================

rm(list = ls())

library(survival)
library(pseudo)
library(rootSolve)

target_dir <- "/Users/JasonWei/Desktop/Spring lab/529"
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
setwd(target_dir)

# ============================================================
# Embedded fun-ahreg.R functions
# ============================================================

kmcens2 <- function(time, status, tau = NULL) {
  if (is.null(tau)) tau <- max(time)
  
  distinct <- unique(sort(c(time, tau)))
  t <- length(distinct)
  n <- length(time)
  
  surv <- rep(0, t)
  nel.wk <- rep(0, t)
  nelson <- rep(0, t)
  
  yi <- sum(as.numeric(time >= distinct[1]))
  di <- sum(as.numeric(time == distinct[1] & status == 0))
  surv[1] <- 1 * (1 - di / yi)
  nel.wk[1] <- di / yi
  
  for (i in 2:t) {
    yi <- sum(as.numeric(time >= distinct[i]))
    di <- sum(as.numeric(time == distinct[i] & status == 0))
    surv[i] <- surv[i - 1] * (1 - di / yi)
    nel.wk[i] <- di / yi
  }
  
  surv[2:t] <- surv[1:(t - 1)]
  surv[1] <- 1
  nel.wk[2:t] <- nel.wk[1:(t - 1)]
  nel.wk[1] <- 0
  
  nelson <- cumsum(nel.wk)
  
  pi_0 <- rep(0, t)
  pi_X <- rep(0, t)
  pi_T <- rep(0, t)
  Fn <- rep(0, t)
  
  for (i in 1:t) {
    yi <- as.numeric(time >= distinct[i])
    ni <- as.numeric(time <= distinct[i] & status == 1)
    pi_0[i] <- mean(yi)
    pi_X[i] <- mean(yi) * surv[i]
    pi_T[i] <- mean(yi) / surv[i]
    Fn[i] <- mean(ni)
  }
  
  Mi <- matrix(0, n, t)
  wk1 <- matrix(0, n, t)
  wk2 <- matrix(0, n, t)
  
  for (j in 1:t) {
    wk1[, j] <- as.numeric(time <= distinct[j] & status == 0)
    wk2[, j] <- as.numeric(time >= distinct[j])
  }
  
  for (k in 1:n) {
    Mi[k, ] <- wk1[k, ] - wk2[k, ] * nelson
  }
  
  psii <- matrix(0, n, t)
  wk1 <- matrix(0, n, t)
  wk2 <- matrix(0, n, t)
  
  for (j in 1:t) {
    wk1[, j] <- as.numeric(time == distinct[j] & status == 0)
    wk2[, j] <- as.numeric(time >= distinct[j])
  }
  
  for (i in 1:n) {
    psii[i, ] <- (
      cumsum(wk1[i, ] / pi_0) -
        cumsum(nel.wk * wk2[i, ] / pi_0)
    ) * surv
  }
  
  ghat.tau <- surv[distinct == tau]
  
  ghat.x <- rep(0, n)
  for (i in 1:n) {
    ghat.x[i] <- surv[distinct == time[i]]
  }
  ghat.x[is.na(ghat.x)] <- 0
  
  pit <- rep(0, n)
  pit[time >= tau] <- ghat.tau
  pit[time <= tau & status == 1] <- ghat.x[time <= tau & status == 1]
  
  vit <- as.numeric(time <= tau) * status + as.numeric(time >= tau)
  
  ipcw <- rep(0, n)
  ipcw[vit == 1] <- 1 / pit[vit == 1]
  
  list(
    surv = surv,
    nelson = nelson,
    distinct = distinct,
    pi_0 = pi_0,
    pi_X = pi_X,
    pi_T = pi_T,
    Mi = Mi,
    psii = psii,
    Fn = Fn,
    ipcw = ipcw
  )
}

vtm <- function(vc, dm) {
  matrix(vc, ncol = length(vc), nrow = dm, byrow = TRUE)
}

cox1 <- function(time, status, covariates, tau = NULL) {
  ncov <- ncol(as.matrix(covariates))
  
  tmpD <- data.frame(cbind(time, status, covariates))
  vars <- colnames(tmpD)[c(-1, -2)]
  fmla <- as.formula(paste("Surv(time, status) ~", paste(vars, collapse = "+")))
  
  ft <- coxph(fmla, data = tmpD)
  beta <- summary(ft)$coef[, 1]
  
  if (ncov == 1) {
    bz <- as.vector(unlist(covariates)) * rep(beta, length(time))
  } else {
    bz <- covariates %*% beta
  }
  
  ebz <- exp(bz)
  
  if (is.null(tau)) tau <- max(time)
  
  distinct <- unique(sort(c(0, time, tau)))
  t <- length(distinct)
  n <- length(time)
  
  wk1 <- vtm(distinct, n)
  wk2 <- t(vtm(time, t))
  
  dNi <- as.matrix(wk2 == wk1 & t(vtm(status, t)) == 1)
  Ni <- t(sapply(data.frame(t(dNi)), cumsum))
  Yi <- as.matrix(wk2 >= wk1) * 1
  
  Yiebz <- Yi * t(vtm(ebz, t))
  S0bar <- apply(Yiebz, 2, mean)
  
  if (ncov != 1) {
    wk1 <- data.frame(covariates)
    tmp1 <- function(x) {
      apply(Yiebz * t(vtm(x, t)), 2, mean)
    }
    
    S1bar <- t(sapply(wk1, tmp1))
    
    wk2 <- data.frame(t(covariates))
    tmp2 <- function(x) {
      as.vector(x %*% t(x))
    }
    
    wk3 <- data.frame(t(sapply(wk2, tmp2)))
    S2bar <- t(sapply(wk3, tmp1))
    
    S2bar.per.S0bar <- S2bar / vtm(S0bar, ncov * ncov)
    S1bar.per.S0bar <- S1bar / vtm(S0bar, ncov)
    
    wk4 <- data.frame(S1bar.per.S0bar)
    S1bar.per.S0bar.2 <- sapply(wk4, tmp2)
    
    v_beta_t <- S2bar.per.S0bar - S1bar.per.S0bar.2
    
    wk5 <- data.frame(t(v_beta_t))
    tmp5 <- function(x) {
      apply(vtm(x, n) * dNi, 1, sum)
    }
    
    wk6 <- sapply(wk5, tmp5)
    Ibeta1 <- matrix(apply(wk6, 2, mean), nrow = ncov, ncol = ncov)
    InvIbeta1 <- solve(Ibeta1)
  }
  
  if (ncov == 1) {
    S1bar <- apply(
      Yiebz * t(vtm(as.matrix(covariates), t)),
      2,
      mean
    )
    
    S2bar <- apply(
      Yiebz * t(vtm(as.matrix(covariates) * as.matrix(covariates), t)),
      2,
      mean
    )
    
    S2bar.per.S0bar <- S2bar / S0bar
    S1bar.per.S0bar <- S1bar / S0bar
    S1bar.per.S0bar.2 <- S1bar.per.S0bar * S1bar.per.S0bar
    
    v_beta_t <- S2bar.per.S0bar - S1bar.per.S0bar.2
    
    Ibeta1 <- mean(dNi %*% v_beta_t)
    InvIbeta1 <- 1 / Ibeta1
  }
  
  breslow <- cumsum(apply(dNi / vtm(S0bar, n), 2, mean))
  
  dAi <- Yiebz * vtm(diff(c(0, breslow)), n)
  Ai <- cumsum(apply(dAi, 2, mean))
  Mi <- Ni - Ai
  dMi <- dNi - dAi
  
  wk1 <- data.frame(t(vtm(1 / S0bar, n) * dMi))
  wk2 <- sapply(wk1, cumsum)
  eta_lam_i <- t(wk2)
  
  if (ncov != 1) {
    wk7 <- v_beta_t * vtm(S0bar, ncov * ncov) *
      vtm(diff(c(0, breslow)), ncov * ncov)
    
    Ibeta <- matrix(apply(wk7, 1, sum), nrow = ncov, ncol = ncov)
    
    tmp2 <- c()
    for (i in 1:ncov) {
      tmp1 <- apply(
        (t(vtm(covariates[, i], t)) - vtm(S1bar.per.S0bar[i, ], n)) * dMi,
        1,
        sum
      )
      tmp2 <- rbind(tmp2, tmp1)
    }
    
    eta_beta_i <- t(-solve(Ibeta) %*% tmp2)
  }
  
  if (ncov == 1) {
    wk7 <- v_beta_t * S0bar * diff(c(0, breslow))
    Ibeta <- sum(wk7)
    
    tmp1 <- apply(
      (t(vtm(as.matrix(covariates), t)) - vtm(S1bar.per.S0bar, n)) * dMi,
      1,
      sum
    )
    
    eta_beta_i <- -1 / Ibeta * tmp1
  }
  
  event <- 1 - status
  event[time >= tau] <- 2
  
  wk1 <- vtm(distinct, n)
  wk2 <- t(vtm(pmin(time, tau), t))
  
  dNCi <- as.matrix(wk2 == wk1)
  dNCi[dNCi == TRUE] <- 1
  
  wk3 <- vtm(breslow, n)
  Lam0 <- apply(wk3 * dNCi, 1, max)
  
  Ghat <- exp(-Lam0 * ebz)
  
  ipcw <- 1 / Ghat
  ipcw[event == 0] <- 0
  
  Z <- list()
  Z$fit <- ft
  Z$distinct <- distinct
  Z$breslow <- breslow
  Z$ebz <- ebz
  Z$InvIbeta <- InvIbeta1
  Z$ipcw <- ipcw
  Z$Lam0.ebz <- Lam0 * ebz
  Z$eta_lam_i <- eta_lam_i
  Z$eta_beta_i <- eta_beta_i
  
  class(Z) <- "cox1"
  Z
}

ahreg <- function(time,
                  status,
                  covariates,
                  tau = NULL,
                  conf.int = 0.95,
                  strata = NULL,
                  link = "log",
                  covariates4cens = NULL) {
  if (is.null(tau)) print("Please specify tau")
  
  if (!is.null(strata) & !is.null(covariates4cens)) {
    stop("Cannot specify both strata and covariates4cens")
  }
  
  if (is.null(strata)) {
    strata <- rep(1, length(time))
    unique_strata <- 1
    nstrata <- 1
  }
  
  if (!is.null(strata)) {
    unique_strata <- sort(unique(strata))
    nstrata <- length(unique_strata)
  }
  
  ah_rs <- sum(status * as.numeric(time < tau)) / sum(pmin(time, tau))
  b_ini <- c(log(ah_rs), rep(0, ncol(covariates)))
  
  design_mat <- cbind(rep(1, length(time)), covariates)
  
  if (is.null(colnames(design_mat))) {
    colnames(design_mat) <- paste0("X", 1:ncol(design_mat) - 1)
  }
  
  if (!is.null(colnames(design_mat))) {
    colnames(design_mat)[1] <- "Intercept"
  }
  
  if (is.null(covariates4cens)) {
    wt <- rep(0, length(time))
    psii_list <- list()
    
    for (i in 1:nstrata) {
      idx <- strata == unique_strata[i]
      kmc <- kmcens2(time[idx], status[idx], tau = tau)
      wt[idx] <- kmc$ipcw
      psii_list[[i]] <- kmc$psii
    }
  }
  
  if (!is.null(covariates4cens)) {
    censoring <- 1 - status
    censoring[time >= tau] <- 0
    fc.cox <- cox1(time, censoring, covariates4cens, tau = tau)
    wt <- fc.cox$ipcw
  }
  
  hx <- c()
  ii <- 0
  qq <- 99999
  beta <- b_ini
  max_itration <- 20
  criterion <- 1e-15
  convergence <- 9
  
  if (link == "log") {
    while (ii < max_itration & qq > criterion) {
      ii <- ii + 1
      
      gbeta <- exp(design_mat %*% beta)
      yy <- (as.numeric(time < tau) - gbeta * pmin(time, tau)) * wt
      yy_mat <- t(matrix(yy, ncol = length(yy), nrow = ncol(design_mat), byrow = TRUE))
      
      S_beta <- apply(design_mat * yy_mat, 2, mean)
      
      zz <- gbeta * pmin(time, tau) * wt
      zz_mat <- t(matrix(zz, ncol = length(zz), nrow = ncol(design_mat), byrow = TRUE))
      
      A <- -t(design_mat * zz_mat) %*% design_mat / nrow(design_mat)
      A_inv <- solve(A)
      
      difference <- A_inv %*% S_beta
      qq <- sum(difference^2)
      
      hx <- rbind(hx, c(ii, qq, beta))
      beta <- beta - difference
      
      if (qq <= criterion) convergence <- 0
    }
  }
  
  if (link == "identity") {
    while (ii < max_itration & qq > criterion) {
      ii <- ii + 1
      
      gbeta <- design_mat %*% beta
      yy <- (as.numeric(time < tau) - gbeta * pmin(time, tau)) * wt
      yy_mat <- t(matrix(yy, ncol = length(yy), nrow = ncol(design_mat), byrow = TRUE))
      
      S_beta <- apply(design_mat * yy_mat, 2, mean)
      
      zz <- pmin(time, tau) * wt
      zz_mat <- t(matrix(zz, ncol = length(zz), nrow = ncol(design_mat), byrow = TRUE))
      
      A <- -t(design_mat * zz_mat) %*% design_mat / nrow(design_mat)
      A_inv <- solve(A)
      
      difference <- A_inv %*% S_beta
      qq <- sum(difference^2)
      
      hx <- rbind(hx, c(ii, qq, beta))
      beta <- beta - difference
      
      if (qq <= criterion) convergence <- 0
    }
  }
  
  convergence_information <- list()
  convergence_information$convergence <- convergence
  convergence_information$itration <- ii
  convergence_information$qq <- hx[, 2]
  
  history <- data.frame(hx)
  colnames(history) <- c("Iteration", "Q", paste0("beta", 1:(ncol(hx) - 2)))
  convergence_information$history <- history
  
  if (link == "log") {
    gbeta <- exp(design_mat %*% beta)
    yy <- (as.numeric(time < tau) - gbeta * pmin(time, tau)) * wt
    yy_mat <- t(vtm(yy, ncol(design_mat)))
    S_beta <- apply(design_mat * yy_mat, 2, mean)
    zz <- gbeta * pmin(time, tau) * wt
  }
  
  if (link == "identity") {
    gbeta <- design_mat %*% beta
    yy <- (as.numeric(time < tau) - gbeta * pmin(time, tau)) * wt
    yy_mat <- t(vtm(yy, ncol(design_mat)))
    S_beta <- apply(design_mat * yy_mat, 2, mean)
    zz <- pmin(time, tau) * wt
  }
  
  zz_mat <- t(vtm(zz, ncol(design_mat)))
  A <- -t(design_mat * zz_mat) %*% design_mat / nrow(design_mat)
  
  if (is.null(covariates4cens)) {
    psii_ks <- matrix(NA, nrow = length(time), ncol = ncol(design_mat))
    
    for (i in 1:nstrata) {
      idx <- strata == unique_strata[i]
      psii <- psii_list[[i]]
      
      distinct <- unique(sort(c(time[idx], tau)))
      
      wk1 <- t(vtm(time[idx], length(distinct)))
      wk2 <- vtm(distinct, sum(idx))
      wk3 <- wk1 <= wk2
      
      K_mat <- matrix(NA, nrow = ncol(design_mat), ncol = length(distinct))
      dK_mat <- matrix(NA, nrow = ncol(design_mat), ncol = length(distinct))
      
      for (jj in 1:ncol(design_mat)) {
        yy2 <- (
          as.numeric(time[idx] < tau) -
            gbeta[idx] * pmin(time[idx], tau)
        ) * wt[idx] * design_mat[idx, jj]
        
        yy2_mat <- t(vtm(yy2, ncol(wk3)))
        wk4 <- wk3 * yy2_mat
        
        K_mat[jj, ] <- apply(wk4, 2, mean)
        
        wk6 <- K_mat[jj, ]
        wk6[distinct > tau] <- 0
        
        wk7 <- diff(c(0, wk6))
        dK_mat[jj, ] <- wk7
      }
      
      psii_ks[idx, ] <- psii %*% t(dK_mat)
    }
    
    Ui <- design_mat * (yy %*% t(rep(1, length(beta)))) - psii_ks
  }
  
  if (!is.null(covariates4cens)) {
    wk1a <- (as.numeric(time < tau) - gbeta * pmin(time, tau)) * wt
    term1 <- t(vtm(wk1a, ncol(design_mat))) * design_mat
    
    wk2a <- (
      as.numeric(time < tau) -
        gbeta * pmin(time, tau)
    ) * wt * fc.cox$Lam0.ebz
    
    wk2b <- t(vtm(wk2a, ncol(design_mat))) * design_mat
    
    K_gamma <- t(wk2b) %*% covariates4cens / length(time)
    term2 <- t(K_gamma %*% t(fc.cox$eta_beta_i))
    
    distinct <- fc.cox$distinct
    
    wk1 <- t(vtm(time, length(distinct)))
    wk2 <- vtm(distinct, length(time))
    wk3 <- (wk1 >= wk2) * 1
    
    K_mat <- matrix(NA, nrow = ncol(design_mat), ncol = length(distinct))
    dK_mat <- matrix(NA, nrow = ncol(design_mat), ncol = length(distinct))
    
    for (jj in 1:ncol(design_mat)) {
      yy2 <- (
        as.numeric(time < tau) -
          gbeta * pmin(time, tau)
      ) * wt * fc.cox$ebz * design_mat[, jj]
      
      yy2_mat <- t(vtm(yy2, ncol(wk3)))
      wk4 <- wk3 * yy2_mat
      
      K_mat[jj, ] <- apply(wk4, 2, mean)
      
      wk6 <- K_mat[jj, ]
      wk6[distinct > tau] <- 0
      
      wk7 <- diff(c(0, wk6))
      dK_mat[jj, ] <- wk7
    }
    
    term3 <- fc.cox$eta_lam_i %*% t(dK_mat)
    
    Ui <- term1 + term2 + term3
  }
  
  B <- (t(Ui) %*% Ui) / nrow(design_mat)
  V <- solve(A) %*% B %*% solve(A)
  
  se <- sqrt(diag(V) / nrow(design_mat))
  
  low <- beta - se * abs(qnorm((1 - conf.int) / 2))
  upp <- beta + se * abs(qnorm((1 - conf.int) / 2))
  
  zstat <- beta / se
  pval <- (1 - pnorm(abs(zstat))) * 2
  
  result <- data.frame(cbind(beta, se, low, upp, zstat, pval))
  colnames(result) <- c(
    "Est",
    "SE",
    paste0("low_", conf.int),
    paste0("upp_", conf.int),
    "Z",
    "p"
  )
  
  variance_information <- list()
  variance_information$A <- A
  variance_information$B <- B
  variance_information$V <- V
  variance_information$n <- nrow(design_mat)
  
  Z <- list()
  Z$result <- result
  Z$beta <- beta
  Z$beta.var <- V / nrow(design_mat)
  Z$Ainv <- solve(A)
  Z$B <- B
  Z$predicted <- gbeta
  Z$link <- link
  
  if (is.null(covariates4cens)) {
    Z$ipcw_information <- kmc
    Z$nstrata <- nstrata
  }
  
  if (!is.null(covariates4cens)) {
    Z$ipcw_information <- fc.cox
    Z$nstrata <- NA
  }
  
  Z$convergence_information <- convergence_information
  Z$variance_information <- variance_information
  
  class(Z) <- "ahreg"
  Z
}

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
print("Embedded fun-ahreg.R functions: kmcens2(), cox1(), vtm(), ahreg()")

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

methods <- c("PO-only", "AH-Cox-only", "Stacked-GMM-AH-Cox")

cat("\nGlobal settings:\n")
print(list(
  n_sample = n_sample,
  nsim = nsim,
  tau = tau,
  conf_level = conf_level,
  methods = methods
))

# ============================================================
# PBC covariate pool
# ============================================================

prepare_pbc_pool <- function() {
  data(pbc, package = "survival")
  
  pbc_clean <- pbc[!is.na(pbc$trt), ]
  pbc_clean <- pbc_clean[
    complete.cases(pbc_clean[, c("age", "bili", "albumin")]),
  ]
  
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
# Data generation
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
    rate_c <- -log(0.50) / tau
    C_potential <- rexp(n, rate = rate_c)
  }
  
  if (censoring == "group_specific") {
    high_bili <- covariates$bili > 1.35
    scale_C <- ifelse(high_bili, 19.61, 5.83)
    C_potential <- rexp(n, rate = 1 / scale_C)
  }
  
  if (censoring == "covariate_dependent") {
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
# PO and AH-Cox wrappers
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

make_ah_covariates <- function(dat) {
  covariates <- as.matrix(dat$covariates[, c("age", "log_bili", "albumin")])
  colnames(covariates) <- c("age", "log_bilirubin", "albumin")
  covariates
}

make_cens_covariates <- function(dat) {
  covariates4cens <- as.matrix(dat$covariates[, c("age", "log_bili")])
  colnames(covariates4cens) <- c("age", "log_bilirubin")
  covariates4cens
}

safe_ahreg_call <- function(time, status, covariates, tau,
                            covariates4cens = NULL) {
  fit <- try(
    ahreg(
      time = time,
      status = status,
      covariates = covariates,
      tau = tau,
      covariates4cens = covariates4cens
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) return(NULL)
  if (is.null(fit$beta) || is.null(fit$result$SE)) return(NULL)
  
  beta_hat <- as.numeric(fit$beta)
  se_hat <- as.numeric(fit$result$SE)
  
  if (length(beta_hat) != length(beta_true)) return(NULL)
  if (length(se_hat) != length(beta_true)) return(NULL)
  
  if (any(!is.finite(beta_hat)) || any(!is.finite(se_hat))) return(NULL)
  
  list(
    fit = fit,
    beta = beta_hat,
    se = se_hat
  )
}

fit_AHCox_uno_with_ipcw <- function(dat, tau) {
  covariates <- make_ah_covariates(dat)
  covariates4cens <- make_cens_covariates(dat)
  
  no_random_censoring <- sum(dat$Delta == 0 & dat$X < tau) == 0
  
  if (no_random_censoring) {
    ah_out <- safe_ahreg_call(
      time = dat$X,
      status = dat$Delta,
      covariates = covariates,
      tau = tau,
      covariates4cens = NULL
    )
  } else {
    ah_out <- safe_ahreg_call(
      time = dat$X,
      status = dat$Delta,
      covariates = covariates,
      tau = tau,
      covariates4cens = covariates4cens
    )
  }
  
  if (is.null(ah_out)) {
    return(list(
      beta = rep(NA_real_, length(beta_true)),
      se = rep(NA_real_, length(beta_true)),
      ipcw = NULL,
      fit = NULL
    ))
  }
  
  wt <- ah_out$fit$ipcw_information$ipcw
  
  if (is.null(wt) || length(wt) != length(dat$X) || any(!is.finite(wt))) {
    return(list(
      beta = ah_out$beta,
      se = ah_out$se,
      ipcw = NULL,
      fit = ah_out$fit
    ))
  }
  
  ipcw <- list(
    weight = as.numeric(wt),
    event_tau = as.numeric(dat$X < tau & dat$Delta == 1),
    time_tau = pmin(dat$X, tau)
  )
  
  list(
    beta = ah_out$beta,
    se = ah_out$se,
    ipcw = ipcw,
    fit = ah_out$fit
  )
}

# ============================================================
# Moment functions and solvers
# ============================================================

moment_PO <- function(beta, Z, po) {
  mu <- as.numeric(exp(pmin(Z %*% beta, 20)))
  residual <- as.numeric(po$F - mu * po$R)
  Z * matrix(residual, nrow = nrow(Z), ncol = ncol(Z))
}

moment_AHCox <- function(beta, Z, ipcw) {
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
  
  A <- -t(Z) %*% (Z * as.numeric(mu * po$R)) / nrow(Z)
  
  M <- moment_PO(b, Z, po)
  se <- sandwich_just_identified(M, A)
  
  list(beta = b, se = se)
}

# ============================================================
# GMM helpers
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

safe_solve <- function(A, ridge = 1e-8) {
  out <- try(solve(A), silent = TRUE)
  
  if (!inherits(out, "try-error")) return(out)
  
  A2 <- A + diag(ridge, nrow(A))
  out2 <- try(solve(A2), silent = TRUE)
  
  if (!inherits(out2, "try-error")) return(out2)
  
  NULL
}

fit_GMM_AHCox_with_se <- function(Z, po, ipcw, start,
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
    M_ip <- moment_AHCox(b, Zc, ipcw)
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
  
  b1 <- as.numeric(fit1$par)
  names(b1) <- names(start)
  
  if (is_bad_beta(b1, limit = beta_limit)) {
    return(list(beta = rep(NA_real_, p), se = rep(NA_real_, p)))
  }
  
  M1 <- build_M(b1, scale0)$M
  n <- nrow(M1)
  
  S1 <- crossprod(scale(M1, center = TRUE, scale = FALSE)) / n
  W1 <- safe_solve(S1, ridge = 1e-6)
  
  if (is.null(W1)) {
    W1 <- diag(ncol(M1))
  }
  
  Q_opt <- function(b) {
    if (length(b) != length(start_c)) return(1e20)
    if (is_bad_beta(b, limit = beta_limit)) return(1e20)
    
    M <- build_M(b, scale0)$M
    g <- colMeans(M)
    
    if (any(!is.finite(g))) return(1e20)
    
    as.numeric(t(g) %*% W1 %*% g)
  }
  
  fit2 <- try(
    optim(
      par = b1,
      fn = Q_opt,
      method = "BFGS",
      control = list(maxit = 1000, reltol = 1e-10)
    ),
    silent = TRUE
  )
  
  if (inherits(fit2, "try-error")) {
    return(list(beta = rep(NA_real_, p), se = rep(NA_real_, p)))
  }
  
  b_c <- as.numeric(fit2$par)
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
  S_hat <- crossprod(scale(M_hat, center = TRUE, scale = FALSE)) / n
  
  bread <- safe_solve(t(D) %*% W1 %*% D, ridge = 1e-8)
  
  if (is.null(bread)) {
    return(list(beta = b_raw, se = rep(NA_real_, p)))
  }
  
  V_c <- bread %*% t(D) %*% W1 %*% S_hat %*% W1 %*% D %*% bread / n
  V_raw <- centered_var_to_raw(V_c, center_vec)
  
  se_raw <- sqrt(diag(V_raw))
  se_raw[!is.finite(se_raw)] <- NA_real_
  
  list(beta = b_raw, se = se_raw)
}

# ============================================================
# Simulation loop
# ============================================================

one_run <- function(n, beta, tau, pbc_pool,
                    shape = 1.0,
                    censoring = "independent",
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
  
  ahcox <- fit_AHCox_uno_with_ipcw(dat, tau)
  
  b_po <- fit_PO_with_se(dat$Z, po, beta)
  
  b_ahcox <- list(
    beta = ahcox$beta,
    se = ahcox$se
  )
  
  if (run_gmm) {
    b_gmm <- fit_GMM_AHCox_with_se(
      dat$Z,
      po,
      ahcox$ipcw,
      beta,
      center_vec = center_vec
    )
  } else {
    b_gmm <- list(
      beta = rep(NA_real_, length(beta)),
      se = rep(NA_real_, length(beta))
    )
  }
  
  list(
    beta = list(
      PO = b_po$beta,
      AHCox = b_ahcox$beta,
      GMM = b_gmm$beta
    ),
    se = list(
      PO = b_po$se,
      AHCox = b_ahcox$se,
      GMM = b_gmm$se
    ),
    observed_random_censor_rate = dat$observed_random_censor_rate,
    potential_censor_rate_at_tau = dat$potential_censor_rate_at_tau,
    total_censor_rate = dat$total_censor_rate
  )
}

run_sim <- function(nsim, n, beta, tau, pbc_pool,
                    shape = 1.0,
                    censoring = "independent",
                    run_gmm = TRUE) {
  p <- length(beta)
  
  center_vec <- get_center_from_pbc_pool(pbc_pool)
  
  res <- list(
    beta = list(
      PO = matrix(NA_real_, nsim, p),
      AHCox = matrix(NA_real_, nsim, p),
      GMM = matrix(NA_real_, nsim, p)
    ),
    se = list(
      PO = matrix(NA_real_, nsim, p),
      AHCox = matrix(NA_real_, nsim, p),
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
        run_gmm = run_gmm,
        center_vec = center_vec
      ),
      silent = TRUE
    )
    
    if (inherits(out, "try-error")) {
      next
    }
    
    res$beta$PO[s, ] <- out$beta$PO
    res$beta$AHCox[s, ] <- out$beta$AHCox
    res$beta$GMM[s, ] <- out$beta$GMM
    
    res$se$PO[s, ] <- out$se$PO
    res$se$AHCox[s, ] <- out$se$AHCox
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
      ipcw_method = "AH-Cox via embedded fun-ahreg.R",
      nsim = nsim
    )
  )
}

# ============================================================
# Summary functions
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
    summ_one(
      res$res$beta$PO,
      res$res$se$PO,
      res$beta_true,
      "PO-only",
      setting_row
    ),
    summ_one(
      res$res$beta$AHCox,
      res$res$se$AHCox,
      res$beta_true,
      "AH-Cox-only",
      setting_row
    ),
    summ_one(
      res$res$beta$GMM,
      res$res$se$GMM,
      res$beta_true,
      "Stacked-GMM-AH-Cox",
      setting_row
    )
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
    "AH-Cox-only",
    "Stacked-GMM-AH-Cox"
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
  
  wide_sd <- reshape(
    tmp[, c(key_cols, "Estimator", "Empirical_SD")],
    idvar = key_cols,
    timevar = "Estimator",
    direction = "wide"
  )
  names(wide_sd) <- gsub("Empirical_SD\\.", "SD_", names(wide_sd))
  
  wide_cov <- reshape(
    tmp[, c(key_cols, "Estimator", "Coverage_Distance")],
    idvar = key_cols,
    timevar = "Estimator",
    direction = "wide"
  )
  names(wide_cov) <- gsub("Coverage_Distance\\.", "CovDist_", names(wide_cov))
  
  out <- Reduce(
    function(x, y) merge(x, y, by = key_cols),
    list(wide_bias, wide_sd, wide_cov)
  )
  
  out$GMM_better_than_AHCox_bias <-
    out$`AbsRelBias_Stacked-GMM-AH-Cox` < out$`AbsRelBias_AH-Cox-only`
  
  out$GMM_SD_smaller_than_PO <-
    out$`SD_Stacked-GMM-AH-Cox` < out$`SD_PO-only`
  
  out$GMM_SD_smaller_than_AHCox <-
    out$`SD_Stacked-GMM-AH-Cox` < out$`SD_AH-Cox-only`
  
  out$GMM_SD_best <-
    out$GMM_SD_smaller_than_PO & out$GMM_SD_smaller_than_AHCox
  
  out$GMM_better_than_AHCox_coverage <-
    out$`CovDist_Stacked-GMM-AH-Cox` < out$`CovDist_AH-Cox-only`
  
  out
}

make_target_check_table <- function(summary_table,
                                    conf_level = 0.95,
                                    bias_tolerance = 0.03,
                                    coverage_tolerance = 0.03) {
  key_cols <- c("Model", "Weibull_shape", "Censoring_label", "Parameter")
  
  wide_bias <- reshape(
    summary_table[, c(key_cols, "Estimator", "Rel_Bias")],
    idvar = key_cols,
    timevar = "Estimator",
    direction = "wide"
  )
  names(wide_bias) <- gsub("Rel_Bias\\.", "Bias_", names(wide_bias))
  
  wide_sd <- reshape(
    summary_table[, c(key_cols, "Estimator", "Empirical_SD")],
    idvar = key_cols,
    timevar = "Estimator",
    direction = "wide"
  )
  names(wide_sd) <- gsub("Empirical_SD\\.", "SD_", names(wide_sd))
  
  wide_cov <- reshape(
    summary_table[, c(key_cols, "Estimator", "Coverage")],
    idvar = key_cols,
    timevar = "Estimator",
    direction = "wide"
  )
  names(wide_cov) <- gsub("Coverage\\.", "Cov_", names(wide_cov))
  
  out <- Reduce(
    function(x, y) merge(x, y, by = key_cols),
    list(wide_bias, wide_sd, wide_cov)
  )
  
  bias_cols <- c(
    "Bias_PO-only",
    "Bias_AH-Cox-only",
    "Bias_Stacked-GMM-AH-Cox"
  )
  
  out$Bias_range_abs <- apply(abs(out[, bias_cols]), 1, function(x) max(x) - min(x))
  out$Bias_similar_across_methods <- out$Bias_range_abs <= bias_tolerance
  
  out$GMM_SD_smaller_than_PO <-
    out$`SD_Stacked-GMM-AH-Cox` < out$`SD_PO-only`
  
  out$GMM_SD_smaller_than_AHCox <-
    out$`SD_Stacked-GMM-AH-Cox` < out$`SD_AH-Cox-only`
  
  out$GMM_SD_best <-
    out$GMM_SD_smaller_than_PO & out$GMM_SD_smaller_than_AHCox
  
  out$PO_cov_close_095 <-
    abs(out$`Cov_PO-only` - conf_level) <= coverage_tolerance
  
  out$AHCox_cov_close_095 <-
    abs(out$`Cov_AH-Cox-only` - conf_level) <= coverage_tolerance
  
  out$GMM_cov_close_095 <-
    abs(out$`Cov_Stacked-GMM-AH-Cox` - conf_level) <= coverage_tolerance
  
  out$All_coverage_close_095 <-
    out$PO_cov_close_095 &
    out$AHCox_cov_close_095 &
    out$GMM_cov_close_095
  
  out
}

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
cat("Expected: 8 settings x 3 estimators x 4 parameters = 96 rows\n")

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