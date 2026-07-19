# ================= 1. Core Mathematical Formulas =================
# Calculates Rho, Lambda, and Delta based on given parameters
calc_overlap_metrics <- function(p1, p2, mu1, mu2, sig1, sig2) {
  eps <- 1e-12
  
  # Rho calculation
  rho <- min(1 - p1, 1 - p2) + sqrt(p1 * p2) * 
    (sqrt(2 * sqrt(sig1) * sqrt(sig2) / (sig1 + sig2)) * 
       exp(-((mu1 - mu2)^2) / (4 * (sig1 + sig2))))
  
  # Lambda calculation
  lambda <- min(1 - p1, 1 - p2) + 
    2 * sqrt(2) * (p1 * p2) / sqrt(sig1 + sig2) * 
    exp((sig2 * sig1 - (mu1 - mu2)^2 - 2 * sig2 * mu1 - 2 * sig1 * mu2) / 
          (2 * (sig2 + sig1))) / 
    ((p1 * exp((sig1 - 4 * mu1) / 4) / sqrt(sig1)) + 
       (p2 * exp((sig2 - 4 * mu2) / 4) / sqrt(sig2)))
  
  # Delta calculation
  if (abs(sig1 - sig2) < eps) {
    delta <- min(1 - p1, 1 - p2) + (p1 + p2) * pnorm(-abs(mu1 - mu2) / (2 * sqrt(sig1)))
  } else {
    term <- sqrt(2 * (sig2 - sig1) * log(sqrt(sig2 / sig1)) + (mu1 - mu2)^2)
    ha1 <- (sig2 * mu1 - sig1 * mu2) / (sig2 - sig1) - sqrt(sig2 * sig1) / abs(sig2 - sig1) * term
    ha2 <- (sig2 * mu1 - sig1 * mu2) / (sig2 - sig1) + sqrt(sig2 * sig1) / abs(sig2 - sig1) * term
    
    cond <- (ha1 - mu1) / sqrt(sig1) < (ha1 - mu2) / sqrt(sig2)
    if (cond) {
      delta <- min(1 - p1, 1 - p2) + 
        p1 * (1 - (pnorm((ha2 - mu1) / sqrt(sig1)) - pnorm((ha1 - mu1) / sqrt(sig1)))) +
        p2 * (pnorm((ha2 - mu2) / sqrt(sig2)) - pnorm((ha1 - mu2) / sqrt(sig2)))
    } else {
      delta <- min(1 - p1, 1 - p2) + 
        p2 * (1 - (pnorm((ha2 - mu2) / sqrt(sig2)) - pnorm((ha1 - mu2) / sqrt(sig2)))) +
        p1 * (pnorm((ha2 - mu1) / sqrt(sig1)) - pnorm((ha1 - mu1) / sqrt(sig1)))
    }
  }
  return(c(rho = rho, lambda = lambda, delta = delta))
}

# ================= 2. Main Permutation Test Function =================
run_permutation_test <- function(N = 100, n1 = 15, n2 = 150, 
                                 mu1 = 3, mu2 = 3, 
                                 sigma1 = 0.5, sigma2 = 0.5, 
                                 p1 = 0.5, p2 = 0.5, 
                                 n_perm = 3000, seed = 222) {
  set.seed(seed)
  
  # Initialize storage vectors for performance
  hrho <- hlambda <- hdelta <- LR1 <- Wa <- numeric(N)
  reject_rho <- reject_lambda <- reject_delta <- logical(N)
  
  # Loop over N simulations
  for (i in 1:N) {
    # --- Generate Group 1 Data ---
    X1 <- rbinom(n1, 1, p1)
    n1_0 <- sum(X1)
    if (n1_0 >= n1 - 1 || n1_0 <= 2) { next }
    n1_1 <- n1 - n1_0
    S1 <- rnorm(n1_0, mu1, sqrt(sigma1))
    Y1 <- c(rep(0, n1_1), S1)
    
    # --- Generate Group 2 Data ---
    X2 <- rbinom(n2, 1, p2)
    n2_0 <- sum(X2)
    if (n2_0 >= n2 - 1 || n2_0 <= 2) { next }
    n2_1 <- n2 - n2_0
    S2 <- rnorm(n2_0, mu2, sqrt(sigma2))
    Y2 <- c(rep(0, n2_1), S2)
    
    # Combine data for permutation
    Y3 <- c(Y1, Y2)
    
    # --- Calculate Observed Statistics ---
    hp1 <- n1_0 / n1; hmu1 <- mean(S1); hsigma1 <- sum((S1 - hmu1)^2) / n1_0
    hp2 <- n2_0 / n2; hmu2 <- mean(S2); hsigma2 <- sum((S2 - hmu2)^2) / n2_0
    
    obs_metrics <- calc_overlap_metrics(hp1, hp2, hmu1, hmu2, hsigma1, hsigma2)
    hrho[i] <- obs_metrics["rho"]
    hlambda[i] <- obs_metrics["lambda"]
    hdelta[i] <- obs_metrics["delta"]
    
    # --- Permutation Loop ---
    perm_rho <- perm_lambda <- perm_delta <- numeric(n_perm)
    for (K in 1:n_perm) {
      idx1 <- sample.int(length(Y3), n1, replace = FALSE)
      Z1 <- Y3[idx1]; Z2 <- Y3[-idx1]
      
      Z11 <- Z1[Z1 != 0]; Z22 <- Z2[Z2 != 0]
      n11_0 <- length(Z11); n22_0 <- length(Z22)
      
      # Skip if variance cannot be calculated
      if (n11_0 <= 1 || n22_0 <= 1) { next }
      
      hp11 <- n11_0 / n1; hmu11 <- mean(Z11); hsigma11 <- sum((Z11 - hmu11)^2) / n11_0
      hp22 <- n22_0 / n2; hmu22 <- mean(Z22); hsigma22 <- sum((Z22 - hmu22)^2) / n22_0
      
      if (hsigma11 == 0 || hsigma22 == 0) { next }
      
      perm_metrics <- calc_overlap_metrics(hp11, hp22, hmu11, hmu22, hsigma11, hsigma22)
      perm_rho[K] <- perm_metrics["rho"]
      perm_lambda[K] <- perm_metrics["lambda"]
      perm_delta[K] <- perm_metrics["delta"]
    }
    
    # --- Calculate P-values based on permutation distribution ---
    valid_perm_rho <- perm_rho[perm_rho != 0]
    valid_perm_lambda <- perm_lambda[perm_lambda != 0]
    valid_perm_delta <- perm_delta[perm_delta != 0]
    
    if (length(valid_perm_rho) > 0) {
      q_rho <- quantile(valid_perm_rho, 0.05, na.rm = TRUE)
      reject_rho[i] <- hrho[i] < q_rho
    }
    if (length(valid_perm_lambda) > 0) {
      q_lambda <- quantile(valid_perm_lambda, 0.05, na.rm = TRUE)
      reject_lambda[i] <- hlambda[i] < q_lambda
    }
    if (length(valid_perm_delta) > 0) {
      q_delta <- quantile(valid_perm_delta, 0.05, na.rm = TRUE)
      reject_delta[i] <- hdelta[i] < q_delta
    }
    
    # --- Likelihood Ratio Test (LR) ---
    pp <- (n1_0 + n2_0) / (n1 + n2)
    uu <- mean(c(S1, S2))
    ww <- (sum((S1 - uu)^2) + sum((S2 - uu)^2)) / (n1_0 + n2_0)
    
    T2 <- (pp / hp1)^n1_0 * ((1 - pp) / (1 - hp1))^n1_1 * 
      ((1 - pp) / (1 - hp2))^n2_1 * (pp / hp2)^n2_0 * 
      (hsigma1 / ww)^(n1_0 / 2) * (hsigma2 / ww)^(n2_0 / 2)
    
    if (is.finite(T2) && T2 > 0) {
      LR1[i] <- -2 * log(T2)
    } else {
      LR1[i] <- NA
    }
    
    # --- Wald Test ---
    C_mat <- matrix(c(1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1), 3, 6)
    eps_w <- 1e-12
    hp1_w <- min(max(hp1, eps_w), 1 - eps_w)
    hp2_w <- min(max(hp2, eps_w), 1 - eps_w)
    
    IW <- diag(c(
      n1_0 / hp1_w^2 + n1_1 / (1 - hp1_w)^2,
      n1_0 / hsigma1,
      n1_0 / (2 * hsigma1^2),
      n2_0 / hp2_w^2 + n2_1 / (1 - hp2_w)^2,
      n2_0 / hsigma2,
      n2_0 / (2 * hsigma2^2)
    ), 6, 6)
    
    H <- matrix(c(hp1 - hp2, hmu1 - hmu2, hsigma1 - hsigma2), 1, 3)
    Wa[i] <- tryCatch(
      as.numeric(H %*% solve(C_mat %*% solve(IW) %*% t(C_mat)) %*% t(H)),
      error = function(e) NA_real_
    )
  }
  
  # ================= 3. Calculate Type I Error Rates =================
  valid_idx <- !is.na(hrho)
  TYR <- mean(reject_rho[valid_idx], na.rm = TRUE)
  TYL <- mean(reject_lambda[valid_idx], na.rm = TRUE)
  TYD <- mean(reject_delta[valid_idx], na.rm = TRUE)
  
  valid_LR <- LR1[!is.na(LR1)]
  TYLR <- mean(valid_LR > qchisq(0.95, 3), na.rm = TRUE)
  
  valid_Wa <- Wa[!is.na(Wa)]
  TYWa <- mean(valid_Wa > qchisq(0.95, 3), na.rm = TRUE)
  
  # Return results
  return(list(
    TYR = TYR,
    TYL = TYL,
    TYD = TYD,
    TYLR = TYLR,
    TYWa = TYWa
  ))
}

# ================= 4. One-Click Execution =================
test_results <- run_permutation_test(N = 1000, n1 = 30, n2 = 150, 
                                     mu1 = 3, mu2 = 2, 
                                     sigma1 = 1, sigma2 = 1, 
                                     p1 = 0.4, p2 = 0.7)

# Print the Type I Error Rates
print(test_results)
