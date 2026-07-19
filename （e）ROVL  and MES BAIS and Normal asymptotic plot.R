# ================= 1. Core Mathematical Formulas =================
# This function calculates the theoretical or estimated values of Rho, Lambda, and Delta
calc_overlap_metrics <- function(p1, p2, mu1, mu2, sig1, sig2) {
  eps <- 1e-12 # Small epsilon value to avoid division by zero or exact equality issues
  
  # Calculate Rho
  rho <- min(1 - p1, 1 - p2) + sqrt(p1 * p2) * 
    (sqrt(2 * sqrt(sig1) * sqrt(sig2) / (sig1 + sig2)) * 
       exp(-((mu1 - mu2)^2) / (4 * (sig1 + sig2))))
  
  # Calculate Lambda
  lambda <- min(1 - p1, 1 - p2) + 
    2 * sqrt(2) * (p1 * p2) / sqrt(sig1 + sig2) * 
    exp((sig2 * sig1 - (mu1 - mu2)^2 - 2 * sig2 * mu1 - 2 * sig1 * mu2) / 
          (2 * (sig2 + sig1))) / 
    ((p1 * exp((sig1 - 4 * mu1) / 4) / sqrt(sig1)) + 
       (p2 * exp((sig2 - 4 * mu2) / 4) / sqrt(sig2)))
  
  # Calculate Delta based on variance equality condition
  if (abs(sig1 - sig2) < eps) {
    # Case: Variances are approximately equal
    delta <- min(1 - p1, 1 - p2) + (p1 + p2) * pnorm(-abs(mu1 - mu2) / (2 * sqrt(sig1)))
  } else {
    # Case: Variances are unequal
    term <- sqrt(2 * (sig2 - sig1) * log(sqrt(sig2 / sig1)) + (mu1 - mu2)^2)
    ha1 <- (sig2 * mu1 - sig1 * mu2) / (sig2 - sig1) - sqrt(sig2 * sig1) / abs(sig2 - sig1) * term
    ha2 <- (sig2 * mu1 - sig1 * mu2) / (sig2 - sig1) + sqrt(sig2 * sig1) / abs(sig2 - sig1) * term
    
    # Determine integration bounds based on density intersection
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

# ================= 2. Main Simulation and Plotting Function =================
run_simulation <- function(N = 2000, n1 = 100, n2 = 100, 
                           mu1 = 2, mu2 = 2, 
                           sigma1 = 1, sigma2 = 2, 
                           p1 = 0.5, p2 = 0.6, 
                           seed = 223) {
  set.seed(seed) # Set seed for reproducibility
  
  # Step 1: Calculate the theoretical true values (fixed across the entire simulation)
  true_metrics <- calc_overlap_metrics(p1, p2, mu1, mu2, sigma1, sigma2)
  true_rho <- true_metrics["rho"]
  true_lambda <- true_metrics["lambda"]
  true_delta <- true_metrics["delta"]
  
  # Step 2: Initialize storage vectors for the estimates
  est_rho <- numeric(N)
  est_lambda <- numeric(N)
  est_delta <- numeric(N)
  
  # Step 3: Monte Carlo Simulation Loop
  for (i in 1:N) {
    # Generate Group 1 data
    X1 <- rbinom(n1, 1, p1)
    n1_0 <- sum(X1)
    # Skip simulation if the number of non-zero observations is too small or too large
    if (n1_0 >= n1 - 1 || n1_0 <= 2) { next } 
    S1 <- rnorm(n1_0, mu1, sqrt(sigma1))
    
    # Generate Group 2 data
    X2 <- rbinom(n2, 1, p2)
    n2_0 <- sum(X2)
    if (n2_0 >= n2 - 1 || n2_0 <= 2) { next }
    S2 <- rnorm(n2_0, mu2, sqrt(sigma2))
    
    # Calculate sample estimates
    hp1 <- n1_0 / n1
    hp2 <- n2_0 / n2
    hmu1 <- mean(S1)
    hmu2 <- mean(S2)
    hsigma1 <- sum((S1 - hmu1)^2) / n1_0
    hsigma2 <- sum((S2 - hmu2)^2) / n2_0
    
    # Call the core function to calculate estimated metrics
    est <- calc_overlap_metrics(hp1, hp2, hmu1, hmu2, hsigma1, hsigma2)
    est_rho[i] <- est["rho"]
    est_lambda[i] <- est["lambda"]
    est_delta[i] <- est["delta"]
  }
  
  # Step 4: Remove NA values generated from skipped iterations
  valid_rho <- est_rho[!is.na(est_rho)]
  valid_lambda <- est_lambda[!is.na(est_lambda)]
  valid_delta <- est_delta[!is.na(est_delta)]
  
  # Step 5: Calculate evaluation metrics (MSE, Bias, Variance)
  MSE <- c(
    rho = mean((valid_rho - true_rho)^2),
    lambda = mean((valid_lambda - true_lambda)^2),
    delta = mean((valid_delta - true_delta)^2)
  )
  
  Bias <- c(
    rho = mean(valid_rho - true_rho),
    lambda = mean(valid_lambda - true_lambda),
    delta = mean(valid_delta - true_delta)
  )
  
  Var <- c(
    rho = var(valid_rho - true_rho),
    lambda = var(valid_lambda - true_lambda),
    delta = var(valid_delta - true_delta)
  )
  
  # Step 6: Plot standardized histograms
  par(mfrow = c(1, 3))
  
  # Rho histogram
  Trho1 <- (valid_rho - true_rho) / sqrt(var(valid_rho))
  hist(Trho1, xlim = c(-3, 3), ylim = c(0, 0.5), probability = TRUE,
       col = 'gray', xlab = "", main = expression(paste((hat(rho) - rho) / var(hat(rho)))))
  lines(density(Trho1, bw = 1), col = 'red', lwd = 3)
  
  # Lambda histogram
  Tlambda1 <- (valid_lambda - true_lambda) / sqrt(var(valid_lambda))
  hist(Tlambda1, xlim = c(-3, 3), ylim = c(0, 0.5), probability = TRUE,
       col = 'gray', xlab = "", main = expression(paste((hat(lambda) - lambda) / var(hat(lambda)))))
  lines(density(Tlambda1, bw = 1), col = 'blue', lwd = 3)
  
  # Delta histogram
  Tdelta1 <- (valid_delta - true_delta) / sqrt(var(valid_delta))
  hist(Tdelta1, xlim = c(-3, 3), ylim = c(0, 0.6), probability = TRUE,
       col = 'gray', xlab = "", main = expression(paste((hat(Delta) - Delta) / var(hat(Delta)))))
  lines(density(Tdelta1, bw = 1), col = 'green', lwd = 3)
  
  # Step 7: Return the summarized results as a list
  result <- list(
    True_Values = c(Rho = true_rho, Lambda = true_lambda, Delta = true_delta),
    MSE = MSE,
    Bias = Bias,
    Variance = Var
  )
  
  return(result)
}

# ================= 3. One-Click Execution Example =================
sim_results <- run_simulation(N = 2000, n1 = 100, n2 = 100, 
                              mu1 = 2, mu2 = 2, 
                              sigma1 = 1, sigma2 = 2, 
                              p1 = 0.5, p2 = 0.6)

# Print the final evaluation results
print(sim_results)