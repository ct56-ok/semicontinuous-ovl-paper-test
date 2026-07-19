library(readxl)

# Define the distribution fitting and information criteria function
fit_distribution_criteria <- function(data, value_col, scale_col, data_label) {
  
  # 1. Data preprocessing and cleaning
  raw_values <- data[[value_col]] * 1000 / data[[scale_col]]
  x <- raw_values[raw_values > 0]
  n <- length(x)
  
  if (n == 0) {
    stop("No positive values found after filtering for ", value_col)
  }
  
  # Plot histograms for raw data and log-transformed data
  par(mfrow = c(1, 2))
  hist(raw_values, main = paste("Raw:", data_label), xlab = value_col, breaks = 30)
  hist(log(x), main = paste("Log-transformed:", data_label), xlab = paste("log(", value_col, ")"), breaks = 30)
  
  # 2. Parameter estimation for each distribution
  # Log-normal
  hmu <- mean(log(x))
  hsigma <- sqrt(sum((log(x) - hmu)^2) / n)
  
  # Exponential
  hlambda <- 1 / mean(x)
  
  # Normal
  hmu1 <- mean(x)
  hsigma1 <- sqrt(sum((x - hmu1)^2) / n)
  
  # Gamma (using Nelder-Mead optimization)
  h1 <- (mean(x))^2 / var(x)
  k1 <- mean(x) / var(x)
  para0 <- c(h1, k1)
  
  # Gamma log-likelihood function (fixed the typo: log(gamma(beta)) -> log(beta))
  Lngamma <- function(para) {
    alpha <- para[1]
    beta <- para[2]
    # Prevent non-positive values causing log errors during optimization
    if (alpha <= 0 || beta <= 0) return(1e10) 
    lngamma <- -n * log(gamma(alpha)) - n * alpha * log(beta) + 
      (alpha - 1) * sum(log(x)) - sum(x) / beta
    return(-lngamma) # optim minimizes by default, so return negative log-likelihood
  }
  
  opt1 <- tryCatch(
    optim(para0, Lngamma, method = "Nelder-Mead", control = list(maxit = 2000))$par,
    error = function(e) c(NA, NA)
  )
  
  # 3. Define log-likelihood functions for each distribution
  Lnlnormal <- function(sigma, mu) {
    -n / 2 * log(2 * pi * sigma^2) - sum(log(x)) - (1 / (2 * sigma^2)) * sum((log(x) - mu)^2)
  }
  
  Lnexp <- function(lambda) {
    n * log(lambda) - lambda * sum(x)
  }
  
  Lnnormal <- function(sigma1, mu1) {
    -n / 2 * log(2 * pi * sigma1^2) - (1 / (2 * sigma1^2)) * sum((x - mu1)^2)
  }
  
  Lngamma_val <- function(alpha, beta) {
    -n * log(gamma(alpha)) - n * alpha * log(beta) + (alpha - 1) * sum(log(x)) - sum(x) / beta
  }
  
  # 4. Calculate AIC and BIC
  k_vec <- c(2, 1, 2, 2) # Number of parameters: LogNormal, Exp, Normal, Gamma
  
  logL_vec <- c(
    Lnlnormal(hsigma, hmu),
    Lnexp(hlambda),
    Lnnormal(hsigma1, hmu1),
    ifelse(any(is.na(opt1)), -Inf, Lngamma_val(opt1[1], opt1[2]))
  )
  
  AIC <- -2 * logL_vec + 2 * k_vec
  BIC <- -2 * logL_vec + log(n) * k_vec
  
  # 5. Format and print results
  cat("\n========== Distribution Fit for:", data_label, "==========\n")
  res_df <- data.frame(
    Distribution = c("Log-Normal", "Exponential", "Normal", "Gamma"),
    AIC = round(AIC, 4),
    BIC = round(BIC, 4)
  )
  print(res_df)
  cat("==================================================\n\n")
  
  return(res_df)
}

# ================= Function Call Example =================

# Read the data
aa <- read_xlsx("C:\\Users\\tingc\\Desktop\\ovl code\\应用数据\\末期141.xlsx")

# Call the function
result_142 <- fit_distribution_criteria(
  data = aa,
  value_col = "F_WHOLE",
  scale_col = "rtkcal",
  data_label = "Dataset 142 - F_WHOLE"
)