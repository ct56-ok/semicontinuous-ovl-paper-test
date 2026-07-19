# =====================================================
# 1. 准备工作：加载必要的包并设置数据路径
# =====================================================
if (!require(readxl)) install.packages("readxl")
library(readxl)

# 设置数据所在的文件夹路径（请根据实际情况修改）
data_dir <- "C:\\Users\\tingc\\Desktop\\ovl code\\应用数据"

# 读取两组数据
aa <- read_xlsx(file.path(data_dir, "末期141.xlsx"))
cc <- read_xlsx(file.path(data_dir, "末期142.xlsx"))

# 检查数据列是否存在，防止读取失败
if (!("V_DOL" %in% names(aa)) || !("V_DOL" %in% names(cc))) {
  stop("Error: Column 'V_DOL' not found in the Excel files.")
}

# =====================================================
# 2. 核心函数定义（已修复 NA 报错问题）
# =====================================================
permutation_test <- function(data1, data2, n_perm = 3000, seed = 222) {
  set.seed(seed)
  
  # 1. Data preprocessing (Log-normal assumption based on original code)
  bb <- data1 * 1000 / data1 
  dd <- data2 * 1000 / data2
  
  # Extract non-zero components for log-transform
  bb1 <- log(bb[bb != 0])
  dd1 <- log(dd[dd != 0])
  
  # Observed statistics
  n1 <- length(bb)
  n2 <- length(dd)
  n1_0 <- length(bb1)
  n2_0 <- length(dd1)
  hp1 <- n1_0 / n1
  hp2 <- n2_0 / n2
  hmu1 <- mean(bb1)
  hmu2 <- mean(dd1)
  hsigma1 <- sum((bb1 - hmu1)^2) / n1_0
  hsigma2 <- sum((dd1 - hmu2)^2) / n2_0
  
  # Pooled statistics for LR test
  pp <- (n1_0 + n2_0) / (n1 + n2)
  uu <- mean(c(bb1, dd1))
  ww <- (sum((bb1 - uu)^2) + sum((dd1 - uu)^2)) / (n1_0 + n2_0)
  
  # 2. Helper function to calculate rho, lambda, delta
  calc_metrics <- function(p1, p2, mu1, mu2, sig1, sig2) {
    eps <- 1e-12
    
    # [FIX] Safety check: If any parameter is NA/NaN or variance is 0, return NA
    if (any(is.na(c(p1, p2, mu1, mu2, sig1, sig2))) || sig1 == 0 || sig2 == 0) {
      return(c(rho = NA_real_, lambda = NA_real_, delta = NA_real_))
    }
    
    # Rho
    rho <- min(1 - p1, 1 - p2) + sqrt(p1 * p2) * 
      (sqrt(2 * sqrt(sig1) * sqrt(sig2) / (sig1 + sig2)) * 
         exp(-((mu1 - mu2)^2) / (4 * (sig1 + sig2))))
    
    # Lambda
    lambda <- min(1 - p1, 1 - p2) + 
      2 * sqrt(2) * (p1 * p2) / sqrt(sig1 + sig2) * 
      exp((sig2 * sig1 - (mu1 - mu2)^2 - 2 * sig2 * mu1 - 2 * sig1 * mu2) / 
            (2 * (sig2 + sig1))) / 
      ((p1 * exp((sig1 - 4 * mu1) / 4) / sqrt(sig1)) + 
         (p2 * exp((sig2 - 4 * mu2) / 4) / sqrt(sig2)))
    
    # Delta
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
  
  # Observed metrics
  obs_metrics <- calc_metrics(hp1, hp2, hmu1, hmu2, hsigma1, hsigma2)
  
  # 3. Permutation loop
  perm_rho <- numeric(n_perm)
  perm_lambda <- numeric(n_perm)
  perm_delta <- numeric(n_perm)
  LR1 <- numeric(n_perm)
  Wa <- numeric(n_perm)
  
  bd <- c(bb, dd)
  
  for (k in 1:n_perm) {
    idx1 <- sample.int(length(bd), n1, replace = FALSE)
    Z1 <- bd[idx1]
    Z2 <- bd[-idx1]
    
    Z11 <- Z1[Z1 != 0]
    Z22 <- Z2[Z2 != 0]
    n11_0 <- length(Z11)
    n22_0 <- length(Z22)
    
    # [FIX] Ensure there are at least 2 non-zero values in each group to calculate variance
    if (n11_0 < 2 || n22_0 < 2) next
    
    hp11 <- n11_0 / n1
    hp22 <- n22_0 / n2
    hmu11 <- mean(Z11)
    hmu22 <- mean(Z22)
    hsigma11 <- sum((Z11 - hmu11)^2) / n11_0
    hsigma22 <- sum((Z22 - hmu22)^2) / n22_0
    
    # [FIX] Skip this permutation if variance is exactly 0 or NA
    if (is.na(hsigma11) || is.na(hsigma22) || hsigma11 == 0 || hsigma22 == 0) next
    
    # Permutation metrics
    pm <- calc_metrics(hp11, hp22, hmu11, hmu22, hsigma11, hsigma22)
    perm_rho[k] <- pm["rho"]
    perm_lambda[k] <- pm["lambda"]
    perm_delta[k] <- pm["delta"]
    
    # LR Test
    n1_1 <- n1 - n1_0
    n2_1 <- n2 - n2_0
    T2 <- (pp / hp1)^n1_0 * ((1 - pp) / (1 - hp1))^n1_1 *
      ((1 - pp) / (1 - hp2))^n2_1 * (pp / hp2)^n2_0 *
      (hsigma1 / ww)^(n1_0 / 2) * (hsigma2 / ww)^(n2_0 / 2)
    LR1[k] <- if (is.finite(T2) && T2 > 0) -2 * log(T2) else NA
    
    # Wald Test
    C_mat <- matrix(c(1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1), 3, 6)
    eps_w <- 1e-12
    hp1_w <- min(max(hp1, eps_w), 1 - eps_w)
    hp2_w <- min(max(hp2, eps_w), 1 - eps_w)
    IW <- diag(c(
      n1_0 / hp1_w^2 + n1_1 / (1 - hp1_w)^2, n1_0 / hsigma1, n1_0 / (2 * hsigma1^2),
      n2_0 / hp2_w^2 + n2_1 / (1 - hp2_w)^2, n2_0 / hsigma2, n2_0 / (2 * hsigma2^2)
    ), 6, 6)
    H <- matrix(c(hp1 - hp2, hmu1 - hmu2, hsigma1 - hsigma2), 1, 3)
    Wa[k] <- tryCatch(
      as.numeric(H %*% solve(C_mat %*% solve(IW) %*% t(C_mat)) %*% t(H)),
      error = function(e) NA_real_
    )
  }
  
  # 4. Calculate p-values and Type I error rates
  valid_rho <- perm_rho[!is.na(perm_rho)]
  valid_lambda <- perm_lambda[!is.na(perm_lambda)]
  valid_delta <- perm_delta[!is.na(perm_delta)]
  valid_LR <- LR1[!is.na(LR1)]
  valid_Wa <- Wa[!is.na(Wa)]
  
  p_rho <- mean(valid_rho <= obs_metrics["rho"])
  p_lambda <- mean(valid_lambda <= obs_metrics["lambda"])
  p_delta <- mean(valid_delta <= obs_metrics["delta"])
  
  lr_stat <- valid_LR[1]
  wa_stat <- valid_Wa[1]
  p_lr <- if (!is.na(lr_stat)) 1 - pchisq(lr_stat, 3) else NA
  p_wa <- if (!is.na(wa_stat)) 1 - pchisq(wa_stat, 3) else NA
  
  TYLR <- if (length(valid_LR) > 0) mean(valid_LR > qchisq(0.95, 3)) else NA
  TYWa <- if (length(valid_Wa) > 0) mean(valid_Wa > qchisq(0.95, 3)) else NA
  
  # 5. Return structured results
  return(list(
    Observed_Metrics = obs_metrics,
    P_Values = c(rho = p_rho, lambda = p_lambda, delta = p_delta, LR = p_lr, Wald = p_wa),
    Type_I_Error_Rates = c(LR = TYLR, Wald = TYWa),
    Permutation_SD = c(
      rho = ifelse(length(valid_rho) >= 2, sd(valid_rho), NA),
      lambda = ifelse(length(valid_lambda) >= 2, sd(valid_lambda), NA),
      delta = ifelse(length(valid_delta) >= 2, sd(valid_delta), NA)
    )
  ))
}

# =====================================================
# 3. 一键调用函数
# =====================================================
cat("Starting permutation test...\n")
results <- permutation_test(
  data1 = aa$V_DOL*1000/aa$rtkcal, 
  data2 = cc$V_DOL*1000/cc$rtkcal, 
  n_perm = 3000
)
cat("Test completed successfully!\n\n")

# =====================================================
# 4. 格式化输出结果
# =====================================================
# 将结果转换为更易读的数据框
final_table <- data.frame(
  Method = c("Rho", "Lambda", "Delta", "Likelihood Ratio (LR)", "Wald Test"),
  Observed_Statistic = c(results$Observed_Metrics, results$P_Values[4:5]), # 注意：这里仅展示统计量
  P_Value = results$P_Values,
  Permutation_SD = c(results$Permutation_SD, NA, NA)
)

# 打印核心结果
print(final_table[, c("Method", "P_Value", "Permutation_SD")])

cat("\n--- Type I Error Rates ---\n")
print(results$Type_I_Error_Rates)
