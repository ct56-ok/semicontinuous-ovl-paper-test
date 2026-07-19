library(readxl)

# ================= 1. 核心数学公式函数化 =================
calc_overlap_metrics <- function(p1, p2, mu1, mu2, sig1, sig2) {
  eps <- 1e-12
  
  # 安全检查：防止出现 NA 或方差为 0 导致报错
  if (any(is.na(c(p1, p2, mu1, mu2, sig1, sig2))) || sig1 == 0 || sig2 == 0) {
    return(c(rho = NA_real_, lambda = NA_real_, delta = NA_real_))
  }
  
  # Rho 计算
  rho <- min(1 - p1, 1 - p2) + sqrt(p1 * p2) * 
    (sqrt(2 * sqrt(sig1) * sqrt(sig2) / (sig1 + sig2)) * 
       exp(-((mu1 - mu2)^2) / (4 * (sig1 + sig2))))
  
  # Lambda 计算
  lambda <- min(1 - p1, 1 - p2) + 
    2 * sqrt(2) * (p1 * p2) / sqrt(sig1 + sig2) * 
    exp((sig2 * sig1 - (mu1 - mu2)^2 - 2 * sig2 * mu1 - 2 * sig1 * mu2) / 
          (2 * (sig2 + sig1))) / 
    ((p1 * exp((sig1 - 4 * mu1) / 4) / sqrt(sig1)) + 
       (p2 * exp((sig2 - 4 * mu2) / 4) / sqrt(sig2)))
  
  # Delta 计算
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

# ================= 2. 主置换检验函数 =================
permutation_test <- function(data1, data2, n_perm = 3000, seed = 222) {
  set.seed(seed)
  
  # 预处理数据
  bb <- data1 * 1000 / data1
  dd <- data2 * 1000 / data2
  
  bb1 <- log(bb[bb != 0]); bb2 <- bb[bb == 0]; bb3 <- c(bb2, bb1)
  dd1 <- log(dd[dd != 0]); dd2 <- dd[dd == 0]; dd3 <- c(dd2, dd1)
  bd <- c(bb3, dd3)
  
  n1 <- length(bb); n2 <- length(dd)
  n1_0 <- length(bb1); n2_0 <- length(dd1)
  hp1 <- n1_0 / n1; hp2 <- n2_0 / n2
  hmu1 <- mean(bb1); hmu2 <- mean(dd1)
  hsigma1 <- sum((bb1 - hmu1)^2) / n1_0
  hsigma2 <- sum((dd1 - hmu2)^2) / n2_0
  
  # 计算观测值
  obs_metrics <- calc_overlap_metrics(hp1, hp2, hmu1, hmu2, hsigma1, hsigma2)
  
  # 初始化置换结果存储
  perm_rho <- numeric(n_perm)
  perm_lambda <- numeric(n_perm)
  perm_delta <- numeric(n_perm)
  LR1 <- numeric(n_perm)
  Wa <- numeric(n_perm)
  
  # 预计算 LR 和 Wald 需要的全局统计量
  n1_1 <- n1 - n1_0
  n2_1 <- n2 - n2_0
  pp <- (n1_0 + n2_0) / (n1 + n2)
  uu <- mean(c(bb1, dd1))
  ww <- (sum((bb1 - uu)^2) + sum((dd1 - uu)^2)) / (n1_0 + n2_0)
  
  # 置换循环
  for (k in 1:n_perm) {
    idx1 <- sample.int(length(bd), n1, replace = FALSE)
    Z1 <- bd[idx1]
    Z2 <- bd[-idx1]
    
    Z11 <- Z1[Z1 != 0]
    Z22 <- Z2[Z2 != 0]
    n11_0 <- length(Z11)
    n22_0 <- length(Z22)
    
    if (n11_0 < 2 || n22_0 < 2) next
    
    hp11 <- n11_0 / n1
    hp22 <- n22_0 / n2
    hmu11 <- mean(Z11)
    hmu22 <- mean(Z22)
    hsigma11 <- sum((Z11 - hmu11)^2) / n11_0
    hsigma22 <- sum((Z22 - hmu22)^2) / n22_0
    
    # 调用函数计算 rho, lambda, delta
    pm <- calc_overlap_metrics(hp11, hp22, hmu11, hmu22, hsigma11, hsigma22)
    perm_rho[k] <- pm["rho"]
    perm_lambda[k] <- pm["lambda"]
    perm_delta[k] <- pm["delta"]
    
    # LR 检验
    T2 <- (pp / hp1)^n1_0 * ((1 - pp) / (1 - hp1))^n1_1 *
      ((1 - pp) / (1 - hp2))^n2_1 * (pp / hp2)^n2_0 *
      (hsigma1 / ww)^(n1_0 / 2) * (hsigma2 / ww)^(n2_0 / 2)
    LR1[k] <- if (is.finite(T2) && T2 > 0) -2 * log(T2) else NA
    
    # Wald 检验
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
  
  # 汇总结果
  valid_rho <- perm_rho[!is.na(perm_rho)]
  valid_lambda <- perm_lambda[!is.na(perm_lambda)]
  valid_delta <- perm_delta[!is.na(perm_delta)]
  valid_LR <- LR1[!is.na(LR1)]
  valid_Wa <- Wa[!is.na(Wa)]
  
  # 计算 P 值
  p_rho <- mean(valid_rho <= obs_metrics["rho"])
  p_lambda <- mean(valid_lambda <= obs_metrics["lambda"])
  p_delta <- mean(valid_delta <= obs_metrics["delta"])
  
  lr_stat <- if (length(valid_LR) > 0) valid_LR[1] else NA
  wa_stat <- if (length(valid_Wa) > 0) valid_Wa[1] else NA
  p_lr <- if (!is.na(lr_stat)) 1 - pchisq(lr_stat, 3) else NA
  p_wa <- if (!is.na(wa_stat)) 1 - pchisq(wa_stat, 3) else NA
  
  # 计算 Type I Error Rates
  TYLR <- if (length(valid_LR) > 0) mean(valid_LR > qchisq(0.95, 3)) else NA
  TYWa <- if (length(valid_Wa) > 0) mean(valid_Wa > qchisq(0.95, 3)) else NA
  TYR <- if (length(valid_rho) > 0) mean(valid_rho > obs_metrics["rho"]) else NA
  TYL <- if (length(valid_lambda) > 0) mean(valid_lambda > obs_metrics["lambda"]) else NA
  TYD <- if (length(valid_delta) > 0) mean(valid_delta > obs_metrics["delta"]) else NA
  
  # 安全计算标准差
  safe_sd <- function(x) if (length(x) >= 2) sd(x) else NA_real_
  
  return(list(
    Observed_Metrics = obs_metrics,
    P_Values = c(rho = p_rho, lambda = p_lambda, delta = p_delta, LR = p_lr, Wald = p_wa),
    Type_I_Error_Rates = c(Rho = TYR, Lambda = TYL, Delta = TYD, LR = TYLR, Wald = TYWa),
    Permutation_SD = c(
      rho = safe_sd(valid_rho), lambda = safe_sd(valid_lambda), delta = safe_sd(valid_delta),
      LR = safe_sd(valid_LR), Wald = safe_sd(valid_Wa)
    )
  ))
}

# ================= 3. 一键调用示例 =================
# data_dir <- "C:\\Users\\tingc\\Desktop\\ovl code\\应用数据"
# aa <- read_xlsx(file.path(data_dir, "末期141.xlsx"))
# cc <- read_xlsx(file.path(data_dir, "末期142.xlsx"))
# results <- permutation_test(aa$V_DOL * 1000 / aa$rtkcal, cc$V_DOL * 1000 / cc$rtkcal, n_perm = 3000)
# print(results)