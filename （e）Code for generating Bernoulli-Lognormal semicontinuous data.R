# 1. Define a function to generate Bernoulli-Lognormal (BNT) semicontinuous data
generate_bnt_data <- function(n, prob_nonzero, mu_log, sigma_log, seed = NULL) {
  # Set random seed to ensure reproducibility
  if (!is.null(seed)) set.seed(seed)
  
  # Bernoulli process: Generate an indicator variable (1 for non-zero, 0 for zero)
  is_nonzero <- rbinom(n, size = 1, prob = prob_nonzero)
  
  # Lognormal process: Generate continuous values for the non-zero part
  # Note: meanlog and sdlog correspond to the mean and standard deviation of the log-transformed values
  continuous_values <- rlnorm(n, meanlog = mu_log, sdlog = sigma_log)
  
  # Combine both parts: Force the values to 0 where the indicator variable is 0
  semicontinuous_data <- is_nonzero * continuous_values
  
  return(semicontinuous_data)
}

# 2. Set the sample size
sample_size <- 1000

# 3. Generate two groups of dietary data
# Group 1: Daily staple intake (fewer zeros, higher overall intake level)
group1_data <- generate_bnt_data(
  n = sample_size, 
  prob_nonzero = 0.80, 
  mu_log = 3.0, 
  sigma_log = 0.5, 
  seed = 42
)

# Group 2: Occasional snack intake (many zeros, lower overall intake level)
group2_data <- generate_bnt_data(
  n = sample_size, 
  prob_nonzero = 0.40, 
  mu_log = 1.5, 
  sigma_log = 0.8, 
  seed = 123
)

# 4. Combine into a data frame and view summary statistics
df_dietary <- data.frame(
  Group_1_Intake = group1_data,
  Group_2_Intake = group2_data
)

cat("Data generation complete! Preview of the first 10 rows:\n")
print(head(df_dietary, 10))

cat("\nBasic summary statistics:\n")
print(summary(df_dietary))

# 5. Basic visualization: Plot histograms for both groups
par(mfrow = c(1, 2)) # Split the canvas into 1 row and 2 columns
hist(group1_data, main = "Group 1: Daily Staple Intake", 
     xlab = "Intake Level", col = "steelblue", border = "white")
hist(group2_data, main = "Group 2: Occasional Snack Intake", 
     xlab = "Intake Level", col = "coral", border = "white")