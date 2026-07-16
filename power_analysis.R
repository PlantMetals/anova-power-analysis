{
  {
    cat("--- ANOVA Sample Size & Power Dynamic Workstation ---\n\n")
    a <- as.numeric(readline(prompt = "1. How many separate treatments (a) are you planning to have in your ANOVA? "))
    delta <- as.numeric(readline(prompt = "2. What is the minimum difference between the most extreme means that you expect to detect? "))
    ms_within <- as.numeric(readline(prompt = "3. What is the estimated MS within (residual variance/MS error) from your preliminary experiments? "))
    alpha <- as.numeric(readline(prompt = "4. Enter your chosen Type-I error rate threshold (alpha, significance level, e.g. 0.05): "))
    target_power <- 0.80
  }
  v1 <- a - 1
  for (n in 2:1000) {
    v2 <- a * (n - 1)
    phi <- sqrt((n * delta^2) / (2 * a * ms_within))
    nc_param <- a * phi^2  # Noncentrality parameter lambda
    
    f_crit <- qf(1 - alpha, v1, v2)
    power <- 1 - pf(f_crit, v1, v2, nc_param)
    
    if (power >= target_power) {
      print(paste("Required replicates per group (n):", n))
      print(paste("Resulting phi value:", round(phi, 3)))
      print(paste("Exact Power:", round(power, 3)))
      break
    }
  }
}