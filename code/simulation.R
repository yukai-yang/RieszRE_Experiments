# Load required packages
library(MASS)
library(ggplot2)
library(tikzDevice)

# Function to generate block dependencies
generate_blocks <- function(n, d) {
  B <- floor(n^(1 - d))  # number of blocks
  block_ids <- sample(rep(1:B, length.out = n))  # assign each unit to a block
  return(block_ids)
}

# Function to simulate a single dataset
simulate_data <- function(n, d, sigma2 = 1) {
  block_ids <- generate_blocks(n, d)
  B <- max(block_ids)
  eta_b <- rnorm(B, mean = 0, sd = sqrt(sigma2))
  
  alpha <- rnorm(n)
  beta <- rnorm(n)
  delta <- rnorm(n)
  x <- rnorm(n)
  z <- rbinom(n, size = 1, prob = 0.5)
  eps <- eta_b[block_ids]
  
  y <- alpha + beta * z + delta * z * x + eps
  y1 <- alpha + beta * 1 + delta * 1 * x + eta_b[block_ids]
  y0 <- alpha + beta * 0 + delta * 0 * x + eta_b[block_ids]
  tau <- mean(y1 - y0)
  
  list(y = y, z = z, x = x, tau = tau, block_ids = block_ids)
}

# Function to compute Riesz estimator
compute_riesz <- function(y, z) {
  mu1 <- mean(z)
  mu0 <- 1 - mu1
  psi <- (z / mu1) - ((1 - z) / mu0)
  tau_hat <- mean(y * psi)
  se_hat <- sd(y * psi) / sqrt(length(y))
  return(list(tau_hat = tau_hat, se_hat = se_hat))
}

# Simulation loop
run_simulation <- function(n, d, reps = 2000) {
  results <- data.frame(riesz = numeric(reps), riesz_se = numeric(reps), tau = numeric(reps))
  
  for (r in 1:reps) {
    data <- simulate_data(n, d)
    r_est <- compute_riesz(data$y, data$z)
    results[r, ] <- c(r_est$tau_hat, r_est$se_hat, data$tau)
  }
  
  results$riesz_std <- (results$riesz - results$tau) / results$riesz_se
  results$coverage_riesz <- as.numeric(abs(results$riesz - results$tau) <= 1.96 * results$riesz_se)
  return(results)
}

# Run simulations for various d
d_values <- c(0.0, 0.1, 0.2, 0.25, 0.3)
n <- 500
all_results <- list()

# For reproducibility
set.seed(123)

for (d in d_values) {
  cat("Running simulation for d =", d, "\\n")
  res <- run_simulation(n, d)
  res$d <- d
  all_results[[as.character(d)]] <- res
}

combined <- do.call(rbind, all_results)

# Plot with TikZ output
tikz("coverage_vs_d.tex", width = 5.5, height = 3.5)

ggplot(combined, aes(x = factor(d), y = coverage_riesz)) +
  stat_summary(fun = mean, geom = "bar", fill = "gray80", colour = "black", width = 0.6) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "black") +
  coord_cartesian(ylim = c(0.90, 1.00)) +
  labs(x = "Dependency growth rate $d$",
       y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(10, 10, 10, 10)
  )

dev.off()


# Bias and RMSE summary
summary_stats <- aggregate(cbind(riesz, tau) ~ d, data = combined, function(x) c(mean = mean(x), sd = sd(x)))
summary_stats <- do.call(data.frame, summary_stats)
summary_stats$bias_riesz <- summary_stats$riesz.mean - summary_stats$tau.mean
summary_stats$rmse_riesz <- sqrt(summary_stats$bias_riesz^2 + summary_stats$riesz.sd^2)

# Export to TikZ
tikz("bias_rmse_vs_d.tex", width = 5.5, height = 3.5)

ggplot(summary_stats, aes(x = d)) +
  geom_line(aes(y = bias_riesz, linetype = "Bias"), colour = "black", linewidth = 0.8) +
  geom_point(aes(y = bias_riesz), colour = "black", shape = 16, size = 3) +
  geom_line(aes(y = rmse_riesz, linetype = "RMSE"), colour = "gray40", linewidth = 0.8) +
  geom_point(aes(y = rmse_riesz), colour = "gray40", shape = 16, size = 3) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_x_continuous(breaks = summary_stats$d) +  # Ensure all d values are shown
  coord_cartesian(ylim = c(-0.01, 0.06)) +
  labs(x = "Dependency growth rate $d$", y = NULL, linetype = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(10, 10, 10, 10)
  )

dev.off()

