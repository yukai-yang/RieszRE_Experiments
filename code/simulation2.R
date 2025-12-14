# Load required packages
library(MASS)
library(ggplot2)
library(tikzDevice)
library(reshape2)
library(dplyr)

# Function to generate block dependencies
generate_blocks <- function(n, d) {
  B <- floor(n^(1 - d))
  block_ids <- rep(1:B, each = ceiling(n/B))[1:n]
  block_ids <- sample(block_ids)
  return(block_ids)
}

# Function to simulate a single dataset (with signed dependence)
simulate_data <- function(n, d, sigma2 = 1) {
  block_ids <- generate_blocks(n, d)
  B <- max(block_ids)
  eta_b <- rnorm(B, mean = 0, sd = sqrt(sigma2))
  gamma <- sample(c(-1, 1), size = n, replace = TRUE)
  nu <- rnorm(n, mean = 0, sd = 1)
  
  alpha <- rnorm(n)
  beta <- 1 + rnorm(n)
  delta <- rnorm(n)
  x <- rnorm(n)
  z <- sample(c(0, 1), n, replace = TRUE)
  
  eps <- gamma * eta_b[block_ids] + nu
  y1 <- alpha + delta * x + eps + beta
  y0 <- alpha + delta * x + eps
  theta <- y1 - y0
  
  y = y1*z + y0*(1-z)
  
  # power
  #list(y = y, y1 = y1, y0 = y0, z = z, x = x, theta = 0, block_ids = block_ids)
  # size
  list(y = y, y1 = y1, y0 = y0, z = z, x = x, theta = theta, block_ids = block_ids)
}

# Compute Riesz estimator with plug-in variance (Equation 4.3)
compute_riesz_plugin <- function(y, z, theta, block_ids) {
  n <- length(y)
  mu1 <- mean(z)
  mu0 <- 1 - mu1
  psi <- (z / mu1) - ((1 - z) / mu0)
  tau_hat <- mean(y * psi)
  zeta <- y * psi - theta
  E_n <- which(outer(block_ids, block_ids, FUN = "=="), arr.ind = TRUE)
  nu <- rep(1 / n, n)
  sigma2_hat <- sum(nu[E_n[,1]] * nu[E_n[,2]] * zeta[E_n[,1]] * zeta[E_n[,2]])
  se_hat <- sqrt(sigma2_hat)
  return(list(tau_hat = tau_hat, se_hat = se_hat))
}

# Simulation function
run_simulation <- function(n, d, reps = 2000) {
  results <- data.frame(riesz = numeric(reps), riesz_se = numeric(reps), tau = numeric(reps))
  for (r in 1:reps) {
    data <- simulate_data(n, d)
    r_est <- compute_riesz_plugin(data$y, data$z, data$theta, data$block_ids)
    results[r, ] <- c(r_est$tau_hat, r_est$se_hat, mean(data$theta))
  }
  results$riesz_std <- (results$riesz - results$tau) / results$riesz_se
  results$coverage_riesz <- as.numeric(abs(results$riesz - results$tau) <= 1.96 * results$riesz_se)
  return(results)
}

# Parameters
d_values <- c(0.0, 0.1, 0.2, 0.25, 0.3)
n_values <- c(100, 200, 500, 1000)
reps <- 2000

set.seed(123)

all_results <- list()
for (n in n_values) {
  for (d in d_values) {
    cat("Running simulation for n =", n, "and d =", d, "\n")
    res <- run_simulation(n, d, reps)
    res$n <- n
    res$d <- d
    all_results[[paste(n, d, sep = "_")]] <- res
  }
}
combined <- bind_rows(all_results)

# Coverage plot
coverage_data <- combined %>%
  group_by(n, d) %>%
  summarise(coverage = mean(coverage_riesz), .groups = "drop")

# numbers
coverage_table <- combined %>%
  group_by(n, d) %>%
  summarise(Coverage_95pct = mean(coverage_riesz), .groups = "drop") %>%
  arrange(n, d) %>%
  mutate(Coverage_95pct = round(Coverage_95pct, 4))
print(coverage_table)

tikz("coverage_multiN.tex", width = 5.8, height = 3.2)

ggplot(coverage_data, aes(x = factor(d), y = coverage, fill = factor(n))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), colour = "black", width = 0.6) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "black") +
  scale_fill_manual(values = c("gray85", "gray65", "gray45", "gray25"), name = "$n=$") +
  labs(x = "Dependency growth rate $d$", y = NULL) +
  coord_cartesian(ylim = c(0.92, 0.98)) +
  theme_minimal(base_size = 11) +
  theme(axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.position = "top",
        legend.title = element_text(size = 10),
        plot.margin = margin(10, 10, 10, 10))

dev.off()

# Bias and RMSE
bias_rmse <- combined %>%
  group_by(n, d) %>%
  summarise(tau_mean = mean(tau), riesz_mean = mean(riesz), riesz_sd = sd(riesz),
            bias = mean(riesz) - mean(tau),
            rmse = sqrt(mean((riesz - tau)^2)), .groups = "drop") %>%
  select(n, d, bias, rmse)
bias_rmse_long <- melt(bias_rmse, id.vars = c("n", "d"))

tikz("bias_rmse_multiN.tex", width = 5.8, height = 3.2)

ggplot(bias_rmse_long, aes(x = factor(d), y = value, group = interaction(variable, n), linetype = factor(n))) +
  geom_line(aes(colour = variable), linewidth = 0.8) +
  geom_point(aes(colour = variable), shape = 16, size = 2) +
  scale_colour_manual(values = c("black", "gray40"), labels = c("Bias", "RMSE")) +
  labs(x = "Dependency growth rate $d$", y = NULL, linetype = "$n=$", colour = NULL) +
  coord_cartesian(ylim = c(-0.03, 0.38)) +
  theme_minimal(base_size = 11) +
  theme(axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        plot.margin = margin(10, 10, 10, 10))

dev.off()


# Rejection rates
combined$rej_1 <- as.numeric(abs(combined$riesz_std) > qnorm(0.995))
combined$rej_5 <- as.numeric(abs(combined$riesz_std) > qnorm(0.975))
combined$rej_10 <- as.numeric(abs(combined$riesz_std) > qnorm(0.95))
rejection_table <- combined %>%
  group_by(n, d) %>%
  summarise(Rej_1pct = mean(rej_1), Rej_5pct = mean(rej_5), Rej_10pct = mean(rej_10), .groups = "drop") %>%
  arrange(n, d) %>%
  mutate(across(starts_with("Rej_"), ~ round(.x, 4)))
print(rejection_table)
