# Task 3: Fit the Poisson changepoint model using BFGS.
# Inputs: y contains observed counts; lambda controls shrinkage;
# eps smooths the LASSO penalty; start provides initial coefficients.
# Initialise at a constant mean if start is not supplied.
# Return the optim result and the number of fitting attempts.

fit_poisson_lasso <- function(y, lambda, eps = 1e-6, start = NULL) {
  
  # Check that y contains finite, non-negative integer counts.
  # A positive mean is required for the log-scale initial value.
  if (!is.numeric(y) || length(y) < 1 || any(!is.finite(y)) ||
      any(y < 0) || any(y != floor(y)) || mean(y) <= 0) {
    stop("y must contain finite counts with a positive sample mean.")
  }
  
  # Check that lambda is non-negative and eps is strictly positive.
  # Both must be finite scalar values.
  if (length(lambda) != 1 || !is.finite(lambda) || lambda < 0 ||
      length(eps) != 1 || !is.finite(eps) || eps <= 0) {
    stop("Check lambda and eps.")
  }
  
  # Use the sample log-mean as the unpenalised baseline.
  # Set all log-mean increments to zero for a constant initial mean.
  if (is.null(start)) {
    start <- c(log(mean(y)), rep(0, length(y) - 1))
  }
  
  # Each observation requires one coefficient in the model.
  if (length(start) != length(y) || any(!is.finite(start))) {
    stop("start must have one finite value per observation.")
  }
  
  # Fit using BFGS with the analytical gradient from Task 2.
  # Retry a non-converged fit with more iterations if necessary.
  for (attempt in seq_len(2)) {
    
    fit <- optim(
      par = start,
      fn = poission_lasso_objective,
      gr = poission_lasso_grad,
      y = y,
      lambda = lambda,
      eps = eps,
      method = "BFGS",
      control = list(maxit = 2000 * attempt, reltol = 1e-10)
    )
    
    # Accept only a converged fit with finite estimates and cost.
    valid <- fit$convergence == 0 && is.finite(fit$value) &&
      all(is.finite(fit$par))
    
    if (valid) {
      fit$attempts <- attempt
      return(fit)
    }
    
    # Reuse the previous estimates when another attempt is needed.
    if (all(is.finite(fit$par))) {
      start <- fit$par
    }
  }
  
  # Stop rather than using unreliable coefficients after failure.
  stop(sprintf("BFGS failed for lambda = %.6g; code %d.",
               lambda, fit$convergence))
}


# Task 4: Evaluate the fitted model for a given penalty.
# Inputs: beta contains fitted coefficients; y contains counts;
# lambda is the penalty value and eps controls penalty smoothing.
# Compute EDF from the penalised and unpenalised Hessians.
# Return EDF, effective changepoint count, log-likelihood and mBIC.

model_score <- function(beta, y, lambda, eps = 1e-6) {
  
  # Evaluate both Hessians at the same fitted coefficients.
  # H_zero excludes the penalty curvature without refitting.
  H_lambda <- poission_lasso_hess(beta, y, lambda, eps)
  H_zero <- poission_lasso_hess(beta, y, lambda = 0, eps = eps)
  
  # EDF = trace(H_lambda^(-1) %*% H_zero).
  # solve(A, B) computes the matrix solution without explicit inversion.
  edf <- sum(diag(solve(H_lambda, H_zero)))
  
  # The unpenalised baseline contributes at least one degree of freedom.
  # The total EDF cannot exceed the number of observations.
  if (!is.finite(edf) || edf < 1 - 1e-6 ||
      edf > length(y) + 1e-6) {
    stop("EDF is outside its expected range; check the Hessians.")
  }
  
  # Remove the baseline degree of freedom to count effective changes.
  # Keep this value fractional when evaluating the information criterion.
  m_eff <- edf - 1
  
  # Recover the Poisson means from the cumulative log-mean increments.
  # Evaluate the unpenalised log-likelihood at the fitted coefficients.
  mu <- exp(cumsum(beta))
  loglik <- sum(dpois(y, lambda = mu, log = TRUE))
  
  # Modified BIC penalises effective changepoints and poor model fit.
  # The preferred lambda is the one producing the smallest mBIC.
  mbic <- 3 * m_eff * log(length(y)) - 2 * loglik
  
  return(c(edf = edf, m_eff = m_eff,
           loglik = loglik, mbic = mbic))
}


# Check that the data and functions from Tasks 1 and 2 are available.
# Run the preceding group member's code in the same R session.
if (!exists("y") || !exists("time") || !exists("mu_true") ||
    !exists("eps") ||
    !exists("poission_lasso_objective", mode = "function") ||
    !exists("poission_lasso_grad", mode = "function") ||
    !exists("poission_lasso_hess", mode = "function")) {
  stop("Run the preceding code before fitting the model.")
}

# Record elapsed time for fitting and penalty selection.
selection_start <- proc.time()[["elapsed"]]


# Task 3: Fit the model with the default penalty lambda = 1.
# Initialise the baseline at log(mean(y)) with no initial changes.

beta_start <- c(log(mean(y)), rep(0, length(y) - 1))
fit_initial <- fit_poisson_lasso(
  y = y, lambda = 1, eps = eps, start = beta_start
)

# Convert the fitted cumulative log-means to Poisson means.
mu_initial <- exp(cumsum(fit_initial$par))

# Compare the initial fit with the data and the true mean function.
plot(time, y, pch = 16, col = "grey65",
     ylim = range(y, mu_true, mu_initial),
     xlab = "Time", ylab = "Count / mean",
     main = "Initial fit: lambda = 1")

lines(time, mu_true, type = "s", lwd = 2)
lines(time, mu_initial, type = "s", col = "steelblue", lwd = 2)

legend("topright",
       c("Observed counts", "True mean", "Initial fit"),
       col = c("grey65", "black", "steelblue"),
       pch = c(16, NA, NA),
       lty = c(NA, 1, 1),
       bty = "n", cex = 0.8)


# Task 4: Select lambda by minimising the modified BIC.
# Search the specified grid of log-lambda values from -1 to 5.
# Exponentiate the grid because the objective uses lambda itself.

lambda_grid <- exp(seq(-1, 5, by = 0.05))

# Store model scores and convergence information for each lambda.
grid_results <- data.frame(
  lambda = lambda_grid,
  edf = NA_real_,
  m_eff = NA_real_,
  loglik = NA_real_,
  mbic = NA_real_,
  convergence = NA_integer_,
  attempts = NA_integer_
)

# Store the fitted coefficients for each penalty value.
# Each column corresponds to one lambda in the grid.
beta_grid <- matrix(
  NA_real_, nrow = length(y), ncol = length(lambda_grid)
)

# Begin with the constant-mean initial coefficients.
start <- beta_start

# Fit the model and evaluate mBIC at every penalty value.
for (i in seq_along(lambda_grid)) {
  
  # Fit using the current penalty and previous starting values.
  fit <- fit_poisson_lasso(
    y = y, lambda = lambda_grid[i], eps = eps, start = start
  )
  
  # Evaluate EDF, effective changes, log-likelihood and mBIC.
  score <- model_score(fit$par, y, lambda_grid[i], eps)
  
  # Save the model scores and optimisation diagnostics.
  grid_results[i, c("edf", "m_eff", "loglik", "mbic")] <- as.list(score)
  grid_results$convergence[i] <- fit$convergence
  grid_results$attempts[i] <- fit$attempts
  
  # Preserve the fitted coefficients for subsequent model selection.
  beta_grid[, i] <- fit$par
  
  # Warm-start the next fit to improve optimisation efficiency.
  start <- fit$par
}

# Select the penalty giving the lowest mBIC over the grid.
best_index <- which.min(grid_results$mbic)
lambda_best <- lambda_grid[best_index]


# Task 5: Refit the model using the selected penalty.
# Initialise from the corresponding grid fit to aid convergence.

fit_best <- fit_poisson_lasso(
  y = y, lambda = lambda_best, eps = eps,
  start = beta_grid[, best_index]
)

# Retain the fitted coefficients and mean function for bootstrapping.
beta_best <- fit_best$par
mu_best <- exp(cumsum(beta_best))

# Recalculate the scores using the final fitted coefficients.
best_score <- model_score(beta_best, y, lambda_best, eps)
edf_best <- unname(best_score["edf"])
m_eff_best <- unname(best_score["m_eff"])

# Correct tiny negative values caused by numerical roundoff.
if (m_eff_best < 0 && m_eff_best > -1e-6) {
  m_eff_best <- 0
}

# Estimate the number of changes by flooring the effective count.
m_hat <- floor(m_eff_best)

# Rank the log-mean increments by absolute magnitude.
# Exclude beta[1] because it represents the unpenalised baseline.
ranked_changes <- order(abs(beta_best[-1]), decreasing = TRUE)

# Select the largest m_hat changes and recover their time indices.
# Adding one accounts for removing beta[1] before ranking.
changepoints <- sort(ranked_changes[seq_len(m_hat)] + 1L)

# Compare the selected mean function with the observed and true data.
plot(time, y, pch = 16, col = "grey65",
     ylim = range(y, mu_true, mu_best),
     xlab = "Time", ylab = "Count / mean",
     main = "Selected Poisson model")

lines(time, mu_true, type = "s", lwd = 2)
lines(time, mu_best, type = "s", col = "steelblue", lwd = 2)

# Mark the estimated changepoint locations with dashed lines.
if (length(changepoints) > 0) {
  abline(v = changepoints, col = "steelblue", lty = 3)
}

legend("topright",
       c("Observed counts", "True mean",
         "Selected fit", "Estimated changes"),
       col = c("grey65", "black", "steelblue", "steelblue"),
       pch = c(16, NA, NA, NA),
       lty = c(NA, 1, 1, 3),
       bty = "n", cex = 0.8)

# Summarise the selected penalty and effective model complexity.
selected_summary <- data.frame(
  lambda = lambda_best,
  edf = edf_best,
  effective_changes = m_eff_best,
  estimated_changes = m_hat,
  mbic = unname(best_score["mbic"])
)

print(selected_summary, digits = 5, row.names = FALSE)

# Report the estimated changepoints and grid convergence results.
cat("Estimated changepoints:", changepoints, "\n")
cat("Converged grid fits:",
    sum(grid_results$convergence == 0),
    "of", length(lambda_grid), "\n")

# Report total elapsed time for fitting and penalty selection.
selection_elapsed <- proc.time()[["elapsed"]] - selection_start
cat(sprintf("Fitting and selection: %.2f seconds\n",
            selection_elapsed))
