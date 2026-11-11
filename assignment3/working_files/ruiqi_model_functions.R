# Poisson model: cost, derivatives and simulated data.
# Run this part first. Tong's code uses the functions and data below.

poission_lasso_objective <- function(beta, y, lambda = 1, eps = 1e-6) {
  # Returns the negative log-likelihood plus the smooth LASSO penalty.
  # beta contains log-mean increments, one per count in y.
  # lambda controls shrinkage; eps smooths the absolute value near zero.
  # Eq. (6): C(beta) = -log L(beta; y)
  #                   + lambda * sum(sqrt(beta[k]^2 + eps)), k = 2, ..., n.
  # mu[k] = exp(sum(beta[1:k])); beta[1] is the unpenalised baseline.
  
  eta <- cumsum(beta)
  mu <- exp(eta)
  
  # An overflowing trial step should not be accepted by optim.
  if (any(!is.finite(mu))) return(Inf)
  
  negative_loglik <- sum(mu - y * eta + lgamma(y + 1))
  penalty <- lambda * sum(sqrt(beta[-1]^2 + eps))
  return(negative_loglik + penalty)
}

poission_lasso_grad <- function(beta, y, lambda = 1, eps = 1e-6) {
  # Gradient of the cost, with one entry per coefficient.
  # Uses the same beta, counts y, penalty lambda and smoothing eps.
  
  mu <- exp(cumsum(beta))
  
  # d(-log L)/d beta[j] = sum(mu[k] - y[k]) over k >= j.
  # rev(cumsum(rev())) vectorised tail-sum.
  gradient <- rev(cumsum(rev(mu - y)))
  # For j >= 2, add lambda * beta[j] / sqrt(beta[j]^2 + eps).
  penalty_gradient <- lambda * beta / sqrt(beta^2 + eps)
  penalty_gradient[1] <- 0
  return(gradient + penalty_gradient)
}

poission_lasso_hess <- function(beta, y, lambda = 1, eps = 1e-6) {
  # Returns the n by n Hessian, where n is the length of beta.
  # beta, lambda and eps have the same meaning as in the cost function.
  # y is kept as an argument for consistency, but drops out on differentiating.
  
  mu <- exp(cumsum(beta))
  n <- length(beta)
  tail_mu <- rev(cumsum(rev(mu)))
  
  # H0[i, j] = sum(mu[k]) over k >= max(i, j).
  first_shared <- outer(seq_len(n), seq_len(n), pmax)
  hessian <- matrix(tail_mu[first_shared], nrow = n, ncol = n)
  
  # Add lambda * eps / (beta[j]^2 + eps)^(3/2) to diagonal j, for j >= 2.
  curvature <- lambda * eps / (beta^2 + eps)^(3 / 2)
  curvature[1] <- 0
  diag(hessian) <- diag(hessian) + curvature
  return(hessian)
}

check_derivatives <- function(beta, y, lambda, eps = 1e-6, h = 1e-6) {
  # Check the derivatives at beta for counts y, using lambda and eps.
  # g[j] approx (C(beta + h*e_j) - C(beta - h*e_j)) / (2*h),
  # where e_j is one at position j and zero elsewhere.
  # Apply the same difference to the gradient for Hessian column j.
  # Returns scaled errors for both derivatives and the Hessian asymmetry.
  
  n <- length(beta)
  gradient_fd <- numeric(n)
  hessian_fd <- matrix(0, n, n)
  
  for (j in seq_len(n)) {
    beta_plus <- beta_minus <- beta
    beta_plus[j] <- beta_plus[j] + h
    beta_minus[j] <- beta_minus[j] - h
    
    gradient_fd[j] <- (
      poission_lasso_objective(beta_plus, y, lambda, eps) -
        poission_lasso_objective(beta_minus, y, lambda, eps)
    ) / (2 * h)
    
    hessian_fd[, j] <- (
      poission_lasso_grad(beta_plus, y, lambda, eps) -
        poission_lasso_grad(beta_minus, y, lambda, eps)
    ) / (2 * h)
  }
  
  gradient <- poission_lasso_grad(beta, y, lambda, eps)
  hessian <- poission_lasso_hess(beta, y, lambda, eps)
  gradient_error <- max(abs(gradient - gradient_fd) /
                          (1 + abs(gradient) + abs(gradient_fd)))
  hessian_error <- max(abs(hessian - hessian_fd) /
                         (1 + abs(hessian) + abs(hessian_fd)))
  
  return(c(gradient_error = gradient_error,
           hessian_error = hessian_error,
           symmetry_error = max(abs(hessian - t(hessian)))))
}

# Try both signs and values close to zero, where the penalty curves sharply.
# lambda = 0 checks the likelihood alone. Also compare two step sizes.
beta_check <- c(log(5), 0, 0.2, -0.1, 1e-4)
y_check <- c(4, 7, 6, 3, 5)
eps <- 1e-6
derivative_checks <- expand.grid(lambda = c(0, 1, 10), h = c(1e-5, 1e-6))
derivative_checks$gradient_error <- NA_real_
derivative_checks$hessian_error <- NA_real_
derivative_checks$symmetry_error <- NA_real_

for (i in seq_len(nrow(derivative_checks))) {
  errors <- check_derivatives(
    beta_check, y_check, derivative_checks$lambda[i], eps,
    h = derivative_checks$h[i]
  )
  derivative_checks[i, 3:5] <- as.list(errors)
}

# Errors are scaled because Hessian entries differ greatly in size.
stopifnot(
  max(derivative_checks$gradient_error) < 1e-4,
  max(derivative_checks$hessian_error) < 1e-4,
  max(derivative_checks$symmetry_error) < 1e-12
)
print(derivative_checks, digits = 3, row.names = FALSE)

# Use the seed and segment means from the practical.
set.seed(12345)
n <- 100
time <- seq_len(n)
tau_true <- c(1, 30, 70, n + 1)
m_true <- length(tau_true) - 2
theta_true <- c(10, 2, 7)
mu_true <- numeric(n)

for (i in seq_along(theta_true)) {
  mu_true[tau_true[i]:(tau_true[i + 1] - 1)] <- theta_true[i]
}
y <- rpois(n, lambda = mu_true)

plot(time, y, pch = 16, col = "grey60",
     xlab = "Time", ylab = "Count / mean", main = "Simulated Poisson data")
lines(time, mu_true, type = "s", lwd = 2)
abline(v = tau_true[-c(1, length(tau_true))], lty = 3)
legend("topright", c("Observed counts", "True mean", "True changepoints"),
       col = c("grey60", "black", "black"), pch = c(16, NA, NA),
       lty = c(NA, 1, 3), bty = "n", cex = 0.8)

