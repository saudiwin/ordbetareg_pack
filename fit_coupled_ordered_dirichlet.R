
# fit_coupled_ordered_dirichlet.R
#
# Example R script for fitting the coupled hurdle-inclusion ordered Dirichlet
# model with cmdstanr.
#
# Expected input:
#   y: N x K matrix of nonnegative values.
#      Rows may sum to 0, or may be positive compositions/count-like masses.
#      Nonzero rows are internally normalized in Stan.
#
#   X: N x P model matrix.
#
# Example:
#   source("fit_coupled_ordered_dirichlet.R")
#   fit <- fit_coupled_ordered_dirichlet(y, X)

library(cmdstanr)
library(posterior)
library(bayesplot)

fit_coupled_ordered_dirichlet <- function(
  y,
  X,
  stan_file = "coupled_ordered_dirichlet.stan",
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 123,
  adapt_delta = 0.95,
  max_treedepth = 12,
  eps = 1e-10
) {
  stopifnot(is.matrix(y))
  stopifnot(is.matrix(X))
  stopifnot(nrow(y) == nrow(X))
  stopifnot(all(y >= 0))

  N <- nrow(y)
  K <- ncol(y)
  P <- ncol(X)

  stan_data <- list(
    N = N,
    K = K,
    P = P,
    X = X,
    y = split(y, row(y)),
    eps = eps
  )

  # cmdstanr can accept array[N] vector[K] as a list of numeric vectors.
  stan_data$y <- lapply(seq_len(N), function(i) as.numeric(y[i, ]))

  mod <- cmdstan_model(stan_file)

  fit <- mod$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    refresh = 100
  )

  fit
}

# -------------------------------------------------------------------
# Helper: simulate from the coupled ordered Dirichlet model
# -------------------------------------------------------------------

rdirichlet1 <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}

simulate_coupled_ordered_dirichlet <- function(
  X,
  beta,
  delta,
  tau = 0,
  a_rho = 1,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  X <- as.matrix(X)
  beta <- as.matrix(beta)

  N <- nrow(X)
  K <- ncol(beta)

  y <- matrix(0, N, K)

  for (i in seq_len(N)) {
    theta <- as.numeric(X[i, ] %*% beta)
    mu <- exp(theta - max(theta))
    mu <- mu / sum(mu)

    rho <- plogis(a_rho * theta - tau)

    z <- rbinom(K, size = 1, prob = rho)
    A <- which(z == 1)

    if (length(A) == 0) {
      y[i, ] <- rep(0, K)
    } else if (length(A) == 1) {
      y[i, A] <- 1
    } else {
      phi <- exp(sum(X[i, ] * delta))
      mu_A <- mu[A] / sum(mu[A])
      alpha_A <- phi * mu_A
      y[i, A] <- rdirichlet1(alpha_A)
    }
  }

  y
}

# -------------------------------------------------------------------
# Minimal reproducible example
# -------------------------------------------------------------------

example_fit <- function() {
  set.seed(123)

  N <- 300
  K <- 4
  P <- 3

  X <- cbind(
    intercept = 1,
    x1 = rnorm(N),
    x2 = rnorm(N)
  )

  # Baseline category K has zero coefficients.
  beta <- matrix(0, P, K)
  beta[, 1] <- c(0.5,  0.8, -0.3)
  beta[, 2] <- c(-0.2, -0.5, 0.4)
  beta[, 3] <- c(0.1,  0.2, 0.7)
  beta[, 4] <- c(0, 0, 0)

  delta <- c(log(12), 0.2, -0.1)

  y <- simulate_coupled_ordered_dirichlet(
    X = X,
    beta = beta,
    delta = delta,
    tau = 0.3,
    a_rho = 1.2,
    seed = 456
  )

  cat("Share all-zero rows:", mean(rowSums(y) == 0), "\\n")
  cat("Share singleton rows:", mean(rowSums(y > 0) == 1), "\\n")
  cat("Share multi-active rows:", mean(rowSums(y > 0) >= 2), "\\n")

  fit <- fit_coupled_ordered_dirichlet(
    y = y,
    X = X,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 500,
    seed = 789
  )

  print(fit$summary(c("tau", "a_rho", "delta")))

  fit
}

# Uncomment to run:
# fit <- example_fit()
