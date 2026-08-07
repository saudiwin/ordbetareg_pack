# Parameter recovery for Ordered Dirichlet – Single Shared Mean
# File: ordered_dirichlet_single_mean.stan
#
# For each observation i:
#   theta_i = beta' * X[i]          K-vector of category log-odds (reference K = 0)
#   phi_i   = exp(delta' * X[i])    precision, log-linear in predictors
#   y[i] ~ ordered_dirichlet_shared(theta_i, phi_i, tau)
#
# ordered_dirichlet_shared definition (inferred from model structure):
#   mu   = softmax(theta)           mean direction of Dirichlet
#   eta[k] = theta[k] - tau         per-dim activity logit; tau is the single
#                                   shared sparsity threshold
#   Condition on non-empty active set; active dims drawn Dirichlet(phi * mu_A).
#
# This is more parsimonious than the K-separate-eta Gumbel model: the same
# theta drives both mean direction and activity probabilities, with tau as a
# single global knob for overall sparsity.
#
# Requires: cmdstanr, posterior, ggplot2, dplyr, tidyr

library(cmdstanr)
library(posterior)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(42)

rdirichlet <- function(n, alpha) {
  k <- length(alpha)
  x <- matrix(rgamma(n * k, rep(alpha, each = n)), nrow = n, ncol = k)
  x / rowSums(x)
}

# ── True parameters ───────────────────────────────────────────────────────────

K <- 3L   # categories; category K is reference (beta[:, K] = 0)
P <- 2L   # predictors: intercept + 1 continuous covariate

# beta is P x K; last column fixed to 0 (reference)
beta_raw_true <- matrix(c( 0.5, -0.3,   # intercepts (row 1)
                            1.0, -0.5),  # slopes     (row 2)
                        nrow = P, ncol = K - 1L, byrow = TRUE)
beta_true <- cbind(beta_raw_true, 0)    # P x K, col K = 0

tau_true   <- 0.5
delta_true <- c(1.5, 0.3)   # intercept + slope for log(phi)

N <- 400L

# ── Simulate ──────────────────────────────────────────────────────────────────

X <- cbind(1, rnorm(N))    # N x P: intercept column + predictor

sim_shared <- function(N, X, beta, tau, delta) {
  K <- ncol(beta)
  Y <- matrix(0.0, nrow = N, ncol = K)
  for (i in seq_len(N)) {
    x_i     <- X[i, ]
    theta_i <- as.vector(t(beta) %*% x_i)      # K-vector
    phi_i   <- exp(sum(delta * x_i))            # scalar > 0
    mu_i    <- exp(theta_i) / sum(exp(theta_i)) # softmax → mean direction
    eta_i   <- theta_i - tau                    # single shared threshold
    rho_i   <- plogis(eta_i)
    repeat {
      active <- rbinom(K, 1L, rho_i) == 1L
      if (any(active)) break
    }
    idx <- which(active)
    if (length(idx) == 1L) {
      Y[i, idx] <- 1.0
    } else {
      mu_a      <- mu_i[idx] / sum(mu_i[idx])
      Y[i, idx] <- rdirichlet(1L, phi_i * mu_a)
    }
  }
  Y
}

Y <- sim_shared(N, X, beta_true, tau_true, delta_true)

cat("Simulated data summary:\n")
cat("Proportion of zeros per dimension:\n")
print(round(colMeans(Y == 0), 3))
cat("Proportion of 1s per dimension:\n")
print(round(colMeans(Y == 1), 3))
cat("Mean of interior observations per dimension:\n")
interior <- Y > 0 & Y < 1
for (k in seq_len(K))
  cat(sprintf("  dim %d: %.3f  (n=%d)\n", k,
              mean(Y[interior[, k], k]), sum(interior[, k])))

# ── Stan model ────────────────────────────────────────────────────────────────

stan_code <- '
functions {
  real ordered_dirichlet_shared_lpdf(
      vector y, vector theta, real phi, real tau) {
    int K = num_elements(y);
    vector[K] mu  = softmax(theta);    // mean direction
    vector[K] eta = theta - tau;       // single shared threshold
    int A_size = 0;
    array[K] int active;
    real lp_A     = 0;
    real log_empty = 0;
    for (k in 1:K) {
      real log_rho   = -log1p_exp(-eta[k]);
      real log1m_rho = -log1p_exp(eta[k]);
      log_empty += log1m_rho;
      if (y[k] > 0) {
        A_size       += 1;
        active[A_size] = k;
        lp_A         += log_rho;
      } else {
        lp_A += log1m_rho;
      }
    }
    lp_A -= log1m_exp(log_empty);     // condition on non-empty active set
    if (A_size == 1) {
      return lp_A;
    } else {
      vector[A_size] y_A;
      vector[A_size] alpha_A;
      real mu_sum = 0;
      for (a in 1:A_size) mu_sum += mu[active[a]];
      for (a in 1:A_size) {
        y_A[a]     = y[active[a]];
        alpha_A[a] = phi * mu[active[a]] / mu_sum;
      }
      return lp_A + dirichlet_lpdf(y_A | alpha_A);
    }
  }
}
data {
  int<lower=1> N;
  int<lower=2> K;
  int<lower=1> P;
  matrix[N, P]        X;
  array[N] vector[K]  y;
}
parameters {
  matrix[P, K - 1] beta_raw;
  real tau;
  vector[P] delta;
}
transformed parameters {
  matrix[P, K] beta;
  for (k in 1:(K - 1))
    beta[, k] = beta_raw[, k];
  beta[, K] = rep_vector(0, P);
}
model {
  to_vector(beta_raw) ~ normal(0, 1);
  tau   ~ normal(0, 2);
  delta ~ normal(0, 1);
  for (i in 1:N) {
    vector[K] theta_i = beta\' * to_vector(X[i]\');
    real      phi_i   = exp(dot_product(to_vector(X[i]\'), delta));
    y[i] ~ ordered_dirichlet_shared(theta_i, phi_i, tau);
  }
}
'

stan_file <- tempfile(fileext = ".stan")
writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

fit <- mod$sample(
  data = list(N = N, K = K, P = P, X = X, y = Y),
  chains          = 4L,
  parallel_chains = 4L,
  iter_warmup     = 1000L,
  iter_sampling   = 1000L,
  seed            = 42L,
  refresh         = 500L
)

# ── Parameter recovery ────────────────────────────────────────────────────────

# Build truth table matching Stan's variable naming convention
# beta_raw[p, k] in Stan → "beta_raw[p,k]" in posterior draws
truth <- bind_rows(
  tibble(
    param = sprintf("beta_raw[%d,%d]", rep(seq_len(P), K - 1L),
                    rep(seq_len(K - 1L), each = P)),
    truth = as.vector(beta_raw_true)    # column-major order matches Stan
  ),
  tibble(param = "tau",   truth = tau_true),
  tibble(param = sprintf("delta[%d]", seq_len(P)),
         truth = delta_true)
)

cat("\n── Recovery summary ─────────────────────────────────────────────\n")
fit$summary(c("beta_raw", "tau", "delta")) |>
  left_join(truth, by = c("variable" = "param")) |>
  select(variable, truth, mean, sd, q5, q95, rhat, ess_bulk) |>
  print(n = Inf)

# Histogram + true-value overlay
post <- fit$draws(c("beta_raw", "tau", "delta"), format = "df") |>
  select(-starts_with(".")) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  left_join(truth, by = "param")

p_hist <- ggplot(post, aes(x = value)) +
  geom_histogram(bins = 60, fill = "#1b7837", alpha = 0.75, colour = NA) +
  geom_vline(aes(xintercept = truth), colour = "red", linewidth = 0.9) +
  facet_wrap(~param, scales = "free") +
  labs(
    title    = "Ordered Dirichlet (Single Shared Mean) — Parameter Recovery",
    subtitle = "Red line = true value",
    x = NULL, y = "Count"
  ) +
  theme_minimal(base_size = 12)

print(p_hist)

# Interval plot
intervals <- fit$summary(c("beta_raw", "tau", "delta"),
                          ~quantile(.x, probs = c(0.05, 0.5, 0.95))) |>
  rename(param = variable, q05 = `5%`, med = `50%`, q95 = `95%`) |>
  left_join(truth, by = "param")

p_int <- ggplot(intervals, aes(y = param)) +
  geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 1.2,
                 colour = "#1b7837", alpha = 0.7) +
  geom_point(aes(x = med), size = 3, colour = "#1b7837") +
  geom_point(aes(x = truth), shape = 18, size = 4, colour = "red") +
  labs(
    title    = "Ordered Dirichlet (Single Shared Mean) — 90% Posterior Intervals",
    subtitle = "Green = posterior median + 90% CI  |  Red diamond = truth",
    x = "Value", y = NULL
  ) +
  theme_minimal(base_size = 12)

print(p_int)


# ── Posterior predictive check ────────────────────────────────────────────────
#
# For S posterior draws, simulate Y_rep using the same X as the observed data.
# Compare three summary statistics per dimension:
#   mean proportion, P(zero), P(one)

S_ppc    <- 100L
draws_ppc <- fit$draws(c("beta_raw", "tau", "delta"), format = "df") |>
  select(-starts_with(".")) |>
  slice_sample(n = S_ppc)

# Helper: extract (beta, tau, delta) from one draw row
extract_params <- function(row) {
  br_cols <- paste0("beta_raw[", rep(seq_len(P), K - 1L), ",",
                    rep(seq_len(K - 1L), each = P), "]")
  beta_s  <- cbind(
    matrix(as.numeric(row[br_cols]), nrow = P, ncol = K - 1L),
    0
  )
  list(
    beta  = beta_s,
    tau   = as.numeric(row["tau"]),
    delta = as.numeric(row[paste0("delta[", seq_len(P), "]")])
  )
}

ppc_rows <- vector("list", S_ppc)
for (s in seq_len(S_ppc)) {
  p      <- extract_params(draws_ppc[s, ])
  Y_rep  <- sim_shared(N, X, p$beta, p$tau, p$delta)
  ppc_rows[[s]] <- data.frame(
    draw  = s,
    dim   = paste0("dim ", rep(seq_len(K), 3L)),
    stat  = rep(c("mean proportion", "P(zero)", "P(one)"), each = K),
    value = c(colMeans(Y_rep),
              colMeans(Y_rep == 0),
              colMeans(Y_rep == 1))
  )
}
ppc_df <- do.call(rbind, ppc_rows)

obs_stats <- data.frame(
  dim   = paste0("dim ", rep(seq_len(K), 3L)),
  stat  = rep(c("mean proportion", "P(zero)", "P(one)"), each = K),
  value = c(colMeans(Y), colMeans(Y == 0), colMeans(Y == 1))
)

p_ppc <- ggplot(ppc_df, aes(x = value)) +
  geom_histogram(bins = 40, fill = "#762a83", alpha = 0.65, colour = NA) +
  geom_vline(data = obs_stats, aes(xintercept = value),
             colour = "orange", linewidth = 1) +
  facet_grid(stat ~ dim, scales = "free_x") +
  labs(
    title    = "Posterior Predictive Check",
    subtitle = "Purple = replicated statistic distribution  |  Orange = observed",
    x = NULL, y = "Count"
  ) +
  theme_minimal(base_size = 11)

print(p_ppc)


# ── Marginal effects on expected proportions ──────────────────────────────────
#
# E[y_k | x2] is computed analytically by enumerating all 2^K - 1 non-empty
# active subsets.  phi does not affect E[y_k] (only the Dirichlet variance).
#
# For each subset A:
#   P(A active | non-empty) = prod(rho[j], j in A) * prod(1-rho[j], j not in A)
#                             / (1 - prod(1-rho[j], all j))
#   E[y_k | A active, |A|=1] = 1
#   E[y_k | A active, |A|>=2] = mu[k] / sum(mu[A])   (Dirichlet mean)

expected_prop <- function(theta, tau) {
  K         <- length(theta)
  th_shift  <- theta - max(theta)                      # stable softmax
  mu        <- exp(th_shift) / sum(exp(th_shift))
  rho       <- plogis(theta - tau)
  log_norm  <- log1p(-exp(sum(log1p(-rho))))           # log(1 - prod(1-rho))
  E_y       <- numeric(K)
  for (mask in seq_len(2L^K - 1L)) {
    bits  <- bitwAnd(mask, 2L^(seq_len(K) - 1L)) > 0L
    A     <- which(bits)
    not_A <- which(!bits)
    log_pA <- sum(log(rho[A])) +
              (if (length(not_A)) sum(log1p(-rho[not_A])) else 0) -
              log_norm
    pA <- exp(log_pA)
    if (length(A) == 1L) {
      E_y[A] <- E_y[A] + pA
    } else {
      mu_A   <- mu[A] / sum(mu[A])
      E_y[A] <- E_y[A] + pA * mu_A
    }
  }
  E_y
}

x2_grid  <- seq(-2.5, 2.5, length.out = 60L)
S_me     <- 200L
draws_me <- fit$draws(c("beta_raw", "tau"), format = "df") |>
  select(-starts_with(".")) |>
  slice_sample(n = S_me)

me_rows <- vector("list", S_me)
for (s in seq_len(S_me)) {
  br_cols <- paste0("beta_raw[", rep(seq_len(P), K - 1L), ",",
                    rep(seq_len(K - 1L), each = P), "]")
  beta_s  <- cbind(
    matrix(as.numeric(draws_me[s, br_cols]), nrow = P, ncol = K - 1L),
    0
  )
  tau_s <- as.numeric(draws_me$tau[s])

  props <- t(vapply(x2_grid, function(x2) {
    expected_prop(as.vector(t(beta_s) %*% c(1, x2)), tau_s)
  }, numeric(K)))                                      # n_grid x K

  me_rows[[s]] <- data.frame(
    draw = s,
    x2   = rep(x2_grid, K),
    dim  = paste0("dim ", rep(seq_len(K), each = length(x2_grid))),
    E_y  = as.vector(props)
  )
}
me_df <- do.call(rbind, me_rows)

me_summary <- me_df |>
  group_by(dim, x2) |>
  summarise(
    med = median(E_y),
    q05 = quantile(E_y, 0.05),
    q95 = quantile(E_y, 0.95),
    .groups = "drop"
  )

# True expected proportions for reference
true_me <- data.frame(
  x2  = rep(x2_grid, K),
  dim = paste0("dim ", rep(seq_len(K), each = length(x2_grid))),
  E_y = as.vector(t(vapply(x2_grid, function(x2) {
    expected_prop(as.vector(t(beta_true) %*% c(1, x2)), tau_true)
  }, numeric(K))))
)

p_me <- ggplot(me_summary, aes(x = x2, colour = dim, fill = dim)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.15, colour = NA) +
  geom_line(aes(y = med), linewidth = 1) +
  geom_line(data = true_me, aes(y = E_y),
            linetype = "dashed", linewidth = 0.75, colour = "black",
            inherit.aes = FALSE) +
  facet_wrap(~dim) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title    = "Marginal Effects: E[proportion] vs Predictor X[,2]",
    subtitle = "Coloured = posterior median + 90% CI  |  Black dashed = true value",
    x        = "X[,2] (continuous predictor)",
    y        = "Expected proportion",
    colour   = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

print(p_me)
