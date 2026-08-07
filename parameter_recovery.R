# Parameter recovery for two ordered Dirichlet distributions
#
# Distribution 1: Symmetric Ordered Dirichlet (Gumbel parameterization)
#   K=3 dimensions; parameters: mu (simplex), phi (precision), eta (K-vector
#   of logit-scale marginal activity probabilities)
#
# Distribution 2: Ordered Dirichlet (cutpoint parameterization)
#   K=2 dimensions; parameters: mu (simplex), phi (precision),
#   cut (K x 2 matrix of ordered cutpoints per dimension)
#   eta = 0 here (intercept-only); in a regression context eta would be a
#   linear predictor added per observation.
#
# Requires: cmdstanr, posterior, ggplot2, dplyr, tidyr

library(cmdstanr)
library(posterior)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(42)

# Minimal rdirichlet using gamma draws (avoids extra dependency)
rdirichlet <- function(n, alpha) {
  k <- length(alpha)
  x <- matrix(rgamma(n * k, rep(alpha, each = n)), nrow = n, ncol = k)
  x / rowSums(x)
}


# ══════════════════════════════════════════════════════════════════════════
# PART 1  ·  Symmetric Ordered Dirichlet  (ordered_dirichlet_gumbel_lpdf)
# ══════════════════════════════════════════════════════════════════════════
#
# Each dimension k is independently "active" (y[k] > 0) with Bernoulli
# probability rho[k] = sigmoid(eta[k]), conditioned on at least one active.
# When |active| >= 2 the active components are jointly Dirichlet with
# concentration  alpha_k = phi * mu[k] / sum(mu[active]).

## Simulate ----------------------------------------------------------------

sim_sym_ord_dirichlet <- function(N, mu, phi, eta) {
  K   <- length(mu)
  rho <- plogis(eta)
  Y   <- matrix(0.0, nrow = N, ncol = K)
  for (i in seq_len(N)) {
    repeat {
      a <- rbinom(K, 1L, rho) == 1L
      if (any(a)) break          # condition on non-empty active set
    }
    idx <- which(a)
    if (length(idx) == 1L) {
      Y[i, idx] <- 1.0
    } else {
      mu_a      <- mu[idx] / sum(mu[idx])
      Y[i, idx] <- rdirichlet(1L, phi * mu_a)
    }
  }
  Y
}

K1        <- 3L
mu_true1  <- c(0.50, 0.30, 0.20)
phi_true1 <- 8.0
eta_true1 <- c(0.8, 0.0, -0.8)
N1        <- 300L

Y1 <- sim_sym_ord_dirichlet(N1, mu_true1, phi_true1, eta_true1)

cat("Part 1 data – proportion zeros per dimension:\n")
print(round(colMeans(Y1 == 0), 3))

## Stan model ---------------------------------------------------------------

sym_stan <- '
functions {
  real symmetric_ordered_dirichlet_lpdf(
      vector y, vector mu, real phi, vector eta) {
    int K = num_elements(y);
    int A_size = 0;
    array[K] int active;
    real lp_subset      = 0;
    real log_prob_empty = 0;
    for (k in 1:K) {
      real log_rho   = -log1p_exp(-eta[k]);
      real log1m_rho = -log1p_exp(eta[k]);
      log_prob_empty += log1m_rho;
      if (y[k] > 0) {
        A_size       += 1;
        active[A_size] = k;
        lp_subset    += log_rho;
      } else {
        lp_subset += log1m_rho;
      }
    }
    lp_subset -= log1m_exp(log_prob_empty);   // condition on non-empty
    if (A_size == 1) {
      return lp_subset;
    } else {
      vector[A_size] y_A;
      vector[A_size] alpha_A;
      real mu_sum = 0;
      for (a in 1:A_size) mu_sum += mu[active[a]];
      for (a in 1:A_size) {
        y_A[a]     = y[active[a]];
        alpha_A[a] = phi * mu[active[a]] / mu_sum;
      }
      return lp_subset + dirichlet_lpdf(y_A | alpha_A);
    }
  }
}
data {
  int<lower=1> N;
  int<lower=2> K;
  array[N] vector[K] y;
}
parameters {
  simplex[K] mu;
  real<lower=0> phi;
  vector[K] eta;
}
model {
  mu  ~ dirichlet(rep_vector(1.0, K));
  phi ~ exponential(0.2);
  eta ~ normal(0, 2);
  for (i in 1:N)
    y[i] ~ symmetric_ordered_dirichlet(mu, phi, eta);
}
'

sym_file <- tempfile(fileext = ".stan")
writeLines(sym_stan, sym_file)
sym_mod  <- cmdstan_model(sym_file)

sym_fit <- sym_mod$sample(
  data            = list(N = N1, K = K1, y = Y1),
  chains          = 4L,
  parallel_chains = 4L,
  iter_warmup     = 1000L,
  iter_sampling   = 1000L,
  seed            = 42L,
  refresh         = 500L
)

## Parameter recovery -------------------------------------------------------

sym_truth <- tibble(
  param = c(sprintf("mu[%d]", seq_len(K1)),
            "phi",
            sprintf("eta[%d]", seq_len(K1))),
  truth = c(mu_true1, phi_true1, eta_true1)
)

sym_post <- sym_fit$draws(c("mu", "phi", "eta"), format = "df") |>
  select(-starts_with(".")) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  left_join(sym_truth, by = "param")

cat("\n── Part 1 recovery summary ──────────────────────────────────────\n")
sym_fit$summary(c("mu", "phi", "eta")) |>
  left_join(sym_truth, by = c("variable" = "param")) |>
  select(variable, truth, mean, sd, q5, q95, rhat, ess_bulk) |>
  print(n = Inf)

p1_hist <- ggplot(sym_post, aes(x = value)) +
  geom_histogram(bins = 60, fill = "#2c7bb6", alpha = 0.75, colour = NA) +
  geom_vline(aes(xintercept = truth), colour = "red", linewidth = 0.9) +
  facet_wrap(~param, scales = "free") +
  labs(
    title    = "Symmetric Ordered Dirichlet — Parameter Recovery",
    subtitle = "Red line = true value",
    x = NULL, y = "Count"
  ) +
  theme_minimal(base_size = 12)

print(p1_hist)

# Interval plot: posterior median + 90 % CI vs. truth
sym_intervals <- sym_fit$summary(c("mu", "phi", "eta"),
                                  ~quantile(.x, probs = c(0.05, 0.5, 0.95))) |>
  rename(param = variable, q05 = `5%`, med = `50%`, q95 = `95%`) |>
  left_join(sym_truth, by = "param")

p1_int <- ggplot(sym_intervals, aes(y = param)) +
  geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 1.2,
                 colour = "#2c7bb6", alpha = 0.7) +
  geom_point(aes(x = med), size = 3, colour = "#2c7bb6") +
  geom_point(aes(x = truth), shape = 18, size = 4, colour = "red") +
  labs(
    title    = "Symmetric Ordered Dirichlet — 90% Posterior Intervals",
    subtitle = "Blue = posterior median + 90% CI  |  Red diamond = truth",
    x = "Value", y = NULL
  ) +
  theme_minimal(base_size = 12)

print(p1_int)


# ══════════════════════════════════════════════════════════════════════════
# PART 2  ·  Ordered Dirichlet  (cutpoint parameterization, K = 2)
# ══════════════════════════════════════════════════════════════════════════
#
# Per-dimension ordered cutpoints cut[k] = (c1, c2) define three states:
#   p0[k]  = P(y[k] = 0)      = inv_logit(c1 - eta[k])
#   pm[k]  = P(y[k] ∈ (0,1)) = inv_logit(c2 - eta[k]) - p0[k]
#   p1[k]  = P(y[k] = 1)      = 1 - inv_logit(c2 - eta[k])
#
# Feasible outcomes (conditions for a valid simplex):
#   · exactly one dimension = 1, all others = 0  (vertex)
#   · two or more dimensions ∈ (0,1), Dirichlet-distributed  (interior)
#
# For K=2, three feasible subsets:  {1}, {2}, {1,2}
# log_norm = log_sum_exp of their face weights (implemented below).

## Simulate ----------------------------------------------------------------

sim_ord_dirichlet <- function(N, mu, phi, cut_mat, eta = NULL) {
  K     <- nrow(cut_mat)
  if (is.null(eta)) eta <- rep(0.0, K)
  p0_v  <- plogis(cut_mat[, 1L] - eta)
  pcum1 <- plogis(cut_mat[, 2L] - eta)
  pm_v  <- pcum1 - p0_v
  p1_v  <- 1.0 - pcum1
  Y     <- matrix(0.0, nrow = N, ncol = K)
  for (i in seq_len(N)) {
    repeat {
      # 1 = zero, 2 = middle, 3 = exactly-one-candidate
      states <- vapply(seq_len(K), function(k) {
        sample.int(3L, 1L, prob = c(p0_v[k], pm_v[k], p1_v[k]))
      }, integer(1L))
      n_one <- sum(states == 3L)
      n_mid <- sum(states == 2L)
      if ((n_one == 1L && n_mid == 0L) ||   # vertex case
          (n_one == 0L && n_mid >= 2L))       # interior case
        break
    }
    if (any(states == 3L)) {
      Y[i, which(states == 3L)] <- 1.0
    } else {
      idx      <- which(states == 2L)
      mu_a     <- mu[idx] / sum(mu[idx])
      Y[i, idx] <- rdirichlet(1L, phi * mu_a)
    }
  }
  Y
}

K2        <- 2L
mu_true2  <- c(0.60, 0.40)
phi_true2 <- 6.0
# cut_true2[k, ] = c(lower_cutpoint, upper_cutpoint) for dimension k
# (eta = 0 for both dims in this intercept-only demo)
cut_true2 <- matrix(c(-1.5,  0.5,   # dim 1: p0≈0.18, pm≈0.44, p1≈0.38
                       -1.0,  1.0), # dim 2: p0≈0.27, pm≈0.46, p1≈0.27
                    nrow = K2, byrow = TRUE)
N2 <- 300L

Y2 <- sim_ord_dirichlet(N2, mu_true2, phi_true2, cut_true2)

cat("\nPart 2 data – state proportions per dimension:\n")
for (k in seq_len(K2)) {
  cat(sprintf("  dim %d  zero=%.2f  mid=%.2f  one=%.2f\n", k,
              mean(Y2[, k] == 0),
              mean(Y2[, k] > 0 & Y2[, k] < 1),
              mean(Y2[, k] == 1)))
}

## Stan model ---------------------------------------------------------------

ord_stan <- '
functions {
  // K=2 variant: log_norm enumerates all 3 feasible subsets.
  // For K>2 replace with general subset enumeration (2^K - 1 terms).
  real ordered_dirichlet_lpdf(
      vector y, vector mu, real phi, vector eta, array[] vector cut) {
    int K = num_elements(y);
    array[K] int active;
    int A_size = 0;
    real log_w_obs = 0;
    real log_norm;
    vector[K] log_p0;
    vector[K] log_pm;
    vector[K] log_p1;
    for (k in 1:K) {
      real p0    = inv_logit(cut[k][1] - eta[k]);
      real p1cum = inv_logit(cut[k][2] - eta[k]);
      log_p0[k] = log(p0);
      log_pm[k] = log(p1cum - p0);
      log_p1[k] = log1m(p1cum);
      if (y[k] > 0) { A_size += 1; active[A_size] = k; }
    }
    // K=2 normalisation: w{1} + w{2} + w{1,2}
    log_norm = log_sum_exp(log_p1[1] + log_p0[2],
               log_sum_exp(log_p0[1] + log_p1[2],
                           log_pm[1]  + log_pm[2]));
    if (A_size == 1) {
      int j = active[1];
      log_w_obs = log_p1[j];
      for (k in 1:K) if (k != j) log_w_obs += log_p0[k];
    } else {
      for (k in 1:K) {
        if (y[k] > 0) log_w_obs += log_pm[k];
        else          log_w_obs += log_p0[k];
      }
    }
    if (A_size >= 2) {
      vector[A_size] y_A;
      vector[A_size] alpha_A;
      real mu_sum = 0;
      for (a in 1:A_size) mu_sum += mu[active[a]];
      for (a in 1:A_size) {
        y_A[a]     = y[active[a]];
        alpha_A[a] = phi * mu[active[a]] / mu_sum;
      }
      return log_w_obs - log_norm + dirichlet_lpdf(y_A | alpha_A);
    } else {
      return log_w_obs - log_norm;
    }
  }
}
data {
  int<lower=1> N;
  array[N] vector[2] y;
}
parameters {
  simplex[2] mu;
  real<lower=0> phi;
  array[2] ordered[2] cut;   // ordered[2] enforces cut[k][1] < cut[k][2]
}
model {
  mu  ~ dirichlet(rep_vector(1.0, 2));
  phi ~ exponential(0.2);
  for (k in 1:2) cut[k] ~ normal(0, 2);
  for (i in 1:N)
    y[i] ~ ordered_dirichlet(mu, phi, rep_vector(0.0, 2), cut);
}
'

ord_file <- tempfile(fileext = ".stan")
writeLines(ord_stan, ord_file)
ord_mod  <- cmdstan_model(ord_file)

ord_fit <- ord_mod$sample(
  data            = list(N = N2, y = Y2),
  chains          = 4L,
  parallel_chains = 4L,
  iter_warmup     = 1000L,
  iter_sampling   = 1000L,
  seed            = 42L,
  refresh         = 500L
)

## Parameter recovery -------------------------------------------------------

ord_truth <- tibble(
  param = c("mu[1]", "mu[2]", "phi",
            "cut[1,1]", "cut[1,2]", "cut[2,1]", "cut[2,2]"),
  truth = c(mu_true2, phi_true2,
            cut_true2[1, 1], cut_true2[1, 2],
            cut_true2[2, 1], cut_true2[2, 2])
)

ord_post <- ord_fit$draws(c("mu", "phi", "cut"), format = "df") |>
  select(-starts_with(".")) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  left_join(ord_truth, by = "param")

cat("\n── Part 2 recovery summary ──────────────────────────────────────\n")
ord_fit$summary(c("mu", "phi", "cut")) |>
  left_join(ord_truth, by = c("variable" = "param")) |>
  select(variable, truth, mean, sd, q5, q95, rhat, ess_bulk) |>
  print(n = Inf)

p2_hist <- ggplot(ord_post, aes(x = value)) +
  geom_histogram(bins = 60, fill = "#d7191c", alpha = 0.75, colour = NA) +
  geom_vline(aes(xintercept = truth), colour = "darkblue", linewidth = 0.9) +
  facet_wrap(~param, scales = "free") +
  labs(
    title    = "Ordered Dirichlet — Parameter Recovery",
    subtitle = "Blue line = true value",
    x = NULL, y = "Count"
  ) +
  theme_minimal(base_size = 12)

print(p2_hist)

ord_intervals <- ord_fit$summary(c("mu", "phi", "cut"),
                                   ~quantile(.x, probs = c(0.05, 0.5, 0.95))) |>
  rename(param = variable, q05 = `5%`, med = `50%`, q95 = `95%`) |>
  left_join(ord_truth, by = "param")

p2_int <- ggplot(ord_intervals, aes(y = param)) +
  geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 1.2,
                 colour = "#d7191c", alpha = 0.7) +
  geom_point(aes(x = med), size = 3, colour = "#d7191c") +
  geom_point(aes(x = truth), shape = 18, size = 4, colour = "darkblue") +
  labs(
    title    = "Ordered Dirichlet — 90% Posterior Intervals",
    subtitle = "Red = posterior median + 90% CI  |  Blue diamond = truth",
    x = "Value", y = NULL
  ) +
  theme_minimal(base_size = 12)

print(p2_int)
