functions {
  // Zero-Adjusted Dirichlet (ZADR) log-density — Tsagris & Stewart (2018).
  // Identical to zadr_lpdf.stan but receives log_theta_b as a scalar so the
  // caller can pass either a fixed or a stochastic (drawn from posterior) value.
  real zadr_lpdf(vector y, vector mu, real phi, real log_theta) {
    int K = num_elements(y);
    real mu_active_sum = 0;
    real lgamma_sum   = 0;
    real y_terms      = 0;

    for (k in 1:K) {
      if (y[k] > 0)
        mu_active_sum += mu[k];
    }

    for (k in 1:K) {
      if (y[k] > 0) {
        real mu_k    = mu[k] / mu_active_sum;
        real alpha_k = phi * mu_k;
        lgamma_sum  += lgamma(alpha_k);
        y_terms     += (alpha_k - 1) * log(y[k]);
      }
    }

    return log_theta + lgamma(phi) - lgamma_sum + y_terms;
  }
}

data {
  int<lower=1> N;
  int<lower=2> K;
  int<lower=1> P;
  int<lower=1> B;                      // number of distinct zero patterns

  matrix[N, P] X;
  array[N] vector[K] y;

  array[N] int<lower=1, upper=B> pattern_id;  // which pattern each obs belongs to

  // Symmetric Dirichlet prior concentration; alpha_0 = 1 => uniform on simplex,
  // alpha_0 < 1 => sparse (concentrates mass on few patterns).
  real<lower=0> alpha_0;
}

parameters {
  // Multinomial logistic coefficients; component 1 is reference (beta[,1] = 0)
  matrix[P, K - 1] beta_raw;

  real<lower=0> phi;      // shared Dirichlet precision

  simplex[B] theta;       // stochastic zero-pattern probabilities
}

transformed parameters {
  matrix[P, K] beta;

  beta[, 1] = rep_vector(0, P);
  for (k in 2:K)
    beta[, k] = beta_raw[, k - 1];
}

model {
  to_vector(beta_raw) ~ normal(0, 1);
  phi ~ gamma(2, 0.1);
  theta ~ dirichlet(rep_vector(alpha_0, B));

  for (i in 1:N) {
    vector[K] mu = softmax(beta' * to_vector(X[i]'));
    y[i] ~ zadr(mu, phi, log(theta[pattern_id[i]]));
  }
}

generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    vector[K] mu = softmax(beta' * to_vector(X[i]'));
    log_lik[i] = zadr_lpdf(y[i] | mu, phi, log(theta[pattern_id[i]]));
  }
}
