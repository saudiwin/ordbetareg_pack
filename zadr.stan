functions {
  // Zero-Adjusted Dirichlet (ZADR) log-density — Tsagris & Stewart (2018).
  // Non-zero components of y are modelled as Dirichlet with renormalized mean
  // mu and shared precision phi. log_theta accounts for the probability of the
  // observed zero pattern (log proportion of obs with the same pattern).
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

  matrix[N, P] X;
  array[N] vector[K] y;

  // For each observation: log(n_b / N) where n_b = number of observations
  // sharing the same zero pattern. Computed in R before fitting.
  vector[N] log_theta;
}

parameters {
  // Multinomial logistic coefficients; component 1 is reference (beta[,1] = 0)
  matrix[P, K - 1] beta_raw;

  real<lower=0> phi;     // shared Dirichlet precision
}

transformed parameters {
  matrix[P, K] beta;

  beta[, 1] = rep_vector(0, P);   // reference category
  for (k in 2:K)
    beta[, k] = beta_raw[, k - 1];
}

model {
  to_vector(beta_raw) ~ normal(0, 1);
  phi ~ gamma(2, 0.1);    // mean 20, weakly informative

  for (i in 1:N) {
    vector[K] mu = softmax(beta' * to_vector(X[i]'));
    y[i] ~ zadr(mu, phi, log_theta[i]);
  }
}

generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    vector[K] mu = softmax(beta' * to_vector(X[i]'));
    log_lik[i] = zadr_lpdf(y[i] | mu, phi, log_theta[i]);
  }
}
