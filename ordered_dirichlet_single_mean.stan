data {
  int<lower=1> N;
  int<lower=2> K;
  int<lower=1> P;

  matrix[N, P] X;
  array[N] vector[K] y;
}

parameters {
  // baseline category K has coefficients fixed to zero
  matrix[P, K - 1] beta_raw;

  real tau;

  vector[P] delta;          // precision model
}

transformed parameters {
  matrix[P, K] beta;

  for (k in 1:(K - 1)) {
    beta[, k] = beta_raw[, k];
  }

  beta[, K] = rep_vector(0, P);
}

model {
  // priors
  to_vector(beta_raw) ~ normal(0, 1);
  tau ~ normal(0, 2);
  delta ~ normal(0, 1);

  for (i in 1:N) {
    vector[K] theta;
    real phi_i;

    theta = beta' * to_vector(X[i]');
    phi_i = exp(dot_product(to_vector(X[i]'), delta));

    y[i] ~ ordered_dirichlet_shared(theta, phi_i, tau);
  }
}