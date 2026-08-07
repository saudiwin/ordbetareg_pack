
// coupled_ordered_dirichlet.stan
//
// Coupled hurdle-inclusion ordered Dirichlet regression.
//
// Allows:
//   y_i = all zeros
//   component-level zeros
//   component-level ones
//   Dirichlet allocation over active components
//
// Key idea:
//   z_ik ~ Bernoulli(rho_ik)
//   A_i = {k : z_ik = 1}
//   A_i = empty  -> y_i = all zeros
//   A_i nonempty -> y_i lies on simplex over A_i
//
// No conditioning constant is needed because the empty active set is allowed.

functions {
  real coupled_ordered_dirichlet_lpdf(vector y,
                                      vector theta,
                                      real phi,
                                      real tau,
                                      real a_rho,
                                      real eps) {
    int K = num_elements(y);
    real total = sum(y);

    vector[K] mu = softmax(theta);
    real lp_A = 0;

    // Case 1: all-zero row.
    if (total <= eps) {
      for (k in 1:K) {
        real eta_k = a_rho * theta[k] - tau;
        lp_A += bernoulli_logit_lpmf(0 | eta_k);
      }
      return lp_A;
    }

    // Case 2: nonzero row.
    {
      vector[K] y_norm = y / total;
      int A_size = 0;
      array[K] int active;

      for (k in 1:K) {
        real eta_k = a_rho * theta[k] - tau;

        if (y_norm[k] > eps) {
          A_size += 1;
          active[A_size] = k;
          lp_A += bernoulli_logit_lpmf(1 | eta_k);
        } else {
          lp_A += bernoulli_logit_lpmf(0 | eta_k);
        }
      }

      // This should not happen because total > eps,
      // but keep as a guard.
      if (A_size == 0) {
        return negative_infinity();
      }

      // Singleton active set -> vertex of simplex.
      // Dirichlet part is degenerate and contributes 0.
      if (A_size == 1) {
        return lp_A;
      }

      // Multiple active components -> Dirichlet over active face.
      {
        vector[A_size] y_A;
        vector[A_size] alpha_A;
        real mu_sum = 0;

        for (a in 1:A_size) {
          mu_sum += mu[active[a]];
        }

        for (a in 1:A_size) {
          int k = active[a];
          y_A[a] = y_norm[k];
          alpha_A[a] = phi * mu[k] / mu_sum;
        }

        return lp_A + dirichlet_lpdf(y_A | alpha_A);
      }
    }
  }
}

data {
  int<lower=1> N;
  int<lower=2> K;
  int<lower=1> P;

  matrix[N, P] X;

  // y can include all-zero rows.
  array[N] vector<lower=0>[K] y;

  real<lower=0> eps;
}

parameters {
  // Baseline category K fixed to zero for softmax identification.
  matrix[P, K - 1] beta_raw;

  // Shared component-inclusion threshold.
  real tau;

  // Inclusion scale. Positive version recommended for interpretation.
  real<lower=0> a_rho;

  // Precision model.
  vector[P] delta;
}

transformed parameters {
  matrix[P, K] beta;

  for (k in 1:(K - 1)) {
    beta[, k] = beta_raw[, k];
  }
  beta[, K] = rep_vector(0, P);
}

model {
  // Priors.
  to_vector(beta_raw) ~ normal(0, 1);
  tau ~ normal(0, 2);
  a_rho ~ lognormal(0, 0.5);
  delta ~ normal(0, 1);

  for (i in 1:N) {
    vector[K] theta_i;
    real phi_i;

    theta_i = beta' * to_vector(X[i]');
    phi_i = exp(dot_product(to_vector(X[i]'), delta));

    target += coupled_ordered_dirichlet_lpdf(
      y[i] | theta_i, phi_i, tau, a_rho, eps
    );
  }
}

generated quantities {
  array[N] real log_lik;

  for (i in 1:N) {
    vector[K] theta_i;
    real phi_i;

    theta_i = beta' * to_vector(X[i]');
    phi_i = exp(dot_product(to_vector(X[i]'), delta));

    log_lik[i] = coupled_ordered_dirichlet_lpdf(
      y[i] | theta_i, phi_i, tau, a_rho, eps
    );
  }
}
