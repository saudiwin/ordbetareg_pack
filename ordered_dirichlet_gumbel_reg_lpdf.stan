functions {
  real ordered_dirichlet_lpdf(vector y,
                              vector mu,
                              real phi,
                              vector eta) {
    int K = num_elements(y);
    int A_size = 0;
    array[K] int active;

    real lp_A = 0;
    real log_empty = 0;

    for (k in 1:K) {
      real log_rho = -log1p_exp(-eta[k]);
      real log1m_rho = -log1p_exp(eta[k]);

      log_empty += log1m_rho;

      if (y[k] > 0) {
        A_size += 1;
        active[A_size] = k;
        lp_A += log_rho;
      } else {
        lp_A += log1m_rho;
      }
    }

    // condition on at least one active component
    lp_A -= log1m_exp(log_empty);

    if (A_size == 1) {
      return lp_A;
    } else {
      vector[A_size] y_A;
      vector[A_size] alpha_A;
      real mu_sum = 0;

      for (a in 1:A_size) {
        mu_sum += mu[active[a]];
      }

      for (a in 1:A_size) {
        int k = active[a];
        y_A[a] = y[k];
        alpha_A[a] = phi * mu[k] / mu_sum;
      }

      return lp_A + dirichlet_lpdf(y_A | alpha_A);
    }
  }
}