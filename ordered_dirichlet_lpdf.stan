functions {
  real ordered_dirichlet_lpdf(vector y,
                              vector mu,
                              real phi,
                              vector eta,
                              array[] vector cut) {
    int K = num_elements(y);
    array[K] int active;
    int A_size = 0;
    real log_w_obs = 0;
    real log_norm;
    
    vector[K] log_p0;
    vector[K] log_pm;
    vector[K] log_p1;

    for (k in 1:K) {
      real p0 = inv_logit(cut[k][1] - eta[k]);
      real p1_cum = inv_logit(cut[k][2] - eta[k]);

      log_p0[k] = log(p0);
      log_pm[k] = log(p1_cum - p0);
      log_p1[k] = log1m(p1_cum);

      if (y[k] > 0) {
        A_size += 1;
        active[A_size] = k;
      }
    }

    // observed face weight
    if (A_size == 1) {
      int j = active[1];
      log_w_obs += log_p1[j];
      for (k in 1:K)
        if (k != j) log_w_obs += log_p0[k];
    } else {
      for (k in 1:K) {
        if (y[k] > 0) log_w_obs += log_pm[k];
        else log_w_obs += log_p0[k];
      }
    }

    // NOTE: for small K, enumerate all feasible subsets B here
    // and compute log_norm = log_sum_exp(log_w_B).
    // For clarity, placeholder:
    log_norm = 0; // replace with subset enumeration

    if (A_size >= 2) {
      vector[A_size] y_A;
      vector[A_size] alpha_A;
      real mu_sum = 0;

      for (a in 1:A_size)
        mu_sum += mu[active[a]];

      for (a in 1:A_size) {
        y_A[a] = y[active[a]];
        alpha_A[a] = phi * mu[active[a]] / mu_sum;
      }

      return log_w_obs - log_norm + dirichlet_lpdf(y_A | alpha_A);
    } else {
      return log_w_obs - log_norm;
    }
  }
}