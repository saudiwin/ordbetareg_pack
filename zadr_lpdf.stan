functions {
  // Zero-Adjusted Dirichlet (ZADR) log-density for a single observation.
  //
  // Zeros in y mark structural/observed zero components. The non-zero
  // components are renormalized so their mu sums to 1, keeping the shared
  // precision phi constant across all zero patterns (Tsagris & Stewart 2018).
  //
  // y        : K-vector composition, may have zeros, non-zero entries sum to 1
  // mu       : K-vector softmax mean (sum = 1); renormalized inside to active set
  // phi      : shared precision (> 0); equals sum of Dirichlet concentrations
  // log_theta: log(proportion of observations sharing this zero pattern)
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
        real mu_k   = mu[k] / mu_active_sum;   // renormalize to active set
        real alpha_k = phi * mu_k;
        lgamma_sum  += lgamma(alpha_k);
        y_terms     += (alpha_k - 1) * log(y[k]);
      }
    }

    // lgamma(phi) = lgamma(sum alpha_k) because renormalized mu sums to 1
    return log_theta + lgamma(phi) - lgamma_sum + y_terms;
  }
}
