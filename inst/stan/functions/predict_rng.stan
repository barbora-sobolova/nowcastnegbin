array[] int predict_rng(array[] real exp_obs, vector nb_size,
                          vector reporting_delay,
                          int model_obs, array[] int P, array[] int p, int d,
                          array[] int D) {
  // Prepare the expected values and the neg. binom sizes
  int n = num_elements(p);
  int n_predict = d * n - sum(p);
  array[n*d] real nb_size_expanded;
  if (model_obs == 1) {
    // NegBinX
    nb_size_expanded = rep_array(nb_size[1], d*n);
  } else if (model_obs == 2) {
    // NegBin2D
    nb_size_expanded = observe_onsets_with_delay(rep_array(nb_size[1], n), reporting_delay, D, rep_array(d, n));
  } else if (model_obs == 3) {
    // NegBin1D
    nb_size_expanded = multiply_array(nb_size[1], exp_obs);
  }

  // Sample the predictions
  array[n_predict] int predicted_counts;
  int predicted_index = 0;
  for (i in 1:n) {
    int missing_reports = d - p[i];
    int P_index = 0;
    int D_index = 0;
    if (i != 1) {
      P_index = P[i - 1];
      D_index = D[i - 1];
    }
    if (missing_reports > 0) {
      for (j in 1:missing_reports) {
        if (model_obs == 0 || model_obs == 4 || model_obs == 5) {
          predicted_counts[predicted_index + j] = poisson_rng(
              exp_obs[D_index + p[i] + j]
            );
        } else {
          predicted_counts[predicted_index + j] = neg_binomial_2_rng(
              exp_obs[D_index + p[i] + j],
              nb_size_expanded[D_index + p[i] + j]
            );
        }
      }
      predicted_index += missing_reports;
    }
  }
  return predicted_counts;
}
