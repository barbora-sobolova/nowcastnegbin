/**
 * Predict the missing partial reports
 *
 * Predicts the entries of the reporting triangle.
 *
 * @param exp_obs Expectations of the partial counts, an array of reals. For the
 * Poisson, NegBin2M and NegBin1M model, it includes the random effect.
 *
 * @param nb_size Size = (1 / Dispersion parameter) for the negative binomial
 * model, an array of size 0 for Poisson (being ignored), same size `exp_obs` as
 * the `exp_obs` paramater for the negative binomial models. Used only for
 * NegBinX, NegBin2D and NegBin1D models.
 *
 * @param model_obs Indicator of the model used (0 for Poisson, 1 for NegBinX,
 * 2 for NegBin2D, 3 for NegBin1D, 4 for NegBin2M and 5 for NegBin1M).
 *
 * @param d maximum delay
 *
 * @param P,p,D arrays of lookup variables
 *
 * @return The predicted partial reports as a flat array concatenated rowwise
 *
 * @note The function selects between the negative binomial models and the
 * Poisson model based on the `model_obs` flag. It uses
 * `neg_binomial_2_rng` for NegBinX, NegBin2D and NegBin1D models and
 * `poisson_rng` for the Poisson, NegBin2M and NegBin1M models.
 **/
array[] int predict_rng(array[] real exp_obs, array[] real nb_size,
                          int model_obs, array[] int P, array[] int p, int d,
                          array[] int D) {
  // Prepare the expected values and the neg. binom sizes
  int n = num_elements(p);
  int n_predict = d * n - sum(p);

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
              nb_size[D_index + p[i] + j]
            );
        }
      }
      predicted_index += missing_reports;
    }
  }
  return predicted_counts;
}
