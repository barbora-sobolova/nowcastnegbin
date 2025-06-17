/**
 * Expand the size parameter of the negative binomial distribution
 *
 * Expands the size parameter of the negative binomial distribution into an
 * array of the same length as the expected value of the partial counts.
 * Required for the NegBinX, NegBin2D and NegBin1D models only. Nothing happens
 * for the Poisson, NegBin2M and NegBin1M models.
 *
 * @param exp_obs Expectations of the partial counts, an array of reals.
 *
 * @param nb_size Size = (1 / Dispersion parameter) for the negative binomial
 * model, a vector of size 0 for Poisson (being ignored), size 1 for the
 * negative binomial models
 *
 * @param reporting_delay reporting delay probability, a vector
 *
 * @param model_obs Indicator of the model used (0 for Poisson, 1 for NegBinX,
 * 2 for NegBin2D, 3 for NegBin1D, 4 for NegBin2M and 5 for NegBin1M).
 *
 * @param P,p arrays of lookup variables
 *
 * @return The expectation of the total
 **/
array[] real expand_nb_size(array[] real exp_obs, vector nb_size,
                                vector reporting_delay, int model_obs,
                                array[] int P, array[] int p) {
  int n = num_elements(p);
  int m = sum(p);
  array[m] real nb_size_expanded;
  if (model_obs == 1) {
     // NegBinX
     nb_size_expanded = rep_array(nb_size[1], m);
  } else if (model_obs == 2) {
    // NegBin2D
    nb_size_expanded = observe_onsets_with_delay(rep_array(nb_size[1], n), reporting_delay, P, p);
  } else if (model_obs == 3) {
    // NegBin1D
    nb_size_expanded = multiply_array(nb_size[1], exp_obs);
  }
  return nb_size_expanded;
}
