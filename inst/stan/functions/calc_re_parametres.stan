/**
 * Calculate the parameters of the random effect distribution
 *
 * Calculates the parameters, which enter as the shape and rate of the gamma
 * distribution of the random effects needed for the NegBin2M and NegBin1M
 * observation model.
 *
 * @param lambda Expected value of the total observations, an array of reals.
 *
 * @param nb_size Size = (1 / Dispersion parameter) for the negative binomial
 * model, a vector of size 0 for Poisson (being ignored), size 1 for the
 * negative binomial models
 *
 * @param model_obs Indicator of the model used (0 for Poisson, 1 for NegBinX,
 * 2 for NegBin2D, 3 for NegBin1D, 4 for NegBin2M and 5 for NegBin1M).
 *
 * @return Array of the random effect distribution parameters, same length as
 * the expected value of the total observations.
 **/
array[] real calc_re_parameters(array[] real lambda, vector nb_size, int model_obs) {
  int n = num_elements(lambda);
  array[n] real re_params;
  if (model_obs == 4) {
    re_params = rep_array(nb_size[1], n);
  } else if (model_obs == 5) {
    re_params = multiply_array(nb_size[1], lambda);
  }
  return(re_params);
}
