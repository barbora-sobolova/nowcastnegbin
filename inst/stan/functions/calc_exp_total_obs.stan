/**
 * Mean value of the total counts
 *
 * Multiplies the mean value of the total observation by the random effect.
 * Only for the NegBin2M and NegBin2M models. For the rest of the models, the
 * mean value remains unchanged.
 *
 * @param lambda Expected value of the total observations without the random
 * effect, an array of reals.
 *
 * @param random_effect multiplicative random effect, used only for the NegBin2M
 * and NegBin1M models, a vector of the same length as `lambda`
 *
 * @param model_obs Indicator of the model used (0 for Poisson, 1 for NegBinX,
 * 2 for NegBin2D, 3 for NegBin1D, 4 for NegBin2M and 5 for NegBin1M).
 *
 * @return The expectation of the total
 **/
array[] real calc_exp_total_obs(array[] real lambda, vector random_effect, int model_obs) {
  int n = num_elements(lambda);
  array[n] real exp_total_obs;
  if (model_obs == 4 || model_obs == 5) {
    exp_total_obs = multiply_array(random_effect, lambda);
  } else {
    exp_total_obs = lambda;
  }
  return(exp_total_obs);
}
