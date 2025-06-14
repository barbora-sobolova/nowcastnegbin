/**
 * Observation likelihood
 *
 * Computes the log probability mass function (LPMF) for observations according
 * to the 6 observational models.
 *
 * @param obs Observations, an array of integers.
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
 * @return The log probability mass of the observations under the specified
 * model.
 *
 * @note The function selects between the negative binomial models and the
 * Poisson model based on the `model_obs` flag. It uses
 * `neg_binomial_2_lpmf` for NegBinX, NegBin2D and NegBin1D models and
 * `poisson_lpmf` for the Poisson, NegBin2M and NegBin1M models.
 */
real obs_lpmf(array[] int obs, array[] real exp_obs, vector nb_size,
                  vector reporting_delay, int model_obs, array[] int P,
                  array[] int p) {
   real tar = 0;
   int m = num_elements(obs);
   int n = num_elements(p);

   if (model_obs == 0 || model_obs == 4 || model_obs == 5) {
     // Poisson, NegBin2M and NegBin1M
     tar = poisson_lpmf(obs | exp_obs);
   } else if (model_obs == 1) {
     // NegBinX
     tar = neg_binomial_2_lpmf(obs | exp_obs, nb_size[1]);
   } else if (model_obs == 2) {
     // NegBin2D
     array[m] real nb_size_expanded = observe_onsets_with_delay(rep_array(nb_size[1], n), reporting_delay, P, p);
     tar = neg_binomial_2_lpmf(obs | exp_obs, nb_size_expanded);
   } else if (model_obs == 3) {
     // NegBin1D
     array[m] real nb_size_expanded = multiply_array(nb_size[1], exp_obs);
     tar = neg_binomial_2_lpmf(obs | exp_obs, nb_size_expanded);
   }
   return(tar);
}
