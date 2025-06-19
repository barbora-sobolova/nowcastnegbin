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
 * model, an array of size 0 for Poisson (being ignored), same size `exp_obs` as
 * the `exp_obs` paramater for the negative binomial models. Used only for
 * NegBinX, NegBin2D and NegBin1D models.
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
real obs_lpmf(array[] int obs, array[] real exp_obs, array[] real nb_size,
                  int model_obs, array[] int P, array[] int p) {
   real tar = 0;

   if (model_obs == 0 || model_obs == 4 || model_obs == 5) {
     // Poisson, NegBin2M and NegBin1M
     tar = poisson_lpmf(obs | exp_obs);
   } else {
     // NegBinX, NegBin2D and NegBin1D
     tar = neg_binomial_2_lpmf(obs | exp_obs, nb_size);
   }
   return(tar);
}
