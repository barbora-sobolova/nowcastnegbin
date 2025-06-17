functions {
  #include "functions/geometric_random_walk.stan"
  #include "functions/observe_onsets_with_delay.stan"
  #include "functions/predict_rng.stan"
  #include "functions/combine_obs_with_predicted_obs.stan"
  #include "functions/multiply_array.stan"
  #include "functions/calc_exp_total_obs.stan"
  #include "functions/calc_re_parametres.stan"
  #include "functions/expand_nb_size.stan"
  #include "functions/obs_lpmf.stan"
}

data {
  int n;                // number of days
  int m;                // number of reports
  array[n] int p;       // number of observations per day
  array[m] int obs;     // observed symptom onsets
  int d;                // number of reporting delays
  int model_obs;        // observation model number
}

transformed data{
  array[n] int P = to_int(cumulative_sum(p));
  array[n] int D = to_int(cumulative_sum(rep_array(d, n)));
}

parameters {
  vector<lower=0>[model_obs == 0 ? 0 : 1] nb_size; // For the NegBin models only
  real<lower=0> init_onsets;
  array[n-1] real rw_noise;
  real<lower=0> rw_sd;
  simplex[d] reporting_delay; // reporting delay distribution
  // random effect needed for the NegBin1M and NegBin2M models only
  vector<lower=0>[model_obs == 4 || model_obs == 5 ? n : 0] random_effect;
}

transformed parameters {
  // Mean value of the total count without the random effect
  array[n] real lambda = geometric_random_walk(init_onsets, rw_noise, rw_sd);
  // Parameters of the random effect distribution
  array[n] real re_params = calc_re_parameters(lambda, nb_size, model_obs);
  // Mean value of the total counts with the random effect
  array[n] real exp_total_obs = calc_exp_total_obs(lambda, random_effect, model_obs);
  // Mean value of the counts split by the delay with the random effect
  // Right truncated reports
  array[m] real exp_obs = observe_onsets_with_delay(exp_total_obs, reporting_delay, P, p);
  // Mean value of the counts split by the delay with the random effect
  // Complete reports
  array[d*n] real exp_obs_complete = observe_onsets_with_delay(exp_total_obs, reporting_delay, D, rep_array(d, n));
  // Size parameter of the NegBin distribution expanded according to the
  // observation modelto match the length of the mean value
  // Right truncated reports
  array[m] real nb_size_expanded = expand_nb_size(exp_obs, nb_size, reporting_delay, model_obs, P, p);
  // Size parameter of the NegBin distribution expanded according to the
  // observation modelto match the length of the mean value
  // Complete reports
  array[d*n] real nb_size_expanded_complete = expand_nb_size(exp_obs_complete, nb_size, reporting_delay, model_obs, D, rep_array(d, n));
}

model {
  // Prior
  init_onsets ~ normal(3, 2) T[0,];
  rw_noise ~ std_normal();
  rw_sd ~ normal(0, 0.2) T[0,];
  reporting_delay ~ dirichlet(rep_vector(1, d));
  nb_size ~ normal(1, 3) T[0, ];
  // Random effect for the NegBin1M and NegBin2M models
  if (model_obs == 4 || model_obs == 5) {
    random_effect ~ gamma(re_params, re_params);
  }
  // Likelihood
  obs ~ obs(exp_obs, nb_size_expanded, model_obs, P, p);
}

generated quantities {
  array[d * n - m] int predicted_counts = predict_rng(
      exp_obs_complete, nb_size_expanded_complete, model_obs, P, p, d, D
     );
  array[n] int nowcast = combine_obs_with_predicted_obs(
      obs, predicted_counts, P, p, d, D
    );
}
