functions {
  #include "functions/geometric_random_walk.stan"
  #include "functions/observe_onsets_with_delay.stan"
  #include "functions/predict_rng.stan"
  #include "functions/combine_obs_with_predicted_obs.stan"
  #include "functions/multiply_array.stan"
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
  array[n] real lambda = geometric_random_walk(init_onsets, rw_noise, rw_sd);
  array[n] real re_params;
  if (model_obs == 4) {
    re_params = rep_array(nb_size[1], n);
  } else if (model_obs == 5) {
    re_params = multiply_array(nb_size[1], lambda);
  }
  // multiply lambda by the random effect, when needed
  array[n] real exp_total_obs;
  if (model_obs == 4 || model_obs == 5) {
    exp_total_obs = multiply_array(random_effect, lambda);
  } else {
    exp_total_obs = lambda;
  }
  array[m] real exp_obs = observe_onsets_with_delay(exp_total_obs, reporting_delay, P, p);
  array[d*n] real exp_obs_complete = observe_onsets_with_delay(exp_total_obs, reporting_delay, D, rep_array(d, n));
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
  obs ~ obs(exp_obs, nb_size, reporting_delay, model_obs, P, p);
}

generated quantities {
  array[d * n - m] int predicted_counts = predict_rng(
      exp_obs_complete, nb_size, reporting_delay, model_obs, P, p, d, D
     );
  array[n] int nowcast = combine_obs_with_predicted_obs(
      obs, predicted_counts, P, p, d, D
    );
}
