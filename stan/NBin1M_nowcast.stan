functions {
  #include "functions/geometric_random_walk.stan"
  #include "functions/observe_onsets_with_delay.stan"
  #include "functions/predict_pois_rng.stan"
  #include "functions/combine_obs_with_predicted_obs.stan"
  #include "functions/multiply_array.stan"
}

data {
  int n;                // number of days
  int m;                // number of reports
  array[n] int p;       // number of observations per day
  array[m] int obs;     // observed symptom onsets
  int d;                // number of reporting delays
}

transformed data{
  array[n] int P = to_int(cumulative_sum(p));
  array[n] int D = to_int(cumulative_sum(rep_array(d, n)));
}

parameters {
  real<lower=0> nb_size;
  real<lower=0> init_onsets;
  array[n-1] real rw_noise;
  vector<lower=0>[n] random_effect;
  real<lower=0> rw_sd;
  simplex[d] reporting_delay; // reporting delay distribution
}

transformed parameters {
  array[n] real lambda = geometric_random_walk(init_onsets, rw_noise, rw_sd);
  array[n] real lambda_times_nb_size = multiply_array(nb_size, lambda);
  array[n] real lambda_with_re = multiply_array(random_effect, lambda);
  array[m] real lambda_with_re_times_pis = observe_onsets_with_delay(lambda_with_re, reporting_delay, P, p);
}

model {
  // Prior
  init_onsets ~ normal(3, 2) T[0,];
  rw_noise ~ std_normal();
  rw_sd ~ normal(0, 0.2) T[0,];
  reporting_delay ~ dirichlet(rep_vector(1, d));
  nb_size ~ normal(1, 3) T[0, ];
  // Model
  random_effect ~ gamma(lambda_times_nb_size, lambda_times_nb_size);
  // Likelihood
  obs ~ poisson(lambda_with_re_times_pis);
}

generated quantities {
  array[d*n] real complete_lambdas_with_re_times_pis = observe_onsets_with_delay(lambda_with_re, reporting_delay, D, rep_array(d, n));
  array[d * n - m] int predicted_counts = predict_pois_rng(complete_lambdas_with_re_times_pis, P, p, d, D);
  array[n] int nowcast = combine_obs_with_predicted_obs(obs, predicted_counts, P, p, d, D);
}
