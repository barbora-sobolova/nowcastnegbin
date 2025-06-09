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
