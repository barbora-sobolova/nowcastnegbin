array[] int predict_nbin_rng(array[] real complete_lambdas_times_pis, array[] real complete_nb_sizes, array[] int P, array[] int p, int d, array[] int D) {
  int n = num_elements(p);
  int n_predict = d * n - sum(p);
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
        predicted_counts[predicted_index + j] = neg_binomial_2_rng(complete_lambdas_times_pis[D_index + p[i] + j], complete_nb_sizes[D_index + p[i] + j]);
      }
      predicted_index += missing_reports;
    }
  }
  return predicted_counts;
}
