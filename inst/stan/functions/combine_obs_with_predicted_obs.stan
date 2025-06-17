/**
 * Combine the observed and predicted partial counts
 *
 * Combines the flat array of observations with the flat array of predicted
 * counts.
 *
 * @param obs Observations, an array of integers
 *
 * @param predicted_counts, an array of the predicted counts
 *
 * @param d maximum delay
 *
 * @param P,p,D arrays of lookup variables
 *
 * @return The completed reporting triangle as a flat array concatenated
 * rowwise
 **/
array[] int combine_obs_with_predicted_obs(array[] int obs,
                                             array[] int predicted_counts,
                                             array[] int P, array[] int p,
                                             int d, array[] int D) {
    int n = num_elements(p);
    array[n] int reported_cases = rep_array(0, n);
    int predicted_index = 0;
    for (i in 1:n) {
      int missing_reports = d - p[i];
      int P_index = 0;
      int D_index = 0;
      int p_increment = d;
      if (i != 1) {
        P_index = P[i - 1];
        D_index = D[i - 1];
        p_increment = p[i];
      }
      if (missing_reports != d) {
        for (j in 1:p[i]) {
          reported_cases[i] += obs[P_index + j];
        }
      }
      if (missing_reports > 0) {
        for (j in 1:missing_reports) {
          reported_cases[i] += predicted_counts[predicted_index + j];
        }
        predicted_index += missing_reports;
      }
    }
    return reported_cases;
  }
