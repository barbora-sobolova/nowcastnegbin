fit_stan_model <- function(
  compiled_model,
  stan_data,
  model_obs,
  stan_settings
) {
  # Fit the model
  fitted_model <- do.call(
    compiled_model,
    c(data = list(c(stan_data, model_obs = model_obs)), stan_settings)
  )
  # Extract the nowcasts
  df_nowcast <- fitted_model |>
    tidybayes::gather_draws(nowcast[week]) |>
    mutate(Distribution = model_obs) |>
    # Keep only the counts that needed correction, which are those located at
    # the last `max_lag - 1`. we can calculate the last positions using the
    # maximum lag `d` and the number of reporting triangle rows `n` from the
    # STAN data
    dplyr::filter(week > stan_data$n - stan_data$d + 1) |>
    dplyr::select(-".variable") |>
    ungroup()
  # Extract the estimates of the expected counts
  df_lambda <- fitted_model |>
    tidybayes::gather_draws(lambda[week]) |>
    mutate(Distribution = model_obs) |>
    dplyr::select(-".variable") |>
    ungroup()
  # Extract the delay probabilities
  df_delay_prob <- fitted_model |>
    tidybayes::gather_draws(reporting_delay[delay]) |>
    mutate(Distribution = model_obs) |>
    dplyr::select(-".variable") |>
    ungroup()
  # Extract the diagnostic summary
  diagnostics <- fitted_model$diagnostic_summary() |>
    as.data.frame() |>
    mutate(Distribution = model_obs)
  # Return the draws as a list
  ret_list <- list(
    nowcast = df_nowcast,
    lambda = df_lambda,
    delay_prob = df_delay_prob,
    diagnostics = diagnostics
  )
  # Extract the draws of the negative binomial size parameter if we don't fit
  # the Poisson model
  if (model_obs != 0) {
    df_nb_size <- fitted_model |>
      tidybayes::gather_draws(nb_size[1]) |>
      mutate(Distribution = model_obs) |>
      dplyr::select(-".variable") |>
      ungroup()
  } else {
    df_nb_size <- NULL
  }
  # return the draws as a list
  ret_list <- c(ret_list, list(nb_size = df_nb_size))
  ret_list
}
