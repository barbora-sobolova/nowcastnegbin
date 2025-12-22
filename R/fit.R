#' Fit the nowcasting model in STAN
#'
#' @description Fit the compiled STAN model and return the draws of selected
#' parameters in a list.
#'
#' @param compiled_model a compiled STAN model.
#' @param stan_data a list of data and parameters accepted by the STAN model
#' returned by the `get_stan_data()` function
#' @param model_obs an integer indicating the observation model. 0 - Poisson,
#' 1 - NegBinX, 2 - NegBin2D, 3 - NegBin1D, 4 - NegBin2M, 5 - NegBin1M.
#' @param stan_settings a list of STAN settings
#'
#' @return list of the data frames with the MCMC draws of different quantities:
#' \describe{
#'   \item{\code{nowcast}}{samples from the nowcasting distribution,}
#'   \item{\code{lambda}}{samples of the mean incidence trajectory,}
#'   \item{\code{delay_prob}}{samples of the delay probability vector,}
#'   \item{\code{diagnostics}}{a diagnostic summary of the Markov chains,}
#'   \item{\code{nb_size}}{The draws of the size parameter of the negative
#'   binomial distribution. Not applicable for the Poisson model.}
#'  }
#'
#' @import dplyr
#' @importFrom tidybayes gather_draws
#'
#' @export
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
    tidybayes::gather_draws(nowcast[week]) |> # nolint
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
    tidybayes::gather_draws(lambda[week]) |> # nolint
    mutate(Distribution = model_obs) |>
    dplyr::select(-".variable") |>
    ungroup()
  # Extract the delay probabilities
  df_delay_prob <- fitted_model |>
    tidybayes::gather_draws(reporting_delay[delay]) |> # nolint
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
      tidybayes::gather_draws(nb_size[1]) |> # nolint
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
