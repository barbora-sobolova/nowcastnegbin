#' Fit multiple MCMC nowcasting models
#'
#' @description This is a wrapper around `fit_stan_model()` that allows fitting
#' for multiple observation models and also possibly multiple time points. The
#' purpose of this function is to return the fits from different observation
#' models for one date in a single bundle to simplify grouping structures in the
#' targets pipeline.
#'
#' @param compiled_model a compiled STAN model.
#' @param stan_data a list of data and parameters accepted by the STAN model
#' returned by the `get_stan_data()` function
#' @param date_of_the_nowcast a vector of dates, when the nowcast is made
#' @param model_obs an integer vector indicating the observation model. Must be
#' of the same length as `date_of_the_nowcast`. 0 - Poisson,
#' 1 - NegBinX, 2 - NegBin2D, 3 - NegBin1D, 4 - NegBin2M, 5 - NegBin1M.
#' @param stan_settings a list of STAN settings
#'
#' @return list of the data frames with the MCMC draws of different quantities:
#' \describe{
#'   \item{\code{nowcast}}{samples from the nowcasting distribution,}
#'   \item{\code{lambda}}{samples of the mean incidence trajectory,}
#'   \item{\code{delay_prob}}{samples of the delay probability vector,}
#'   \item{\code{rw_sd}}{samples of the standard deviation of the random walk,}
#'   \item{\code{diagnostics}}{a diagnostic summary of the Markov chains,}
#'   \item{\code{nb_size}}{The draws of the size parameter of the negative
#'   binomial distribution. Not applicable for the Poisson model.}
#'  }
#'
#' @importFrom dplyr bind_rows
#' @importFrom purrr map
#'
#' @export
fit_all_stan_models <- function(
  compiled_model,
  stan_data,
  model_obs,
  date_of_the_nowcast,
  stan_settings
) {
  fits <- vector("list", length(model_obs))
  for (k in seq_along(model_obs)) {
    fits[[k]] <- fit_stan_model(
      compiled_model,
      stan_data,
      model_obs[k],
      # The grouping structure ensures that date is identical for all
      # observation models, so we could also pass just `date_of_the_nowcast[1]`
      # each time
      date_of_the_nowcast[k],
      stan_settings
    )
  }
  ret_list <- list(
    nowcast = bind_rows(purrr::map(fits, "nowcast")),
    lambda = bind_rows(purrr::map(fits, "lambda")),
    delay_prob = bind_rows(purrr::map(fits, "delay_prob")),
    rw_sd = bind_rows(purrr::map(fits, "rw_sd")),
    diagnostics = bind_rows(purrr::map(fits, "diagnostics")),
    nb_size = bind_rows(purrr::map(fits, "nb_size"))
  )
  ret_list
}

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
#' @param date_of_the_nowcast a date, when the nowcast is made to add as a
#' column to the data frame with results
#' @param stan_settings a list of STAN settings
#'
#' @return list of the data frames with the MCMC draws of different quantities:
#' \describe{
#'   \item{\code{nowcast}}{samples from the nowcasting distribution,}
#'   \item{\code{lambda}}{samples of the mean incidence trajectory,}
#'   \item{\code{delay_prob}}{samples of the delay probability vector,}
#'   \item{\code{rw_sd}}{samples of the standard deviation of the random walk,}
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
  date_of_the_nowcast,
  stan_settings
) {
  # A helper function that adds the model names and date of the nowcast into
  # a data frame with parameter samples
  add_meta <- function(df) {
    mutate(
      df,
      # Observation models are numbered from zero
      Distribution = get_model_names()[model_obs + 1],
      nowcast_date = date_of_the_nowcast
    )
  }
  # Fit the model
  fitted_model <- do.call(
    compiled_model,
    c(data = list(c(stan_data, model_obs = model_obs)), stan_settings)
  )
  # Extract the diagnostic summary
  diagnostics <- fitted_model$diagnostic_summary() |>
    suppressMessages() |>
    as.data.frame() |>
    mutate(seed = stan_settings$seed, nowcast_date = date_of_the_nowcast) |>
    # Add the information about the sampling duration
    cbind(fitted_model$time()$chains)
  # Refit the model, if we get too many divergent transitions, or the ebfmi is
  # low in at least one chain.
  refit <- 0
  while (
    (any(diagnostics$num_divergent >= 100) || any(diagnostics$ebfmi < 0.3)) &&
      refit < 3
  ) {
    # Shift the seed and refit
    stan_settings$seed <- stan_settings$seed + 1
    refitted_model <- do.call(
      compiled_model,
      c(data = list(c(stan_data, model_obs = model_obs)), stan_settings)
    )
    diagnostics_refit <- refitted_model$diagnostic_summary() |>
      suppressMessages() |>
      as.data.frame()
    refit <- refit + 1
    # If the diagnostics of the refitted model is better than for the first fit,
    # store the refit
    store_refit <- (
      # Accept if divergences improve and ebfmi doesn't get critically worse
      max(diagnostics$num_divergent) > max(diagnostics_refit$num_divergent) &&
        min(diagnostics$ebfmi) - min(diagnostics_refit$ebfmi) < 0.1
    ) ||
      (
        # Or if ebfmi improves and divergences don't get critically worse
        min(diagnostics$ebfmi) < min(diagnostics_refit$ebfmi) &&
          max(diagnostics$num_divergent) -
            max(diagnostics_refit$num_divergent) > -50
      )
    if (store_refit) {
      fitted_model <- refitted_model
      diagnostics <- diagnostics_refit |>
        mutate(seed = stan_settings$seed, nowcast_date = date_of_the_nowcast) |>
        cbind(fitted_model$time()$chains)
    }
  }

  # Extract the nowcasts
  df_nowcast <- fitted_model |>
    tidybayes::gather_draws(nowcast[week]) |> # nolint
    ungroup() |>
    add_meta() |>
    # Keep only the counts that needed correction, which are those located at
    # the last `max_lag - 1`. we can calculate the last positions using the
    # maximum lag `d` and the number of reporting triangle rows `n` from the
    # STAN data
    dplyr::filter(week > stan_data$n - stan_data$d + 1) |>
    dplyr::select(-".variable")
  # Extract the estimates of the expected counts
  df_lambda <- fitted_model |>
    tidybayes::gather_draws(lambda[week]) |> # nolint
    ungroup() |>
    add_meta() |>
    dplyr::select(-".variable")
  # Extract the delay probabilities
  df_delay_prob <- fitted_model |>
    tidybayes::gather_draws(reporting_delay[delay]) |> # nolint
    ungroup() |>
    add_meta() |>
    dplyr::select(-".variable")
  # Extract the standard error of the random walk
  df_rw_sd <- fitted_model |>
    tidybayes::gather_draws(rw_sd) |> # nolint
    ungroup() |>
    add_meta() |>
    dplyr::select(-".variable")
  # Add the seed and the model number to the diagnostic summary
  diagnostics <- diagnostics |> add_meta()
  # Return the draws as a list
  ret_list <- list(
    nowcast = df_nowcast,
    lambda = df_lambda,
    delay_prob = df_delay_prob,
    rw_sd = df_rw_sd,
    diagnostics = diagnostics
  )
  # Extract the draws of the negative binomial size parameter if we don't fit
  # the Poisson model
  if (model_obs != 0) {
    df_nb_size <- fitted_model |>
      tidybayes::gather_draws(nb_size[1]) |> # nolint
      ungroup() |>
      add_meta() |>
      tidyr::unnest(.data$.value) |>
      dplyr::select(-".variable")
  } else {
    df_nb_size <- NULL
  }
  # return the draws as a list
  ret_list <- c(ret_list, list(nb_size = df_nb_size))
  ret_list
}

#' Select the distributional family of a GAMLSS model
#'
#' @param model_name string, name of the observational model to be fit using the
#' GAMLSS framework
#' @return list containing the GAMLSS family and the formula for its scale
#' parameter sigma
select_gamlss_model <- function(
  model_name = c("Poisson", "NegBinX", "NegBin2D", "NegBin1D")
) {
  # Family as defined in the gamlss package. Note that NBII() denotes what we
  # call NegBin1 - a negative binomial distribution with linear
  # mean-variance relationship. Conversely NBI() denotes what is NegBin2 for us,
  # i.e. a negative binomial distribution with quadratic mean-variance
  # relationship.
  family <- switch(
    model_name,
    Poisson = gamlss.dist::PO(),
    NegBinX = gamlss.dist::NBI(),
    NegBin2D = gamlss.dist::NBI(),
    NegBin1D = gamlss.dist::NBII()
  )
  # Sigma formula is different for NegBin2D, where we have to multiply the
  # overdispersion parameter by the delay probability, in order for the
  # negative binomial distribution to be preserved.
  sigma_formula <- switch(
    model_name,
    Poisson = NULL,
    NegBinX = ~ 1,
    NegBin2D = ~ delay,
    NegBin1D = ~ 1
  )
  list(family = family, sigma_formula = sigma_formula)
}

#' Generate the nowcasts based on a fitted GLM
#'
#' @description This function implements the nowcasting method by
#' van de Kassteele (2019). It generates the nowcasts from the predictive
#' distribution and returns them together with the model parameters sampled from
#' the multivariate normal distribution. The smooth term and the categorical
#' terms are sampled independently due to the limitations of \code{gamlss2()},
#' which does not return the cross-correlation between the smooth terms and the
#' fixed terms.
#'
#' @param fitted_gamlss_obj a fitted GLM model using the \code{gamlss2()}
#' function
#' @param spline_basis a matrix with \code{t_len} rows
#' @param glm_data_all the data frame used for fitting the GLM model including
#' the entries dropped from the likelihood. The dropped, but already observed
#' entries are needed to be added to the predicted part, when calculating the
#' nowcast.
#' @param t_len an integer indicating how many time units are spanned by the
#' data
#' @param max_lag an integer, the maximum reporting delay. Here, they are
#' numbered from 1, even though the initial delay is often regarded as the
#' zeroth delay.
#' @param model_name model_name a string indicating the observation model. One
#' of "Poisson", "NegBinX", "NegBin2D" and "NegBin1D".
#' @param date_of_the_nowcast a date, when the nowcast is made to add as a
#' column to the data frame with results
#' @param n_samples how many samples from the nowcasting distribution we draw
#'
#' @return List of the data frames with the draws of different model parameters:
#' \describe{
#'   \item{\code{nowcast}}{samples from the nowcasting distribution,}
#'   \item{\code{lambda}}{samples of the mean incidence trajectory,}
#'   \item{\code{delay_prob}}{samples of the delay probability vector,}
#'   \item{\code{nb_size}}{The draws of the size parameter of the negative
#'   binomial distribution. Not applicable for the Poisson model.}
#'  }
#' Additionally, the number of iterations \code{iter} necessary for fitting the
#' model is returned in the list.
#'
#' @references
#' van de Kassteele J (2019).
#' “Nowcasting the number of new symptomatic cases during infectious disease
#' outbreaks using constrained p-spline smoothing.”
#' \emph{Epidemiology}, 30(5), 737–745.
#' \href{
#' https://doi.org/10.1097/EDE.0000000000001050
#' }{doi:10.1097/EDE.0000000000001050}
#'
#' @import dplyr
#' @importFrom MASS mvrnorm
#'
#' @export
generate_glm_nowcasts <- function(
  fitted_gamlss_obj,
  spline_basis,
  glm_data_all,
  t_len,
  max_lag,
  model_name,
  date_of_the_nowcast,
  n_samples = 4000
) {
  # Create nowcasts using the sampling following van de Kasstelee 2019:
  # (1) Draw the parameter values from the multivariate normal distribution
  #     using the estimates and its standard errors.
  # (2) Draw a sample from the observation model based on the parameter values.
  # (3) Repeat (1)-(2) until we have a sample of the same size as from the STAN
  #     model.
  # (4) Calculate the median and the quantiles.
  # For simplicity, we sample the parameter values separately for the fixed
  # effects and for the smooth term. This approximation is valid, since
  # empirically, the cross-correlation between these coefficients is low in our
  # case studies.

  # Smooth term
  vcov_smooth <- fitted_gamlss_obj$fitted.specials$mu$`s(week)`$vcov
  coeffs_smooth <- fitted_gamlss_obj$fitted.specials$mu$`s(week)`$coefficients
  sampled_pars_smooth <- MASS::mvrnorm(
    n_samples,
    mu = coeffs_smooth,
    Sigma = vcov_smooth
  )

  # Delays and also the overdispersion parameter in the case of a negbin model
  vcov_fixed <- vcov(fitted_gamlss_obj)
  coeffs_fixed <- unlist(fitted_gamlss_obj$coefficients)
  sampled_pars_fixed <- MASS::mvrnorm(
    n_samples,
    mu = coeffs_fixed,
    Sigma = vcov_fixed
  )
  # Create indices for selecting the columns corresponding to the mean value
  # and to the scale parameter
  mu_ind <- grep("mu", colnames(sampled_pars_fixed))
  sigma_ind <- grep("sigma", colnames(sampled_pars_fixed))

  # We fit the model with an intercept, so we have to calculate the linear
  # predictor using the right "contrasts". The first factor level is the
  # reference, so the linear predictor for the fixed part is constructed in the
  # standard way, that is, the effect for the first level is beta_0, for the
  # second level, it's beta_0 + beta_1 and so on.
  contrasts_mu <- diag(max_lag)
  contrasts_mu[1, ] <- 1

  # To do the prediction, we need to calculate the mean of each reporting
  # triangle cell. We expand the time with the delay first and then do an inner
  # join with the corresponding values (there will be `n_samples` of them). This
  # way we ensure that the smooth and fixed terms are assigned to the correct
  # time and delay. We could include only the last part of the reporting table,
  # where we want to do the prediction, but will calculate everything to be able
  # to return the whole mean process (lambda_t).
  df_skeleton_grid <- tidyr::expand_grid(
    week = seq_len(t_len),
    delay = seq_len(max_lag)
  )

  # Calculate the smooth curve using the sampled parameters. This will not
  # equal the estimate of the mean process, which we call lambda_t, since there
  # is still a multiplicative constant absorbed by the fixed terms.
  df_smooth_sampled <- data.frame(
    # The resulting matrix has dimensions n_samples x t_len. Concatenation
    # is done by column, therefore the `week` column is defined using the `each`
    # parameter of rep().
    smooth_lpred = c(sampled_pars_smooth %*% t(spline_basis)),
    week = rep(seq_len(t_len), each = n_samples),
    # What sample/draw number this is. We use the ".draw" name to match the
    # column names of the output from the `fit_stan_model()` function.
    .draw = rep(seq_len(n_samples), times = t_len)
  )

  # Calculate the fixed terms using the sampled parameters.
  df_fixed_sampled <- data.frame(
    # The resulting matrix has dimensions n_samples x max_lag. Concatenation
    # is done by column, therefore the `delay` column is defined using the
    # `each` parameter of rep().
    fixed_lpred = c(sampled_pars_fixed[, mu_ind] %*% contrasts_mu),
    delay = rep(seq_len(max_lag), each = n_samples),
    # What sample/draw number this is.
    .draw = rep(seq_len(n_samples), times = max_lag)
  ) |>
    # Extract the probabilities from the sampled fixed parameters by
    # exponentiating and standardizing
    group_by(.data$.draw) |>
    mutate(
      probs_sampled = exp(.data$fixed_lpred) / sum(exp(.data$fixed_lpred))
    ) |>
    ungroup()

  # Join the fixed and smooth terms together to calculate the mean of each cell
  # of the reporting triangle
  df_mu_sampled <- dplyr::inner_join(
    df_skeleton_grid,
    df_smooth_sampled,
    by = "week",
    # `df_all_samped` has `t_len` x `max_lag` rows. `df_smooth_sampled` has
    # `t_len` x `n_samples` rows and the joint data frame is supposed to have
    # `t_len` x `max_lag` x `n_samples` rows.
    relationship = "many-to-many"
  ) |>
    inner_join(
      df_fixed_sampled,
      by = c("delay", ".draw"),
      # `df_fixed_sampled` has `max_lag` x `n_samples` and the joint data frame
      # is supposed to have `t_len` x `max_lag` x `n_samples` rows.
    ) |>
    mutate(
      # Calculate the mean of the counts in each cell of the reporting triangle.
      mu = exp(.data$smooth_lpred + .data$fixed_lpred),
      # Add the code of the model to match the `fit_stan_model()` output
      Distribution = model_name
    ) |>
    # Calculate the value of lambda for each time point. We need to multiply the
    # sum of `mu` across delays by a constant, that got absorbed by the factor
    # terms
    group_by(.data$.draw, .data$week) |>
    mutate(
      lambda_sampled = exp(.data$smooth_lpred) * sum(exp(.data$fixed_lpred))
    ) |>
    ungroup()

  # Calculate the overdispersion parameter for different models. We calculate
  # two different values. Firstly, we calculate the sampled value of the size
  # parameter which is used in the theoretical parametrization of the negative
  # binomial distribution. This value does not depend on time or the delay.
  # Secondly, we prepare the value, which is used in `rnbinom()` to sample the
  # counts. Here we need to adjust it accordingly based on time and delay and
  # it will not be a part of returned results.
  if (model_name %in% c("NegBinX", "NegBin1D")) {
    df_nb_size <- data.frame(
      Distribution = model_name,
      nowcast_date = date_of_the_nowcast,
      # The theoretical overdispersion parameter is identical for all times and
      # delays.
      nb_size = rep(
        exp(-sampled_pars_fixed[, sigma_ind]),
        times = t_len * max_lag
      ),
      # Construct the rest of the columns so that the data frame has the same
      # amount of rows as `df_mu_sampled`.
      week = rep(seq_len(t_len), each = max_lag * n_samples),
      delay = rep(rep(seq_len(max_lag), each = n_samples), times = t_len),
      .draw = rep(seq_len(n_samples), times = max_lag * t_len)
    ) |>
      mutate(
        # The overdispersion parameter to pass to `rnbinom()`. For NegBin1D we
        # will have to multiply it by the mean value later.
        nb_size_internal = .data$nb_size
      )
  } else if (model_name == "NegBin2D") {
    # For NegBin2D, we have multiple columns of the sampled parameters for the
    # scale parameter. The model with intercept is constructed in the same way
    # as the fixed terms for the mean value, so we can use the same contrasts.
    df_nb_size <- data.frame(
      Distribution = model_name,
      nowcast_date = date_of_the_nowcast,
      # Matrix of size `n_samples` x `max_delay`
      nb_size_internal = rep(
        c(exp(-sampled_pars_fixed[, sigma_ind] %*% contrasts_mu)),
        times = t_len
      ),
      # Construct the rest of the columns so that the data frame has the same
      # amount of rows as `df_mu_sampled`.
      week = rep(seq_len(t_len), each =  max_lag * n_samples),
      delay = rep(rep(seq_len(max_lag), each = n_samples), times = t_len),
      .draw = rep(seq_len(n_samples), times = max_lag * t_len)
    ) |>
      group_by(.data$week, .data$.draw) |>
      mutate(
        nb_size = sum(.data$nb_size_internal)
      ) |>
      ungroup()
  }

  # How many count samples we generate. This is the number of parameter samples
  # times the number of missing cells of the reporting triangle.
  n_cells_to_predict <-  n_samples * sum(seq_len(max_lag - 1))
  # Sample from the count distribution. For Poisson, we need only the mean and
  # can we sample directly from the Poisson distribution. For the negative
  # binomial counts, we have to join the mean part and the size part.
  if (model_name == "Poisson") {
    df_predicted <- df_mu_sampled |>
      filter(.data$week + .data$delay > t_len + 1) |>
      mutate(counts = rpois(n_cells_to_predict, .data$mu))
    # Create the return list, which is empty for Poisson at this point.
    ret_list <- list()
  } else {
    # Join the mu part and the scale part.
    df_predicted <- inner_join(
      df_mu_sampled,
      df_nb_size,
      by = c("week", "delay", ".draw", "Distribution"),
      relationship = "one-to-one"
    ) |>
      # For NegBin1D, we have to multiply the size parameter by the expectation
      # before we pass it to `rnbinom()`.
      mutate(
        nb_size_internal = if (model_name == "NegBin1D") {
          .data$nb_size_internal * .data$mu
        } else {
          .data$nb_size_internal
        }
      ) |>
      # Sample the counts only for the days we are predicting. We predict for
      # all days, where week + delay > last week + 1, since the delay
      filter(.data$week + .data$delay > t_len + 1) |>
      mutate(
        counts = rnbinom(
          n_cells_to_predict,
          mu = .data$mu,
          size = .data$nb_size_internal
        )
      )
    # Create the return list, which contains the data frame with the samples of
    # the size of the negative binomial distribution
    ret_list <- list(
      nb_size = df_nb_size |>
        select("nb_size", ".draw", "Distribution") |>
        # Remove the duplicated values of the overdispersion parameter
        unique() |>
        rename(".value" = "nb_size") |>
        mutate(".variable" = "nb_size", nowcast_date = date_of_the_nowcast)
    )
  }

  # Aggregate the counts by delay for both the predicted counts and the
  # observed counts. The observed counts have to contain the observations that
  # were skipped due to the reporting anomaly around Christmas, in order to
  # include them in the nowcast.
  df_obs <- glm_data_all |>
    group_by(.data$week) |>
    summarize(
      obs_counts = sum(.data$obs)
    ) |>
    # Keep only the time points, where we do the nowcasting
    tail(max_lag - 1)
  df_nowcast <- df_predicted |>
    group_by(.data$week, .data$.draw, .data$Distribution) |>
    summarize(
      predicted_counts = sum(.data$counts),
      .groups = "drop"
    ) |>
    # Join with the observed counts
    inner_join(
      df_obs,
      by = "week",
      relationship = "many-to-one"
    ) |>
    mutate(
      # Make the column names match those produced by the STAN
      # procedure in `fit_stan_model()`
      .value = .data$obs_counts + .data$predicted_counts,
      .variable = "nowcast",
      nowcast_date = date_of_the_nowcast
    ) |>
    # Drop redundant columns
    select(-c("obs_counts", "predicted_counts"))

  # Arrange all sampled parameters into the return list
  ret_list <- c(
    ret_list,
    list(
      nowcast = df_nowcast,
      lambda = df_mu_sampled |>
        select("week", "lambda_sampled", ".draw", "Distribution") |>
        # Remove duplicate lambda samples
        unique() |>
        rename(".value" = "lambda_sampled") |>
        mutate(".variable" = "lambda", nowcast_date = date_of_the_nowcast),
      delay_prob = df_mu_sampled |>
        select("delay", "probs_sampled", ".draw", "Distribution") |>
        # Remove duplicate samples of the delay probabilities
        unique() |>
        rename(".value" = "probs_sampled") |>
        mutate(
          ".variable" = "reporting_delay",
          nowcast_date = date_of_the_nowcast
        ),
      # Store the information about the number of iterations of the `gamlss2()`
      # optimizing function.
      iter = fitted_gamlss_obj$iterations
    )
  )
  ret_list
}

#' Fit the nowcasting model using GLMs
#'
#' @description This function fits a generalized linear model to the counts and
#' returns draws of selected model parameters in a list. The parameters here are
#' NOT sampled using MCMC. Instead, they are generated from the multivariate
#' normal distribution based on the regression model output - parameter
#' estimates and their standard errors - to account for uncertainty in the
#' estimation. This method is taken from van de Kassteele (2019). The model is
#' fit using the \href{https://github.com/gamlss-dev/gamlss2}{\code{gamlss2}}
#' package that implements the GAMLSS (GAM for location scale and shape)
#' framework by Stasinopoulos and Rigby (2007).
#'
#' @param stan_data a list of data and parameters accepted by the STAN model
#' returned by the \code{get_stan_data()} function. This list is reused here to
#' create a data frame for the regression model.
#' @param date_of_the_nowcast a date, when the nowcast is made to add as a
#' column to the data frame with results
#' @param model_name a string indicating the observation model. One of
#' "Poisson", "NegBinX", "NegBin2D" and "NegBin1D".
#' @param n_samples how many samples from the nowcasting distribution we draw
#'
#' @details We fit a simple generalized regression model, that has one smooth
#' term for the mean process of the total counts and that takes the delay as a
#' categorical covariate. The smooth term is constructed using \eqn{N} number
#' of P-spline bases. The model can be written as:
#'
#' \deqn{\log(\mu_{t,d}) = \sum_{i = 1}^N\alpha_i s_i(t) + \beta_0 +
#' \sum_{j = 1}^D \beta_j \mathbb{I}(j = d),}
#'
#' The smooth term can be used to reconstruct the mean process of the total
#' counts up to a multiplicative constant, while the categorical terms can
#' be used to calculate the delay probabilities (again up to a multiplicative
#' constant). The model is fitted with the intercept, since this is required by
#' \code{gamlss2()} function.
#'
#' The Poisson, NegBinX and NegBin1D models are fitted using the \code{family}
#' parameter of the \code{gamlss2()} function. For NegBin2D, we fit
#' an additional regression model for the scale parameter:
#'
#' \deqn{\log(\sigma_{d}) = \gamma_0 +
#' \sum_{j = 1}^D \gamma_j \mathbb{I}(j = d),}
#'
#' The resulting model will not strictly be the NegBin2D model, since the
#' \eqn{\beta} and \eqn{\gamma} parameters can't be enforced to be identical.
#' However, this is an acceptable approximation of the model.
#'
#' @return List of the data frames with the draws of different model parameters:
#' \describe{
#'   \item{\code{nowcast}}{samples from the nowcasting distribution,}
#'   \item{\code{lambda}}{samples of the mean incidence trajectory,}
#'   \item{\code{delay_prob}}{samples of the delay probability vector,}
#'   \item{\code{nb_size}}{The draws of the size parameter of the negative
#'   binomial distribution. Not applicable for the Poisson model.}
#'  }
#' Additionally, the number of iterations \code{iter} necessary for fitting the
#' model is returned in the list. The columns of the data frames are a subset of
#' the columns of the corresponding data frame returned by
#' \code{fit_stan_model}. Columns `.chain` and `.iteration`, that are relevant
#' only for the MCMC sampling, are not present here. The return list is designed
#' similarly as the output of \code{fit_stan_model()} to facilitate the
#' processing and plotting of the results.
#'
#' @references
#' Stasinopoulos DM, Rigby RA (2007).
#' “Generalized Additive Models for Location Scale and Shape (GAMLSS) in R.”
#' \emph{Journal of Statistical Software}, 23(7), 1–46.
#' \href{https://doi.org/10.18637/jss.v023.i07}{doi:10.18637/jss.v023.i07}
#'
#' van de Kassteele J (2019).
#' “Nowcasting the number of new symptomatic cases during infectious disease
#' outbreaks using constrained p-spline smoothing.”
#' \emph{Epidemiology}, 30(5), 737–745.
#' \href{
#' https://doi.org/10.1097/EDE.0000000000001050
#' }{doi:10.1097/EDE.0000000000001050}
#'
#' @import dplyr
#' @importFrom gamlss2 gamlss2
#' @importFrom mgcv gam
#' @importFrom mgcv s
#'
#' @export
fit_glm_model <- function(
  stan_data,
  date_of_the_nowcast,
  model_name = c("Poisson", "NegBinX", "NegBin2D", "NegBin1D"),
  n_samples = 4000
) {
  model_name <- match.arg(model_name)

  # Select the gamlss.family and sigma.formula based on the model
  gamlss_specs <- select_gamlss_model(model_name)  # nolint

  # The list of data for the STAN model contains all the information to
  # construct the data frame for the GLM. We can drop some observations from the
  # likelihood by filtering according to the `idx_include` vector.
  glm_data_all <- data.frame(
    obs = stan_data$obs,
    delay = factor(sequence(stan_data$p)),
    week = rep(seq_len(stan_data$n), times = stan_data$p)
  )
  glm_data <- glm_data_all[stan_data$idx_include, ]

  # We use the length of the time period over 5 to determine the number of
  # spline spline basis functions based on van de Kasstelee 2019. The s() way
  # of defining the smooth term in the model is lifted from the mgcv package,
  # where the number of basis functions is k - 1.
  n_basis_functions <- round(stan_data$n / 5) + 1

  # Fit the GAMLSS model
  fit <- substitute(
    gamlss2(
      # We need to include the intercept, otherwise the design matrix will be
      # inconsistent between the Poisson and NegBin models.
      obs ~ s(week, k = n_basis) + delay,
      sigma.formula = gamlss_specs$sigma_formula,
      family = gamlss_specs$family,
      data = glm_data,
      # For some models (NegBin2D), it takes some time to converge
      maxit = 900,
      trace = FALSE
    ),
    list(n_basis = n_basis_functions)
  ) |> eval()

  # Fit the Poisson model using the gam() function from the mgcv package. This
  # is necessary to extract the spline basis, which is not returned by the
  # gamlss2() function. Since gamlss2() uses mgcv under the hood, the bases
  # are identical and we can use it to reconstruct the spline curve.
  mod_mgcv <- gam(
    obs ~ s(week, k = n_basis_functions) + delay,
    data = glm_data,
    family = poisson
  )
  smooth_coeffs_inds <- grep(pattern = "s()", names(mod_mgcv$coefficients))
  # Extract the basis using `predict.gam(type = "lpmatrix", ...)`. Often,
  # we can just extract the design matrix as is, but when we skip certain
  # observations, some weeks might not be represented in the data. For these
  # cases, we need to create a data frame with no gaps.
  basis <- mgcv::predict.gam(
    mod_mgcv,
    newdata = data.frame(
      week = seq_len(max((glm_data$week))),
      delay = 1
    ),
    type = "lpmatrix"
  )[, smooth_coeffs_inds]

  # Generate nowcasts by the van de Kasstelee 2019 method, using sampling from
  # the multivariate normal distribution for the parameters.
  ret_list <- generate_glm_nowcasts(
    fit,
    basis,
    glm_data_all,
    stan_data$n,
    stan_data$d,
    model_name,
    date_of_the_nowcast,
    n_samples
  )
  ret_list
}

#' Fit multiple GLM nowcasting models
#'
#' @description This is a wrapper around `fit_glm_model()` that allows fitting
#' for multiple observation models and also possibly multiple time points. The
#' purpose of this function is to return the fits from different observation
#' models for one date in a single bundle to simplify grouping structures in the
#' targets pipeline.
#'
#' @param stan_data a list of data and parameters accepted by the STAN model
#' returned by the \code{get_stan_data()} function. This list is reused here to
#' create a data frame for the regression model.
#' @param date_of_the_nowcast a vector of dates, when the nowcast is made
#' @param model_name a character vector indicating observation models to be
#' fit. Must be of the same length as `date_of_the_nowcast`
#' @param n_samples how many samples from the nowcasting distribution we draw
#'
#' @return List of the data frames with the draws of different model parameters
#' and additional quantities:
#' \describe{
#'   \item{\code{nowcast}}{samples from the nowcasting distribution,}
#'   \item{\code{lambda}}{samples of the mean incidence trajectory,}
#'   \item{\code{delay_prob}}{samples of the delay probability vector,}
#'   \item{\code{nb_size}}{The draws of the size parameter of the negative
#'   binomial distribution. Not applicable for the Poisson model,}
#'   \item{\code{iter}}{a scalar, the number of iterations necessary for fitting
#'   the gamlss model.}
#'  }
#' The columns of the data frames are a subset of the columns of the
#' corresponding data frame returned by \code{fit_stan_model}. Columns `.chain`
#' and `.iteration`, that are relevant only for the MCMC sampling, are not
#' present here. The return list is designed similarly as the output of
#' \code{fit_stan_model()} to facilitate the processing and plotting of the
#' results.
#'
#' @importFrom dplyr bind_rows
#' @importFrom purrr map
#'
#' @export
fit_all_glm_models <- function(
  stan_data,
  date_of_the_nowcast,
  model_name,
  n_samples = 4000
) {
  fits <- vector("list", length(model_name))
  for (k in seq_along(model_name)) {
    fits[[k]] <- fit_glm_model(
      stan_data,
      # The grouping structure ensures that date is identical for all
      # observation models, so we could also pass just `date_of_the_nowcast[1]`
      # each time
      date_of_the_nowcast[k],
      model_name[k],
      n_samples
    )
  }
  ret_list <- list(
    nowcast = bind_rows(purrr::map(fits, "nowcast")),
    lambda = bind_rows(purrr::map(fits, "lambda")),
    delay_prob = bind_rows(purrr::map(fits, "delay_prob")),
    nb_size = bind_rows(purrr::map(fits, "nb_size")),
    iter = unlist(purrr::map(fits, "iter"))
  )
  ret_list
}
