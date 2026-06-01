test_that("Model output from the MCMC method is stable and plausible", {
  # Compile the model
  mod <- cmdstanr::cmdstan_model(system.file(
    "stan",
    "nowcast.stan",
    package = "nowcastnegbin"
  ))

  # Example parameters
  params <- list()
  params$log_lambda0 <- log(100)
  params$rw_sd <- 0.01
  lgt <- 100
  params$max_lag <- 3
  params$probs <- c(0.5, 0.3, 0.2)
  params$nb_size <- 1.5
  # Include skipping one "Christmas" week
  params$skip_rows <- 50
  # Prior parameters of the reporting delay distribution
  params$prior_delay_param <- c(1, 1, 1)

  stan_settings <- list(
    seed = 123,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    show_messages = FALSE,
    show_exceptions = FALSE,
    refresh = 0
  )

  # Loop over the model types
  model_names <- c(
    "Poisson",
    "NegBinX",
    "NegBin2D",
    "NegBin1D",
    "NegBin2M",
    "NegBin1M"
  )
  for (model_obs in 0:5) {
    # Generate data
    set.seed(123456)
    obs_full <- with(
      params,
      generate_reports(
        lgt,
        max_lag,
        probs,
        log_lambda0,
        rw_sd,
        nb_size,
        model = model_names[model_obs + 1]
      )
    )
    data_list <- get_stan_data(
      obs_full$reports,
      prior_delay_param = params$prior_delay_param,
      skip_rows = params$skip_rows
    )

    # Run sampling with fixed seed
    fit <- fit_stan_model(
      mod$sample,
      data_list,
      model_obs,
      date_of_the_nowcast = as.Date("2024-12-09"),  # Arbitrary date
      mean_log = 0,
      sd_log = 1.5,
      stan_settings
    )

    # Extract the quantiles
    probs_sampled <- fit$delay_prob |>
      group_by(delay) |>
      summarise(
        quantile_2.5 = quantile(.value, 0.025),
        quantile_97.5 = quantile(.value, 0.975)
      )
    lambda_sampled <- fit$lambda |>
      dplyr::filter(week >= lgt - 1) |>
      group_by(week) |>
      summarise(
        quantile_2.5 = quantile(.value, 0.025),
        quantile_97.5 = quantile(.value, 0.975)
      )

    # Compare, whether the true value is inside the 95% CI
    if (model_obs != 0) {
      nb_size_sampled <- fit$nb_size |>
        # The `.value` column contains a list with one element per row
        mutate(.value = unlist(.value)) |>
        summarise(
          quantile_2.5 = quantile(.value, 0.025),
          quantile_97.5 = quantile(.value, 0.975)
        )
      expect_lt(nb_size_sampled$quantile_2.5, params$nb_size)
      expect_gt(nb_size_sampled$quantile_97.5, params$nb_size)
    }
    expect_true(all(probs_sampled$quantile_2.5 < params$probs))
    expect_true(all(probs_sampled$quantile_97.5 > params$probs))
    expect_true(
      all(lambda_sampled$quantile_2.5 < obs_full$exp_obs_total[lgt - 1:0])
    )
    expect_true(
      all(lambda_sampled$quantile_97.5 > obs_full$exp_obs_total[lgt - 1:0])
    )
  }
})

test_that("Model output from the GLM method is stable and plausible", {
  # Example parameters
  params <- list()
  params$log_lambda0 <- log(100)
  params$rw_sd <- 0.01
  lgt <- 100
  params$max_lag <- 3
  params$probs <- c(0.5, 0.3, 0.2)
  params$nb_size <- 1.5
  # Include skipping one "Christmas" week
  params$skip_rows <- 50

  # Loop over the model types
  model_names <- get_model_names()
  for (model_obs in 0:3) {
    # Generate data
    set.seed(123456)
    obs_full <- with(
      params,
      generate_reports(
        lgt,
        max_lag,
        probs,
        log_lambda0,
        rw_sd,
        nb_size,
        model = model_names[model_obs + 1]
      )
    )

    # Fit the model
    fit <- fit_glm_model(
      stan_data = get_stan_data(obs_full$reports, skip_rows = params$skip_rows),
      date_of_the_nowcast = as.Date("2024-12-09"),  # Arbitrary date
      model_name = model_names[model_obs + 1]
    )

    # Extract the quantiles
    probs_sampled <- fit$delay_prob |>
      group_by(delay) |>
      summarize(
        quantile_2.5 = quantile(.data$.value, probs = 0.025),
        quantile_97.5 = quantile(.data$.value, probs = 0.975)
      )
    lambda_sampled <- fit$lambda |>
      dplyr::filter(week > lgt - params$max_lag + 1) |>
      group_by(week) |>
      summarize(
        quantile_2.5 = quantile(.data$.value, probs = 0.025),
        quantile_97.5 = quantile(.data$.value, probs = 0.975)
      )

    # Compare, whether the true value is inside the 95% CI.
    if (model_obs != 0) {
      nb_size_sampled <- c(
        quantile(fit$nb_size$.value, probs = 0.025),
        quantile(fit$nb_size$.value, probs = 0.975)
      )
      expect_lt(nb_size_sampled[1], params$nb_size)
      expect_gt(nb_size_sampled[2], params$nb_size)
    }
    expect_true(all(probs_sampled$quantile_2.5 < params$probs))
    expect_true(all(probs_sampled$quantile_97.5 > params$probs))
    expect_true(
      all(lambda_sampled$quantile_2.5 < obs_full$exp_obs_total[lgt - 1:0])
    )
    expect_true(
      all(lambda_sampled$quantile_97.5 > obs_full$exp_obs_total[lgt - 1:0])
    )
  }
})
