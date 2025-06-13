test_that("Model sampling output is stable and plausible", {
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

  index_mat <- lower.tri(
    matrix(nrow = lgt, ncol = params$max_lag),
    diag = TRUE
  )[lgt:1, ]
  data_list <- list(
    n = nrow(index_mat),
    m = sum(index_mat),
    p = apply(index_mat, 1, sum),
    d = params$max_lag
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
        log_lambda0,
        rw_sd,
        probs,
        nb_size,
        model = model_names[model_obs + 1]
      )
    )
    obs_flat <- t(obs_full$reports)[t(index_mat)]

    # Run sampling with fixed seed
    fit <- mod$sample(
      data = c(data_list, model_obs = model_obs, obs = list(obs_flat)),
      seed = 123,
      parallel_chains = 4,
      iter_sampling = 1000,
      iter_warmup = 1000,
      refresh = 0,
      show_messages = FALSE,
      show_exceptions = FALSE
    )

    # Extract the quantiles
    probs_sampled <- fit$draws(
      variables = paste0("reporting_delay[", 1:params$max_lag, "]")
    ) |>
      apply(3, quantile, probs = c(0.025, 0.975))
    lambda_sampled <- fit$draws(variables = paste0("lambda[", lgt - 1:0, "]")) |>
      apply(3, quantile, probs = c(0.025, 0.975))

    # Compare, whether the true value is inside the 95% CI
    if (model_obs != 0) {
      nb_size_sampled <- fit$draws(variables = paste0("nb_size[1]")) |>
        apply(3, quantile, probs = c(0.025, 0.975))
      expect_lt(nb_size_sampled[1], params$nb_size)
      expect_gt(nb_size_sampled[2], params$nb_size)
    }
    expect_true(all(probs_sampled[1, ] < params$probs))
    expect_true(all(probs_sampled[2, ] > params$probs))
    expect_true(all(lambda_sampled[1, ] < obs_full$exp_obs_total[lgt - 1:0]))
    expect_true(all(lambda_sampled[2, ] > obs_full$exp_obs_total[lgt - 1:0]))
    # expect_lt(lambda_sampled[1], obs_full$exp_obs_total[lgt])
    # expect_gt(lambda_sampled[2], obs_full$exp_obs_total[lgt])
  }
})
