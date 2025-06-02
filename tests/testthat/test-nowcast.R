test_that("Model sampling output is stable (snapshot)", {

  # Example data
  lgt <- 10
  max_lag <- 3
  obs_full <- c(5, 2, 0, 3, 3, 0, 6, 0, 1, 9, 0, 0, 11, 1, 1, 5, 1, 0, 2, 3, 1,
                9, 10, 0, 2, 6, 1, 1, 13, 0) |>
    matrix(nrow = lgt, ncol = max_lag, byrow = TRUE)
  index_mat <- lower.tri(obs_full, diag = TRUE)[10:1, ]
  obs_flat <- t(obs_full)[t(index_mat)]

  # Compile the model
  mod <- cmdstanr::cmdstan_model(system.file(
    "stan",
    "nowcast.stan",
    package = "nowcastnegbin"
  ))

  data_list <- list(
    n = nrow(index_mat),
    m = length(obs_flat),
    p = apply(index_mat, 1, sum),
    obs = obs_flat,
    d = max_lag
  )

  # Loop over the model types
  for (model_obs in 0:5) {
    # Run sampling with fixed seed
    fit <- mod$sample(
      data = c(data_list, model_obs = model_obs),
      seed = 123,
      chains = 2,
      iter_sampling = 1000,
      iter_warmup = 1000,
      refresh = 0,
      show_messages = FALSE,
      show_exceptions = FALSE
    )

    # Extract summary as output to snapshot
    summary_df <- fit$summary()

    # You can snapshot a subset if output is too large
    snapshot_output <- summary_df[
      summary_df$variable %in% paste0("nowcast[", 1:lgt, "]"),
      c("mean", "median", "sd", "mad", "q5", "q95")
    ]

    # Snapshot the output
    expect_snapshot_output(snapshot_output)
  }
})
