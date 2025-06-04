test_that("moments as expected", {
  models <- c("Poisson", "NegBinX", "NegBin2D", "NegBin1D", "NegBin2M",
              "NegBin1M")
  lgt <- 1e5
  d <- 3
  # No noise => constant mean over time
  rw_noise_sd <- 0
  # initial mean number of cases
  log_lambda0 <- log(10)
  probs <- c(0.55, 0.15, 0.3)
  # size of the NegBin distribution
  r <- 0.8
  exp_obs <- exp(log_lambda0)

  # Theoretical means and covariances
  means_theo <- exp_obs * probs
  cov_theo <- list(
    Poisson = diag(means_theo),
    NegBinX = diag(means_theo * (1 + means_theo / r)),
    NegBin2D = diag(means_theo * (1 + exp_obs / r)),
    NegBin1D = diag(means_theo * (1 + 1 / r)),
    NegBin2M = probs %*% t(probs) * exp_obs^2 / r + diag(probs) * exp_obs,
    NegBin1M = probs %*% t(probs) * exp_obs / r + diag(probs) * exp_obs
  )

  for (k in seq_along(models)) {
    trajectory <- generate_reports(
      lgt,
      length(probs),
      log_lambda0,
      rw_noise_sd,
      probs,
      nb_size = r,
      model = models[k]
    )
    means_observed <- apply(trajectory$reports, 2, mean)
    cov_observed <- cov(trajectory$reports)
    expect_equal(means_observed, means_theo, tolerance = 0.1)
    expect_equal(cov_observed, cov_theo[[k]], tolerance = 1.5)
  }
})
