create_totals_data_frame <- function(
  obs_counts,
  max_lag,
  start_date,
  length_of_train_data
) {
  obs_mat <- as.matrix(obs_counts)
  # Data frame to plot the observations - final counts and counts available at
  # the time of creating the nowcast
  data.frame(
    counts = c(
      rowSums(obs_mat),
      rowSums(mock_unobserved(obs_mat), na.rm = TRUE)
    ),
    date = rep(start_date + (seq_len(length_of_train_data) - 1) * 7, 2),
    data = rep(c("Final", "Preliminary"), each = nrow(obs_counts))
  )
}
summarize_nowcast <- function(
  df_nowcast,
  df_total
) {
  # Recover the beginning and the length of the timeline from the total counts
  start_date <- min(df_total$date)
  length_of_train_data <- nrow(df_total)

  # Reformat the `week` column of the dat aframe with the nowcasts, so that it's
  # aligned with he actual date
  df_nowcast <- df_nowcast |>
    mutate(date = start_date + (week - 1) * 7) |>
    dplyr::select(-"week")
  # Quantiles of the MCMC sample to calculate - the median and quantiles for
  # constructing the 95% and 50% prediction interval
  quantiles_to_get <- c(50, 2.5, 25, 75, 97.5)
  # Summarize the sample of the nowcasts
  df_nowcast_plot <- df_total |>
    # The `df_total` data frame contains the final counts and the preliminary.
    # We need the final observed counts to evaluate the nowcasts via CRPS.
    filter(data == "Final") |>
    inner_join(df_nowcast, by = "date") |>
    group_by(date, Distribution) |>
    summarise(
      # Calculate the CRPS from the MCMC sample. We pass counts[1] as the true
      # observed value, since  this is the same value for each date.
      CRPS = scoringutils::crps_sample(
        observed = counts[1],
        predicted = .value
      ),
      # Calculate the mean and the quantiles
      mean = mean(.value),
      quantiles = list(
        as_tibble(
          as.list(
            quantile(.value, probs = quantiles_to_get / 100)
          )
        )
      )
    ) |>
    unnest(quantiles)
  # Rename the quantile columns to have nicer names
  cols_to_rename <- ncol(df_nowcast_plot) + 1 -
    rev(seq_along(quantiles_to_get))
  colnames(df_nowcast_plot)[cols_to_rename] <- paste(
    "quantile",
    quantiles_to_get,
    sep = "_"
  )
  df_nowcast_plot
}
