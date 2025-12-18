create_totals_data_frame <- function(
  train_data,
  start_date
) {
  # Data frame to plot the observations - final counts and counts available at
  # the time of creating the nowcast
  data.frame(
    counts = c(
      rowSums(train_data),
      rowSums(mock_unobserved(train_data), na.rm = TRUE)
    ),
    date = rep(start_date + (seq_len(nrow(train_data)) - 1) * 7, 2),
    data = rep(c("Final", "Preliminary"), each = nrow(train_data))
  )
}

summarize_nowcast <- function(
  df_nowcast,
  df_total,
  date_of_the_nowcast
) {
  # Recover the beginning of the estimation window from the total
  # counts
  start_date <- min(df_total$date)
  # Recover the length of the estimation window from the total counts.
  # We have to add 1 to the difference of the two dates, since
  # `date_of_the_nowcast` is included as the last point of the estimation
  # window.
  length_of_train_data <- (date_of_the_nowcast - start_date + 1) / 7

  # Reformat the `week` column of the data frame with the nowcasts, so that it's
  # aligned with the actual date
  df_nowcast <- df_nowcast |>
    mutate(
      date = start_date + (week - 1) * 7,
      nowcast_date = date_of_the_nowcast
    ) |>
    # Remove the, now redundant, week column
    dplyr::select(-"week")
  # Quantiles of the MCMC sample to calculate - the median and quantiles for
  # constructing the 95% and 50% prediction interval
  quantiles_to_get <- c(50, 2.5, 25, 75, 97.5)
  # Summarize the sample of the nowcasts
  df_nowcast_plot <- df_total |>
    # The `df_total` data frame contains the final and the preliminary counts.
    # We need the final observed counts to evaluate the nowcasts via CRPS.
    filter(data == "Final") |>
    # Join the two data frames to put the predicted and the true values together
    inner_join(df_nowcast, by = "date") |>
    group_by(date, nowcast_date, Distribution) |>
    summarize(
      # Calculate the CRPS from the MCMC sample. We pass counts[1] as the true
      # observed value, since this is the same value for each date.
      CRPS = scoringutils::crps_sample(
        observed = counts[1],
        predicted = .value
      ),
      # Keep the true value
      true_val = counts[1],
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
    unnest(quantiles) |>
    ungroup()
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
