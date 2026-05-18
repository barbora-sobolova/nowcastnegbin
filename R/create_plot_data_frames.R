#' Create a data frame from the reporting table
#'
#' @description This function creates a ready-to-be-plotted data frame of the
#' total incidence from the reporting table
#'
#' @param train_data the reporting table in a matrix format
#' @param start_date date (in a date format) when the incidence begins
#'
#' @return a data frame with columns `counts` and `date`. The data frame is in a
#' long format. Its first half contains the complete total incidence. The second
#' half contains the incomplete sums. The data versions are indicated by a
#' string in the `data` column.
#'
#' @export
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

#' Summarize the sample from the predictive distribution of the nowcasts
#'
#' @description This function creates a data frame containing summarized
#' nowcast draws from the MCMC, or GLM model-fitting procedure. The summaries
#' calculated are the mean, the median, the CRPS and selected quantiles. We
#' calculate the 2.5%, 25%, 75% and the 97.5% quantiles to prepare for the
#' plotting of the 50% and 95% prediction intervals.
#'
#' @param df_nowcast a data frame with columns `week`, `.value`, `Distribution`
#' @param df_total a data frame containing columns `date`, `counts` and `data`.
#' The last column `data` is an indicator, whether the values in the `counts`
#' column are the final sums of the counts, or the preliminary data version.
#'
#' @return a data frame with columns
#' \describe{
#'   \item{\code{date}}{date of the nowcasting target,}
#'   \item{\code{nowcast_date}}{date when the nowcast was calculated,}
#'   \item{\code{Distribution}}{factor, label of the observation model,}
#'   \item{\code{crps}}{the CRPS calculated from the sample of the nowcasts,}
#'   \item{\code{dispersion}}{the dispersion component of the CRPS,}
#'   \item{\code{overprediction}}{the overprediction component of the CRPS,}
#'   \item{\code{underprediction}}{the underprediction component of the CRPS,}
#'   \item{\code{true_val}}{the final value  of the incidence to compare the
#'   nowcast to,}
#'   \item{\code{mean}}{mean of the sampled nowcasts,}
#'   \item{\code{quantile_*}}{quantiles of the sampled nowcasts ranging from 5 %
#'    to 95 % wit a 5 % step and also the 2.5 % and 97.5 % quantiles}
#' }
#'
#' @import dplyr
#' @importFrom scoringutils crps_sample
#' @importFrom stats quantile
#' @importFrom tibble as_tibble
#'
#' @export
summarize_nowcast <- function(
  df_nowcast,
  df_total
) {
  # Recover the beginning of the estimation window from the total
  # counts
  start_date <- min(df_total$date)

  # Reformat the `week` column of the data frame with the nowcasts, so that it's
  # aligned with the actual date
  df_nowcast <- df_nowcast |>
    mutate(date = start_date + (.data$week - 1) * 7) |>
    # Remove the, now redundant, week column
    dplyr::select(-"week")
  # Quantiles of the sampled nowcasts to calculate - 5 % to 95 % quantiles with
  # a 5 % step in between and the 2.5 % and 97.5 % quantiles for constructing
  # the 95% and 50% prediction interval
  quantiles_to_get <- c(2.5, seq(5, 95, by = 5), 97.5)
  # Summarize the sample of the nowcasts
  df_nowcast_plot <- df_total |>
    # The `df_total` data frame contains the final and the preliminary counts.
    # We need the final observed counts to evaluate the nowcasts via CRPS.
    dplyr::filter(data == "Final") |>
    # Join the two data frames to put the predicted and the true values together
    dplyr::inner_join(df_nowcast, by = "date") |>
    dplyr::group_by(.data$date, .data$nowcast_date, .data$Distribution) |>
    dplyr::summarize(
      # Calculate the CRPS from the sample. We pass counts[1] as the true
      # observed value, since this is the same value for each date.
      CRPS = as_tibble(scoringutils::crps_sample(
        observed = .data$counts[1],
        predicted = .data$.value,
        separate_results = TRUE
      )),
      # Keep the true value
      true_val = .data$counts[1],
      # Calculate the mean and the quantiles
      mean = mean(.data$.value),
      quantiles = list(
        tibble::as_tibble(
          as.list(
            quantile(.data$.value, probs = quantiles_to_get / 100)
          )
        )
      ),
      .groups = "drop"
    ) |>
    tidyr::unnest(c("quantiles", "CRPS")) |>
    mutate(
      # Calculate the nowcasting horizon and save it as a factor for easier
      # plotting
      delay = factor(as.numeric(date - .data$nowcast_date) / 7),
      # Convert to factor to make sure, the plotting order of the models is
      # consistent.
      Distribution = factor(.data$Distribution, levels = get_model_names())
    )
  # Rename the quantile columns to have nicer names
  cols_to_rename <- colnames(df_nowcast_plot) %in% paste0(quantiles_to_get, "%")
  colnames(df_nowcast_plot)[cols_to_rename] <- paste(
    "quantile",
    quantiles_to_get,
    sep = "_"
  )
  df_nowcast_plot
}
