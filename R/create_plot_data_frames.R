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
#' @param fitting_method a method used for fitting the nowcasting model, either
#' "mcmc", or "glm"
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
  df_total,
  fitting_method = c("mcmc", "glm")
) {
  fitting_method <- match.arg(fitting_method)

  # Recover the beginning of the estimation window from the total
  # counts
  start_date <- min(df_total$date)

  # If we fit the GLM model, we don't do the sensitivity analysis, as it
  # concerns only some priors. To allow for the data frame grouping used below,
  # we add the "" string as the sensitivity analysis scenario, which indicates
  # the main analysis when we use the MCMC method.
  if (!("sensitivity_sc" %in% colnames(df_nowcast))) {
    df_nowcast <- df_nowcast |> mutate(sensitivity_sc = "")
  }

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
    dplyr::group_by(
      .data$date,
      .data$nowcast_date,
      .data$Distribution,
      .data$sensitivity_sc
    ) |>
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
      Distribution = factor(.data$Distribution, levels = get_model_names()),
      method = fitting_method
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

#' Combine data frames of results from different methods
#'
#' @description This function combines data frames containing results from the
#' MCMC and GLM method while also filtering only results corresponding to
#' selected dates and observation model.
#'
#' @param df_mcmc a data frame with results obtained via the MCMC method
#' @param df_glm a data frame with results obtained via the GLM method
#' @param dates_to_show a selection of 4-6 consecutive dates for which we want
#' to show the estimates.
#' @param model_to_show a string indicating an observation model, from which we
#' want to show the estimates
#'
#' @return a filtered data frame combining the GLM and MCMC results
#'
#' @import dplyr
#'
#' @export
filter_and_combine_methods <- function(
  df_mcmc,
  df_glm,
  dates_to_show,
  model_to_show
) {
  df_combined <- bind_rows(
    mutate(
      dplyr::filter(
        df_mcmc,
        # Keep only rows with the desired date and observation model
        .data$nowcast_date %in% dates_to_show &
          .data$Distribution %in% model_to_show
      ),
      method = "mcmc"
    ),
    mutate(
      dplyr::filter(
        df_glm,
        # Keep only rows with the desired date and observation model
        .data$nowcast_date %in% dates_to_show &
          .data$Distribution %in% model_to_show
      ),
      method = "glm"
    )
  )
  df_combined
}

#' Extract estimates of the mean process and nowcasts for selected time windows
#'
#' @description This function filters data frame rows corresponding to selected
#' dates and observation model.
#'
#' @param df_nowcast_mcmc a data frame with columns `week`, `.value`,
#' `Distribution`, `nowcast_date` and `sensitivity_sc` that contains the
#' distribution of the nowcast obtained from the MCMC method in a sample format
#' @param df_nowcast_glm same as \code{df_nowcast_mcmc} but the results are
#' obtained from the GLM method
#' @param df_total a data frame containing columns `date`, `counts` and `data`.
#' The last column `data` is an indicator, whether the values in the `counts`
#' column are the final sums of the counts, or the preliminary data version.
#' Needed to plot the observations alongside the nowcasts.
#' @param df_lambda_mcmc a data frame with columns `week`, `.value`,
#' `Distribution` and `nowcast_date`, that contains the distribution of the mean
#' process obtained from the MCMC method in a sample format. If NULL, the
#' distribution of the mean process is ignored.
#' @param df_lambda_glm same as \code{df_lambda_mcmc} but the results are
#' obtained from the GLM method
#' @param dates_to_show a selection of 4-6 consecutive dates for which we want
#' to show the estimates.
#' @param model_to_show a character vector indicating an observation models,
#' from which we want to show the estimates
#'
#' @return a list containing 3 data frames
#' \describe{
#'   \item{\code{lambda}}{with columns `week` (the time as an integer starting
#'   from 1), `.value` (sampled value of lambda), `Distribution` (observation
#'   model name), `nowcast_date` (date when the nowcast was calculated),
#'   `sensitivity_sc` (name of the sensitivity analysis scenario, applicable
#'   only for the MCMC method) and `method` (either "mcmc" or "glm"),}
#'   \item{\code{nowcast}}{with columns `date` (date of the nowcasting target),
#'   `delay` (reporting delay in weeks), `nowcast_date`, `Distribution`,
#'   `quantile_2.5`,`quantile_25`, `quantile_50`, `quantile_75`, `quantile_97.5`
#'   and `method`,}
#'   \item{\code{total}}{a data frame containing columns `date`, `counts` and
#'   `data` (indicator, whether the counts are preliminary or final).}
#' }
#' The returned data frames are empty, if the input data contained no entries
#' from the selected dates.
#'
#' @import dplyr
#'
#' @export
filter_nowcast_example_dates <- function(
  df_nowcast_mcmc,
  df_nowcast_glm,
  df_total,
  df_lambda_mcmc = NULL,
  df_lambda_glm = NULL,
  dates_to_show = c("2019-03-24", "2019-03-31", "2019-04-07", "2019-04-14"),
  model_to_show = c(
    "Poisson",
    "NegBinX",
    "NegBin2D",
    "NegBin1D",
    "NegBin2M",
    "NegBin1M"
  )
) {
  # We filter out all rolling windows we don't want to show. Since the function
  # is called for each dynamic branch, in most cases, the filtered data frame
  # will have 0 rows.
  df_nowcast <- filter_and_combine_methods(
    df_nowcast_mcmc,
    df_nowcast_glm,
    dates_to_show,
    model_to_show
  )
  if (!is.null(df_lambda_mcmc) && !is.null(df_lambda_glm)) {
    df_lambda <- filter_and_combine_methods(
      df_lambda_mcmc,
      df_lambda_glm,
      dates_to_show,
      model_to_show
    )
  } else {
    df_lambda <- NULL
  }

  ret_list <- list(lambda = df_lambda, nowcast = df_nowcast)
  # We also select the data frame with the total counts only for the selected
  # dates. Otherwise we return NULL.
  if (max(df_total$date) %in% dates_to_show) {
    # Add the date of the nowcast for easier faceting during plotting
    df_total <- df_total |> mutate(nowcast_date = max(df_total$date))
    ret_list <- c(ret_list, list(total = df_total))
  } else {
    ret_list <- c(ret_list, list(total = NULL))
  }
  ret_list
}
