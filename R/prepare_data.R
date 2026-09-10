#' Create the data frame defining the estimation windows
#'
#' @description This function creates a data frame containing the beginning and
#' the end point of estimation windows. The time points are then used to extract
#' the training data.
#'
#' @param start_date a date in the date format, the starting date for a single
#' estimation window.
#' @param timesteps_to_fit the number of rolling estimation windows.
#' @param length_of_train_data the size of one estimation window.
#' @param skip_dates a date vector containing the dates, where we don't wish to
#' calculate the nowcast
#' @param time_step an integer indicating the time resolution of the data. The
#' default \code{time_step = 7} imply weekly data resolution
#'
#' @return a data frame containing columns `train_data_begin`, with the
#' estimation window starts in the date format, and `nowcast_date`, where the
#' estimation window ends. Both dates will be included in the training data
#'
#' @importFrom dplyr mutate
#'
#' @export
get_time_horizons <- function(
  start_date,
  timesteps_to_fit,
  length_of_train_data,
  time_step = 7,
  skip_dates = NULL
) {
  data.frame(
    train_data_begin = start_date + (seq_len(timesteps_to_fit) - 1) * time_step
  ) |>
    mutate(
      nowcast_date = .data$train_data_begin +
        (length_of_train_data - 1) * time_step
    ) |>
    filter(!(.data$nowcast_date %in% skip_dates))
}

#' Download and save the ILI data
#'
#' @description This function downloads the ILI data from the Fluview project
#' using the DELPHI epidemiological data API and converts it to the same format
#' as the SARI data.
#'
#' @param start_date a date in the date format, the starting date of the first
#' rolling window of the main analysis. This date is used to pivot around when
#' determining the data points used for the main and the auxiliary analysis.
#' @param timesteps_to_fit the number of rolling estimation windows in the main
#' analysis.
#' @param aux_timesteps_to_fit the number of rolling estimation windows in the
#' auxiliary analysis.
#' @param length_of_train_data the size of one estimation window, same for the
#' main and auxiliary analysis
#' @param max_lag integer, the number of columns of the reporting triangle
#'
#' @return a data frame containing columns `date` and columns `value_0w`,
#' `value_1w`, etc. until `max_lag - 1`. The same data frame is saved onto the
#' disc in the CSV format.
#'
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_wider
#' @importFrom epidatr pub_fluview epirange
#' @importFrom readr write_csv
#'
#' @export
process_ili_data <- function(
  start_date,
  timesteps_to_fit,
  aux_timesteps_to_fit,
  length_of_train_data,
  max_lag = 6
) {
  # Last point of the last rolling window
  analysis_end_date <- start_date +
    (timesteps_to_fit + length_of_train_data - 2) * 7
  # Find the start of the data to be downloaded. Before the first window begins,
  # we need enough data to determine priors.
  download_start <- start_date -
    (aux_timesteps_to_fit + length_of_train_data - 1) * 7
  # Download snapshots from some weeks after the analysis ends to consolidate
  # the final data.
  download_end <- analysis_end_date + (max_lag - 1) * 7

  # Download the data using the Delphi API. The data are in a snapshot format.
  ili_all <- epidatr::pub_fluview(
    regions = "nat",
    # Weeks we want to download
    epiweeks = epidatr::epirange(
      # Convert the date to the Year-week format readable by `epirange()`
      format(download_start, "%Y%W"),
      format(analysis_end_date, "%Y%W")
    ),
    # Which data versions to download.
    issues = epidatr::epirange(
      format(download_start, "%Y%W"),
      format(download_end, "%Y%W")
    )
  ) |>
    # The dates are shifted by one day in comparison to the McGough analysis.
    # We shift them back to be compatible with the paper and with the SARI data,
    # which assume Monday to be the reference day.
    mutate(nowcast_date = .data$issue + 1, date = .data$epiweek + 1) |>
    select(c("nowcast_date", "date", "num_ili"))
  # Process the snapshots to the reporting triangle format
  ili_processed <- ili_all |>
    mutate(
      delay = as.numeric((.data$nowcast_date - .data$date) / 7),
      dummy_col_name = "value"
    ) |>
    # Truncate the data at the maximum lag. The truncated data becomes our
    # underlying truth.
    filter(.data$delay < max_lag) |>
    # Arrange by date to calculate the increments from the cumulative data.
    arrange(.data$date, .data$nowcast_date) |>
    group_by(date) |>
    # Calculate the weekly increments. They occasionally become negative, which
    # we replace by 0.
    mutate(counts = diff(c(0, .data$num_ili))) |>
    ungroup() |>
    select(c("date", "counts", "delay", "dummy_col_name")) |>
    # Pivot the data to obtain a reporting triangle same as in the SARI case
    # study
    tidyr::pivot_wider(
      id_cols = "date",
      names_from = c("dummy_col_name", "delay"),
      values_from = "counts"
    ) |>
    # Add a dummy age group column to match the SARI data
    mutate(age_group = "00+") |>
    # Adjust the column names with the reported counts by adding the "w" letter
    rename_with(~paste0(.x, "w"), starts_with("value"))
  # Some negative values occur in the reporting table. When that happens, we
  # subtract the corresponding number of counts from the last positive column to
  # the left.
  subtract <- rep(0, nrow(ili_processed))
  for (k in rev(seq_len(max_lag) - 1)) {
    col_name <- paste0("value_", k, "w")
    ili_processed[, col_name] <- ili_processed[, col_name] + subtract
    subtract <- pmin(0, ili_processed[, col_name][[1]])
    ili_processed[, col_name] <- pmax(0, ili_processed[, col_name][[1]])
  }

  # Save as a CSV file
  readr::write_csv(ili_processed, "inst/extdata/fluview_ili.csv")
  # Return the NULL value. To retrieve the data, they must be loaded from the
  # disc
  NULL
}

#' Load the data in the reporting triangle format
#'
#' @description This function loads the full preprocessed dataset, that is
#' already in the form of a reporting triangle, and restricts it to the period
#' used in the case study
#'
#' @param path the path of the data file
#' @param start_date a date in the date format, the starting date for the whole
#' case study.
#' @param num_of_weeks number of weeks starting from the `start_date`
#' (included), we want to load.
#'
#' @return a data frame containing columns `date` and columns
#' `value_0w`, `value_1w`, etc. until `max_lag - 1`.
#'
#' @export
load_preprocessed_data <- function(path, start_date, num_of_weeks) {
  # Set the end date. It will be excluded from the dataset
  analysis_end_date <- start_date + num_of_weeks * 7
  # Load the full dataset
  readr::read_csv(
    path,
    show_col_types = FALSE
  ) |>
    dplyr::filter(
      # No stratification, we work with the aggregate numbers only
      .data$age_group == "00+",
      # Filter only the desired time period
      date >= start_date & date < analysis_end_date
    )
}

#' Simulate the reporting table of a nowcasting problem based on existing data
#'
#' @description Generate the full reporting table, i.e. with no right censoring.
#' The mean process is based on a real data set, where we smooth the observed
#' counts using a moving average process.
#'
#' @param df_series Data frame with a column `value`, where observed counts are
#' stored. These are smoothed by a moving average and taken as the mean process
#' to be passed to \code{generate_reports}.
#' @param ma_degree Integer, the order of the moving average process.
#' @param max_lag Integer, the maximum reporting delay, the width of the table.
#' @param probs Numeric, a numeric vector specifying the delay distribution.
#' @param nb_size Numeric, a positive real value specifying the size of the
#' negbin distribution. The lower, the more dispersed
#' @param model name of the observation model
#' @param seed An integer for seeding the simulation
#'
#' @return The reporting table in a data frame format
#'
#' @importFrom dplyr mutate
simulate_full_data <- function(
  df_series,
  ma_degree,
  max_lag,
  probs,
  nb_size = NULL,
  model = c(
    "Poisson",
    "NegBinX",
    "NegBin2D",
    "NegBin1D",
    "NegBin2M",
    "NegBin1M"
  ),
  seed = 123456
) {
  model <- match.arg(
    model,
    c("Poisson", "NegBinX", "NegBin2D", "NegBin1D", "NegBin2M", "NegBin1M")
  )

  lgt <- nrow(df_series) - ma_degree + 1

  # Use tail to skip the initial `ma_degree` - 1 observations that are set to
  # NA.
  mean_proc <- tail(
    rowSums(
      sapply(seq_len(ma_degree) - 1, dplyr::lag, x = df_series$value)
    ) / ma_degree,
    lgt
  )
  # Generate the reporting table in the format of a matrix
  reporting_table <- generate_reports(
    lgt,
    max_lag,
    probs,
    nb_size = nb_size,
    model = model,
    fixed_lambda = mean_proc,
    seed = seed
  )$reports
  # Coerce the matrix to a data frame and name the columns appropriately
  colnames(reporting_table) <- paste0("value_", seq_len(max_lag) - 1, "w")
  reporting_table <- reporting_table |>
    as.data.frame() |>
    dplyr::mutate(
      date = tail(df_series$date, lgt),
      mean_proc = mean_proc,
      Distribution = model
    )
  reporting_table
}

#' Extract data from one rolling window
#'
#' @description This function extracts the data of a single rolling window.
#' The data are returned in the form of a list containing a matrix and indices
#' of skipped rows.
#'
#' @param full_data a data frame containing columns `date` and columns
#' `value_0w`, `value_1w`, etc. until `max_lag - 1`.
#' @param start_date a date (indeed in the date format), where the training data
#' start. The starting point will be included.
#' @param end_date a date (indeed in the date format), where the training data
#' ends. The ending point will be included.
#' @param max_lag maximum reporting delay represented by the number of columns
#' of the reporting table. In this way, the 0-th lag counts as the first, 1-st
#' lag as the second and so on.
#' @param skip_dates a vector of dates, on which we don't calculate the nowcast.
#' This is typically one or two weeks during the Christmas period.
#'
#' @return a list of two elements:
#' \describe{
#'   \item{\code{train_data}}{matrix with `max_lag` columns containing the
#'   partial counts of the reporting table. The bottom-right part, which is
#'   usually unobserved, still contains the partial count values, which will be
#'   hidden later,}
#'   \item{\code{skip_rows}}{indices of rows, corresponding to the dates, on
#'   which we don't calculate the nowcast due to Christmas. These are used by
#'   \code{get_stan_data} to calculate the specific reporting pattern of the
#'   Christmas period.}
#' }
#'
#' @export
filter_train_period <- function(
  full_data,
  start_date,
  end_date,
  max_lag,
  skip_dates = NULL
) {
  full_data_filtered <- full_data |> dplyr::filter(
    # Filter only the desired time period including the last date
    date >= start_date & date <= end_date
  )
  # Select only the columns containing the values
  train_data <- full_data_filtered |>
    dplyr::select(paste0("value_", seq_len(max_lag) - 1, "w")) |>
    # Convert to matrix for simpler calculations
    as.matrix()
  # Find indices of rows, where we don't want to do the nowcasting due to the
  # Christmas break. These will be used in `get_stan_data()` to calculate the
  # indices of cells, we want to skip. The Christmas break has a very specific
  # pattern. Change accordingly here and in `get_stan_data()`, if the pattern
  # becomes more general.
  skip_rows <- which(full_data_filtered$date %in% skip_dates)
  list(train_data = train_data, skip_rows = skip_rows)
}

#' Mock the unobserved counts
#'
#' @description This function replaces the bottom-right part of the reporting
#' table to mimic the real-time observation process and create the triangular
#' shape.
#'
#' @param obs_counts a matrix of the partial counts with the complete delayed
#' counts in its columns.
#'
#' @return a matrix with the same number of columns, with the bottom-right part
#' of the reporting triangle filled by `NA` values.
#'
#' @export
mock_unobserved <- function(obs_counts) {
  max_lag <- ncol(obs_counts)
  # Coerce to a matrix in case a data frame was supplied
  obs_mat <- as.matrix(obs_counts)
  # Create an indexing matrix to flag the entries considered as unobserved.
  index_mat <- upper.tri(
    matrix(nrow = nrow(obs_mat), ncol = max_lag),
    diag = FALSE
  )[rev(seq_len(nrow(obs_mat))), ]  # reverse the rows to get the right triangle
  # Replace the last data points by NA to imitate the unknown entries in the
  # reporting triangle.
  obs_mat[index_mat] <- NA
  obs_mat
}

#' Get the list of STAN data
#'
#' @description This function prepares the list of data and parameters to pass
#' to the STAN model.
#'
#' @param train_data a matrix of the partial counts with the complete delayed
#' counts in its columns.
#' @param skip_rows indices of rows, corresponding to the  dates, on which we
#' don't calculate the nowcast due to Christmas. These are used to calculate the
#' specific reporting pattern of the Christmas period.
#'
#' @return a list of parameters for the STAN model:
#' \describe{
#'   \item{\code{n}}{number of the timepoints - rows of the reporting triangle,}
#'   \item{\code{m}}{number of the partial counts in total - the number of cells
#'   with an observed entry}
#'   \item{\code{p}}{a vector with the number of cells with an observed entry
#'   per row}
#'   \item{\code{d}}{the maximum lag - the number of columns of the reporting
#'   triangle}
#'   \item{\code{obs}}{the flattened counts from the reporting triangle.
#'   The flattening is done by row with unobserved entries (`NA`s) skipped.}
#'   \item{\code{idx_include}}{the indices of observations in the flat
#'   observation matrix, we want to include in the likelihood.}
#'   \item{\code{n_idx_include}}{number of observations in the flat observation
#'   matrix, we want to include in the likelihood. Equals to \code{m} if all
#'   shall be included.}
#' }
#'
#' @importFrom tidyr replace_na
#'
#' @export
get_stan_data <- function(
  train_data,
  skip_rows = NULL
) {
  # Grab the maximum lag
  max_lag <- ncol(train_data)
  # Replace the known counts by NAs to create the reporting triangle
  obs_mat_truncated <- mock_unobserved(train_data)
  # Flatten the observation matrix by row
  obs_flat <- obs_mat_truncated |> t() |> c()
  stan_data <- list(
    # Total number of days/weeks/time units. This is the number of rows of the
    # reporting triangle
    n = nrow(obs_mat_truncated),
    # Total number of non-empty cells of the reporting triangle, can be
    # calculated as the total number of non-NA elements of the observation
    # matrix.
    m = sum(!is.na(obs_mat_truncated)),
    # Vector indicating how many observations we have available in each row of
    # the reporting triangle.
    p = apply(obs_mat_truncated, 1, function(x) sum(!is.na(x))),
    # Maximum lag with delay zero counting as the first lag. This is the number
    # of columns of the reporting triangle.
    d = max_lag,
    # Observations in the flat format with the unobserved entries skipped. Must
    # be defined like this, otherwise the look-up indices in the STAN algorithm
    # won't work.
    obs = obs_flat[!is.na(obs_flat)]
  )

  # Find indices of cells, where we drop observations due to the Christmas
  # break. The Christmas break has a very specific pattern. Change accordingly
  # here and in `filter_train_period()`, if the pattern becomes more general.
  # Currently, we remove the diagonal observations corresponding to the days,
  # where nothing was reported. In addition we remove the cells to the right of
  # the last dropped diagonal, because here the reports compensate for the
  # lack of reports in the previous week(s).
  # The `skip_rows` vector is empty, if the skipped dates are not present in the
  # filtered time period.
  if (!purrr::is_empty(skip_rows)) {
    # Calculate the pairs of indices of the cells we skip. They are constructed
    # in a way, that we start in the top-right corner and follow the diagonal to
    # the bottom-left corner. Then we move one diagonal down and go again from
    # its top-right corner to the bottom-left one.
    skip_row_ind <- c(
      # The diagonals we skip completely
      sequence(
        rep(max_lag, length(skip_rows)),
        from = skip_rows - max_lag + 1
      ),
      # The diagonal after the Christmas break, where we keep the bottom-left
      # cell
      seq(from = max(skip_rows) - max_lag + 2, to = max(skip_rows))
    )
    skip_col_ind <- c(
      # The diagonals we skip completely
      rev(sequence(rep(max_lag, length(skip_rows)))),
      # The diagonal after the Christmas break, where we keep the bottom-left
      # cell
      max_lag:2
    )
    # If the skipped dates are at the beginning of the training data, the
    # calculated row indices might be negative.
    skip_ind_pair <- cbind(skip_row_ind, skip_col_ind)[skip_row_ind > 0, ]

    # Create a copy of the observation matrix and fill the skipped cells by an
    # arbitrary value, which can't possibly appear in the data. In our case,
    # anything negative will do.
    obs_mat_truncated_copy <- obs_mat_truncated
    obs_mat_truncated_copy[skip_ind_pair] <- -100
    # Flatten the observation matrix by row
    obs_flat_copy <- obs_mat_truncated_copy |> t() |> c()
    # Remove the unobserved part from the flat observation matrix and locate the
    # negative values. This way we can correctly identify, which observations
    # from the flat matrix we want to include and which shall be skipped.
    stan_data$idx_include <- which(obs_flat_copy[!is.na(obs_flat_copy)] != -100)
    stan_data$n_idx_include <- length(stan_data$idx_include)
  } else {
    # If nothing is skipped, we simply include the whole flat observation
    # matrix.
    stan_data$idx_include <- seq_len(stan_data$m)
    stan_data$n_idx_include <- stan_data$m
  }
  stan_data
}

#' Find the prior for the dispersion parameter
#'
#' @description This function takes the estimates of the dispersion parameter
#' from the auxiliary analysis and calculates the parameters of the log-normal
#' prior distribution that is used for the dispersion parameter in the main
#' analysis. The scale of the prior distribution is calculated based on the
#' standard deviation and point estimates from the auxiliary analysis.
#'
#' @param log_disp_par a data frame with columns `log_disp_hat` (the point
#' estimate from the GLM method), `log_disp_se` (the standard error from the GLM
#' method), `Distribution` (the name of the observation model)
#' @param disp_par_prior_scale_factor a numeric value controlling the spread of
#' the prior distribution. The higher the value, the flatter the prior is.
#' The default value of 3 corresponds to the main scenario.
#' @return a data frame with columns `mean_log` (location parameter of the
#' log-normal distribution),`sd_log` (scale parameter of the log-normal
#' distribution) and `model_name`. The data frame has 6 rows, one for each
#' observation model. For Poisson model, we use placeholder values -1, for the
#' NegBin2M and NegBin1M models, we reuse the prior of other models.
calc_disp_par_prior <- function(log_disp_par, disp_par_prior_scale_factor = 3) {
  prior_pars_from_glm <- log_disp_par |>
    group_by(.data$Distribution) |>
    summarize(
      mean_log = mean(.data$log_disp_hat),
      # Loosely inspired by Rubin's rules. The scale factor is there to make
      # prior distribution even wider and can be subjected to a sensitivity
      # analysis.
      sd_log = sqrt(
        mean(.data$log_disp_se^2) + (1 + 1 / n()) * var(.data$log_disp_hat)
      ) *
        disp_par_prior_scale_factor,
      disp_par_factor = disp_par_prior_scale_factor
    ) |>
    # Set placeholder values for the Poisson model
    tidyr::replace_na(list(mean_log = -1, sd_log = -1))

  # For the NegBin2M we will use the same prior as for NegBinX, as these have
  # identical marginals. For NegBin1M, we will take the parameters of NegBin1D.
  prior_pars_assigned <- bind_rows(
    mutate(
      filter(prior_pars_from_glm, .data$Distribution == "NegBinX"),
      Distribution = "NegBin2M"
    ),
    mutate(
      filter(prior_pars_from_glm, .data$Distribution == "NegBin1D"),
      Distribution = "NegBin1M"
    )
  )
  bind_rows(
    prior_pars_from_glm,
    prior_pars_assigned
  ) |> rename(
    "model_name" = "Distribution"
  )
}

#' Find the prior for the delay probability
#'
#' @description This function takes the data from the auxiliary analysis and
#' calculates the parameters of the Dirichlet prior distribution that is used
#' for the delay probability in the main analysis.
#'
#' @param full_data a data frame with columns `date` and columns
#' `value_0w`, `value_1w`, etc. until `max_lag - 1`. The value of `max_lag` is
#' not checked here and is only derived from the columns of \code{full_data}.
#' @param start_date a date in the date format, the beginning of the auxiliary
#' analysis
#' @param end_date a date in the date format, the endpoint of the auxiliary
#' analysis
#' @param prior_scale_factor a numeric value controlling the spread of
#' the prior distribution. The higher the value, the flatter the prior is.
#' The default value of 4 corresponds to the main scenario.
#' @return a vector of length `max_lag` with the parameters of the Dirichlet
#' distribution.
calc_delay_prob_prior <- function(
  full_data,
  start_date,
  end_date,
  prior_scale_factor = 4
) {
  prob_vec <- full_data |>
    filter(date >= start_date & date < end_date) |>
    select(starts_with("value_")) |>
    as.matrix() |>
    apply(1, function(x) x / sum(x)) |>
    t() |>
    # The scale factor of controls the "flatness" of the prior distribution. May
    # be varied as a part of a sensitivity analysis.
    apply(2, mean)
  names(prob_vec) <- paste("delay", seq_along(prob_vec) - 1, sep = "_")
  ret <- as.data.frame(prior_scale_factor %*% t(prob_vec)) |>
    mutate(delay_prob_factor = prior_scale_factor)
  ret
}

#' Create a data frame with info about the dynamic branching structure
#'
#' @description This function creates a grouped data frame of dates and models
#' to be used to create a dynamic branching structure. Per one nowcasting date,
#' a bundle of observation models should be fitted.
#'
#' @param time_horizons a data frame with columns `train_data_begin` and
#' `nowcast_date`
#' @param disp_par_prior a data frame with columns `mean_log`, `sd_log` and
#' `model_name`, indicating the prior parameters for the negative binomial
#' dispersion parameter. Relevant only for \code{fitting_method = "mcmc"},
#' otherwise NULL
#' @param delay_prob_prior a data frame with columns containing the Dirichlet
#' prior parameters for the reporting delay probability vector, called
#' `delay_0`, `delay_1` until the maximum delay, and column `delay_prob_factor`
#' which indicates the scale of the prior as a sum of the Dirichlet parameters.
#' Relevant only for \code{fitting_method = "mcmc"}, otherwise NULL
#' @param obs_model_glm a vector of model names to be fitted using the GLM
#' method. Only relevant, when \code{fitting_method = "glm"}, otherwise NULL
#' @param fitting_method a string indicating the model fitting procedure. For
#' GLM, we don't have the NegBin2M and NegBin1M models
#' @param sensitivity_scenarios a data frame with 3 columns defining the
#' sensitivity analysis scenarios for the MCMC procedure. To fit the MCMC
#' method, we need to specify the scale parameter of the prior distribution for
#' the dispersion parameter and the delay probability vector. These are
#' contained in columns `disp_par_factor` and `delay_prob_factor`. The last
#' column `scenario_name` indicates the name of the sensitivity analysis
#' scenario. The main scenario is indicated by "". This data frame is ignored
#' for \code{fitting_method = "glm"}.
#'
#' @return a data frame with columns `train_data_begin`, `nowcast_date` (The
#' data frame will be grouped by these 2 columns in the pipeline.),
#' `model_name` and `scenario_name`. For \code{fitting_method = "mcmc"} we also
#' have columns `disp_par_factor`, `mean_log`, `delay_prob_factor`, `sd_log`,
#' `delay_0`, `delay_1`, etc. until the maximum lag
#'
#' @import dplyr
#' @importFrom tibble tibble
#' @importFrom tidyr unnest
#'
#' @export
group_branches <- function(
  time_horizons,
  disp_par_prior = NULL,
  delay_prob_prior = NULL,
  obs_model_glm = NULL,
  fitting_method = c("mcmc", "glm"),
  sensitivity_scenarios = data.frame(
    delay_prob_factor = 4,
    disp_par_factor = 3,
    scenario_name = ""
  )
) {
  fitting_method <- match.arg(fitting_method)
  if (fitting_method == "mcmc") {
    if (is.null(delay_prob_prior)) {
      stop("'delay_prob_prior' must be provided for fitting_method = 'mcmc'")
    }
    if (is.null(disp_par_prior)) {
      stop("'disp_par_prior' must be provided for fitting_method = 'mcmc'")
    }
    prior_pars <- inner_join(
      disp_par_prior,
      sensitivity_scenarios,
      by = "disp_par_factor",
      relationship = "many-to-many"
    ) |>
      unique() |>
      inner_join(
        delay_prob_prior,
        by = "delay_prob_factor",
        relationship = "many-to-many"
      )
    ret <- time_horizons |>
      mutate(
        # For the MCMC method, we will fit all 6 models and we need to store
        # their numbers to call the fitting function.
        nested_col = list(
          tibble::tibble(model_name = get_model_names(), model_code = 0:5)
        )
      ) |>
      tidyr::unnest("nested_col") |>
      # Join with the data frame of prior parameters for the negative binomial
      # dispersion parameter
      inner_join(
        prior_pars,
        by = "model_name",
        relationship = "many-to-many"
      )
  } else {
    ret <- time_horizons |>
      mutate(model_name = list(obs_model_glm)) |>
      tidyr::unnest("model_name")
  }
  ret
}
