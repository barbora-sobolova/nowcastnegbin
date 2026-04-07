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
  skip_dates = NULL
) {
  data.frame(
    train_data_begin = start_date + (seq_len(timesteps_to_fit) - 1) * 7
  ) |>
    mutate(
      nowcast_date = .data$train_data_begin + (length_of_train_data - 1) * 7
    ) |>
    filter(!(.data$nowcast_date %in% skip_dates))
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
#' @param prior_delay_param a vector of positive real values, parameters of the
#' prior Dirichlet distribution of the reporting delay. The higher the sum of
#' its elements is, the more informative the prior distribution of the reporting
#' delay becomes.
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
#' @export
get_stan_data <- function(
  train_data,
  prior_delay_param = NULL,
  skip_rows = NULL
) {
  # Grab the maximum lag
  max_lag <- ncol(train_data)
  #
  if (length(prior_delay_param) != max_lag) {
    stop("The vector of prior parameters of the reporting delay must have the same length as there are columns in the reporting triangle.")  # nolint
  }
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
    # Parameters of the prior Dirichlet distribution of the reporting delay
    prior_delay_param = prior_delay_param,
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
