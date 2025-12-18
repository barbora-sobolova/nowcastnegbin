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
  length_of_train_data
) {
  data.frame(
    train_data_begin = start_date + (seq_len(timesteps_to_fit) - 1) * 7
  ) |>
    mutate(
      nowcast_date = .data$train_data_begin + (length_of_train_data - 1) * 7
    )
}

#' Load the data in the reporting triangle format
#'
#' @description This function extracts the data to fit the nowcasting model to
#' a single rolling window. The data are returned in the form of a matrix.
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
#' The data are returned in the form of a matrix.
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
#'
#' @return a matrix with `max_lag` columns containing the partial counts of the
#' reporting table The bottom-right part, which is usually unobserved,
#' still contains the partial count values, which will be hidden later.
#'
#' @export
filter_train_period <- function(full_data, start_date, end_date, max_lag) {
  full_data |> dplyr::filter(
    # Filter only the desired time period including the last date
    date >= start_date & date <= end_date
  ) |>
    dplyr::select(paste0("value_", 1:max_lag - 1, "w")) |>
    # Convert to matrix for simpler calculations
    as.matrix()
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
#' @return a matrix with the same number of columns, with the bottom-left part
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
#' }
#'
#' @export
get_stan_data <- function(train_data) {
  # Replace the known counts by NAs to create the reporting triangle
  obs_mat_truncated <- mock_unobserved(train_data)

  # Flatten the observation matrix by row
  obs_flat <- obs_mat_truncated |> t() |> c()
  list(
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
    d = ncol(obs_mat_truncated),
    # Observations in the flat format with the unobserved entries skipped. Must
    # be defined like this, otherwise the look-up indices in the STAN algorithm
    # won't work.
    obs = obs_flat[!is.na(obs_flat)]
  )
}
