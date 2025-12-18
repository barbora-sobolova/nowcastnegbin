get_time_horizons <- function(
  start_date,
  timesteps_to_fit,
  length_of_train_data
) {
  data.frame(
    train_data_begin = start_date + (seq_len(timesteps_to_fit) - 1) * 7
  ) |>
    mutate(
      nowcast_date = train_data_begin + (length_of_train_data - 1) * 7
    )
}

load_preprocessed_data <- function(path, start_date, num_of_weeks) {
  # Set the end date. It will be excluded from the dataset
  analysis_end_date <- start_date + num_of_weeks * 7
  # Load the full dataset
  readr::read_csv(
    path,
    show_col_types = FALSE
  ) |>
    filter(
      # No stratification, we work with the aggregate numbers only
      age_group == "00+",
      # Filter only the desired time period
      date >= start_date & date < analysis_end_date
    )
}

filter_train_period <- function(full_data, start_date, end_date, max_lag) {
  full_data |> dplyr::filter(
    # Filter only the desired time period including the last date
    date >= start_date & date <= end_date
  ) |>
    dplyr::select(paste0("value_", 1:max_lag - 1, "w")) |>
    # Convert to matrix for simpler calculations
    as.matrix()
}

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
