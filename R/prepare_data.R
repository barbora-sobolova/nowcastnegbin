load_triangle <- function(path, start_date, num_of_weeks, max_lag) {
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
      date >= analysis_start_date & date < analysis_end_date
    )
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

get_stan_data <- function(obs_mat) {
  # Assume the latest counts to be unobserved
  obs_mat_truncated <- obs_mat |> mock_unobserved()
  # Flatten the observation matrix by row
  obs_flat <- obs_mat_truncated |> t() |> c()
  list(
    # Total number of days/weeks/time units. This is the number of rows of the
    # reporting triangle
    n = nrow(obs_mat),
    # Total number of non-empty cells of the reporting triangle, can be
    # calculated as the total number of non-NA elements of the observation
    # matrix.
    m = sum(!is.na(obs_mat_truncated)),
    # Vector indicating how many observations we have available in each row of
    # the reporting triangle.
    p = apply(obs_mat_truncated, 1, function(x) sum(!is.na(x))),
    # Maximum lag with delay zero counting as the first lag. This is the number
    # of columns of the reporting triangle.
    d = ncol(obs_mat),
    # Observations in the flat format with the unobserved entries skipped. Must
    # be defined like this, otherwise the look-up indices in the STAN algorithm
    # won't work.
    obs = obs_flat[!is.na(obs_flat)]
  )
}
