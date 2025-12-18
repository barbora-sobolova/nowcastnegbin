library(targets)
library(tarchetypes)
library(qs2)

# Set targets options
tar_option_set(
  packages = c(
    "scoringutils",
    "ggplot2",
    "ggpubr",
    "purrr",
    "here",
    "tidyr",
    "dplyr",
    "qs2",
    "readr"
  ),
  format = "qs", # Use qs format (qs2 is used via repository option)
  memory = "transient", # Free memory after each target completes
  garbage_collection = TRUE, # Run garbage collection
  repository = "local", # Use qs2 backend for storage
  error = "continue" # Continue pipeline when targets fail
)
tar_source(files = "R")

# set the ggplot theme
ggplot2::theme_set(ggplot2::theme_bw())

# Set the global objects =======================================================
model_colors <- c(
  "NegBin2D" = "#009E73",
  "NegBin1D" = "#56B4E9",
  "NegBin2M" = "#004282",
  "NegBin1M" = "#F0E442",
  "Poisson" = "#CC79A7",
  "NegBinX" = "#D55E00"
)

max_lag <- 5

# Where is the beginning of the data used for the case study
analysis_start_date <- as.Date("2024-07-07")
# How many weeks we want to include as "training" data.
# This includes the last `max_lag - 1` weeks for which we calculate the nowcast.
length_of_train_data <- 28
# For how many dates we want to do the fitting. For each time step, we shift the
# window of the train data to include a new week of observations mimicking a
# real-time analysis.
timesteps_to_fit <- 44

# A data frame encoding the observation model
obs_model <- data.frame(
  model_name = c(
    "Poisson",
    "NegBinX",
    "NegBin2D",
    "NegBin1D",
    "NegBin2M",
    "NegBin1M"
  ),
  model_number = 0:5
)

# Define the pipeline ==========================================================
list(
  # Compile the STAN model
  tar_target(compiled_model, {
    cmdstanr::cmdstan_model(here::here("inst", "stan", "nowcast.stan"))
  }),
  # STAN settings
  tar_target(stan_settings, {
    list(
      chains = 4,
      parallel_chains = 1,
      iter_warmup = 1000,
      iter_sampling = 1000,
      show_messages = FALSE,
      show_exceptions = FALSE,
      refresh = 0,
      seed = 12345
    )
  }),
  # Data frame storing the beginning and end points of the training data to
  # keep track of the rolling windows
  tar_target(time_horizons, {
    get_time_horizons(
      analysis_start_date,
      timesteps_to_fit,
      length_of_train_data
    )
  }),
  # Load the preprocessed data with no stratification, restricted to the time
  # period of interest
  tar_target(full_data, {
    load_preprocessed_data(
      here::here(
        "data",
        "SARI",
        "reporting_triangle-icosari-sari-preprocessed.csv"
      ),
      start_date = analysis_start_date,
      num_of_weeks = timesteps_to_fit + length_of_train_data - 1
    )
  }),
  # Create a matrix containing the training data for each date. This one has all
  # observations there, so it's not in the triangular form yet.
  tar_target(
    train_data,
    filter_train_period(
      full_data,
      start_date = time_horizons$train_data_begin,
      end_date = time_horizons$nowcast_date,
      max_lag = max_lag
    ),
    pattern = map(time_horizons),
    iteration = "list"
  ),
  # Create the list of data and parameters to pass to the STAN model
  tar_target(
    stan_data,
    get_stan_data(train_data),
    pattern = map(time_horizons, train_data),
    iteration = "list"
  ),
  # Calculate the reporting table rowsums and partial rowsums for each date.
  # The total sum (final counts) is used for plotting and evaluating the
  # prediction. The partial sums are used only for plotting.
  tar_target(df_total, {
    create_totals_data_frame(
      train_data,
      time_horizons$train_data_begin
    )
  },
  pattern = map(time_horizons, train_data),
  iteration = "list"),
  # Static branching over 6 observational models defined in `obs_model`.
  # By using static branching, we can easily remove some of the models (possibly
  # NegBin1M and NegBin2M) if required.
  tar_map(
    unlist = TRUE,
    values = obs_model,
    names = model_name,
    # Fit each observational model to each rolling window
    tar_target(fitted, {
      fit_stan_model(
        compiled_model$sample,
        stan_data = stan_data,
        model_obs = model_number,
        stan_settings = stan_settings
      )
    },
    pattern = map(time_horizons, stan_data),
    iteration = "list"),
    # Summarize the STAN draws of the nowcasts
    tar_target(df_summarized_nowcast, {
      summarize_nowcast(fitted$nowcast, df_total, time_horizons$nowcast_date)
    },
    pattern = map(time_horizons, df_total, fitted),
    iteration = "list"
    )
  ),
  # Combine the summarized nowcasts into a single table
  tar_target(df_summarized_nowcast, {
    dplyr::bind_rows(
        df_summarized_nowcast_Poisson,
        df_summarized_nowcast_NegBinX,
        df_summarized_nowcast_NegBin2D,
        df_summarized_nowcast_NegBin1D,
        df_summarized_nowcast_NegBin2M,
        df_summarized_nowcast_NegBin1M
    ) |>
      mutate(
        # Calculate the nowcasting horizon and save it as a factor for easier
        # plotting
        delay = factor(as.numeric(date - nowcast_date) / 7),
        # Replace the model number by the text label of the model
        Distribution = factor(Distribution, labels = obs_model$model_name)
      )
  },
  pattern = map(
    df_summarized_nowcast_Poisson,
    df_summarized_nowcast_NegBinX,
    df_summarized_nowcast_NegBin2D,
    df_summarized_nowcast_NegBin1D,
    df_summarized_nowcast_NegBin2M,
    df_summarized_nowcast_NegBin1M
  )
  ),
  # Plot the nowcasts for each estimation window
  tar_target(nowcast_plot, {
    plot_nowcast(
      df_summarized_nowcast,
      df_total,
      model_codes = setNames(obs_model$model_name, obs_model$model_number),
      model_colors = model_colors,
      nowcast_date = time_horizons$nowcast_date
    )
  },
  pattern = map(df_total, df_summarized_nowcast, time_horizons),
  iteration = "list"),
  # Plot the overall coverage of the models
  tar_target(coverage_plot, {
    plot_coverage(
      df_summarized_nowcast,
      model_codes = setNames(obs_model$model_name, obs_model$model_number),
      model_colors = model_colors
    )
  }),
  # Plot the distribution of the CRPS
  tar_target(crps_plot, {
    plot_crps(
      df_summarized_nowcast,
      model_colors = model_colors
    )
  }),
  # Plot the whole incidence trajectory highlighting the first and the last
  # estimation windows
  tar_target(whole_trajectory_plot, {
    plot_trajectory(
      full_data,
      analysis_start_date,
      length_of_train_data,
      max_lag
    )}
  )
)
