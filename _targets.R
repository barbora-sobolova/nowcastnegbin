library(targets)
library(tarchetypes)
library(qs2)

# Set targets options
tar_option_set(
  packages = c(
    "cmdstanr",
    "dplyr",
    "gamlss.dist",
    "gamlss2",
    "ggplot2",
    "ggpubr",
    "here",
    "purrr",
    "qs2",
    "readr",
    "scoringutils",
    "tidyr"
  ),
  format = "qs", # Use qs format
  memory = "transient", # Free memory after each target completes
  garbage_collection = TRUE, # Run garbage collection
  repository = "local", # Use qs2 backend for storage
  error = "continue" # Continue pipeline when targets fail
)
tar_source(files = "R")

# set the ggplot theme
ggplot2::theme_set(ggplot2::theme_bw())

# Set the global objects =======================================================

max_lag <- 5

# Where is the beginning of the data used for the case study
analysis_start_date <- as.Date("2024-07-28")
# How many weeks we want to include as "training" data.
# This includes the last `max_lag - 1` weeks for which we calculate the nowcast.
length_of_train_data <- 20
# For how many dates we want to do the fitting. For each time step, we shift the
# window of the train data to include a new week of observations mimicking a
# real-time analysis.
timesteps_to_fit <- 55
# What dates shall be skipped due to the Christmas break. These dates indicate
# two things:
#  1. No nowcast will be produced on these days
#  2. The diagonal of the reporting triangle corresponding to these dates and
#     most of the one directly following will be dropped from the likelihood.
skip_dates <- as.Date(c("2024-12-22", "2024-12-29", "2025-12-21", "2025-12-28"))
# Parameters of the prior reporting delay distribution.
prior_delay_param <- c(5, 1.5, 0.5, 0.25, 0.25)

# Define the pipeline ==========================================================
list(
  # Select the names of models we want to fit with the GLM method to branch over
  # it.
  tar_target(obs_model_glm, c("Poisson", "NegBinX", "NegBin2D", "NegBin1D")),
  # Compile the STAN model
  tar_target(compiled_model, {
    cmdstanr::cmdstan_model(here::here("inst", "stan", "nowcast.stan"))
  }),
  # STAN settings
  tar_target(stan_settings, {
    list(
      parallel_chains = 4,
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
      length_of_train_data,
      skip_dates = skip_dates
    )
  }),
  # Create grouped data frames to group targets by date. As a result, the models
  # will be stored and subsequently loaded in bundles  of 4 (for GLM), or 6
  # (for MCMC)
  tar_group_by(branches_mcmc, {
    time_horizons |>
      mutate(
        nested_col = list(
          tibble(model_name = get_model_names(), model_code = 0:5)
        )
      ) |>
      unnest(nested_col)
  },
  train_data_begin,
  nowcast_date
  ),
  tar_group_by(branches_glm, {
    time_horizons |>
      mutate(model_obs = list(obs_model_glm)) |>
      unnest(model_obs)
  },
  train_data_begin,
  nowcast_date
  ),
  # Load the preprocessed data with no stratification, restricted to the time
  # period of interest
  tar_target(full_data, {
    load_preprocessed_data(
      here::here(
        "inst",
        "extdata",
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
      max_lag = max_lag,
      skip_dates = skip_dates
    ),
    pattern = map(time_horizons),
    iteration = "list"
  ),
  # Create the list of data and parameters to pass to the STAN model
  tar_target(
    stan_data,
    get_stan_data(
      train_data$train_data,
      prior_delay_param,
      train_data$skip_rows
    ),
    pattern = map(time_horizons, train_data),
    iteration = "list"
  ),
  # Calculate the reporting table rowsums and partial rowsums for each date.
  # The total sum (final counts) is used for plotting and evaluating the
  # prediction. The partial sums are used only for plotting.
  tar_target(df_total, {
    create_totals_data_frame(
      train_data$train_data,
      time_horizons$train_data_begin
    )
  },
  pattern = map(time_horizons, train_data),
  iteration = "list"),
  # Select the names of models we want to fit with the MCMC method. This is all
  # 6 observation models
  tar_target(obs_model, get_model_names()),
  # Fitting of all models using dynamic branching over the rolling windows
  # which are defined as groups of `branches_mcmc`
  tar_target(fitted_mcmc, {
    fit_all_stan_models(
      compiled_model$sample,
      stan_data = stan_data,
      model_obs = branches_mcmc$model_code,
      date_of_the_nowcast = branches_mcmc$nowcast_date,
      stan_settings = stan_settings
    )
  },
  pattern = map(branches_mcmc, stan_data),
  iteration = "list"
  ),
  # Calculate the quantiles and CRPS of the nowcasts obtained by the MCMC method
  tar_target(summarized_nowcast_mcmc, {
    summarize_nowcast(fitted_mcmc$nowcast, df_total = df_total)
  },
  pattern = map(fitted_mcmc, df_total),
  iteration = "list"
  ),
  # Create plots for each rolling window. For the MCMC procedure we plot:
  # - the nowcast,
  # - posterior density of the delay probability,
  # - posterior density of the dispersion parameter on a scale, where 0 means
  #   the Poisson model and higher values indicate more dispersion,
  # - the scatter plot of the dispersion parameter against the standard
  #   deviation of the random walk.
  tar_target(rolling_plots_mcmc, {
    plot_per_window(
      summarized_nowcast_mcmc,
      fitted_mcmc$delay_prob,
      fitted_mcmc$nb_size,
      fitted_mcmc$rw_sd,
      df_total,
      obs_model,
      time_horizons$nowcast_date,
      fitting_method = "mcmc"
    )
  },
  pattern = map(
    fitted_mcmc,
    time_horizons,
    df_total,
    summarized_nowcast_mcmc
  ),
  iteration = "list"),
  # Create plots of aggregated results from the MCMC method. We plot:
  # - the coverage of nowcasts,
  # - the crps decomposition.
  tar_target(aggreg_plots_mcmc, {
    plot_aggregated(
      bind_rows(summarized_nowcast_mcmc),
      obs_model,
      fitting_method = "mcmc"
    )
  }),
  # Plot the diagnostic summaries for the MCMC models
  tar_target(plot_diagnostics, {
    plot_mcmc_diagnostics(
      bind_rows(map(fitted_mcmc, "diagnostics")),
      obs_model
    )
  }),
  # Fit the gamlss models
  tar_target(fitted_glm, {
    fit_all_glm_models(
      stan_data = stan_data,
      date_of_the_nowcast = branches_glm$nowcast_date,
      model_name = branches_glm$model_obs
      )
  },
  pattern = map(branches_glm, stan_data),
  iteration = "list"
  ),
  # Calculate the quantiles and CRPS of the nowcasts obtained by the GLM method
  tar_target(summarized_nowcast_glm, {
    summarize_nowcast(fitted_glm$nowcast, df_total = df_total)
  },
  pattern = map(fitted_glm, df_total),
  iteration = "list"
  ),
  # Create plots for each rolling window. For the GLM procedure we plot:
  # - the nowcast,
  # - posterior density of the delay probability,
  # - posterior density of the dispersion parameter on a scale, where 0 means
  #   the Poisson model and higher values indicate more dispersion.
  tar_target(rolling_plots_glm, {
    plot_per_window(
      summarized_nowcast_glm,
      fitted_glm$delay_prob,
      fitted_glm$nb_size,
      NULL,  # We don't have the random walk parameters
      df_total,
      obs_model_glm,
      time_horizons$nowcast_date,
      fitting_method = "glm"
    )
  },
  pattern = map(fitted_glm, time_horizons, df_total, summarized_nowcast_glm),
  iteration = "list"),
  # Create plots of aggregated results from the GLM method. We plot:
  # - the coverage of nowcasts,
  # - the crps decomposition.
  tar_target(aggreg_plots_glm, {
    plot_aggregated(
      bind_rows(summarized_nowcast_glm),
      obs_model_glm,
      fitting_method = "glm"
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
    )
  })
)
