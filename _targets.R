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
model_colors <- c(
  "Poisson" = "#CC79A7",
  "NegBinX" = "#D55E00",
  "NegBin2D" = "#009E73",
  "NegBin1D" = "#56B4E9",
  "NegBin2M" = "#004282",
  "NegBin1M" = "#F0E442"
)

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
  # A data frame encoding the observation model
  tar_target(obs_model, {
    data.frame(model_name = get_model_names(), model_number = 0:5)
  }),
  # Fitting of all models using dynamic branching over 6 observational models
  # defined in `obs_model` and rolling windows defined in `time_horizons`.
  tar_target(fitted_mcmc, {
    fit_stan_model(
      compiled_model$sample,
      stan_data = stan_data,
      model_obs = obs_model$model_number,
      stan_settings = stan_settings
    )
  },
  pattern = cross(obs_model, map(time_horizons, stan_data)),
  iteration = "list"
  ),
  tar_target(df_summarized_nowcast_mcmc, {
    summarize_nowcast(
      fitted_mcmc$nowcast,
      df_total,
      time_horizons$nowcast_date
    )
  },
  pattern = map(cross(obs_model, map(time_horizons, df_total)), fitted_mcmc)
  ),
  # Select the names of models we want to fit with the GLM method to branch over
  # it.
  tar_target(obs_model_glm, c("Poisson", "NegBinX", "NegBin2D", "NegBin1D")),
  # Fit the gamlss models
  tar_target(fitted_glm, {
    fit_glm_model(stan_data = stan_data, model_name = obs_model_glm)
  },
  pattern = cross(stan_data, obs_model_glm),
  iteration = "list"
  ),
  tar_target(df_summarized_nowcast_glm, {
    summarize_nowcast(
      fitted_glm$nowcast,
      df_total,
      time_horizons$nowcast_date
    )
  },
  pattern = map(fitted_glm, cross(map(time_horizons, df_total), obs_model_glm))
  ),
  # Collect the diagnostic summaries for the MCMC models
  tar_target(diagnostic_summaries, {
    fitted_mcmc$diagnostics |>
      mutate(
        date_of_the_nowcast = time_horizons$nowcast_date
      )
  },
  pattern = map(cross(obs_model, time_horizons), fitted_mcmc)
  ),
  # Plot the diagnostics of the MCMC procedure
  tar_target(plot_diagnostics, {
    plot_mcmc_diagnostics(diagnostic_summaries, model_colors)
  }),
  # Plot the nowcasts from the STAN model for each estimation window
  tar_target(nowcast_plot_mcmc, {
    plot_nowcast(
      df_summarized_nowcast_mcmc,
      df_total,
      model_codes = setNames(obs_model$model_name, obs_model$model_number),
      model_colors = model_colors,
      date_of_the_nowcast = time_horizons$nowcast_date,
      fitting_method = "mcmc"
    )
  },
  pattern = map(df_total, time_horizons),
  iteration = "list"),
  # Plot the nowcasts from the GLM model for each estimation window
  tar_target(nowcast_plot_glm, {
    plot_nowcast(
      df_summarized_nowcast_glm,
      df_total,
      # Select only the codes and colors of the first 4 models (that is
      # excluding NegBin2M and NegBin1M)
      model_codes = setNames(
        obs_model_glm,
        obs_model$model_number[obs_model$model_name %in% obs_model_glm]
      ),
      model_colors = model_colors[obs_model_glm],
      date_of_the_nowcast = time_horizons$nowcast_date,
      fitting_method = "glm"
    )
  },
  pattern = map(df_total, time_horizons),
  iteration = "list"),
  # Plot the overall coverage of the models
  tar_target(coverage_plot_mcmc, {
    plot_coverage(
      df_summarized_nowcast_mcmc,
      model_codes = setNames(obs_model$model_name, obs_model$model_number),
      model_colors = model_colors,
      fitting_method = "mcmc"
    )
  }),
  tar_target(coverage_plot_glm, {
    plot_coverage(
      df_summarized_nowcast_glm,
      model_codes = setNames(
        obs_model_glm,
        obs_model$model_number[obs_model$model_name %in% obs_model_glm]
      ),
      model_colors = model_colors[obs_model_glm],
      fitting_method = "glm"
    )
  }),
  # Plot the distribution of the CRPS
  tar_target(crps_plot_mcmc, {
    plot_crps(
      df_summarized_nowcast_mcmc,
      model_colors = model_colors,
      fitting_method = "mcmc"
    )
  }),
  tar_target(crps_plot_glm, {
    plot_crps(
      df_summarized_nowcast_glm,
      model_colors = model_colors[obs_model_glm],
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
