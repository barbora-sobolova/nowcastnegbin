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

# Where the beginning of the data used for the case study is
analysis_start_date <- as.Date("2024-06-23")
# How many weeks we want to include as "training" data.
# This includes the last `max_lag - 1` weeks for which we calculate the nowcast.
length_of_train_data <- 20
# We run an auxiliary case study one year before the actual one to determine
# the prior distributions.
aux_analysis_start_date <- analysis_start_date - (52 + length_of_train_data) * 7
# For how many dates we want to do the fitting. For each time step, we shift the
# window of the train data to include a new week of observations mimicking a
# real-time analysis.
timesteps_to_fit <- 75
# What dates shall be skipped due to the Christmas break. These dates indicate
# two things:
#  1. No nowcast will be produced on these days
#  2. The diagonal of the reporting triangle corresponding to these dates and
#     most of the one directly following will be dropped from the likelihood.
skip_dates <- as.Date(c("2024-12-22", "2024-12-29", "2025-12-21", "2025-12-28"))

# Where the beginning of the data used for the simulation study is. For the
# simulation study, we take the total SARI counts from several years back,
# smooth them to obtain a mean process and then simulate the counts according to
# one of our models.
sim_start_date <- as.Date("2014-10-05")
# We smooth the data using moving average of degree 3.
ma_degree <- 3
# For how many dates we want to do the fitting.
sim_timesteps_to_fit <- 512
# Delay probabilities used in the simulation.
sim_delay_prob <- c(0.5, 0.3, 0.2, 0.1)
# Size of the negative binomial distribution used in the simulation.
sim_nb_size <- 0.5
# Selected models for the simulation study
sim_obs_model <- data.frame(
  model_name = c("NegBinX", "NegBin2D", "NegBin1D"),
  model_number = c(1, 2, 3)
)
# Seed used to simulate the counts
sim_seed <- 2436

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


  # Load the case study data ---------------------------------------------------

  # Load the preprocessed data with no stratification. This data spans the
  # whole period from the beginning of the auxiliary case study to the end of
  # the actual case study.
  tar_target(full_data, {
    load_preprocessed_data(
      here::here(
        "inst",
        "extdata",
        "reporting_triangle-icosari-sari-preprocessed.csv"
      ),
      start_date = aux_analysis_start_date,
      # How many weeks of data (rows of the reporting triangle) we want to load.
      # This is the length of the auxiliary case study (52 weeks + train data)
      # and the length of the actual case study (train data + desired number of
      # rolling windows). The -1 part is included to get the exact number of
      # rolling windows, since we count the "zeroth" window as a first one.
      num_of_weeks = 52 + 2 * length_of_train_data + timesteps_to_fit - 1
    )
  }),

  # Fit the GLM models to the previous year to obtain the priors ---------------

  # Data frame storing the beginning and end points of the training data for the
  # auxiliary analysis to keep track of the rolling windows
  tar_target(time_horizons_prev_year, {
    get_time_horizons(
      aux_analysis_start_date,
      52,
      length_of_train_data,
      skip_dates = as.Date("2023-12-24")
    )
  }),
  # Create a matrix containing the training data for each date in the auxiliary
  # analysis. This matrix contains all observations. To obtain the triangular
  # form, latest observations will be masked by the `get_stan_data()` function
  # further downstream.
  tar_target(
    train_data_prev_year,
    filter_train_period(
      full_data,
      start_date = time_horizons_prev_year$train_data_begin,
      end_date = time_horizons_prev_year$nowcast_date,
      max_lag = max_lag,
      skip_dates = as.Date("2023-12-24")
    ),
    pattern = map(time_horizons_prev_year),
    iteration = "list"
  ),
  # Create the list of data and parameters to pass to the STAN model for the
  # auxiliary analysis
  tar_target(
    stan_data_prev_year,
    get_stan_data(
      train_data_prev_year$train_data,
      prior_delay_param,
      train_data_prev_year$skip_rows
    ),
    pattern = map(time_horizons_prev_year, train_data_prev_year),
    iteration = "list"
  ),
  # Fit the GLM models to the auxiliary data.
  tar_target(fitted_glm_prev_year, {
    fit_glm_model(
      stan_data = stan_data_prev_year,
      date_of_the_nowcast = time_horizons_prev_year$nowcast_date,
      model_name = obs_model_glm
    )
  },
  pattern = cross(
    map(stan_data_prev_year, time_horizons_prev_year),
    obs_model_glm
  ),
  iteration = "list"
  ),
  # Extract the dispersion parameter estimates from the fit
  tar_target(
    glm_log_disp_par_prev_year,
    fitted_glm_prev_year$log_disp_coeff,
    pattern = map(fitted_glm_prev_year)
  ),
  # Calculate the prior parameters based on the estimates of the dispersion
  # parameter
  tar_target(
    disp_par_prior,
    calc_disp_par_prior(glm_log_disp_par_prev_year)
  ),
  # Calculate the parameters of the Dirichlet prior from the auxiliary data
  # only, without looking at the GLM estimates.
  tar_target(prior_delay_param, {
    calc_delay_prob_prior(
      full_data,
      aux_analysis_start_date,
      aux_analysis_start_date + (length_of_train_data + 52) * 7
      )
  }),

  # Case study -----------------------------------------------------------------

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
  # will be stored and subsequently loaded in bundles of 4 (for GLM), or 6
  # (for MCMC)
  tar_group_by(branches_mcmc, {
    time_horizons |>
      mutate(
        nested_col = list(
          tibble(model_name = get_model_names(), model_code = 0:5)
        )
      ) |>
      unnest(nested_col) |>
      inner_join(disp_par_prior, relationship = "many-to-one")
  },
  train_data_begin,
  nowcast_date
  ),
  tar_group_by(branches_glm, {
    time_horizons |>
      mutate(model_name = list(obs_model_glm)) |>
      unnest(model_name)
  },
  train_data_begin,
  nowcast_date
  ),
  # Create a matrix containing the training data for each date. This matrix
  # contains all observations. To obtain the triangular form, latest
  # observations will be masked by the `get_stan_data()` function further
  # downstream.
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
      mean_log = branches_mcmc$mean_log,
      sd_log = branches_mcmc$sd_log,
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
      fitting_method = "mcmc",
      prior_delay_param,
      disp_par_prior
    )
  },
  pattern = map(
    fitted_mcmc,
    time_horizons,
    df_total,
    summarized_nowcast_mcmc,
    branches_mcmc
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
      model_name = branches_glm$model_name
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
      max_lag,
      aux_analysis_start_date
    )
  })
)
