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
analysis_start_date <- as.Date("2023-12-24")
# How many weeks we want to include as "training" data.
# This includes the last `max_lag - 1` weeks for which we calculate the nowcast.
length_of_train_data <- 20
# We run an auxiliary case study of 32 time windows before the actual one to
# determine the prior distributions.
aux_timesteps_to_fit <- 32
aux_analysis_start_date <- analysis_start_date -
  (aux_timesteps_to_fit + length_of_train_data - 1) * 7
# For how many dates we want to do the fitting. For each time step, we shift the
# window of the train data to include a new week of observations mimicking a
# real-time analysis.
timesteps_to_fit <- 104
# What dates shall be skipped due to the Christmas break. These dates indicate
# two things:
#  1. No nowcast will be produced on these days
#  2. The diagonal of the reporting triangle corresponding to these dates and
#     most of the one directly following will be dropped from the likelihood.
skip_dates <- as.Date(
  c("2023-12-24", "2024-12-22", "2024-12-29", "2025-12-21", "2025-12-28")
)

# Where the beginning of the data used for the simulation study is. For the
# simulation study, we take the total SARI counts from several years back,
# smooth them to obtain a mean process and then simulate the counts according to
# one of our models.
sim_start_date <- as.Date("2015-10-11")
# Like in the case study, we run an auxiliary simulation study on the first
# "year" of the simulated data to determine the prior distributions.
# This date should fall to 19. October 2014, which is the first date, where
# simulated data are available.
aux_sim_start_date <- sim_start_date -
  (aux_timesteps_to_fit + length_of_train_data - 1) * 7
# We smooth the data using moving average of degree 3.
ma_degree <- 3
# For how many rolling windows we want to do the fitting.
sim_timesteps_to_fit <- 500
# Delay probabilities used in the simulation.
sim_delay_prob <- c(0.5, 0.3, 0.2, 0.1)
# Dispersion parameter of the negative binomial distribution used in the
# simulation. No single value can be used, as the dispersion parameters are on a
# different scale for each model.
sim_disp_par <- c("NegBinX" = 0.04, "NegBin2D" = 0.01, "NegBin1D" = 100)
# Selected models for the simulation study
sim_obs_model <- data.frame(
  model_obs = c("NegBinX", "NegBin2D", "NegBin1D"),
  model_number = c(1, 2, 3)
)
# Seed used to simulate the counts
sim_seed <- 2436

# Define the pipeline
list(
  # Names of GLM models we use in the auxiliary case/simulation study to
  # determine the priors. NegBin2D can't be fitted exactly using the GLM method,
  # but we consider it to be close enough to produce plausible estimates for
  # determining the prior distribution for the negative binomial dispersion
  # parameter in the MCMC method.
  tar_target(
    obs_model_glm_for_auxiliary,
    c("Poisson", "NegBinX", "NegBin2D", "NegBin1D")
  ),
  # Names of GLM models we use in the main case/simulation study. We don't fit
  # NegBin2D in the main part anymore
  tar_target(obs_model_glm, c("Poisson", "NegBinX", "NegBin1D")),
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

  # Simulation study ===========================================================

  # Generate data for the simulation study -------------------------------------

  # Load the data in the snapshot format. The historical time series in the
  # snapshot format is available all the way back to 2014. We smooth the curve
  # using a MA-process to obtain the mean process of the total counts.
  tar_target(sim_data_series, {
    load_preprocessed_data(
      here::here(
        "inst",
        "extdata",
        "latest_data-SARI-sari.csv"
      ),
      # We need to load MA-degree - 1 weeks of extra data in order to be able to
      # smooth the data using a MA-process.
      start_date = aux_sim_start_date - (ma_degree - 1) * 7,
      # How many weeks of data (rows of the reporting triangle) we want to load.
      # This is the length of the auxiliary simulation study (auxiliary windows
      # + train data) and the length of the actual simulation study
      # (train data + desired number of rolling windows). The -1 part is
      # included to get the exact number of rolling windows, since we count the
      # "zeroth" window as the first one. The length of training data is
      # identical here and in the case study.
      num_of_weeks = aux_timesteps_to_fit +
        2 * length_of_train_data + sim_timesteps_to_fit - 1
    )
  }),
  # Data frame storing the beginning and end points of the training data for the
  # auxiliary simulation study to keep track of the rolling windows
  tar_target(
    sim_time_horizons_prev_year,
    get_time_horizons(
      aux_sim_start_date,
      aux_timesteps_to_fit,
      length_of_train_data,
      skip_dates = NULL
    )
  ),
  # Define the rolling windows for the main part of the simulation study
  tar_target(
    sim_time_horizons,
    get_time_horizons(
      sim_start_date,
      sim_timesteps_to_fit,
      length_of_train_data,
      skip_dates = NULL
    )
  ),
  # We simulate from 3 observation models: NegBinX, NegBin2D and NegBin1D.
  # For each model we find priors, do the fitting using MCMC and GLM and plot
  # the results.
  tar_map(
    unlist = TRUE,
    values = sim_obs_model,
    names = model_obs,
    # Simulate the reporting triangle using the smoothed version of the
    # historical data as the mean process. It is identical for all 3 observation
    # models.
    tar_target(sim_full_data, {
      simulate_full_data(
        sim_data_series,
        ma_degree,
        max_lag = length(sim_delay_prob),
        probs = sim_delay_prob,
        nb_size = 1 / sim_disp_par[model_obs],
        model = model_obs,
        seed = sim_seed
      )
    }),
    # Extract the last part of the mean process in each rolling window that will
    # be plotted alongside the estimates
    tar_target(
      sim_mean_process_tail, {
        filter(
          sim_full_data,
          date > sim_time_horizons$nowcast_date -
            (length(sim_delay_prob) - 1) * 7 &
            date <= sim_time_horizons$nowcast_date
        ) |>
          select(c("date", "mean_proc"))
      },
      pattern = map(sim_time_horizons),
      iteration = "list"
    ),

    # Fit the GLM models to the beginning of the simulated data to obtain the
    # priors -------------------------------------------------------------------

    # Create a matrix containing the training data for each date in the
    # auxiliary simulation study.
    tar_target(
      sim_train_data_prev_year,
      filter_train_period(
        sim_full_data,
        start_date = sim_time_horizons_prev_year$train_data_begin,
        end_date = sim_time_horizons_prev_year$nowcast_date,
        max_lag = length(sim_delay_prob),
        skip_dates = NULL
      ),
      pattern = map(sim_time_horizons_prev_year),
      iteration = "list"
    ),
    # Create the list of data and parameters that we would pass to the STAN
    # model. In the auxiliary simulation study, the list will be passed to the
    # GLM model only.
    tar_target(
      sim_stan_data_prev_year,
      get_stan_data(sim_train_data_prev_year$train_data),
      pattern = map(sim_time_horizons_prev_year, sim_train_data_prev_year),
      iteration = "list"
    ),
    # Fit the GLM model to the first part of the simulated data to determine
    # the prior parameters of the negative binomial overdispersion parameter
    tar_target(
      sim_fitted_glm_prev_year,
      fit_glm_model(
        stan_data = sim_stan_data_prev_year,
        date_of_the_nowcast = sim_time_horizons_prev_year$nowcast_date,
        model_name = obs_model_glm_for_auxiliary
      ),
      pattern = cross(
        map(sim_stan_data_prev_year, sim_time_horizons_prev_year),
        obs_model_glm_for_auxiliary
      ),
      iteration = "list"
    ),
    # Extract the dispersion parameter estimates from the fit
    tar_target(
      sim_glm_log_disp_par_prev_year,
      sim_fitted_glm_prev_year$log_disp_coeff,
      pattern = map(sim_fitted_glm_prev_year)
    ),
    # Calculate the prior parameters based on the GLM estimates of the
    # dispersion parameter
    tar_target(
      sim_disp_par_prior,
      calc_disp_par_prior(sim_glm_log_disp_par_prev_year)
    ),
    # Calculate the parameters of the Dirichlet prior from the auxiliary
    # simulated data, without looking at the GLM estimates.
    tar_target(
      sim_prior_delay_param,
      calc_delay_prob_prior(
        sim_full_data,
        aux_sim_start_date,
        aux_sim_start_date + (length_of_train_data + aux_timesteps_to_fit) * 7
      )
    ),

    # Main part of the simulation study ----------------------------------------

    # Create grouped data frames to group targets by date. As a result, the
    # models will be stored and subsequently loaded in bundles of 4 (for GLM),
    # or 6 (for MCMC)
    tar_group_by(
      sim_branches_mcmc,
      group_branches(
        sim_time_horizons,
        sim_disp_par_prior,
        sim_prior_delay_param,
        fitting_method = "mcmc"
      ),
      train_data_begin,
      nowcast_date
    ),
    tar_group_by(
      sim_branches_glm,
      group_branches(
        sim_time_horizons,
        obs_model_glm = obs_model_glm,
        fitting_method = "glm"
      ),
      train_data_begin,
      nowcast_date
    ),
    # Create a matrix containing the training data for each date of the
    # simulation study.
    tar_target(
      sim_train_data,
      filter_train_period(
        sim_full_data,
        start_date = sim_time_horizons$train_data_begin,
        end_date = sim_time_horizons$nowcast_date,
        max_lag = length(sim_delay_prob),
        skip_dates = NULL
      ),
      pattern = map(sim_time_horizons),
      iteration = "list"
    ),
    # Prepare the STAN data for each date of the simulation study.
    tar_target(
      sim_stan_data,
      get_stan_data(sim_train_data$train_data),
      pattern = map(sim_time_horizons, sim_train_data),
      iteration = "list"
    ),
    # Calculate the reporting table rowsums and partial rowsums of the simulated
    # dataset for each date. The total sum (final counts) is used for plotting
    # and evaluating the prediction. The partial sums are used only for
    # plotting.
    tar_target(
      sim_df_total,
      create_totals_data_frame(
        sim_train_data$train_data,
        sim_time_horizons$train_data_begin
      ),
      pattern = map(sim_time_horizons, sim_train_data),
      iteration = "list"
    ),
    # Fit each observational model to each rolling window of the simulation
    # study using the MCMC method.
    tar_target(
      sim_fitted_mcmc,
      fit_all_stan_models(
        compiled_model$sample,
        stan_data = sim_stan_data,
        model_obs = sim_branches_mcmc$model_code,
        date_of_the_nowcast = sim_branches_mcmc$nowcast_date,
        prior_delay_param = select(sim_branches_mcmc, starts_with("delay_")),
        mean_log = sim_branches_mcmc$mean_log,
        sd_log = sim_branches_mcmc$sd_log,
        stan_settings = stan_settings,
        sensitivity_scenario_name = branches_mcmc$scenario_name
      ),
      pattern = map(sim_branches_mcmc, sim_stan_data),
      iteration = "list"
    ),
    # Calculate the quantiles and CRPS of the nowcasts in the simulation study
    # obtained by the MCMC method
    tar_target(
      sim_summarized_nowcast_mcmc,
      summarize_nowcast(sim_fitted_mcmc$nowcast, df_total = sim_df_total),
      pattern = map(sim_fitted_mcmc, sim_df_total),
      iteration = "list"
    ),
    # Create plots for a random sample of the rolling windows from the
    # simulation study. For the MCMC procedure we plot:
    # - the nowcast,
    # - posterior density of the delay probability,
    # - posterior density of the dispersion parameter on a scale, where 0 means
    #   the Poisson model and higher values indicate more dispersion,
    # - posterior density of the mean process for the weeks, where we perform
    #   nowcasting
    # - the scatter plot of the dispersion parameter against the standard
    #   deviation of the random walk.
    tar_target(
      sim_rolling_plots_mcmc,
      plot_per_window(
        sim_summarized_nowcast_mcmc,
        sim_fitted_mcmc$delay_prob,
        sim_fitted_mcmc$nb_size,
        sim_fitted_mcmc$lambda,
        sim_fitted_mcmc$rw_sd,
        sim_df_total,
        # What observation models we fitted
        obs_model,
        sim_time_horizons$nowcast_date,
        fitting_method = "mcmc",
        prob_prior_pars = sim_prior_delay_param,
        disp_prior_pars = sim_disp_par_prior,
        # From which observation model we simulated the data
        data_origin = model_obs,
        # True values of the model parameters used to generate the data
        prob_true_val = sim_delay_prob,
        disp_true_val = sim_disp_par[model_obs],
        lambda_true_val = sim_mean_process_tail$mean_proc
      ),
      pattern = sample(
        map(
          sim_fitted_mcmc,
          sim_time_horizons,
          sim_df_total,
          sim_summarized_nowcast_mcmc,
          sim_mean_process_tail
        ),
        n = 15
      ),
      iteration = "list"
    ),
    # Create plots of aggregated results from the simulation study for the MCMC
    # method. We plot:
    # - the coverage of nowcasts,
    # - the crps decomposition.
    tar_target(sim_aggreg_plots_mcmc, {
      plot_aggregated(
        bind_rows(sim_summarized_nowcast_mcmc),
        # What observation models we fitted
        obs_model,
        fitting_method = "mcmc",
        # From which observation model we simulated the data
        data_origin = model_obs
      )
    }),
    # Extract the diagnostic summaries for the MCMC models in the simulation
    # study. We do it per branch to avoid loading all fits at once when we want
    # to plot the diagnostics into a single plot.
    tar_target(
      sim_diagnostics,
      sim_fitted_mcmc$diagnostics,
      pattern = map(sim_fitted_mcmc)
    ),
    # Plot the diagnostic summaries for the MCMC models in the simulation study
    tar_target(
      sim_plot_diagnostics,
      plot_mcmc_diagnostics(sim_diagnostics, obs_model, data_origin = model_obs)
    ),
    # Plot the whole incidence trajectory highlighting the first and the last
    # estimation windows
    tar_target(
      sim_whole_trajectory_plot,
      plot_trajectory(
        sim_full_data,
        sim_start_date,
        length_of_train_data,
        length(sim_delay_prob),
        aux_sim_start_date,
        data_origin = model_obs
      )
    ),
    # Fit each observational model to each rolling window of the simulation
    # study using the GLM method.
    tar_target(
      sim_fitted_glm,
      fit_all_glm_models(
        stan_data = sim_stan_data,
        date_of_the_nowcast = sim_branches_glm$nowcast_date,
        model_name = sim_branches_glm$model_name
      ),
      pattern = map(sim_branches_glm, sim_stan_data),
      iteration = "list"
    ),
    # Calculate the quantiles and CRPS of the nowcasts in the simulation study
    # obtained by the GLM method
    tar_target(
      sim_summarized_nowcast_glm,
      summarize_nowcast(sim_fitted_glm$nowcast, df_total = sim_df_total),
      pattern = map(sim_fitted_glm, sim_df_total),
      iteration = "list"
    ),
    # Create plots for each rolling window of the simulation study. For the GLM
    # procedure we plot:
    # - density of estimates of the delay probability,
    # - density of estimates of the mean process for the weeks, where we
    #   perform nowcasting
    # - density of estimates of the dispersion parameter on a scale, where 0
    #   means the Poisson model and higher values indicate more dispersion.
    tar_target(
      sim_rolling_plots_glm,
      plot_per_window(
        sim_summarized_nowcast_glm,
        sim_fitted_glm$delay_prob,
        sim_fitted_glm$nb_size,
        sim_fitted_glm$lambda,
        NULL,  # We don't have the random walk parameters
        sim_df_total,
        obs_model_glm,
        sim_time_horizons$nowcast_date,
        fitting_method = "glm",
        # From which observation model we simulated the data
        data_origin = model_obs,
        # True values of the model parameters used to generate the data
        prob_true_val = sim_delay_prob,
        disp_true_val = sim_disp_par[model_obs],
        lambda_true_val = sim_mean_process_tail$mean_proc
      ),
      pattern = sample(
        map(
          sim_fitted_glm,
          sim_time_horizons,
          sim_df_total,
          sim_summarized_nowcast_glm,
          sim_mean_process_tail
        ),
        n = 15
      ),
      iteration = "list"
    ),
    # Create plots of aggregated results from the simulation study for the GLM
    # method. We plot:
    # - the coverage of nowcasts,
    # - the crps decomposition.
    tar_target(sim_aggreg_plots_glm, {
      plot_aggregated(
        bind_rows(sim_summarized_nowcast_glm),
        # What observation models we fitted
        obs_model_glm,
        fitting_method = "glm",
        # From which observation model we simulated the data
        data_origin = model_obs
      )
    })
  ),

  # Case study =================================================================

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
      # This is the length of the auxiliary case study (auxiliary windows +
      # train data) and the length of the actual case study (train data +
      # desired number of rolling windows). The -1 part is included to get the
      # exact number of rolling windows, since we count the "zeroth" window as
      # the first one.
      num_of_weeks = aux_timesteps_to_fit +
        2 * length_of_train_data + timesteps_to_fit - 1
    )
  }),

  # Fit the GLM models to the previous year to obtain the priors ---------------

  # Data frame storing the beginning and end points of the training data for the
  # auxiliary analysis to keep track of the rolling windows
  tar_target(time_horizons_prev_year, {
    get_time_horizons(
      aux_analysis_start_date,
      aux_timesteps_to_fit,
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
  # Create the list of data and parameters that we would pass to the STAN
  # model. In the auxiliary analysis, the list will be passed to the GLM model
  # only.
  tar_target(
    stan_data_prev_year,
    get_stan_data(
      train_data_prev_year$train_data,
      skip_rows = train_data_prev_year$skip_rows
    ),
    pattern = map(time_horizons_prev_year, train_data_prev_year),
    iteration = "list"
  ),
  # Fit the GLM models to the auxiliary data.
  tar_target(fitted_glm_prev_year, {
    fit_glm_model(
      stan_data = stan_data_prev_year,
      date_of_the_nowcast = time_horizons_prev_year$nowcast_date,
      model_name = obs_model_glm_for_auxiliary
    )
  },
  pattern = cross(
    map(stan_data_prev_year, time_horizons_prev_year),
    obs_model_glm_for_auxiliary
  ),
  iteration = "list"
  ),
  # Extract the dispersion parameter estimates from the fit
  tar_target(
    glm_log_disp_par_prev_year,
    fitted_glm_prev_year$log_disp_coeff,
    pattern = map(fitted_glm_prev_year)
  ),
  # Determine the sensitivity analysis scenario. We examine the robustness of
  # the MCMC model with respect to the flatness of the delay probability prior
  # and the dispersion parameter prior. All other parameters are either fixed as
  # some standard values, or estimated from the auxiliary data. We change only
  # one parameter at a time to avoid fitting the whole grid of parameter
  # combinations.
  tar_target(sensitivity_scenarios, {
    data.frame(
      scenario_name = c(
        "",  # Main analysis has no name
        "disp_low",  # Low multiplicative factor ~ more informative prior
        "disp_high",  # High multiplicative factor ~ less informative prior
        "prob_low",  # Low multiplicative factor ~ more informative prior
        "prob_high"  # High multiplicative factor ~ less informative prior
      ),
      delay_prob_factor = c(4, 4, 4, 1, 16),
      disp_par_factor = c(3, 1, 9, 3, 3)
    )
  }),
  # Calculate the prior parameters based on the estimates of the dispersion
  # parameter
  tar_target(
    disp_par_prior,
    calc_disp_par_prior(
      glm_log_disp_par_prev_year,
      sensitivity_scenarios$disp_par_factor
    ),
    pattern = map(sensitivity_scenarios)
  ),
  # Calculate the parameters of the Dirichlet prior from the auxiliary data
  # only, without looking at the GLM estimates.
  tar_target(prior_delay_param, {
    calc_delay_prob_prior(
      full_data,
      aux_analysis_start_date,
      aux_analysis_start_date +
        (length_of_train_data + aux_timesteps_to_fit) * 7,
      unique(sensitivity_scenarios$delay_prob_factor)
    )
  }),

  # Main part of the case study ------------------------------------------------

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
  tar_group_by(
    branches_mcmc,
    group_branches(
      time_horizons,
      disp_par_prior,
      prior_delay_param,
      fitting_method = "mcmc",
      sensitivity_scenarios = sensitivity_scenarios
      ),
    train_data_begin,
    nowcast_date
  ),
  tar_group_by(
    branches_glm,
    group_branches(
      time_horizons,
      obs_model_glm = obs_model_glm,
      fitting_method = "glm"
    ),
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
      train_data$skip_rows
    ),
    pattern = map(branches_mcmc, train_data),
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
      prior_delay_param = select(branches_mcmc, starts_with("delay_")),
      mean_log = branches_mcmc$mean_log,
      sd_log = branches_mcmc$sd_log,
      stan_settings = stan_settings,
      sensitivity_scenario_name = branches_mcmc$scenario_name
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
  # - posterior density of the mean process for the weeks, where we perform
  #   nowcasting
  # - the scatter plot of the dispersion parameter against the standard
  #   deviation of the random walk.
  tar_target(rolling_plots_mcmc, {
    plot_per_window(
      summarized_nowcast_mcmc,
      fitted_mcmc$delay_prob,
      fitted_mcmc$nb_size,
      fitted_mcmc$lambda,
      fitted_mcmc$rw_sd,
      df_total,
      obs_model,
      time_horizons$nowcast_date,
      fitting_method = "mcmc",
      prob_prior_pars = select(
        branches_mcmc,
        c("model_name", "scenario_name", paste0("delay_", seq_len(max_lag) - 1))
      ),
      disp_prior_pars = select(
        branches_mcmc,
        c("model_name", "mean_log", "sd_log", "scenario_name")
      ),
      data_origin = "case_study"
    )
  },
  pattern = sample(
    map(
      fitted_mcmc,
      time_horizons,
      df_total,
      summarized_nowcast_mcmc,
      branches_mcmc
    ),
    n = 15
  ),
  iteration = "list"),
  # Create plots of aggregated results from the MCMC method. We plot:
  # - the coverage of nowcasts,
  # - the crps decomposition.
  tar_target(aggreg_plots_mcmc, {
    plot_aggregated(
      bind_rows(summarized_nowcast_mcmc),
      obs_model,
      fitting_method = "mcmc",
      data_origin = "case_study"
    )
  }),
  # Extract the diagnostic summaries for the MCMC models. We do it per branch to
  # avoid loading all fits at once when we want to plot the diagnostics into a
  # single plot.
  tar_target(diagnostics, fitted_mcmc$diagnostics, pattern = map(fitted_mcmc)),
  # Plot the diagnostic summaries for the MCMC models
  tar_target(plot_diagnostics, {
    plot_mcmc_diagnostics(
      diagnostics,
      obs_model,
      data_origin = "case_study"
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
  # - density of estimates of the delay probability,
  # - density of estimates of the mean process for the weeks, where we
  #   perform nowcasting
  # - density of estimates of the dispersion parameter on a scale, where 0 means
  #   the Poisson model and higher values indicate more dispersion.
  tar_target(rolling_plots_glm, {
    plot_per_window(
      summarized_nowcast_glm,
      fitted_glm$delay_prob,
      fitted_glm$nb_size,
      fitted_glm$lambda,
      NULL,  # We don't have the random walk parameters
      df_total,
      obs_model_glm,
      time_horizons$nowcast_date,
      fitting_method = "glm",
      data_origin = "case_study"
    )
  },
  pattern = sample(
    map(fitted_glm, time_horizons, df_total, summarized_nowcast_glm),
    n = 15
  ),
  iteration = "list"),
  # Create plots of aggregated results from the GLM method. We plot:
  # - the coverage of nowcasts,
  # - the crps decomposition.
  tar_target(aggreg_plots_glm, {
    plot_aggregated(
      bind_rows(summarized_nowcast_glm),
      obs_model_glm,
      fitting_method = "glm",
      data_origin = "case_study"
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
      aux_analysis_start_date,
      data_origin = "case_study"
    )
  })
)
