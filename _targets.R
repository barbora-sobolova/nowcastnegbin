library(targets)
library(tarchetypes)
library(qs2)

# Set targets options
tar_option_set(
  packages = c(
    "scoringutils",
    "ggplot2",
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
analysis_start_date <- as.Date("2023-10-08")
# How many weeks we want to include as "training" data.
# This includes the last `max_lag - 1` weeks for which we calculate the nowcast.
length_of_train_data <- 56

obs_model <- data.frame(
  model_name = c(
    "Poisson",
    "NegBinX",
    "NegBin2D",
    "NegBin1D"#,
    #  "NegBin2M",
    #  "NegBin1M"
  ),
  model_number = 0:3# 0:5
)

# Define the pipeline ==========================================================
list(
  tar_target(compiled_model, {
    cmdstanr::cmdstan_model(here::here("inst", "stan", "nowcast.stan"))
  }),
  tar_target(
    df_data,
    load_triangle(
      here::here(
        "data",
        "SARI",
        "reporting_triangle-icosari-sari-preprocessed.csv"
      ),
      start_date = as.Date("2023-10-08"),
      num_of_weeks = length_of_train_data,
      max_lag = max_lag
    )
  ),
  tar_target(df_total, {
    create_totals_data_frame(
      dplyr::select(df_data, paste0("value_", 1:max_lag - 1, "w")),
      max_lag,
      analysis_start_date,
      length_of_train_data
    )
  }),
  tar_target(
    stan_data,
    # Pass only the columns with the counts that are named as "value_0w",
    # "value_1w" and so on
    get_stan_data(dplyr::select(df_data, paste0("value_", 1:max_lag - 1, "w")))
  ),
  tar_target(stan_settings, {
    list(
      chains = 2,
      parallel_chains = 1,
      iter_warmup = 1000,
      iter_sampling = 1000,
      show_messages = FALSE,
      show_exceptions = FALSE,
      refresh = 0,
      seed = 12345
    )
  }),
  # Static branching over 6 observational models defined in `obs_model`.
  tar_map(
    unlist = TRUE,
    values = obs_model,
    names = model_name,
    # Fit for each observational model
    tar_target(fitted, {
      fit_stan_model(
        compiled_model$sample,
        stan_data = stan_data,
        model_obs = model_number,
        stan_settings = stan_settings
      )
    }),
    tar_target(df_summarized_nowcast, {
      summarize_nowcast(fitted$nowcast, df_total)
    })
  ),
  tar_target(nowcast_plot, {
    plot_nowcast(
      dplyr::bind_rows(
        df_summarized_nowcast_Poisson,
        df_summarized_nowcast_NegBinX,
        df_summarized_nowcast_NegBin2D,
        df_summarized_nowcast_NegBin1D
      ),
      df_total,
      model_codes = setNames(obs_model$model_name, obs_model$model_number),
      model_colors = model_colors
    )
  }),
  tar_target(saved_plot, {
    ggsave("nowcast_plot.png", nowcast_plot)
  })
)
