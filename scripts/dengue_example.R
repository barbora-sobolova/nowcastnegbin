library("tidyverse")
library("cmdstanr")
library("tidybayes")
theme_set(theme_bw())

# Script running through 6 relatively simple nowcasting models, demonstrating
# different results. 

######################################
## Load and preprocess the dengue data
######################################

load("data/denguedat.RData")

# Create the reporting table
df_dengue <- denguedat |> group_by(onset_week, report_week) |>
  summarise(
    reports = n()
  ) |>
  mutate(
    delay = as.numeric((report_week - onset_week)) / 7
  )
dengue_observed <- df_dengue |> group_by(onset_week) |>
  summarise(observed = sum(reports))

# Way too long
ggplot(dengue_observed, aes(x = onset_week, y = observed)) +
  geom_line()

# Select a subset
begin_date <- as.Date("2005-05-16", format = "%Y-%m-%d")
end_date <- as.Date("2007-09-24", format = "%Y-%m-%d")  # inclusive
dates_seq <- seq(begin_date, end_date, by = 7)
df_dengue_subset <- filter(
  df_dengue,
  onset_week >= begin_date &
    report_week <= end_date
)
dengue_observed_subset <- df_dengue_subset |> group_by(onset_week) |>
  summarise(observed = sum(reports))

ggplot(dengue_observed_subset, aes(x = onset_week, y = observed)) +
  geom_line()

# What is the maximum lag? d = maximum lag + 1, because minimum lag is 0
max_lag_df <- df_dengue_subset |> group_by(onset_week) |>
  summarise(max_lag = max(delay))

# Create the indexing variables
d <- max(max_lag_df$max_lag) + 1
n <- as.numeric((end_date - begin_date) / 7 + 1)  # including the endpoints
p <- c(rep(d, n - d + 1), (d - 1):1)

# We need to pad the week-delay pairs, where no reports occurred, by zeros
df_all_dates <- tibble(
  onset_week = rep(dates_seq, times = p),
  delay = unlist(sapply(p, function (x) {0:(x - 1)}))
) |> mutate(report_week = onset_week + delay * 7)
df_dengue_subset_padded <- full_join(
  df_all_dates, 
  df_dengue_subset,
  by = c("onset_week", "report_week", "delay")
)
df_dengue_subset_padded$reports <- replace_na(df_dengue_subset_padded$reports, 0)

# List to pass to the stan models
dengue_stan_data <- list(
  n = n,
  m = nrow(df_dengue_subset_padded),
  p = p,
  obs = df_dengue_subset_padded$reports,
  d = d
)

#################
## Fitting models
#################

# Declare models
mod <- cmdstan_model("stan/nowcast.stan")

# Sample from the models
pois_nowcast_fit <- mod$sample(
  data = c(dengue_stan_data, model_obs = 0),
  parallel_chains = 4
)
nbinX_nowcast_fit <- mod$sample(
  data = c(dengue_stan_data, model_obs = 1), 
  parallel_chains = 4
)
nbin2D_nowcast_fit <- mod$sample(
  data = c(dengue_stan_data, model_obs = 2),
  parallel_chains = 4
)
nbin1D_nowcast_fit <- mod$sample(
  data = c(dengue_stan_data, model_obs = 3), 
  parallel_chains = 4
)
nbin2M_nowcast_fit <- mod$sample(
  data = c(dengue_stan_data, model_obs = 4), 
  parallel_chains = 4
)
nbin1M_nowcast_fit <- mod$sample(
  data = c(dengue_stan_data, model_obs = 5),
  parallel_chains = 4
)


#####################
## Gather the results
#####################

# Nowcasts
nbinX_nowcast <- nbinX_nowcast_fit |>
  gather_draws(nowcast[week]) |>
  ungroup() |>
  mutate(Distribution = "NegBinX", week_date = dates_seq[week])
nbin2D_nowcast <- nbin2D_nowcast_fit |>
  gather_draws(nowcast[week]) |>
  ungroup() |>
  mutate(Distribution = "NegBin2D", week_date = dates_seq[week])
nbin1D_nowcast <- nbin1D_nowcast_fit |>
  gather_draws(nowcast[week]) |>
  ungroup() |>
  mutate(Distribution = "NegBin1D", week_date = dates_seq[week])
nbin2M_nowcast <- nbin2M_nowcast_fit |>
  gather_draws(nowcast[week]) |>
  ungroup() |>
  mutate(Distribution = "NegBin2M", week_date = dates_seq[week])
nbin1M_nowcast <- nbin1M_nowcast_fit |>
  gather_draws(nowcast[week]) |>
  ungroup() |>
  mutate(Distribution = "NegBin1M", week_date = dates_seq[week])
pois_nowcast <- pois_nowcast_fit |>
  gather_draws(nowcast[week]) |>
  ungroup() |>
  mutate(Distribution = "Poisson", week_date = dates_seq[week])
nowcast <- rbind(
  nbinX_nowcast,
  nbin2D_nowcast,
  nbin1D_nowcast,
  nbin2M_nowcast,
  nbin1M_nowcast,
  pois_nowcast
) |> filter(
  week %in% (n - d + 1):n
)
nowcast_plot_tib <- nowcast  |> 
  group_by(week, week_date, Distribution) |>
  summarise(
    mean = mean(.value),
    quantiles = list(
      as_tibble(
        as.list(
          quantile(.value, probs = c(0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 
                                     0.75, 0.9, 0.95, 0.975, 0.99, 0.995))
        )
      )
    )
  ) |>
  unnest(quantiles)

# Calculate the variances of the nowcasts
tib_variances <- nowcast |> group_by(Distribution, week, week_date) |>
  summarise(variance = var(.value)) |> 
  filter(variance > 0)

# Delays
nbinX_delays <- nbinX_nowcast_fit |>
  gather_draws(reporting_delay[week]) |>
  ungroup() |>
  mutate(Distribution = "NegBinX", week_date = dates_seq[week])
nbin2D_delays <- nbin2D_nowcast_fit |>
  gather_draws(reporting_delay[week]) |>
  ungroup() |>
  mutate(Distribution = "NegBin2D", week_date = dates_seq[week])
nbin1D_delays <- nbin1D_nowcast_fit |>
  gather_draws(reporting_delay[week]) |>
  ungroup() |>
  mutate(Distribution = "NegBin1D", week_date = dates_seq[week])
nbin2M_delays <- nbin2M_nowcast_fit |>
  gather_draws(reporting_delay[week]) |>
  ungroup() |>
  mutate(Distribution = "NegBin2M", week_date = dates_seq[week])
nbin1M_delays <- nbin1M_nowcast_fit |>
  gather_draws(reporting_delay[week]) |>
  ungroup() |>
  mutate(Distribution = "NegBin1M", week_date = dates_seq[week])
pois_delays <- pois_nowcast_fit |>
  gather_draws(reporting_delay[week]) |>
  ungroup() |>
  mutate(Distribution = "Poisson", week_date = dates_seq[week])
delays <- rbind(nbinX_delays, nbin2D_delays, nbin1D_delays, nbin1M_delays,
                nbin2M_delays, pois_delays)

# Dispersion parameter
nbinX_disp <- nbinX_nowcast_fit |>
  gather_draws(nb_size) |>
  mutate(Distribution = "NegBinX")
nbin2D_disp <- nbin2D_nowcast_fit |>
  gather_draws(nb_size) |>
  mutate(Distribution = "NegBin2D")
nbin1D_disp <- nbin1D_nowcast_fit |>
  gather_draws(nb_size) |>
  mutate(Distribution = "NegBin1D")
nbin2M_disp <- nbin2M_nowcast_fit |>
  gather_draws(nb_size) |>
  mutate(Distribution = "NegBin2M")
nbin1M_disp <- nbin1M_nowcast_fit |>
  gather_draws(nb_size) |>
  mutate(Distribution = "NegBin1M")
df_disp <- rbind(nbinX_disp, nbin2D_disp, nbin1D_disp, nbin1M_disp, nbin2M_disp)

# Final data and the observed data
timepoint_to_plot <- d - 1 + 10
df_true <- df_dengue |> 
  rename("week_date" = "onset_week") |> 
  group_by(week_date) |>
  summarise(true_val = sum(reports)) |>
  filter(week_date %in% dates_seq) |>
  mutate(week = 1:n, Data = "Final") |> 
  filter(week %in% (n - timepoint_to_plot):n)
df_available <- df_dengue_subset_padded |> 
  rename("week_date" = "onset_week") |> 
  group_by(week_date) |>
  summarise(true_val = sum(reports)) |>
  filter(week_date %in% dates_seq) |>
  mutate(week = 1:n, Data = "Available") |> 
  filter(week %in% (n - timepoint_to_plot):n)
df_obs <- rbind(df_available, df_true) |> 
  mutate( Data = factor(Data, levels = c("Available", "Final"), ordered = TRUE))

###################
## Plot the results
###################

model_colors <- c("NegBin2D" = "darkgreen", "NegBin1D" = "dodgerblue",
                  "NegBin1M" = "darkgoldenrod1", "NegBinX" = "#E31A1C",
                  "NegBin2M" = "navyblue", "Poisson" = "magenta")

# Plot the nowcasts
ggplot() +
  geom_line(
    data = df_obs,
    aes(x = week_date, y = true_val, color = Data)
    ) +
  geom_line(
    data = nowcast_plot_tib,
    aes(x = week_date, y = `50%`, color = Distribution)
  ) +
  geom_ribbon(
    data = nowcast_plot_tib,
    aes(x = week_date, ymin = `2.5%`, ymax = `97.5%`, color = Distribution, 
        fill = Distribution),
    alpha = 0.2,
    linetype = 0
  ) +
  scale_color_manual(
    values = c(model_colors, "Available" = "gray", "Final" = "black"),
    breaks = c("Available", "Final")
    ) +
  scale_fill_manual(values = model_colors) +
  scale_x_date(date_labels = "%Y %b") +
  labs(x = "Week", y = "Cases") +
  facet_wrap(~Distribution)
ggsave("figure/Nowcasts_example.pdf", width = 7, height = 5)

# Plot the nowcast variances in time
ggplot(tib_variances, aes(x = week, y = variance, color = Distribution)) +
  scale_color_manual(values = model_colors) +
  geom_line(linewidth = 0.5)

# Plot the posterior of the dispersion parameter
ggplot(df_disp, aes(x = .value, color = Distribution)) +
  geom_density(linewidth = 0.8) +
  scale_color_manual(
    values = model_colors[-6],
  ) +
  labs(x = "size")

# Plot the delays
ggplot(delays, aes(x = .value, color = Distribution)) +
  geom_density(linewidth = 0.8) +
  scale_color_manual(
    values = model_colors,
  ) +
  labs(x = "Delay") +
  facet_wrap(~week, scales = "free")
