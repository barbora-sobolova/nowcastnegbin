#' Plot and save the nowcasts
#'
#' @description This function plots and possibly saves the nowcasts from all
#' observation models for a single date. The nowcasts are plotted along the
#' data for comparison.
#'
#' @param df_summarized_nowcast a data frame containing columns `Distribution`
#' (containing the name of the observation model), `quantile_50`
#' (the point nowcasts), `quantile_2.5`, `quantile_25`, `quantile_75`,
#' `quantile_97.5` (bounds of the prediction intervals) and `date` (x-axis
#' dates)
#' @param df_total a data frame with columns `date`, `counts` and `data`
#' returned by the function \code{create_totals_data_frame()}
#' @param model_names a vector of names of the observation models, we wish to
#' plot.
#' @param date_of_the_nowcast a date, when the nowcast is made
#' @param fitting_method a method used for fitting the nowcasting model, either
#' "mcmc", or "glm"
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggplot2::ggsave()}
#'
#' @return a ggplot object with one facet per observation model, or NULL if
#' \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#'
#' @export
plot_nowcast <- function(
  df_summarized_nowcast,
  df_total,
  model_names,
  date_of_the_nowcast,
  fitting_method = c("mcmc", "glm"),
  save_plot = TRUE
) {
  fitting_method <- match.arg(fitting_method)
  df_total <- df_total |>
    mutate(data = factor(data, levels = c("Preliminary", "Final")))
  # Plot the nowcasts
  nowcasts_plot <- ggplot() +
    # Point prediction
    geom_line(
      data = df_summarized_nowcast,
      mapping = aes(
        x = .data$date,
        y = .data$quantile_50,
        color = .data$Distribution,
        linetype = "Nowcast"
      )
    ) +
    # Different data versions - preliminary and final
    geom_line(
      data = df_total,
      mapping = aes(
        x = .data$date,
        y = .data$counts,
        color = .data$data,
        linetype = .data$data
      )
    ) +
    # 95% prediction intervals
    geom_ribbon(
      data = df_summarized_nowcast,
      mapping = aes(
        x = .data$date,
        ymin = .data$quantile_2.5,
        ymax = .data$quantile_97.5,
        fill = .data$Distribution,
        alpha = "PI_95"
      )
    ) +
    # 50% prediction intervals
    geom_ribbon(
      data = df_summarized_nowcast,
      mapping = aes(
        x = .data$date,
        ymin = .data$quantile_25,
        ymax = .data$quantile_75,
        fill = .data$Distribution,
        alpha = "PI_50"
      )
    ) +
    # Set the color of the models and the data versions
    scale_color_manual(
      values = c(
        get_model_colors()[model_names],
        "Final" = "black",
        "Preliminary" = "gray60"
      ),
      guide = "none"
    ) +
    scale_linetype_manual(
      name = "Type of data",
      values = c(
        "Preliminary" = "solid",
        "Final" = "solid",
        "Nowcast" = "dashed"
      ),
    ) +
    # Set the transparency of the prediction intervals. The transparency is the
    # same value for both prediction intervals, since the 50% interval is inside
    # the 95% one. Hence, the transparency adds up.
    scale_alpha_manual(
      values = c("PI_50" = 0.2, "PI_95" = 0.2),
      labels = c("50%", "95%"),
      name = "Prediction interval"
    ) +
    guides(
      # Customize the legend of the data versions
      linetype = guide_legend(
        override.aes = list(color = c("gray60", "black", "black"))
      ),
      # Customize the legend of the confidence interval transparency, which for
      # the reader appear to be 0.4 for the 95% interval due to the "stacking"
      # of the layers.
      alpha = guide_legend(override.aes = list(alpha = c(0.4, 0.2)))
    ) +
    scale_fill_manual(values = get_model_colors()[model_names]) +
    labs(x = "Date", y = "Incidence") +
    facet_wrap(~Distribution)
  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(
      nowcasts_plot,
      path = paste(
        "inst/figure/nowcast_plots/nowcast",
        fitting_method,
        date_of_the_nowcast,
        sep = "_"
      ),
      width = 9,
      height = 7
    )
    ret <- NULL
  } else {
    ret <- nowcasts_plot
  }
  ret
}

#' Plot and save the coverage
#'
#' @description This function plots and possibly saves the chart of empirical
#' vs. nominal coverage. Only 50% and 95% coverage is considered.
#'
#' @param df_summarized_nowcast a data frame containing columns `Distribution`
#' (containing the name of the observation model), `quantile_2.5`,
#' `quantile_25`, `quantile_75`, `quantile_97.5` (bounds of the prediction
#' intervals), `delay` (the nowcasting horizon) and the true value of the
#' prediction target `true_val`
#' @param model_names a vector of names of the observation models, we wish to
#' plot.
#' @param fitting_method a method used for fitting the nowcasting model, either
#' "mcmc", or "glm"
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggsave()}
#'
#' @return a ggplot object with one facet per nowcasting horizon, or NULL if
#' \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#'
#' @export
plot_coverage <- function(
  df_summarized_nowcast,
  model_names,
  fitting_method = c("mcmc", "glm"),
  save_plot = TRUE
) {
  # Calculate the empirical coverage
  df_coverage <- df_summarized_nowcast |>
    group_by(.data$delay, .data$Distribution) |>
    summarize(
      coverage_50 = sum(
        .data$true_val >= .data$quantile_25 &
          .data$true_val <= .data$quantile_75
      ) / n(),
      coverage_95 = sum(
        .data$true_val >= .data$quantile_2.5 &
          .data$true_val <= .data$quantile_97.5
      ) / n(),
      .groups = "drop"
    ) |>
    # Pivot for easier definition of the alpha aesthetic
    tidyr::pivot_longer(
      cols = starts_with("coverage"),
      names_to = "nominal_coverage",
      values_to = "empirical_coverage"
    )
  # Plot the empirical coverage as horizontal bars
  coverage_plot <- ggplot(
    df_coverage,
    aes(
      x = .data$empirical_coverage,
      y = .data$Distribution,
      fill = .data$Distribution,
      alpha = .data$nominal_coverage
    )
  ) +
    geom_col(position = "identity") +
    # Highlight the 50% and 95% nominal coverage
    geom_vline(xintercept = c(0.5, 0.95), linetype = "dotted") +
    scale_alpha_manual(
      values = c("coverage_50" = 1, "coverage_95" = 0.5),
      labels = c("50% coverage", "95% coverage"),
      name = ""
    ) +
    scale_fill_manual(values = get_model_colors()[model_names]) +
    labs(x = "Empirical coverage", title = "Empirical coverage by horizon") +
    facet_wrap(~delay)
  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(
      coverage_plot,
      paste("inst/figure/coverage_plot", fitting_method, sep = "_"),
      width = 9,
      height = 7
    )
    ret <- NULL
  } else {
    ret <- coverage_plot
  }
  ret
}

#' Plot and save the decomposition of the CRPS
#'
#' @description This function plots and possibly saves the decomposition of
#' the average CRPS decomposed according to the spread, underprediction and
#' overprediction.
#'
#' @param df_summarized_nowcast a data frame containing columns `Distribution`
#' (containing the name of the observation model), `dispersion`,
#' `underprediction`, `overprediction` and `delay` (the nowcasting horizon).
#' @param model_names a vector of names of the observation models, we wish to
#' plot.
#' @param fitting_method a method used for fitting the nowcasting model, either
#' "mcmc", or "glm"
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggsave()}
#'
#' @return a ggplot object with one facet per nowcasting horizon, or NULL if
#' \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#'
#' @export
plot_crps_decomp <- function(
  df_summarized_nowcast,
  model_names,
  fitting_method = c("mcmc", "glm"),
  save_plot = TRUE
) {
  # Calculate the decomposition of the average CRPS
  df_crps_decomp <- df_summarized_nowcast |>
    group_by(.data$delay, .data$Distribution) |>
    summarize(
      Spread = mean(.data$dispersion),
      Underprediction = mean(.data$underprediction),
      Overprediction = mean(.data$overprediction),
      .groups = "drop"
    ) |>
    # Pivot for easier definition of the alpha aesthetic
    tidyr::pivot_longer(
      cols = c("Spread", "Overprediction", "Underprediction"),
      names_to = "Component",
      values_to = "CRPS"
    )
  # Plot the empirical coverage as horizontal bars
  crps_decomp_plot <- ggplot(
    df_crps_decomp,
    aes(
      x = .data$CRPS,
      y = .data$Distribution,
      fill = .data$Distribution,
      alpha = .data$Component
    )
  ) +
    geom_col(position = "stack") +
    scale_alpha_manual(
      values = c("Underprediction" = 1, "Spread" = 0.4, "Overprediction" = 0.7),
      name = ""
    ) +
    scale_fill_manual(values = get_model_colors()[model_names]) +
    labs(x = "Mean CRPS", title = "CRPS decomposition by horizon") +
    facet_wrap(~delay, scales = "free_x")
  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(
      crps_decomp_plot,
      paste("inst/figure/crps_decomposition_plot", fitting_method, sep = "_"),
      width = 9,
      height = 7
    )
    ret <- NULL
  } else {
    ret <- crps_decomp_plot
  }
  ret
}

#' Plot and save the density plot of the dispersion parameter estimates
#'
#' @description This function plots and possibly saves the densities of
#' dispersion parameter estimates for all fitted negative binomial models.
#'
#' @param df_nb_size a data frame containing columns `Distribution`
#' (containing the name of the observation model), `.value` (the empirical
#' distribution of the dispersion parameter estimates) and `nowcast_date` (the
#' date when the nowcast is calculated)
#' @param model_names a vector of names of the observation models, we wish to
#' plot.
#' @param date_of_the_nowcast a date, when the nowcast is made to filter the
#' \code{df_nb_size} table
#' @param fitting_method a method used for fitting the nowcasting model, either
#' "mcmc", or "glm"
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggsave()}
#'
#' @return a ggplot object with one facet showing the density of the dispersion
#' parameter estimates, or NULL if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#'
#' @export
plot_disp_par <- function(
  df_nb_size,
  model_names,
  date_of_the_nowcast,
  fitting_method = c("mcmc", "glm"),
  save_plot = TRUE
) {
  # Filter only values from the corresponding time window
  df_nb_size <- df_nb_size |>
    mutate(
      # Plot the dispersion parameter on the inverted scale, where higher values
      # indicate more dispersion and 0 corresponds to the Poisson model in the
      # limit.
      phi = 1 / .data$.value
    )

  disp_par_plot <- ggplot() +
    # Plot the density of the dispersion parameter estimates
    geom_line(
      df_nb_size,
      mapping = aes(x = .data$phi, color = .data$Distribution),
      stat = "density",
      alpha = 0.7
    )
  if (fitting_method == "mcmc") {
    # Draw a line representing the prior distribution. This must be checked
    # manually to correspond to the prior we are using in STAN.

    # Currently we use the inverse gamma distribution with parameters 0.41 and
    # 0.29. Since we plot on the inverted scale, we will plot the density of the
    # gamma distribution with the same parameters.
    df_prior <- data.frame(
      phi = seq(0, max(df_nb_size$phi), length = 500)
    ) |>
      mutate(
        dens = dgamma(.data$phi, 0.41, 0.29)
      )

    disp_par_plot <- disp_par_plot +
      geom_line(
        data = df_prior,
        aes(x = .data$phi, y = .data$dens, color = "Prior")
      ) +
      geom_segment(
        aes(x = 0, y = -0.2, xend = 54, yend = -0.2),
        arrow = arrow()
      ) +
      geom_text(aes(x = 27, y = -0.3, label = "more dispersion")) +
      scale_color_manual(values = c(get_model_colors(), "Prior" = "black")) +
      labs(
        x = "dispersion parameter",
        title = "Dispersion parameter posterior"
      ) +
      coord_cartesian(ylim = c(-0.4, 2))
  } else {
    disp_par_plot <- disp_par_plot +
      geom_segment(
        aes(x = 0, y = -0.05, xend = 130, yend = -0.05),
        arrow = arrow()
      ) +
      geom_text(aes(x = 65, y = -0.1, label = "more dispersion")) +
      scale_color_manual(
        values = c(get_model_colors()[model_names], "Prior" = "black")
      ) +
      labs(
        x = "dispersion parameter",
        title = "Posterior of the dispersion parameter"
      ) +
      coord_cartesian(ylim = c(-0.15, 0.5))
  }

  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(
      disp_par_plot,
      paste(
        "inst/figure/disp_plots/disp_par_plot",
        fitting_method,
        date_of_the_nowcast,
        sep = "_"
      ),
      width = 7,
      height = 5.5
    )
    ret <- NULL
  } else {
    ret <- disp_par_plot
  }
  ret
}

#' Plot and save the density plot of the delay probability estimates
#'
#' @description This function plots and possibly saves the densities of
#' delay probability estimates for all fitted models.
#'
#' @param df_delay_prob a data frame containing columns `Distribution`
#' (containing the name of the observation model), `.value` (the empirical
#' distribution of the delay probability estimates), `nowcast_date` (the
#' date when the nowcast is calculated) and `delay` (the discrete delay time)
#' @param model_names a vector of names of the observation models, we wish to
#' plot.
#' @param date_of_the_nowcast a date, when the nowcast is made to filter the
#' \code{df_delay_prob} table
#' @param fitting_method a method used for fitting the nowcasting model, either
#' "mcmc", or "glm"
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggsave()}
#'
#' @return a ggplot object with one facet per delay showing the density of the
#' delay probability estimates per delay, or NULL if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#'
#' @export
plot_delay_prob <- function(
  df_delay_prob,
  model_names,
  date_of_the_nowcast,
  fitting_method = c("mcmc", "glm"),
  save_plot = TRUE
) {
  # Filter only values from the corresponding time window
  df_delay_prob <- df_delay_prob |>
    # Turn the delay into a factor to allow for easier faceting
    mutate(delay = factor(.data$delay))

  delay_prob_plot <- ggplot(
    df_delay_prob,
    aes(x = .data$.value, color = .data$Distribution)
  ) +
    # Plot the density of the dispersion parameter estimates
    geom_line(stat = "density", alpha = 0.6, bounds = c(0, 1)) +
    scale_color_manual(values = get_model_colors()[model_names]) +
    labs(
      x = "delay probability",
      title = "Posterior of the delay probability"
    ) +
    coord_cartesian(ylim = c(0, 80)) +
    facet_wrap(~delay, scales = "free")
  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(
      delay_prob_plot,
      paste(
        "inst/figure/delay_prob_plots/delay_prob_plot",
        fitting_method,
        date_of_the_nowcast,
        sep = "_"
      ),
      width = 7,
      height = 5.5
    )
    ret <- NULL
  } else {
    ret <- delay_prob_plot
  }
  ret
}

#' Plot and save the scatter plot of the dispersion parameter against the random
#' walk standard deviation
#'
#' @description This function plots and possibly saves the scatter plot of the
#' dispersion parameter estimates against the estimates of the standard
#' deviation of the random walk.
#'
#' @param df_rw_sd a data frame containing columns `Distribution`
#' (containing the name of the observation model), `.value` (the empirical
#' distribution of the estimates of the random walk standard deviation) and
#' `nowcast_date` (the date when the nowcast is calculated)
#' @param df_nb_size a data frame containing columns `Distribution`
#' (containing the name of the observation model), `.value` (the empirical
#' distribution of the dispersion parameter estimates) and `nowcast_date` (the
#' date when the nowcast is calculated)
#' @param model_names a vector of names of the observation models, we wish to
#' plot.
#' @param date_of_the_nowcast a date, when the nowcast is made to filter the
#' \code{df_delay_prob} table
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggsave()}
#'
#' @return a ggplot object with one facet per distribution showing the scatter
#' plot of the dispersion parameter vs. random walk standard deviation sampled
#' values, or NULL if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#'
#' @export
plot_rw_sd <- function(
  df_rw_sd,
  df_nb_size,
  model_names,
  date_of_the_nowcast,
  save_plot = TRUE
) {
  # Rename the columns with the parameter values to avoid two columns with the
  # name ".value"
  df_rw_sd <- df_rw_sd |>
    rename("rw_sd" = ".value") |>
    filter(.data$Distribution != "Poisson")
  df_nb_size <- df_nb_size |>
    rename("nb_size" = ".value") |>
    # Invert the scale, so that larger values mean more overdispersion
    mutate(phi = 1 / .data$nb_size)
  df_join <- inner_join(
    df_nb_size,
    df_rw_sd,
    by = c(".draw", "Distribution", ".chain", ".iteration", "nowcast_date"),
    relationship = "one-to-one"
  )

  rw_sd_scatter <- ggplot(
    df_join,
    aes(x = .data$phi, y = .data$rw_sd, color = .data$Distribution)
  ) +
    # Plot the density of the dispersion parameter estimates
    geom_point(alpha = 0.5, shape = 1) +
    scale_color_manual(values = get_model_colors()[model_names]) +
    labs(
      x = "dispersion parameter",
      y = "sd of the random walk increments"
    ) +
    facet_wrap(~Distribution, scales = "free")
  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(
      rw_sd_scatter,
      paste(
        "inst/figure/rw_sd_plots/rw_sd_scatter_plot",
        date_of_the_nowcast,
        sep = "_"
      ),
      width = 7,
      height = 5.5
    )
    ret <- NULL
  } else {
    ret <- rw_sd_scatter
  }
  ret
}

#' A wrapper around plotting functions creating all relevant per-window plots
#'
#' @description This is a wrapper around functions that plot the posterior
#' distributions of parameters in each rolling-window: \code{plot_nowcast()},
#' \code{plot_delay_prob()}, \code{plot_disp_par} and in the case of the MCMC
#' method also \code{plot_rw_sd()}.
#'
#' @param df_nowcast a data frame containing columns `Distribution`
#' (containing the name of the observation model), `quantile_50`
#' (the point nowcasts), `quantile_2.5`, `quantile_25`, `quantile_75`,
#' `quantile_97.5` (bounds of the prediction intervals) and `date` (x-axis
#' dates)
#' @param df_delay_prob a data frame containing columns `Distribution`
#' (containing the name of the observation model), `.value` (the empirical
#' distribution of the delay probability estimates), `nowcast_date` (the
#' date when the nowcast is calculated) and `delay` (the discrete delay time)
#' @param df_disp_par a data frame containing columns `Distribution`
#' (containing the name of the observation model), `.value` (the empirical
#' distribution of the dispersion parameter estimates) and `nowcast_date` (the
#' date when the nowcast is calculated)
#' @param df_rw_sd a data frame containing columns `Distribution`
#' (containing the name of the observation model), `.value` (the empirical
#' distribution of the estimates of the random walk standard deviation) and
#' `nowcast_date` (the date when the nowcast is calculated)
#' @param df_total a data frame with columns `date`, `counts` and `data`
#' returned by the function \code{create_totals_data_frame()}
#' @param model_names a vector of names of the observation models, we wish to
#' plot.
#' @param date_of_the_nowcast a date, when the nowcast is made
#' @param fitting_method a method used for fitting the nowcasting model, either
#' "mcmc", or "glm"
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggsave()}
#'
#' @return a list of ggplot objects or list of NULLs if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#'
#' @export
plot_per_window <- function(
  df_nowcast,
  df_delay_prob,
  df_disp_par,
  df_rw_sd,
  df_total,
  model_names,
  date_of_the_nowcast,
  fitting_method = c("mcmc", "glm"),
  save_plot = TRUE
) {
  p_nowcast <- plot_nowcast(
    df_nowcast,
    df_total,
    model_names,
    date_of_the_nowcast,
    fitting_method,
    save_plot
  )
  p_disp <- plot_disp_par(
    df_disp_par,
    model_names,
    date_of_the_nowcast,
    fitting_method,
    save_plot
  )
  p_prob <- plot_delay_prob(
    df_delay_prob,
    model_names,
    date_of_the_nowcast,
    fitting_method,
    save_plot
  )
  ret_list <- list(nowcast = p_nowcast, delay_prob = p_prob, disp = p_disp)
  if (fitting_method == "mcmc") {
    p_rw_sd <- plot_rw_sd(
      df_rw_sd,
      df_disp_par,
      model_names,
      date_of_the_nowcast,
      save_plot
    )
    ret_list <- c(ret_list, p_rw_sd)
  }
  ret_list
}

#' A wrapper around plotting functions creating all relevant aggregated plots
#'
#' @description This is a wrapper around functions that plot the aggregated
#' results: \code{plot_coverage()} and \code{plot_crps_decomp()}.
#'
#' @param df_nowcast a data frame containing columns `Distribution`
#' (containing the name of the observation model), `quantile_50`
#' (the point nowcasts), `quantile_2.5`, `quantile_25`, `quantile_75`,
#' `quantile_97.5` (bounds of the prediction intervals), `dispersion` (spread
#' component of the CRPS), `underprediction` (CRPS component) and
#' `overprediction` (CRPS component). This data frame contains the whole period,
#' where we nowcasting has been done.
#' @param model_names a vector of names of the observation models, we wish to
#' plot.
#' @param fitting_method a method used for fitting the nowcasting model, either
#' "mcmc", or "glm"
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggsave()}
#'
#' @return a list of ggplot objects or list of NULLs if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#'
#' @export
plot_aggregated <- function(
  df_nowcast,
  model_names,
  fitting_method = c("mcmc", "glm"),
  save_plot = TRUE
) {
  p_coverage <- plot_coverage(
    df_nowcast,
    model_names,
    fitting_method,
    save_plot
  )
  p_crps_decomp <- plot_crps_decomp(
    df_nowcast,
    model_names,
    fitting_method,
    save_plot
  )
  ret_list <- list(coverage = p_coverage, crps_decomp = p_crps_decomp)
  ret_list
}

#' Plot the whole incidence trajectory
#'
#' @description This function plots and possibly saves the incidence trajectory
#' used for the case study. The first and the last estimation windows will be
#' highlighted to see the chunk of the data we use for model training.
#'
#' @param full_data a data frame of the whole trajectory containing columns
#' `date` and columns `value_0w`, `value_1w`, etc. until `max_lag - 1`.
#' @param start_date a date (indeed in the date format), where the training data
#' start. The starting point will be included.
#' @param length_of_train_data a number, the length of the estimation window
#' (endpoints included)
#' @param max_lag maximum reporting delay represented by the number of columns
#' of the reporting table. In this way, the 0-th lag counts as the first, 1-st
#' lag as the second and so on.
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggplot2::ggsave()}
#'
#' @return a ggplot object, or NULL if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#' @importFrom ggpubr geom_bracket
#'
#' @export
plot_trajectory <- function(
  full_data,
  start_date,
  length_of_train_data,
  max_lag,
  save_plot = TRUE
) {
  # Arrange the whole trajectory into a data frame for plotting
  totals <- full_data |>
    dplyr::select(paste0("value_", 1:max_lag - 1, "w")) |>
    create_totals_data_frame(start_date) |>
    # `create_totals_data_frame()` returns a long data frame containing the
    # final and the preliminary state of the data. For plotting the whole
    # trajectory we are interested only in the final values.
    dplyr::filter(data == "Final")
  # The end point of the first estimation window and the beginning of the last
  # estimation window to be highlighted in the plot.
  # The estimation windows will be highlighted by braces drawn by
  # `ggpubr::geom_bracket()`.
  first_window_end <- start_date + (length_of_train_data - 1) * 7
  last_window_beg <- start_date + (nrow(totals) - length_of_train_data - 1) * 7
  # We need to find the maximum number of cases in the first and last estimation
  # window in order to place the brace correctly above them.
  first_window_max_cases <- totals |>
    filter(date <= first_window_end) |>
    pull(.data$counts) |>
    # na.rm = TRUE is usually not needed, but it prevents the plot element to
    # disappear in the case of missing values
    max(na.rm = TRUE)
  last_window_max_cases <- totals |>
    filter(date >= last_window_beg) |>
    pull(.data$counts) |>
    max(na.rm = TRUE)
  # 5% offset of the braces to avoid overplotting the trajectory
  bracket_offset <- first_window_max_cases * 0.05
  trajectory_plot <- ggplot(totals, aes(x = .data$date, y = .data$counts)) +
    geom_line() +
    # Highlight the first window of training data excluding the nowcasting part
    ggpubr::geom_bracket(
      xmin = start_date,
      xmax = first_window_end - (max_lag - 2) * 7 - 1,
      y.position = first_window_max_cases + bracket_offset,
      label = "First chunk of\ntraining data",
      label.size = 3
    ) +
    # Highlight the first nowcasting target
    ggpubr::geom_bracket(
      xmin = first_window_end - (max_lag - 2) * 7 + 1,
      xmax = first_window_end,
      y.position = first_window_max_cases + bracket_offset,
      label = "First\nnowcasting\ntarget",
      label.size = 3
    ) +
    # Highlight the last window of training data excluding the nowcasting part
    ggpubr::geom_bracket(
      xmin = last_window_beg,
      xmax = last_window_beg + (length_of_train_data - max_lag + 2) * 7 - 1,
      y.position = last_window_max_cases + bracket_offset,
      label = "Last chunk of\ntraining data",
      label.size = 3
    ) +
    # Highlight the last nowcasting target
    ggpubr::geom_bracket(
      xmin = last_window_beg + (length_of_train_data - max_lag + 2) * 7 + 1,
      xmax = last_window_beg + length_of_train_data * 7,
      y.position = last_window_max_cases + bracket_offset,
      label = "Last\nnowcasting\ntarget",
      label.size = 3
    ) +
    labs(title = "SARI incidence", y = "Incidence")
  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(
      trajectory_plot,
      "inst/figure/SARI_trajectory",
      width = 9,
      height = 7
    )
    ret <- NULL
  } else {
    ret <- trajectory_plot
  }
  ret
}

#' Plot the summary of MCMC diagnostics
#'
#' @description This function plots and possibly saves the summary of the
#' MCMC fitting diagnostics for each time step.
#'
#' @param df_diagnostics a data frame containing the diagnostic summaries for
#' each model and each run (timesteps). It contains columns `num_divergent`,
#' `num_max_treedepth`, `ebfmi`, `Distribution`, `date_of_the_nowcast`.
#' @param model_names a vector of names of the observation models, we wish to
#' plot.
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggplot2::ggsave()}
#'
#' @return a ggplot object, or NULL if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#' @importFrom tidyr pivot_longer
#'
#' @export
plot_mcmc_diagnostics <- function(
  df_diagnostics,
  model_names,
  save_plot = TRUE
) {
  df_diagnostics_long <- df_diagnostics |>
    group_by(.data$Distribution, .data$nowcast_date) |>
    summarise(
      # Sum the numbers of problematic transitions from different chains
      num_max_treedepth = sum(.data$num_max_treedepth),
      num_divergent = sum(.data$num_divergent),
      # Find the chain with the lowest ebmfi value
      min_ebfmi = min(.data$ebfmi),
      .groups = "drop"
    ) |>
    # Create a long data frame to plot the number of divergent transitions,
    # the maximum tree depth and the lowest ebfmi in different facets
    pivot_longer(
      cols = c("num_divergent", "num_max_treedepth", "min_ebfmi"),
      names_to = "Quantity",
      values_to = "Value"
    ) |>
    group_by(.data$Distribution, .data$nowcast_date, .data$Quantity)

  date_breaks <- seq(
    min(df_diagnostics$nowcast_date),
    max(df_diagnostics$nowcast_date),
    by = 4 * 7
  )
  diag_plot <- ggplot(
    df_diagnostics_long,
    aes(
      x = .data$nowcast_date,
      y = .data$Value,
      color = .data$Distribution
    )
  ) +
    geom_line() +
    labs(y = NULL, x = "Date") +
    scale_color_manual(values = get_model_colors()[model_names]) +
    scale_x_date(breaks = date_breaks, date_labels = "%d %b") +
    facet_wrap(~Quantity, nrow = 3, scales = "free_y")

  # Save the plot if required, the width, height and path iare hard-coded here
  if (save_plot) {
    save_figure(
      diag_plot,
      "inst/figure/diagnostics_plot",
      width = 9,
      height = 7
    )
    ret <- NULL
  } else {
    ret <- diag_plot
  }
  ret
}

#' Save a figure in the PDF and the PNG format
#'
#' @param figure a ggplot chart to be saved
#' @param path a file path indicating where to save the plot, typically starting
#'  with "inst/figure". No file extension should be included.
#' @param width,height plot size as accepted by the \code{ggplot2::ggsave()}
#' function
#'
#' @return the path to the saved file resulting from the last
#' \code{ggplot2::ggsave()} will be returned as a string
#'
#' @import ggplot2
#'
#' @export
save_figure <- function(figure, path, width, height) {
  # Save the figure in PDF for a LaTeX manuscript
  ggsave(
    paste0(path, ".pdf"),
    figure,
    width = width,
    height = height,
    create.dir = TRUE
  )
  # Save the figure in PNG for non-LaTeX documents and version comparison tools
  ggsave(
    paste0(path, ".png"),
    figure,
    width = width,
    height = height,
    create.dir = TRUE
  )
}
