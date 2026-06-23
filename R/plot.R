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
#' @param date_of_the_nowcast a date, when the nowcast is made to name the
#' saved file correctly
#' @param fitting_method a method used for fitting the nowcasting model, either
#' "mcmc", or "glm"
#' @param data_origin a string indicating the data generating process of
#' simulated data, or whether the data correspond to the case study. Possible
#' values are "case_study", "NegBinX", "NegBin2D" and "NegBin1D"
#' @param sensitivity_sc a string indicating the sensitivity analysis scenario
#' of the MCMC method. Empty string "" indicates the main analysis.
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
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  sensitivity_sc = "",
  save_plot = TRUE
) {
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
    facet_wrap(~Distribution, nrow = 2)
  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(
      nowcasts_plot,
      path = paste0(
        paste(
          paste0("inst/figure/nowcast_plots/", data_origin, "/nowcast"),
          fitting_method,
          date_of_the_nowcast,
          sep = "_"
        ),
        sensitivity_sc
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
#' @param data_origin a string indicating the data generating process of
#' simulated data, or whether the data correspond to the case study. Possible
#' values are "case_study", "NegBinX", "NegBin2D" and "NegBin1D"
#' @param sensitivity_sc a string indicating the sensitivity analysis scenario
#' of the MCMC method. Empty string "" indicates the main analysis.
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
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  sensitivity_sc = "",
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
    facet_wrap(~delay, nrow = 2)
  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(
      coverage_plot,
      paste0(
        paste(
          "inst/figure/coverage_plot",
          data_origin,
          fitting_method,
          sep = "_"
        ),
        sensitivity_sc
      ),
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
#' @param data_origin a string indicating the data generating process of
#' simulated data, or whether the data correspond to the case study. Possible
#' values are "case_study", "NegBinX", "NegBin2D" and "NegBin1D"
#' @param sensitivity_sc a string indicating the sensitivity analysis scenario
#' of the MCMC method. Empty string "" indicates the main analysis.
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
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  sensitivity_sc = "",
  save_plot = TRUE
) {
  # Calculate the decomposition of the average CRPS
  df_crps_decomp <- df_summarized_nowcast |>
    # Calculate the absolute error
    mutate(AE = abs(.data$true_val - .data$quantile_50)) |>
    group_by(.data$delay, .data$Distribution) |>
    summarize(
      # Mean absolute error
      MAE = mean(.data$AE),
      # CRPS components
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
  crps_decomp_plot <- ggplot() +
    geom_col(
      df_crps_decomp,
      mapping = aes(
        x = .data$CRPS,
        y = .data$Distribution,
        fill = .data$Distribution,
        alpha = .data$Component
      ),
      position = "stack"
    ) +
    scale_alpha_manual(
      values = c("Underprediction" = 1, "Spread" = 0.4, "Overprediction" = 0.7),
      name = ""
    ) +
    geom_point(
      df_crps_decomp,
      mapping = aes(x = .data$MAE, y = .data$Distribution)
    ) +
    scale_fill_manual(values = get_model_colors()[model_names]) +
    labs(x = "Mean CRPS/AE", title = "CRPS decomposition by horizon") +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~delay, scales = "free_x", nrow = 2)
  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(
      crps_decomp_plot,
      paste0(
        paste(
          "inst/figure/crps_decomposition_plot",
          data_origin,
          fitting_method,
          sep = "_"
        ),
        sensitivity_sc
      ),
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
#' @param date_of_the_nowcast a date, when the nowcast is made to name the
#' saved file correctly
#' @param fitting_method a method used for fitting the nowcasting model, either
#' "mcmc", or "glm"
#' @param disp_prior_pars a data frame with columns `model_name`, `mean_log` and
#' `sd_log`, which contain the parameters for the prior log-normal distribution
#' of the dispersion parameter. The data frame should have 6 rows, one for each
#' model, although only the NegBinX, NegBin2D and NegBin1D rows will be used.
#' NegBin2M shares the prior with NegBinX and NegBin1M has the same prior as
#' NegBin1D.
#' @param data_origin a string indicating the data generating process of
#' simulated data, or whether the data correspond to the case study. Possible
#' values are "case_study", "NegBinX", "NegBin2D" and "NegBin1D"
#' @param sensitivity_sc a string indicating the sensitivity analysis scenario
#' of the MCMC method. Empty string "" indicates the main analysis.
#' @param true_value NULL for \code{data_origin = "case_study"}, otherwise the
#' true value of the dispersion parameter used to generate the data
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggsave()}
#'
#' @return a ggplot object with one facet showing the density of the dispersion
#' parameter estimates, or NULL if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#' @importFrom tidyr expand_grid
#'
#' @export
plot_disp_par <- function(
  df_nb_size,
  model_names,
  date_of_the_nowcast,
  fitting_method = c("mcmc", "glm"),
  disp_prior_pars = NULL,
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  sensitivity_sc = "",
  true_value = NULL,
  save_plot = TRUE
) {
  df_nb_size <- df_nb_size |>
    mutate(
      # Plot the dispersion parameter on the inverted scale, where higher values
      # indicate more dispersion and 0 corresponds to the Poisson model in the
      # limit.
      phi = 1 / .data$.value
    )
  # The right limit of the x-axis to allow us to zoom-in to the more interesting
  # part of the plot
  x_max <- quantile(df_nb_size$phi, 0.98)

  disp_par_plot <- ggplot() +
    # Plot the density of the dispersion parameter estimates
    geom_line(
      df_nb_size,
      mapping = aes(
        x = .data$phi,
        color = .data$Distribution,
        linetype = "Posterior"
      ),
      stat = "density",
      alpha = 0.7
    ) +
    geom_segment(
      aes(x = 0, y = -0.05, xend = 150, yend = -0.05),
      arrow = arrow()
    ) +
    geom_text(aes(x = 75, y = -0.1, label = "more dispersion")) +
    scale_color_manual(values = get_model_colors()[model_names]) +
    labs(
      x = "dispersion parameter",
      y = "density"
    ) +
    coord_cartesian(
      ylim = c(-0.15, 0.5),
      xlim = c(0, x_max)
    )
  if (fitting_method == "mcmc") {
    if (is.null(disp_prior_pars) || nrow(disp_prior_pars) == 0) {
      stop("`disp_prior_pars` must be provided when fitting_method = 'mcmc'.")
    }
    # Draw a line representing the prior distribution.
    df_prior <- tidyr::expand_grid(
      model_name = disp_prior_pars$model_name,
      phi = seq(0, x_max, length = 500)
    ) |>
      # Only 3 models have distinct priors. NegBin2M shares the prior with
      # NegBinX and NegBin1M has the same prior as NegBin1D.
      filter(.data$model_name %in% c("NegBinX", "NegBin2D", "NegBin1D")) |>
      inner_join(
        disp_prior_pars,
        by = "model_name",
        relationship = "many-to-one"
      ) |>
      rename("Distribution" = "model_name") |>
      mutate(
        dens = dlnorm(.data$phi, meanlog = .data$mean_log, sdlog = .data$sd_log)
      )

    disp_par_plot <- disp_par_plot +
      geom_line(
        data = df_prior,
        aes(
          x = .data$phi,
          y = .data$dens,
          color = .data$Distribution,
          linetype = "Prior"
        )
      ) +
      scale_linetype_manual(
        values = c("Prior" = "dotted", "Posterior" = "solid")
      ) +
      labs(title = "Dispersion parameter posterior")
  } else {
    # For the GLM method, remove the linetype aesthetics distinguishing between
    # the prior and posterior distribution from the legend
    disp_par_plot <- disp_par_plot +
      scale_linetype(guide = "none") +
      labs(title = "Dispersion parameter asymptotic distribution")
  }

  # If the data comes from a simulation, we will also plot the true parameter
  # value
  if (data_origin != "case_study") {
    disp_par_plot <- disp_par_plot +
      geom_vline(aes(xintercept = true_value), linetype = "dashed")
  }

  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(
      disp_par_plot,
      paste0(
        paste(
          paste0("inst/figure/disp_plots/", data_origin, "/disp_par_plot"),
          fitting_method,
          date_of_the_nowcast,
          sep = "_"
        ),
        sensitivity_sc
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
#' @param date_of_the_nowcast a date, when the nowcast is made to name the
#' saved file correctly
#' @param fitting_method a method used for fitting the nowcasting model, either
#' "mcmc", or "glm"
#' @param prob_prior_pars a vector of the prior parameters of the Dirichlet
#' delay probability distribution
#' @param data_origin a string indicating the data generating process of
#' simulated data, or whether the data correspond to the case study. Possible
#' values are "case_study", "NegBinX", "NegBin2D" and "NegBin1D"
#' @param sensitivity_sc a string indicating the sensitivity analysis scenario
#' of the MCMC method. Empty string "" indicates the main analysis.
#' @param true_value NULL for \code{data_origin = "case_study"}, otherwise the
#' true value of the delay probability vector used to generate the data
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggsave()}
#'
#' @return a ggplot object with one facet per delay showing the density of the
#' delay probability estimates, or NULL if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#'
#' @export
plot_delay_prob <- function(
  df_delay_prob,
  model_names,
  date_of_the_nowcast,
  fitting_method = c("mcmc", "glm"),
  prob_prior_pars = NULL,
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  sensitivity_sc = "",
  true_value = NULL,
  save_plot = TRUE
) {
  df_delay_prob <- df_delay_prob |>
    # Turn the delay into a factor to allow for easier faceting
    mutate(delay = factor(.data$delay))

  delay_prob_plot <- ggplot() +
    # Plot the density of the delay probability estimates
    geom_line(
      data = df_delay_prob,
      mapping = aes(
        x = .data$.value,
        color = .data$Distribution,
        linetype = "Posterior"
      ),
      stat = "density", alpha = 0.6
    ) +
    scale_color_manual(values = get_model_colors()[model_names]) +
    labs(x = "delay probability", y = "density") +
    coord_cartesian(ylim = c(0, 80)) +
    facet_wrap(~delay, scales = "free_x", nrow = 2)
  if (fitting_method == "mcmc") {
    if (is.null(prob_prior_pars) || length(prob_prior_pars) == 0) {
      stop("`prob_prior_pars` must be provided when fitting_method = 'mcmc'.")
    }
    # Draw a line representing the prior distribution.
    max_lag <- length(unique(df_delay_prob$delay))
    # Plot the prior distribution only within the region, where we have a
    # non-zero density
    x_axis_bounds <- df_delay_prob |>
      group_by(.data$delay) |>
      summarize(min_val = min(.data$.value), max_val = max(.data$.value)) |>
      mutate(
        max_val = ifelse(.data$max_val > 0.5, 1, .data$max_val),
        min_val = ifelse(.data$min_val < 0.5, 0, .data$min_val)
      )
    df_prior <- expand.grid(
      p = seq(0, 1, length = 800),
      delay = factor(seq_len(max_lag))
    ) |>
      inner_join(x_axis_bounds, by = "delay") |>
      filter(.data$p > .data$min_val, .data$p < .data$max_val) |>
      # The prior distribution is Dirichlet, so each marginal is beta
      # distributed
      mutate(
        beta_par1 = prob_prior_pars[.data$delay],
        beta_par2 = sum(prob_prior_pars) - prob_prior_pars[.data$delay],
        dens = dbeta(.data$p, .data$beta_par1, .data$beta_par2)
      )

    delay_prob_plot <- delay_prob_plot +
      geom_line(
        data = df_prior,
        mapping = aes(x = .data$p, y = .data$dens, linetype = "Prior")
      ) +
      scale_linetype_manual(
        values = c("Prior" = "dotted", "Posterior" = "solid")
      ) +
      labs(title = "Posterior of the delay probability")
  } else {
    # For the GLM method, remove the linetype aesthetics distinguishing between
    # the prior and posterior distribution from the legend
    delay_prob_plot <- delay_prob_plot +
      scale_linetype(guide = "none") +
      labs(title = "Delay probability asymptotic distribution")
  }

  if (data_origin != "case_study") {
    df_true_value <- data.frame(
      true_value = true_value,
      delay = factor(seq_along(true_value))
    )
    delay_prob_plot <- delay_prob_plot +
      geom_vline(
        data = df_true_value,
        mapping = aes(xintercept = true_value),
        linetype = "dashed"
      )
  }

  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(
      delay_prob_plot,
      paste0(
        paste(
          paste0(
            "inst/figure/delay_prob_plots/",
            data_origin,
            "/delay_prob_plot"
          ),
          fitting_method,
          date_of_the_nowcast,
          sep = "_"
        ),
        sensitivity_sc
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

#' Plot and save the density plot of the mean process estimates
#'
#' @description This function plots and possibly saves the densities of
#' estimates of the mean value of the total counts for all fitted models.
#'
#' @param df_lambda a data frame containing columns `Distribution`
#' (containing the name of the observation model), `.value` (the empirical
#' distribution of the delay probability estimates), `nowcast_date` (the
#' date when the nowcast is calculated) and `week` (the week number from the
#' beginning of the rolling window)
#' @param model_names a vector of names of the observation models, we wish to
#' plot.
#' @param date_of_the_nowcast a date, when the nowcast is made to name the
#' saved file correctly
#' @param max_lag an integer indicating the number of columns of the reporting
#' triangle
#' @param fitting_method a method used for fitting the nowcasting model, either
#' "mcmc", or "glm"
#' @param data_origin a string indicating the data generating process of
#' simulated data, or whether the data correspond to the case study. Possible
#' values are "case_study", "NegBinX", "NegBin2D" and "NegBin1D"
#' @param sensitivity_sc a string indicating the sensitivity analysis scenario
#' of the MCMC method. Empty string "" indicates the main analysis.
#' @param true_value NULL for \code{data_origin = "case_study"}, otherwise the
#' vector of true values of the mean of the total counts used to generate the
#' data. The length of the vector must be \code{max_lag - 1}
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggsave()}
#'
#' @return a ggplot object with one facet per week showing the density of the
#' estimates of the mean process, or NULL if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#'
#' @export
plot_mean_proc <- function(
  df_lambda,
  model_names,
  date_of_the_nowcast,
  max_lag,
  fitting_method = c("mcmc", "glm"),
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  sensitivity_sc = "",
  true_value = NULL,
  save_plot = TRUE
) {
  train_data_lgt <- max(df_lambda$week)
  df_lambda <- df_lambda |>
    # We plot only the weeks, where we perform nowcasting
    filter(.data$week > train_data_lgt - max_lag + 1) |>
    # Turn the delay into a factor to allow for easier faceting
    mutate(week = factor(.data$week))

  lambda_plot <- ggplot() +
    # Plot the density of the dispersion parameter estimates
    geom_line(
      data = df_lambda,
      mapping = aes(
        x = .data$.value,
        color = .data$Distribution
      ),
      stat = "density",
      alpha = 0.6
    ) +
    scale_color_manual(values = get_model_colors()[model_names]) +
    coord_cartesian(ylim = c(0, 0.001)) +
    labs(x = expression(lambda), y = "density") +
    facet_wrap(~week, nrow = 2)

  if (fitting_method == "mcmc") {
    lambda_plot <- lambda_plot +
      labs(title = expression(Posterior~of~lambda[t]))  # nolint
  } else {
    lambda_plot <- lambda_plot +
      scale_linetype(guide = "none") +
      labs(title = expression(Asymptotic~distribution~of~lambda[t]))  # nolint
  }

  if (data_origin != "case_study") {
    df_true_value <- data.frame(
      true_value = true_value,
      week = factor(seq_along(true_value) + train_data_lgt - max_lag + 1)
    )
    lambda_plot <- lambda_plot +
      geom_vline(
        data = df_true_value,
        mapping = aes(xintercept = true_value),
        linetype = "dashed"
      )
  }

  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(
      lambda_plot,
      paste0(
        paste(
          paste0(
            "inst/figure/mean_proc_plots/",
            data_origin,
            "/mean_proc_plot"
          ),
          fitting_method,
          date_of_the_nowcast,
          sep = "_"
        ),
        sensitivity_sc
      ),
      width = 7,
      height = 5.5
    )
    ret <- NULL
  } else {
    ret <- lambda_plot
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
#' @param date_of_the_nowcast a date, when the nowcast is made to name the
#' saved file correctly
#' @param data_origin a string indicating the data generating process of
#' simulated data, or whether the data correspond to the case study. Possible
#' values are "case_study", "NegBinX", "NegBin2D" and "NegBin1D"
#' @param sensitivity_sc a string indicating the sensitivity analysis scenario
#' of the MCMC method. Empty string "" indicates the main analysis.
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
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  sensitivity_sc = "",
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
    facet_wrap(~Distribution, scales = "free", nrow = 2)
  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(
      rw_sd_scatter,
      paste0(
        paste(
          paste0(
            "inst/figure/rw_sd_plots/",
            data_origin,
            "/rw_sd_scatter_plot"
          ),
          date_of_the_nowcast,
          sep = "_"
        ),
        sensitivity_sc
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
#' @param df_lambda a data frame containing columns `Distribution`
#' (containing the name of the observation model), `.value` (the empirical
#' distribution of the delay probability estimates), `nowcast_date` (the
#' date when the nowcast is calculated) and `week` (the week number from the
#' beginning of the rolling window)
#' @param df_rw_sd a data frame containing columns `Distribution`
#' (containing the name of the observation model), `.value` (the empirical
#' distribution of the estimates of the random walk standard deviation) and
#' `nowcast_date` (the date when the nowcast is calculated)
#' @param df_total a data frame with columns `date`, `counts` and `data`
#' returned by the function \code{create_totals_data_frame()}
#' @param model_names a vector of names of the observation models, we wish to
#' plot.
#' @param date_of_the_nowcast a date, when the nowcast is made to name saved
#' files correctly
#' @param fitting_method a method used for fitting the nowcasting model, either
#' "mcmc", or "glm"
#' @param prob_prior_pars a data frame of the prior parameters of the Dirichlet
#' delay probability distribution. The data frame should have columns `delay_0`,
#' `delay_1` until the maximum delay. The number of rows should be the number of
#' models time the number of sensitivity analysis scenarios. The values for each
#' model should be identical within each scenario.
#' @param disp_prior_pars a data frame with columns `model_name`, `mean_log` and
#' `sd_log`, which contain the parameters for the prior log-normal distribution
#' of the dispersion parametr. The data frame should have 6 rows, one for each
#' model
#' @param data_origin a string indicating the data generating process of
#' simulated data, or whether the data correspond to the case study. Possible
#' values are "case_study", "NegBinX", "NegBin2D" and "NegBin1D"
#' @param prob_true_val  NULL for \code{data_origin = "case_study"}, otherwise
#' the true value of the delay probability vector used to generate the data
#' @param disp_true_val NULL for \code{data_origin = "case_study"}, otherwise
#' the true value of the dispersion parameter used to generate the data
#' @param lambda_true_val NULL for \code{data_origin = "case_study"}, otherwise
#' the vector of true values of the mean of the total counts used to generate
#' the data. The length of the vector must be \code{max_lag - 1}
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
  df_lambda,
  df_rw_sd,
  df_total,
  model_names,
  date_of_the_nowcast,
  fitting_method = c("mcmc", "glm"),
  prob_prior_pars = NULL,
  disp_prior_pars = NULL,
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  prob_true_val = NULL,
  disp_true_val = NULL,
  lambda_true_val = NULL,
  save_plot = TRUE
) {
  data_origin <- match.arg(data_origin)
  fitting_method <- match.arg(fitting_method)

  scenario <- unique(df_nowcast$sensitivity_sc)
  # If we fitted the model using the GLM method, there are no scenarios of the
  # sensitivity analysis. To make data frame filtering based on these scenarios
  # work, we add the "" string to the results, which signifies the main
  # analysis.
  if (fitting_method == "glm") {
    df_disp_par <- df_disp_par |> mutate(sensitivity_sc = "")
    df_delay_prob <- df_delay_prob |> mutate(sensitivity_sc = "")
    df_lambda <- df_lambda |> mutate(sensitivity_sc = "")
  }

  # Loop over the sensitivity analysis scenarios. The data frame is always
  # filtered to contain only values from the corresponding scenario. For the GLM
  # method, we perform no sensitivity analysis and there will be only a single
  # loop to be executed.
  ret_list <- vector("list", length(scenario))
  for (k in seq_along(scenario)) {
    p_nowcast <- plot_nowcast(
      filter(df_nowcast, .data$sensitivity_sc == scenario[k]),
      df_total,
      model_names,
      date_of_the_nowcast,
      fitting_method,
      data_origin,
      scenario[k],
      save_plot
    )
    p_disp <- plot_disp_par(
      filter(df_disp_par, .data$sensitivity_sc == scenario[k]),
      model_names,
      date_of_the_nowcast,
      fitting_method,
      filter(disp_prior_pars, .data$scenario_name == scenario[k]),
      data_origin,
      scenario[k],
      disp_true_val,
      save_plot
    )
    # If fitting was done using the MCMC method, create the data frame encoding
    # the prior distribution for the reporting delay.
    if (fitting_method == "mcmc") {
      # Prior parameters for the delay probability are supplied as a data frame.
      # The parameters are typically different for each sensitivity analysis
      # scenario, but they are identical for each model, so we can use `slice`
      # to extract the parameters as a vector from the first data frame row.
      prior_prob_pars <- prob_prior_pars |>
        filter(.data$scenario_name == scenario[k]) |>
        select(starts_with("delay")) |>
        slice(1) |>
        c(recursive = TRUE)
    } else {
      prior_prob_pars <- NULL
    }
    p_prob <- plot_delay_prob(
      filter(df_delay_prob, .data$sensitivity_sc == scenario[k]),
      model_names,
      date_of_the_nowcast,
      fitting_method,
      # Coerce a data frame row to a vector. `prior_prob_pars` has length
      # number of models x number of scenarios
      prior_prob_pars,
      data_origin,
      scenario[k],
      prob_true_val,
      save_plot
    )
    max_lag <- max(df_delay_prob$delay)
    p_lambda <- plot_mean_proc(
      filter(df_lambda, .data$sensitivity_sc == scenario[k]),
      model_names,
      date_of_the_nowcast,
      max_lag,
      fitting_method,
      data_origin,
      scenario[k],
      lambda_true_val,
      save_plot
    )
    ret_list[[k]] <- list(
      nowcast = p_nowcast,
      lambda = p_lambda,
      delay_prob = p_prob,
      disp = p_disp
    )
    if (fitting_method == "mcmc") {
      p_rw_sd <- plot_rw_sd(
        filter(df_rw_sd, .data$sensitivity_sc == scenario[k]),
        filter(df_disp_par, .data$sensitivity_sc == scenario[k]),
        model_names,
        date_of_the_nowcast,
        data_origin,
        scenario[k],
        save_plot
      )
      ret_list[[k]] <- c(ret_list[[k]], rw_sd = p_rw_sd)
    }
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
#' `quantile_97.5` (bounds of the prediction intervals), `date` (x-axis
#' dates), `delay` (nowcast horizon)
#' @param df_nowcast a data frame containing columns `Distribution`
#' (containing the name of the observation model), `quantile_50`
#' (the point nowcasts), `quantile_2.5`, `quantile_25`, `quantile_75`,
#' `quantile_97.5` (bounds of the prediction intervals), `dispersion` (spread
#' component of the CRPS), `underprediction` (CRPS component) and
#' `overprediction` (CRPS component), `date` (x-axis dates), `delay`
#' (nowcast horizon) and `sensitivity_sc` (the name of a sensitivity analysis
#' scenario). This data frame contains the whole period, where nowcasting has
#' been done.
#' @param model_names a vector of names of the observation models, we wish to
#' plot.
#' @param fitting_method a method used for fitting the nowcasting model, either
#' "mcmc", or "glm"
#' @param data_origin a string indicating the data generating process of
#' simulated data, or whether the data correspond to the case study. Possible
#' values are "case_study", "NegBinX", "NegBin2D" and "NegBin1D"
#' @param save_plot logical indicator, whether to save the plots using
#' \code{ggsave()}
#'
#' @return a list of ggplot objects or list of NULLs if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#'
#' @export
plot_aggregated <- function(
  df_nowcast,
  full_data,
  model_names,
  skip_dates,
  fitting_method = c("mcmc", "glm"),
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  save_plot = TRUE
) {
  data_origin <- match.arg(data_origin)
  fitting_method <- match.arg(fitting_method)

  # Loop over the sensitivity analysis scenarios. The data frame is always
  # filtered to contain only values from the corresponding scenario. For the GLM
  # method, we perform no sensitivity analysis and there will be only a single
  # loop to be executed.
  scenario <- unique(df_nowcast$sensitivity_sc)
  ret_list <- vector("list", length(scenario))
  for (k in seq_along(scenario)) {
    p_coverage <- plot_coverage(
      filter(df_nowcast, .data$sensitivity_sc == scenario[k]),
      model_names,
      fitting_method,
      data_origin,
      scenario[k],
      save_plot
    )
    p_crps_decomp <- plot_crps_decomp(
      filter(df_nowcast, .data$sensitivity_sc == scenario[k]),
      model_names,
      fitting_method,
      data_origin,
      scenario[k],
      save_plot
    )
    p_nowcast_bands <- plot_nowcast_bands(
      full_data,
      filter(df_nowcast, .data$sensitivity_sc == scenario[k]),
      model_names,
      skip_dates,
      fitting_method,
      data_origin,
      scenario[k],
      save_plot
    )
    ret_list[[k]] <- list(
      coverage = p_coverage,
      crps_decomp = p_crps_decomp,
      nowcast_bands = p_nowcast_bands
    )
  }
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
#' @param aux_study_start a date (indeed in the date format), where the
#' auxiliary case study period used for determining the priors starts.
#' @param data_origin a string indicating the data generating process of
#' simulated data, or whether the data correspond to the case study. Possible
#' values are "case_study", "NegBinX", "NegBin2D" and "NegBin1D"
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
  aux_study_start,
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  save_plot = TRUE
) {
  data_origin <- match.arg(data_origin)
  # The auxiliary analysis ends exactly one week before the main analysis
  aux_study_end <- start_date - 7
  # Arrange the whole trajectory into a data frame for plotting
  totals <- full_data |>
    dplyr::select(paste0("value_", 1:max_lag - 1, "w")) |>
    create_totals_data_frame(aux_study_start) |>
    # `create_totals_data_frame()` returns a long data frame containing the
    # final and the preliminary state of the data. For plotting the whole
    # trajectory we are interested only in the final values.
    dplyr::filter(data == "Final")
  # The end point of the first estimation window and the beginning of the last
  # estimation window to be highlighted in the plot.
  # The estimation windows will be highlighted by braces drawn by
  # `ggpubr::geom_bracket()`.
  first_window_end <- start_date + (length_of_train_data - 1) * 7
  last_window_beg <- aux_study_start +
    (nrow(totals) - length_of_train_data - 1) * 7
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
    # Highlight the period used for determining the priors
    ggpubr::geom_bracket(
      xmin = aux_study_start,
      xmax = aux_study_end,
      y.position = first_window_max_cases + bracket_offset,
      label = "Data used to\ndetermine priors",
      label.size = 3
    ) +
    ylim(c(0, NA))
  if (data_origin == "case_study") {
    figure_path <- "inst/figure/SARI_trajectory"
    # If we plot the case study data, we will highlight the rolling window and
    # also divide it into the purely training data and the part, where
    # nowcasting is being done
    trajectory_plot <- trajectory_plot +
      # Highlight the first window of training data excluding the nowcasting
      # part
      ggpubr::geom_bracket(
        xmin = start_date,
        xmax = first_window_end - (max_lag - 2) * 7 - 1,
        y.position = first_window_max_cases + bracket_offset,
        label = "First\ntraining\ndata",
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
        label = "Last\ntraining\ndata",
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
  } else {
    figure_path <- paste0("inst/figure/", data_origin, "_simulation_trajectory")
    # If we plot the simulation study data, we will highlight the rolling window
    # without further differentiation, which part of data is complete and what
    # the nowcasting target is
    trajectory_plot <- trajectory_plot +
      # Highlight the first window of training data including the nowcasting
      # part
      ggpubr::geom_bracket(
        xmin = start_date,
        xmax = first_window_end,
        y.position = first_window_max_cases + bracket_offset,
        label = "First\nwindow",
        label.size = 3
      ) +
      # Highlight the last window of training data including the nowcasting part
      ggpubr::geom_bracket(
        xmin = last_window_beg,
        xmax = last_window_beg + length_of_train_data * 7,
        y.position = last_window_max_cases + bracket_offset,
        label = "Last\nwindow",
        label.size = 3
      ) +
      labs(
        title = paste0(data_origin, " simulation incidence"),
        y = "Incidence"
      )
  }

  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(trajectory_plot, figure_path, width = 11, height = 7)
    ret <- NULL
  } else {
    ret <- trajectory_plot
  }
  ret
}

#' Plot the prediction bands for all horizons
#'
#' @description This function creates a patchwork picture composed of the
#' incidence trajectory with the prediction intervals as bands around the
#' observed data for all time horizons.
#'
#' @param full_data a data frame of the whole trajectory containing columns
#' `date` and columns `value_0w`, `value_1w`, etc. until `max_lag - 1`.
#' @param df_nowcast a data frame containing columns `Distribution`
#' (containing the name of the observation model), `quantile_50`
#' (the point nowcasts), `quantile_2.5`, `quantile_25`, `quantile_75`,
#' `quantile_97.5` (bounds of the prediction intervals), `date` (x-axis
#' dates) and `delay` (nowcast horizon)
#' @param start_date a date (in the date format), where the nowcasting ends.
#' The starting point will be excluded.
#' @param sensitivity_sc a string indicating the sensitivity analysis scenario
#' of the MCMC method. Empty string "" indicates the main analysis.
#'
#' @return a ggplot object, or NULL if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#' @importFrom tidyselect starts_with any_of
#' @importFrom patchwork wrap_plots
#'
#' @export
plot_nowcast_bands <- function(
  full_data,
  df_nowcast,
  model_names,
  skip_dates = NULL,
  fitting_method = c("mcmc", "glm"),
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  sensitivity_sc = "",
  save_plot = TRUE
) {
  data_origin <- match.arg(data_origin)
  fitting_method <- match.arg(fitting_method)

  # End points of the trajectory is recovered from the data frame containing the
  # nowcasting results. In case we begin or end with skipped dates (Christmas),
  # we might need to adjust this.
  start_date <- min(df_nowcast$date)
  # In case we display a long trajectory like in the simulation study,
  # we will show only the first 150 weeks
  end_date <- min(max(df_nowcast$date), start_date + 150 * 7)

  full_filtered <- full_data |>
    filter(.data$date >= start_date & .data$date < end_date) |>
    select(starts_with("value_"))
  true_data <- full_filtered |> as.matrix() |> rowSums()
  # Loop over the nowcasting horizons
  horizons <- unique(df_nowcast$delay)
  patches <- vector("list", length(horizons))
  for (k in seq_along(horizons)) {
    # Which columns of the full data to sum
    value_colnames <- paste0(
      "value_",
      # The nowcsting horizon is a negative number, but the delay starts from
      # zero to the maximum delay, which we need to account for
      seq_len(length(horizons) - k + 1) - 1,
      "w"
    )
    prelim_data <- full_filtered |>
      select(any_of(value_colnames)) |>
      as.matrix() |>
      rowSums()
    df_nowcast_filtered <- df_nowcast |>
      filter(.data$delay == horizons[k] & .data$date < end_date)
    patches[[k]] <- plot_nowcast_bands_per_horizon(
      true_data,
      prelim_data,
      df_nowcast_filtered,
      start_date,
      end_date,
      skip_dates,
      horizons[k]
    )
  }
  # Arrange the patches
  arranged <- patchwork::wrap_plots(
    patches,
    nrow = length(horizons),
    ncol = 1,
    guides = "collect",
    axes = "collect"
  )
  if (save_plot) {
    plot_height <- if (fitting_method == "mcmc") {
      17.5
    } else {
      13
    }
    save_figure(
      arranged,
      paste0(
        paste(
          "inst/figure/nowcast_bands",

          data_origin,
          fitting_method,
          sep = "_"
        ),
        sensitivity_sc
      ),
      width = 10,
      height = plot_height
    )
    ret <- NULL
  } else {
    ret <- arranged
  }
  ret
}

#' Plot the prediction bands around the trajectory per horizon
#'
#' @description This function plots the incidence trajectory with the prediction
#' intervals as bands around the observed data. The plot is created for one
#' specific nowcasting horizon.
#'
#' @param true_data a vector of the final state of the incidence time series
#' @param prelim_data a vector of the preliminary state of the incidence time
#' series, which is the sum of partial counts until \code{delay}
#' @param start_date a date (in the date format), where the nowcasting starts.
#' The starting point will be included.
#' @param start_date a date (in the date format), where the nowcasting ends.
#' The starting point will be included.
#' @param horizon integer, the data up to this reporting delay are included in
#' @param model_names a vector of names of the observation models, we wish to
#' plot.
#'
#' @return a ggplot object
#'
#' @import dplyr ggplot2
#' @importFrom tidyr pivot_longer
#'
#' @export
plot_nowcast_bands_per_horizon <- function(
  true_data,
  prelim_data,
  df_nowcast,
  start_date,
  end_date,
  skip_dates,
  horizon,
  model_names
) {
  if (length(true_data) != length(prelim_data)) {
    stop("The length of preliminary data 'prelim_data' to plot must be the same as the length of 'true_data'.")  # nolint
  }

  # Arrange the whole trajectory into a data frame for plotting
  totals <- data.frame(
    date = start_date + (seq_along(true_data) - 1) * 7,
    true_data = true_data,
    prelim_data = prelim_data
  ) |>
    pivot_longer(
      c("true_data", "prelim_data"),
      values_to = "counts",
      names_to = "type"
    ) |>
    mutate(
      type = factor(
        .data$type,
        levels = c("true_data", "prelim_data"),
        labels = c("Final", "Preliminary")
      )
    )

  # In case there are dates, where we skip nowcasting, due to the reporting
  # irregularities around Christmas, we would like to break the nowcast band
  # into segments. For this reason, we split the data frame with the nowcasts
  # into several parts, which will be plotted separately.
  lower_bound <- c(start_date, skip_dates)
  upper_bound <- c(skip_dates, end_date)
  df_nowcast_splitted <- sapply(
    seq_along(lower_bound),
    function(ind) {
      filter(
        df_nowcast,
        df_nowcast$nowcast_date >= lower_bound[ind] &
          df_nowcast$nowcast_date <= upper_bound[ind]
      )
    },
    simplify = FALSE
  )
  splitted_lengths <- lapply(df_nowcast_splitted, nrow) |> unlist()
  df_nowcast_splitted <- df_nowcast_splitted[splitted_lengths > 0]

  # Set the x-axis breaks
  date_breaks <- seq(
    start_date,
    end_date,
    length = 3
  )

  nowcast_band_plot <- ggplot() +
    # True and preliminary data as black and gray solid lines
    geom_line(
      totals,
      mapping = aes(
        x = .data$date,
        y = .data$counts,
        color = .data$type,
        linetype = .data$type
      ),
      linewidth = 0.15
    )

  # Plot the nowcast segments
  for (k in seq_along(df_nowcast_splitted)) {
    nowcast_band_plot <- nowcast_band_plot +
      # Nowcast as a colored, dashed line
      geom_line(
        df_nowcast_splitted[[k]],
        mapping = aes(
          x = .data$date,
          y = .data$quantile_50,
          color = .data$Distribution,
          linetype = "Nowcast"
        ),
        linewidth = 0.15
      ) +
      # 95 % prediction interval
      geom_ribbon(
        df_nowcast_splitted[[k]],
        mapping = aes(
          x = .data$date,
          ymin = .data$quantile_2.5,
          ymax = .data$quantile_97.5,
          fill = .data$Distribution,
          alpha = "PI_95"
        )
      ) +
      # 50 % prediction interval
      geom_ribbon(
        df_nowcast_splitted[[k]],
        mapping = aes(
          x = .data$date,
          ymin = .data$quantile_25,
          ymax = .data$quantile_75,
          fill = .data$Distribution,
          alpha = "PI_50"
        )
      )
  }
  # Finish the plot
  nowcast_band_plot <- nowcast_band_plot +
    scale_x_date(breaks = date_breaks) +
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
    scale_fill_manual(values = get_model_colors()[model_names]) +
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
    labs(
      x = "Date",
      y = "Incidence",
      title = paste("Horizon:", horizon, "weeks", sep = " ")
    ) +
    get_plot_theme() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    facet_wrap(~Distribution, nrow = 2)
  nowcast_band_plot
}

#' Plot the summary of MCMC diagnostics
#'
#' @description This function plots and possibly saves the summary of the
#' MCMC fitting diagnostics for each time step.
#'
#' @param df_diagnostics a data frame containing the diagnostic summaries for
#' each model and each run (timesteps). It contains columns `num_divergent`,
#' `num_max_treedepth`, `ebfmi`, `Distribution`, `nowcast_date`.
#' @param model_names a vector of names of the observation models, we wish to
#' plot.
#' @param data_origin a string indicating the data generating process of
#' simulated data, or whether the data correspond to the case study. Possible
#' values are "case_study", "NegBinX", "NegBin2D" and "NegBin1D"
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
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  save_plot = TRUE
) {
  data_origin <- match.arg(data_origin)
  df_diagnostics_long <- df_diagnostics |>
    group_by(.data$sensitivity_sc, .data$Distribution, .data$nowcast_date) |>
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
    )

  # Set the x-axis breaks
  date_breaks <- seq(
    min(df_diagnostics$nowcast_date),
    max(df_diagnostics$nowcast_date),
    length = 13
  )

  # Loop over the sensitivity analysis scenarios. The data frame is always
  # filtered to contain only values from the corresponding scenario. For the GLM
  # method, we perform no sensitivity analysis and there will be only a single
  # loop to be executed.
  scenario <- unique(df_diagnostics$sensitivity_sc)
  ret <- vector("list", length(scenario))
  for (k in seq_along(scenario)) {
    diag_plot <- ggplot(
      filter(df_diagnostics_long, .data$sensitivity_sc == scenario[k]),
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

    # Save the plot if required, the width, height and path are hard-coded here
    if (save_plot) {
      save_figure(
        diag_plot,
        paste0(
          paste("inst/figure/diagnostics_plot", data_origin, sep = "_"),
          scenario[k]
        ),
        width = 9,
        height = 7
      )
      ret <- NULL
    } else {
      ret[[k]] <- diag_plot
    }
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
