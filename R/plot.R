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
#' dates). If we want to plot models fitted by both, the MCMC andthe  GLM method
#' \cofe{fitting_method = "both"}, the data frame must contain an additional
#' column `model_method_interact`.
#' @param df_total a data frame with columns `date`, `counts` and `data`
#' returned by the function \code{create_totals_data_frame()}
#' @param model_names a vector of names of the observation models, we wish to
#' plot.
#' @param date_of_the_nowcast a date, when the nowcast is made to name the
#' saved file correctly
#' @param fitting_method a method used for fitting the nowcasting model: "mcmc",
#' "glm", or "both", when we want to plot both methods next to each other
#' @param data_origin a string indicating the data generating process of
#' simulated data, or whether the data correspond to the case study. Possible
#' values are "case_study", "NegBinX", "NegBin2D" and "NegBin1D"
#' @param sensitivity_sc a string indicating the sensitivity analysis scenario
#' of the MCMC method. Empty string "" indicates the main analysis.
#' @param axis_limits a list with 2 elements named `x` and `y`. Each element
#' is a vector of length 2 indicating the lower and upper limit of the
#' corresponding plot axis. `x` should be a date, or a string convertible to
#' date. `y` should be numeric. For no automatic axis limits, supply a
#' \code{c(NA, NA)} vector.
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
  fitting_method = c("mcmc", "glm", "both"),
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  sensitivity_sc = "",
  axis_limits = list(x = c(NA, NA), y = c(NA, NA)),
  save_plot = TRUE
) {
  fitting_method <- match.arg(fitting_method)
  if (fitting_method == "both") {
    model_colors <- get_interaction_colors()
    color_variable <- "model_method_interact"
  } else {
    model_colors <- get_model_colors()[model_names]
    color_variable <- "Distribution"
  }

  df_total <- df_total |>
    mutate(data = factor(data, levels = c("Preliminary", "Final")))
  # Plot the nowcasts
  nowcasts_plot <- ggplot() +
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
    # Point prediction
    geom_line(
      data = df_summarized_nowcast,
      mapping = aes(
        x = .data$date,
        y = .data$quantile_50,
        color = .data[[color_variable]],
        linetype = "Nowcast"
      )
    ) +
    # 95% prediction intervals
    geom_ribbon(
      data = df_summarized_nowcast,
      mapping = aes(
        x = .data$date,
        ymin = .data$quantile_2.5,
        ymax = .data$quantile_97.5,
        fill = .data[[color_variable]],
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
        fill = .data[[color_variable]],
        alpha = "PI_50"
      )
    ) +
    # Set the color of the models and the data versions
    scale_color_manual(
      values = c(model_colors, "Final" = "black", "Preliminary" = "gray60"),
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
    # Set the transparency of the prediction intervals. The step between the
    # transparency values is rather moderate, since the 50% interval is inside
    # the 95% one. Hence, the transparency adds up.
    scale_alpha_manual(
      values = c("PI_50" = 0.5, "PI_95" = 0.2),
      labels = c("50%", "95%"),
      name = "Prediction interval"
    ) +
    guides(
      # Customize the legend of the data versions
      linetype = guide_legend(
        override.aes = list(color = c("gray60", "black", "black"))
      ),
      # Customize the legend of the confidence interval transparency, which for
      # the reader appear to be 0.6 for the 50% interval due to the "stacking"
      # of the layers.
      alpha = guide_legend(override.aes = list(alpha = c(0.7, 0.2)))
    ) +
    scale_fill_manual(values = model_colors, name = "Model") +
    labs(
      x = "Date",
      y = "Incidence",
      title = paste0(
        "Nowcasts on ",
        format(as.Date(date_of_the_nowcast), "%d %b %Y")
      )
    ) +
    xlim(as.Date(axis_limits$x)) +
    ylim(as.numeric(axis_limits$y)) +
    get_plot_theme() +
    # The final plot will have 3 facets (GLM), 6 facets (MCMC), or 9 facets
    # (both methods together). Dividing the length of the color vector by 3
    # will result in an integer number of rows.
    facet_wrap(~.data[[color_variable]], nrow = length(model_colors) / 3)
  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    # For the MCMC, we plot 6 models in 2 rows, for GLM only 3 models in a
    # single row, for both methods we plot 9 models in 3 rows.
    plot_height <- if (fitting_method == "mcmc") {
      7
    } else if (fitting_method == "glm") {
      # The height for a single line plot is higher to fit all legend items
      4.5
    } else {
      10.5
    }
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
      # Shorten the height based on the number of lines we have
      height = plot_height
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
#' @importFrom ggpubr get_legend
#' @importFrom patchwork wrap_elements
#'
#' @export
plot_coverage <- function(
  df_summarized_nowcast,
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  sensitivity_sc = "",
  save_plot = TRUE
) {
  data_origin <- match.arg(data_origin)

  # Highlight the model aligned with the true data generating process if known.
  model_y_labels <- get_y_axis_model_labels(data_origin)

  # Calculate the empirical coverage
  df_coverage <- df_summarized_nowcast |>
    group_by(.data$delay, .data$Distribution, .data$method) |>
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
    ) |>
    mutate(
      model_method_interact = factor(
        interaction(.data$Distribution, .data$method),
        levels = names(get_interaction_names()),
        labels = get_interaction_names()
      ),
      # Reverse the levels of the delay factor to show the nowcasting horizon 0
      # on the left and the past horizons more to the right
      delay = factor(.data$delay, levels = rev(levels(.data$delay)))
    )

  # Grab the maximum delay in order to label the plot facets according to the
  # corresponding delay
  max_lag <- length(unique(df_coverage$delay))
  # Set the number of rows for `facet_wrap`. For the case study, we always have
  # a 2x2 grid. For the simulation study, we can decide, whether we want to plot
  # the coverage for each data generating process individually, or patchwork
  # them into a single plot. In the latter case, the individual plots must have
  # only one row to fit.
  n_rows <- if (data_origin == "case_study") {
    2
  } else {
    1
  }
  # Grab the number of models in order to set the height of the plot
  # accordingly. We plot only the 6 MCMC-based models for the scenarios of the
  # sensitivity analysis. Otherwise we plot 9 models (6  MCMC, and 3 GLM).
  n_models <- length(unique(df_coverage$model_method_interact))
  if (!(n_models %in% c(6, 9))) {
    stop("Number of models must be 6 or 9.")
  }

  # Plot the empirical coverage as horizontal bars
  coverage_plot <- ggplot(
    df_coverage,
    aes(
      x = .data$empirical_coverage,
      y = .data$model_method_interact,
      fill = .data$model_method_interact,
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
    scale_fill_manual(values = get_interaction_colors(), name = "Model") +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.25),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    scale_y_discrete(labels = model_y_labels) +
    labs(x = "Empirical coverage", y = "Model") +
    # By default, the colors in the legend show up in the reverse order
    # compared to the barplot
    guides(fill = guide_legend(reverse = TRUE)) +
    get_plot_theme() +
    theme(legend.background = element_blank()) +
    facet_wrap(
      ~delay,
      nrow = n_rows,
      labeller = as_labeller(label_horizon_facet(max_lag))
    )

  # If we intend to plot results from a simulation, we display all simulations
  # together. For this reason we have to add the name of the data generating
  # process and place the legend to the bottom in a way it fits.
  if (data_origin != "case_study") {
    coverage_plot <- coverage_plot +
      labs(title = get_dgp_title(data_origin)) +
      guides(alpha = guide_legend(nrow = 2), fill = guide_legend(nrow = 3)) +
      theme(legend.position = "bottom")
  }
  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    if (data_origin != "case_study") {
      warning("Trying to save the coverage plot for a simulation scenario. The plot size might not be optimal.")  # nolint
    }
    save_figure(
      coverage_plot,
      paste0("inst/figure/coverage_plot_", data_origin, sensitivity_sc),
      width = 10,
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
#' @description This function plots and possibly saves the average CRPS
#' decomposed according to the spread, underprediction and overprediction. Also
#' the mean absolute error is displayed in this figure.
#'
#' @param df_summarized_nowcast a data frame containing columns `Distribution`
#' (containing the name of the observation model), `dispersion`,
#' `underprediction`, `overprediction`, `delay` (the nowcasting horizon),
#' `quantile_50` and `true_val` (to calculate the mean absolute error).
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
#' @importFrom ggpubr get_legend
#' @importFrom patchwork wrap_elements
#'
#' @export
plot_crps_decomp <- function(
  df_summarized_nowcast,
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  sensitivity_sc = "",
  save_plot = TRUE
) {
  data_origin <- match.arg(data_origin)

  # Highlight the model aligned with the true data generating process if known.
  model_y_labels <- get_y_axis_model_labels(data_origin)

  # Calculate the decomposition of the average CRPS
  df_crps <- df_summarized_nowcast |>
    # Calculate the absolute error
    mutate(AE = abs(.data$true_val - .data$quantile_50)) |>
    group_by(.data$delay, .data$Distribution, .data$method) |>
    summarize(
      # Mean absolute error
      MAE = mean(.data$AE),
      # CRPS components
      Dispersion = mean(.data$dispersion),
      Underprediction = mean(.data$underprediction),
      Overprediction = mean(.data$overprediction),
      Total = mean(.data$crps),
      .groups = "drop"
    ) |>
    # Calculate the x-coordinate of the labels denoting the total CRPS
    mutate(
      lab_position = ifelse(
        .data$delay == "0",
        # For horizon 0, the bar is wide enough to place the label inside it
        max(.data$Total) / 35,
        # For other horizons, we place the label outside of the bar, after the
        # dot denoting the MAE
        .data$MAE * 1.3
      )
    ) |>
    mutate(
      model_method_interact = factor(
        interaction(.data$Distribution, .data$method),
        levels = names(get_interaction_names()),
        labels = get_interaction_names()
      ),
      # Reverse the levels of the delay factor to show the nowcasting horizon 0
      # on the left and the past horizons more to the right
      delay = factor(.data$delay, levels = rev(levels(.data$delay)))
    )

  # Grab the maximum delay in order to label the facets according to the
  # corresponding delay
  max_lag <- length(unique(df_crps$delay))
  # Set the number of rows for `facet_wrap`. For the case study, we always have
  # a 2x2 grid. For the simulation study, we can decide, whether we want to plot
  # the coverage for each data generating process individually, or patchwork
  # them into a single plot. In the latter case, the individual plots must have
  # only one row to fit.
  n_rows <- if (data_origin == "case_study") {
    2
  } else {
    1
  }
  # Grab the number of models in order to set the height of the plot
  # accordingly. We plot only the 6 MCMC-based models for the scenarios of the
  # sensitivity analysis. Otherwise we plot 9 models (6  MCMC, and 3 GLM).
  n_models <- length(unique(df_crps$model_method_interact))
  if (!(n_models %in% c(6, 9))) {
    stop("Number of models must be 6 or 9.")
  }

  df_crps_decomp <- df_crps |>
    # Pivot for easier definition of the alpha aesthetic
    tidyr::pivot_longer(
      cols = c("Dispersion", "Overprediction", "Underprediction"),
      names_to = "Component",
      values_to = "CRPS"
    ) |>
    mutate(
      Component = factor(
        .data$Component,
        levels = c("Overprediction", "Dispersion", "Underprediction")
      )
    )

  # Plot the CRPS as horizontal bars
  crps_decomp_plot <- ggplot() +
    geom_col(
      df_crps_decomp,
      mapping = aes(
        x = .data$CRPS,
        y = .data$model_method_interact,
        fill = .data$model_method_interact,
        alpha = .data$Component
      ),
      position = "stack"
    ) +
    geom_label(
      df_crps,
      mapping = aes(
        x = .data$lab_position,
        y = .data$model_method_interact,
        label = round(.data$Total, 2)
      ),
      border.color = "black",
      text.color = "black",
      color = "white",
      hjust = 0
    ) +
    scale_alpha_manual(
      values = c(
        "Underprediction" = 1,
        "Dispersion" = 0.4,
        "Overprediction" = 0.7
      ),
      name = "Component"
    ) +
    geom_point(
      df_crps_decomp,
      mapping = aes(x = .data$MAE, y = .data$model_method_interact)
    ) +
    scale_fill_manual(values = get_interaction_colors(), name = "Model") +
    scale_y_discrete(labels = model_y_labels) +
    labs(x = "Mean CRPS/AE", y = "Model") +
    # By default, the colors in the legend show up in the reverse order
    # compared to the barplot
    guides(fill = guide_legend(reverse = TRUE)) +
    get_plot_theme() +
    theme(legend.background = element_blank()) +
    facet_wrap(
      ~delay,
      nrow = n_rows,
      labeller = as_labeller(label_horizon_facet(max_lag))
    )

  # If we intend to plot results from a simulation, we display all simulations
  # together. For this reason we have to add the name of the data generating
  # process and place the legend to the bottom in a way it fits.
  if (data_origin != "case_study") {
    crps_decomp_plot <- crps_decomp_plot +
      labs(title = get_dgp_title(data_origin)) +
      guides(alpha = guide_legend(nrow = 3), fill = guide_legend(nrow = 3)) +
      theme(legend.position = "bottom")
  }
  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    if (data_origin != "case_study") {
      warning("Trying to save the CRPS plot for a simulation scenario. The plot size might not be optimal.")  # nolint
    }
    save_figure(
      crps_decomp_plot,
      paste0(
        "inst/figure/crps_decomposition_plot_",
        data_origin,
        sensitivity_sc
      ),
      width = 10,
      height = 7
    )
    ret <- NULL
  } else {
    ret <- crps_decomp_plot
  }
  ret
}

patchwork_sim_results <- function(summarized_nowcast_list, save_plot = TRUE) {
  data_origin <- names(summarized_nowcast_list)
  plot_list <- replicate(
    2,
    setNames(vector("list", length(data_origin)), data_origin),
    simplify = FALSE
  )
  names(plot_list) <- c("coverage", "crps")
  # Loop over the data generating processes
  for (dgp in data_origin) {
    plot_list$coverage[[dgp]] <- plot_coverage(
      summarized_nowcast_list[[dgp]],
      data_origin = dgp,
      sensitivity_sc = "",
      save_plot = FALSE
    )
    plot_list$crps[[dgp]] <- plot_crps_decomp(
      summarized_nowcast_list[[dgp]],
      data_origin = dgp,
      sensitivity_sc = "",
      save_plot = FALSE
    )
  }
  coverage_plot <- wrap_plots(plot_list$coverage, nrow = 3) +
    plot_layout(axes = "collect", guides = "collect") &
    theme(legend.position = "bottom", legend.justification = "right")
  crps_plot <- wrap_plots(plot_list$crps, nrow = 3) +
    plot_layout(axes = "collect", guides = "collect") &
    theme(legend.position = "bottom", legend.justification = "right")
  if (save_plot) {
    save_figure(coverage_plot, "inst/figure/coverage_plot_simulation", 9, 10)
    save_figure(
      crps_plot,
      "inst/figure/crps_decomposition_plot_simulation",
      width = 10,
      height = 14
    )
    ret <- vector("list", 2)
  } else {
    ret <- list(coverage = coverage_plot, crps = crps_plot)
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
#' @param example logical indicator, whether we want to point out the plot as
#' an example of different posteriors between the NegBinX and NegBin2D/1D
#' models. If \code{example = TRUE}, than we plot/save a version of the plot
#' more focused on the selected models.
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
  example = FALSE,
  save_plot = TRUE
) {
  max_lag <- max(df_delay_prob$delay)

  df_delay_prob <- df_delay_prob |>
    # Turn the delay into a factor to allow for easier faceting
    mutate(delay = factor(.data$delay))

  facet_titles <- rlang::set_names(
    paste0("Delay: ", seq_len(max_lag), " weeks"),
    seq_len(max_lag)
  )

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
    facet_wrap(
      ~delay,
      scales = "free_x",
      nrow = 2,
      labeller = as_labeller(facet_titles)
    )

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
        values = c("Prior" = "dotted", "Posterior" = "solid"),
        name = ""
      ) +
      labs(
        title = "Posterior of the delay probability",
        x = "Delay probability",
        y = "Density"
      )
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

  # If we want to save the plot and show it in the manuscript, we show only
  # models we are interested in and adjust the plot size in order to emphasize
  # the comparison between the selected models.
  if (example) {
    # Show only selected models
    delay_prob_plot$layers[[1]]$data <- filter(
      delay_prob_plot$layers[[1]]$data,
      .data$Distribution %in% model_names
    )
    # Add a new faceting specification to make all scales varying
    delay_prob_plot <- delay_prob_plot +
      facet_wrap(
        ~delay,
        scales = "free",
        nrow = 2,
        labeller = as_labeller(facet_titles)
      )
    # Remove the y-limit
    delay_prob_plot$coordinates$limits$y <- c(0, NA)
    # Remove the title
    delay_prob_plot$labels$title <- NULL
    if (save_plot) {
      save_figure(
        delay_prob_plot,
        paste(
          "inst/figure/delay_posterior_difference",
          date_of_the_nowcast,
          sep = "_"
        ),
        width = 7,
        height = 4
      )
    }
  }

  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot && !example) {
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
      list(x = c(NA, NA), y = c(NA, NA)),
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
      # Within this wrapper, we don't plot the more focused example of the delay
      # probability posterior, thus example = FALSE
      FALSE,
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
#' results: \code{plot_coverage()}, \code{plot_crps_decomp()} and
#' \code{plot_nowcast_bands()}. For the MCMC method, we run 4 additional
#' scenarios as a robustness check. For these scenarios, the coverage and CRPS
#' plots are saved as glued together by pairs. For the simulations, we don't
#' create the coverage and CRPS. This is done in the
#' \code{patchwork_sim_results} function, where we merge together results from
#' all simulation study scenarios.
#'
#' @param df_nowcast a data frame containing columns `Distribution`
#' (containing the name of the observation model), `quantile_50`
#' (the point nowcasts), `quantile_2.5`, `quantile_25`, `quantile_75`,
#' `quantile_97.5` (bounds of the prediction intervals), `dispersion` (spread
#' component of the CRPS), `underprediction` (CRPS component) and
#' `overprediction` (CRPS component), `date` (x-axis dates), `delay`
#' (nowcast horizon) and `sensitivity_sc` (the name of a sensitivity analysis
#' scenario). This data frame contains the whole period, where nowcasting has
#' been done.
#' @param full_data a data frame with columns `date` and columns
#' `value_0w`, `value_1w`, etc. until `max_lag - 1` used to plot the preliminary
#' and final data alongside the nowcasts
#' @param skip_dates a vector of dates, where no nowcasting has been done and
#' where we should leave gaps in the plot of the nowcasts.
#' \code{skip_dates = NULL} if no gaps are to be plotted.
#' @param data_origin a string indicating the data generating process of
#' simulated data, or whether the data correspond to the case study. Possible
#' values are "case_study", "NegBinX", "NegBin2D" and "NegBin1D"
#' @param save_plot logical indicator, whether to save the plots using
#' \code{ggsave()}
#'
#' @return a nested list of ggplot objects or NULLs if \code{save_plot = TRUE}.
#' There is one list element per sensitivity analysis scenario, each list
#' element contains 3 more elements (reserved for the coverage plot, CRPS plot
#' and prediction interval bands). For the simulation study,
#' (i.e. \code{data_origin} other than "case_study"), there is always NULL for
#' the CRPS and coverage plot.
#'
#' @import dplyr ggplot2
#' @importFrom patchwork plot_layout plot_spacer wrap_elements
#' @importFrom ggpubr get_legend
#'
#' @export
plot_aggregated <- function(
  df_nowcast,
  full_data,
  skip_dates,
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  save_plot = TRUE
) {
  data_origin <- match.arg(data_origin)

  # Loop over the sensitivity analysis scenarios. The data frame is always
  # filtered to contain only values from the corresponding scenario. For the GLM
  # method, we perform no sensitivity analysis and there will be only a single
  # loop to be executed.
  scenario <- unique(df_nowcast$sensitivity_sc)
  ret_list <- vector("list", length(scenario))
  names(ret_list) <- scenario
  for (k in seq_along(scenario)) {
    # Logical indicating whether we are plotting results from the main scenario.
    # For the sensitivity analysis scenario, we do not save the coverage and
    # CRPS plots individually. Rather, we glue them together using patchwork to
    # fit better on page of the manuscript.
    main_scenario <- scenario[k] == ""
    if (data_origin == "case_study") {
      p_coverage <- plot_coverage(
        filter(df_nowcast, .data$sensitivity_sc == scenario[k]),
        data_origin,
        scenario[k],
        # Avoid saving the individual plot for other scenarios than the main one
        save_plot && main_scenario
      )
      p_crps_decomp <- plot_crps_decomp(
        filter(df_nowcast, .data$sensitivity_sc == scenario[k]),
        data_origin,
        scenario[k],
        # Avoid saving the individual plot for other scenarios than the main one
        save_plot && main_scenario
      )
    } else {
      p_coverage <- NULL
      p_crps_decomp <- NULL
    }
    p_nowcast_bands <- plot_nowcast_bands(
      full_data,
      filter(df_nowcast, .data$sensitivity_sc == scenario[k]),
      skip_dates,
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
  # If the plots are to be saved, glue together the coverage and CRPS plots from
  # the sensitivity analysis.
  if (save_plot && data_origin == "case_study") {
    save_patchwork_plots(ret_list)
    # If we save a plot, we usually return NULL in place of the individual
    # plots. For consistency, we reconstruct the list of NULLs with a
    # corresponding structure here.
    ret_list <- replicate(length(scenario), vector("list", 3), simplify = FALSE)
    names(ret_list) <- scenario
  }
  ret_list
}

#' Arrange and save figures from the robustness check
#'
#' @description This functions takes the individual coverage and CRPS plots of
#' the results from the sensitivity analysis, groups them together and saves
#' them. For the MCMC method, we run 4 additional
#' scenarios as a robustness check. For these scenarios, the coverage and CRPS
#' plots are saved as glued together by pairs.
#' \itemize{
#'   \item coverage in scenarios with stronger and weaker prior on the delay
#'   probability,
#'   \item CRPS in scenarios with stronger and weaker prior on the delay
#'   probability,
#'   \item coverage in scenarios with stronger and weaker prior on the
#'   dispersion parameter,
#'   \item CRPS in scenarios with stronger and weaker prior on the
#'   dispersion parameter.
#' }
#'
#' @param plot_list a list containing the individual plots from the 4
#' sensitivity analysis scenarios. The outer list is indexed by the scenarios,
#' the inner list by plot type. The names of the sensitivity screnarios are:
#' `prob_high`, `prob_low`, `disp_high`, `disp_low`. The names of the relevant
#' individual plots are `coverage` and `crps_decomp`.
#'
#' @return NULL
#'
#' @import ggplot2
#' @importFrom patchwork plot_layout plot_spacer
#'
#' @export
save_patchwork_plots <- function(plot_list) {
  # Add a plot title to distinguish between more and less informative priors.
  # The same title is used for the delay probability and the dispersion
  # parameter.
  p_theme_chunk_stronger <- list(
    labs(title = "Stronger prior"),
    theme(plot.title = element_text(hjust = 0.5, size = 16))
  )
  p_theme_chunk_weaker <- list(
    labs(title = "Weaker prior"),
    theme(plot.title = element_text(hjust = 0.5, size = 16))
  )
  # Layout of the patchwork plot. The plots will be placed next to each other,
  # so we add a narrow spacer between them.
  p_layout <- patchwork::plot_layout(
    guides = "collect",
    axes = "collect_y",
    widths = c(7.4, 0.2, 7.4)
  )
  # List of the plots to iterate over
  patchworked_list <- vector("list", 4)
  # Names that will be used as file names
  names(patchworked_list) <- c(
    "coverage_plot_case_study_prob",
    "coverage_plot_case_study_disp",
    "crps_decomposition_plot_case_study_prob",
    "crps_decomposition_plot_case_study_disp"
  )
  # Coverage plot for scenarios modifying the dispersion of the delay
  # probability prior
  patchworked_list[[1]] <- (
    (plot_list$prob_high$coverage + p_theme_chunk_weaker) |
      patchwork::plot_spacer() |
      (plot_list$prob_low$coverage + p_theme_chunk_stronger)
  ) + p_layout
  # Coverage plot for scenarios modifying the dispersion of the dispersion
  # parameter prior
  patchworked_list[[2]] <- (
    (plot_list$disp_high$coverage + p_theme_chunk_weaker) |
      patchwork::plot_spacer() |
      (plot_list$disp_low$coverage + p_theme_chunk_stronger)
  ) + p_layout
  # CRPS plot for scenarios modifying the dispersion of the delay probability
  # prior
  patchworked_list[[3]] <- (
    (plot_list$prob_high$crps_decomp + p_theme_chunk_weaker) |
      patchwork::plot_spacer() |
      (plot_list$prob_low$crps_decomp + p_theme_chunk_stronger)
  ) + p_layout
  # CRPS plot for scenarios modifying the scale of the dispersion parameter
  # prior
  patchworked_list[[4]] <- (
    (plot_list$disp_high$crps_decomp + p_theme_chunk_weaker) |
      patchwork::plot_spacer() |
      (plot_list$disp_low$crps_decomp + p_theme_chunk_stronger)
  ) + p_layout
  # Iterate over the plots and save them
  for (k in seq_along(patchworked_list)) {
    save_figure(
      patchworked_list[[k]],
      paste0("inst/figure/", names(patchworked_list)[k]),
      width = 13,
      height = 6
    )
  }
  NULL
}

#' Plot the whole incidence trajectory
#'
#' @description This function plots and possibly saves an incidence trajectory.
#' The first and the last estimation windows will be highlighted to see the
#' chunk of the data we use for model training. The final plot has two
#' panels - the trajectory (top) and the proportion of reports per
#' delay (bottom).
#'
#' @param full_data a data frame of the whole trajectory containing columns
#' `date` and columns `value_0w`, `value_1w`, etc. until `max_lag - 1`.
#' @param start_date a date (indeed in the date format), where the training data
#' start. The starting point will be included.
#' @param end_date a date (in the date format), where the training data end. The
#' endpoint will be included.
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
#' @return a patchwork plot, or NULL if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#' @importFrom ggpubr geom_bracket
#'
#' @export
plot_trajectory <- function(
  full_data,
  start_date,
  end_date,
  length_of_train_data,
  max_lag,
  aux_study_start,
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  save_plot = TRUE
) {
  data_origin <- match.arg(data_origin)
  # The auxiliary analysis ends exactly one week before the main analysis
  aux_study_end <- start_date - 7
  # Filter the full data to contain only the selected time period
  full_data <- full_data |> filter(date >= aux_study_start, date <= end_date)
  # Arrange the whole trajectory into a data frame for plotting
  totals <- full_data |>
    dplyr::select(paste0("value_", 1:max_lag - 1, "w")) |>
    create_totals_data_frame(aux_study_start) |>
    # `create_totals_data_frame()` returns a long data frame containing the
    # final and the preliminary state of the data. For plotting the whole
    # trajectory we are interested only in the final values.
    dplyr::filter(data == "Final")
  # The end point of the first estimation window to be highlighted in the plot.
  # The estimation windows will be highlighted by braces drawn by
  # `ggpubr::geom_bracket()`.
  first_window_end <- start_date + (length_of_train_data - 1) * 7
  # We need to find the maximum number of cases in the first and last estimation
  # window in order to place the brace correctly above them.
  first_window_max_cases <- totals |>
    filter(date <= first_window_end) |>
    pull(.data$counts) |>
    # na.rm = TRUE is usually not needed, but it prevents the plot element to
    # disappear in the case of missing values
    max(na.rm = TRUE)
  overall_max_cases <- totals |>
    pull(.data$counts) |>
    max(na.rm = TRUE)
  # 5% offset of the braces to avoid overplotting the trajectory
  bracket_offset <- first_window_max_cases * 0.1
  # For the simulation study, place the bracket indicating the first window a
  # little bit higher, since it is located near a season peak.
  if (data_origin == "case_study") {
    figure_path <- "inst/figure/SARI_trajectory"
    first_window_bracket_y <- first_window_max_cases + bracket_offset
  } else {
    figure_path <- paste0("inst/figure/", data_origin, "_simulation_trajectory")
    first_window_bracket_y <- first_window_max_cases + 7 * bracket_offset
  }

  p_trajectory <- ggplot(totals, aes(x = .data$date, y = .data$counts)) +
    geom_line() +
    # Highlight the period used for determining the priors
    ggpubr::geom_bracket(
      xmin = aux_study_start,
      xmax = aux_study_end,
      y.position = first_window_max_cases + bracket_offset,
      label = "Data to\ninform priors",
      label.size = 4
    ) +
    # Highlight the first window of training data including the nowcasting
    # part
    ggpubr::geom_bracket(
      xmin = start_date,
      xmax = first_window_end,
      y.position = first_window_bracket_y,
      label = "First\nwindow",
      label.size = 4
    ) +
    # Highlight the last window of training data including the nowcasting part
    ggpubr::geom_bracket(
      xmin = end_date - (length_of_train_data - 1) * 7,
      xmax = end_date,
      y.position = overall_max_cases + bracket_offset,
      label = "Last\nwindow",
      label.size = 4
    ) +
    labs(y = "Incidence", x = "Date") +
    ylim(
      c(0, max(overall_max_cases, first_window_bracket_y) + 4 * bracket_offset)
    ) +
    get_plot_theme()

  if (data_origin != "case_study") {
    p_trajectory <- p_trajectory + labs(title = get_dgp_title(data_origin))
  }

  # We add the plot of proportions of cases reported for each delay
  df_prop_reported <- full_data |>
    mutate(total = rowSums(across(starts_with("value_")))) |>
    pivot_longer(
      starts_with("value_"),
      names_to = "delay",
      values_to = "counts"
    ) |>
    mutate(
      # Extract the reporting delay from the former string of column names
      # and reverse the factor to plot the 0-delay on the bottom.
      delay = factor(
        gsub("[^0-9]", "", .data$delay),
        levels = rev(seq_len(max_lag) - 1)
      )
    ) |>
    # Calculate the proportions for each date
    group_by(.data$date) |>
    mutate(
      prop_reported = ifelse(.data$total > 0, .data$counts / .data$total, 0)
    ) |>
    ungroup()
  # Plot the proportions
  p_prop_reported <- ggplot(
    df_prop_reported,
    aes(x = .data$date, y = .data$prop_reported, fill = .data$delay)
  ) +
    geom_area() +
    labs(
      y = "Fraction of\nfinal reports",
      fill = "Delay in\nweeks",
      x = "Date"
    ) +
    scale_fill_viridis_d(direction = -1) +
    get_plot_theme() +
    theme(axis.title.y = element_text(hjust = 0))
  # Glue together the trajectory and the proportions
  trajectory_plot <- patchwork::wrap_plots(
    p_trajectory,
    p_prop_reported,
    nrow = 2,
    axes = "collect",
    heights = c(5, 3)
  )

  # Save the plot if required, the width, height and path are hard-coded here.
  # If the plot is saved on the disc, we don't return the ggplot object.
  if (save_plot) {
    save_figure(trajectory_plot, figure_path, width = 9, height = 5)
    ret <- NULL
  } else {
    ret <- trajectory_plot
  }
  ret
}

#' Plot the incidence trajectory for all simulation scenarios
#'
#' @description A wrapper around \code{plot_trajectory}. This function plots and
#' possibly saves the incidence trajectory used in the simulation study. The
#' first and the last estimation windows will be highlighted to see the chunk of
#' the data we use for model training. Each trajectory is accompanied by a plot
#' of proportions of cases per delay. The plots corresponding to different
#' data-generating processes are placed on top of each other.
#'
#' @param full_data a data frame of the whole trajectory containing columns
#' `date`, `Distribution` (the data-generating process) and columns `value_0w`,
#' `value_1w`, etc. until `max_lag - 1`.
#' @param start_date a date (indeed in the date format), where the training data
#' start. The starting point will be included.
#' @param end_date a date (in the date format), where the training data end. The
#' endpoint will be included.
#' @param length_of_train_data a number, the length of the estimation window
#' (endpoints included)
#' @param max_lag maximum reporting delay represented by the number of columns
#' of the reporting table. In this way, the 0-th lag counts as the first, 1-st
#' lag as the second and so on.
#' @param aux_study_start a date (indeed in the date format), where the
#' auxiliary case study period used for determining the priors starts.
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggplot2::ggsave()}
#'
#' @return a patchwork, or NULL if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#' @importFrom ggpubr geom_bracket
#' @importFrom patchwork wrap_plots
#'
#' @export
plot_all_sim_trajectories <- function(
  full_data,
  start_date,
  end_date,
  length_of_train_data,
  max_lag,
  aux_study_start,
  save_plot = TRUE
) {
  data_origin <- unique(full_data$Distribution)
  # Loop over the data generating processes
  trajectory_patches <- vector("list", length(data_origin))
  names(trajectory_patches) <- data_origin
  for (dgp in data_origin) {
    full_data_filtered <- full_data |> filter(.data$Distribution == dgp)
    trajectory_patches[[dgp]] <- plot_trajectory(
      full_data_filtered,
      start_date,
      end_date,
      length_of_train_data,
      max_lag,
      aux_study_start,
      data_origin = dgp,
      save_plot = FALSE
    )
    # In the end, we will get a patchwork plot with a nested layou, where it's
    # not possible anymore to collect the x-axis guide. Therefore, we remove it
    # manually for all plots except for the last one.
    if (dgp != tail(data_origin, 1)) {
      trajectory_patches[[dgp]] <- trajectory_patches[[dgp]] &
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
    }
  }
  # Put the patchwork plot together
  arranged <- patchwork::wrap_plots(
    trajectory_patches,
    nrow = 3,
    guides = "collect"
  )
  if (save_plot) {
    ret <- NULL
    save_figure(
      arranged,
      "inst/figure/sim_trajectory",
      width = 9,
      height = 11
    )
  } else {
    ret <- arranged
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
#' `date` and columns `value_0w`, `value_1w`, etc. until `max_lag - 1` used to
#' plot the preliminary
#' and final data alongside the nowcasts
#' @param df_nowcast a data frame containing columns `Distribution`
#' (containing the name of the observation model), `quantile_50`
#' (the point nowcasts), `quantile_2.5`, `quantile_25`, `quantile_75`,
#' `quantile_97.5` (bounds of the prediction intervals), `date` (x-axis
#' dates), `nowcast_date` (when the nowcast was issued) and `delay` (nowcast
#' horizon)
#' @param skip_dates a vector of dates, where no nowcasting has been done and
#' where we should leave gaps in the plot of the nowcasts.
#' \code{skip_dates = NULL} if no gaps are to be plotted.
#' @param data_origin a string indicating the data generating process of
#' simulated data, or whether the data correspond to the case study. Possible
#' values are "case_study", "NegBinX", "NegBin2D" and "NegBin1D"
#' @param sensitivity_sc a string indicating the sensitivity analysis scenario
#' of the MCMC method. Empty string "" indicates the main analysis.
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggplot2::ggsave()}
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
  skip_dates = NULL,
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  sensitivity_sc = "",
  save_plot = TRUE
) {
  data_origin <- match.arg(data_origin)

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
  # Loop over the nowcasting horizons sorted from 0 to -3
  horizons <- unique(df_nowcast$delay)
  horizons <- rev(horizons[order(as.numeric(as.character(horizons)))])
  patches <- vector("list", length(horizons))
  for (k in seq_along(horizons)) {
    # Which columns of the full data to sum. The delay is indexed from zero
    value_colnames <- paste0("value_", seq_len(k) - 1, "w")
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
  # For some settings, we split the patchwork plot into 2 figures. The first one
  # showing only the nowcasts for horizon 0 will go to the main manuscript, the
  # rest into the supplement.
  # We need to split:
  # - the case study fitted using the MCMC method
  # - the case study fitted using the GLM method
  # - the NegBinX simulation study fitted using the MCMC method
  # - the NegBinX simulation study fitted using the GLM method
  split_figures <- data_origin == "case_study" && sensitivity_sc == ""
  if (split_figures) {
    # Arrange all the patches
    arranged <- patchwork::wrap_plots(
      # The first patch corresponds to delay zero, which is plotted separately
      tail(patches, -1),
      nrow = length(horizons) - 1,
      ncol = 1,
      guides = "collect",
      axes = "collect"
    )
    ret <- list(
      arranged_plot = arranged,
      separated_plot = patches[[1]]
    )
  } else {
    # Arrange all the patches
    arranged <- patchwork::wrap_plots(
      patches,
      nrow = length(horizons),
      ncol = 1,
      guides = "collect",
      axes = "collect"
    )
    ret <- list(
      arranged_plot = arranged,
      separated_plot = NULL
    )
  }

  if (save_plot) {
    # For the case study, we have 4 nowcasting horizons, for the simulation
    # study only 3
    plot_height <- if (data_origin == "case_study" && sensitivity_sc == "") {
      20
    } else {
      15
    }
    plot_path <- paste0(
      "inst/figure/nowcast_bands_",
      data_origin,
      sensitivity_sc
    )
    if (split_figures) {
      save_figure(
        arranged,
        plot_path,
        width = 11.5,
        height = plot_height * (length(horizons) - 1) / length(horizons)
      )
      save_figure(
        patches[[1]],
        paste(plot_path, "delay0", sep = "_"),
        width = 11.5,
        # Increase the height to make enough space for the axis labels and
        # legend items
        height = plot_height / length(horizons) + 2
      )
    } else {
      save_figure(
        arranged,
        plot_path,
        width = 11.5,
        height = plot_height
      )
    }
    ret <- NULL
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
#' @param df_nowcast a data frame containing columns `Distribution`
#' (containing the name of the observation model), `quantile_50`
#' (the point nowcasts), `quantile_2.5`, `quantile_25`, `quantile_75`,
#' `quantile_97.5` (bounds of the prediction intervals), `date` (x-axis
#' dates) and `nowcast_date` (when the nowcast was issued)
#' @param start_date a date (in the date format), where the nowcasting starts.
#' The starting point will be included.
#' @param end_date a date (in the date format), where the nowcasting ends.
#' The ending point will be included.
#' @param skip_dates a vector of dates, where no nowcasting has been done and
#' where we should leave gaps in the plot of the nowcasts.
#' \code{skip_dates = NULL} if no gaps are to be plotted.
#' @param horizon integer, the data up to this reporting delay are included in
#' the preliminary data
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
  horizon
) {
  if (length(true_data) != length(prelim_data)) {
    stop("The length of preliminary data 'prelim_data' to plot must be the same as the length of 'true_data'.")  # nolint
  }
  # Set the factor for the different models and methods and its names
  df_nowcast <- df_nowcast |>
    mutate(
      model_method_interact = factor(
        interaction(.data$Distribution, .data$method),
        levels = names(get_interaction_names(nowcast_bands_ordering = TRUE)),
        labels = get_interaction_names(nowcast_bands_ordering = TRUE)
      )
    )

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
        levels = c("prelim_data", "true_data"),
        labels = c("Preliminary", "Final")
      )
    )
  # Grab the number of models in order to set the height of the plot
  # accordingly. We plot only the 6 MCMC-based models for the scenarios of the
  # sensitivity analysis. Otherwise we plot 9 models (6  MCMC, and 3 GLM).
  n_models <- length(unique(df_nowcast$model_method_interact))
  if (!(n_models %in% c(6, 9))) {
    stop("Number of models must be 6 or 9.")
  }

  # In case there are dates, where we skip nowcasting, due to the reporting
  # irregularities around Christmas, we would like to break the nowcast band
  # into segments. For this reason, we split the data frame with the nowcasts
  # into several parts, which will be plotted separately.
  lower_bound <- c(start_date, skip_dates)
  # In case that skip_dates = NULL, the resulting vector will be coerced to
  # numeric. For this reason, we need to keep as.Date. This is not the case
  # for the lower bound, since there we always start with a date.
  upper_bound <- as.Date(c(skip_dates, end_date))
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
    length = 4
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
      linewidth = 0.3
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
          color = .data$model_method_interact,
          linetype = "Nowcast"
        ),
        linewidth = 0.25
      ) +
      # 95 % prediction interval
      geom_ribbon(
        df_nowcast_splitted[[k]],
        mapping = aes(
          x = .data$date,
          ymin = .data$quantile_2.5,
          ymax = .data$quantile_97.5,
          fill = .data$model_method_interact,
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
          fill = .data$model_method_interact,
          alpha = "PI_50"
        )
      )
  }
  # Finish the plot
  nowcast_band_plot <- nowcast_band_plot +
    scale_x_date(breaks = date_breaks, date_labels = "%b %Y") +
    scale_color_manual(
      values = c(
        get_interaction_colors(),
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
    scale_fill_manual(values = get_interaction_colors(), name = "Model") +
    # Set the transparency of the prediction intervals. The step between the
    # transparency values is rather moderate, since the 50% interval is inside
    # the 95% one. Hence, the transparency adds up.
    scale_alpha_manual(
      values = c("PI_50" = 0.5, "PI_95" = 0.2),
      labels = c("50%", "95%"),
      name = "Prediction\ninterval"
    ) +
    guides(
      # Customize the legend of the data versions
      linetype = guide_legend(
        override.aes = list(color = c("gray60", "black", "black"))
      ),
      # Customize the legend of the confidence interval transparency, which for
      # the reader appear to be 0.4 for the 95% interval due to the "stacking"
      # of the layers.
      alpha = guide_legend(override.aes = list(alpha = c(0.7, 0.2)))
    ) +
    labs(
      x = "Date",
      y = "Incidence",
      title = paste("Horizon:", horizon, "weeks", sep = " ")
    ) +
    get_plot_theme() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    # Plot 3 models per row, since `n_models` can be either 6, or 9
    facet_wrap(~model_method_interact, nrow = n_models / 3)
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

#' Plot the nowcasts for several dates
#'
#' @description This function creates a patchwork picture composed of the
#' nowcast plots as calculated in "real time". Individual panels correspond
#' to different dates, when the nowcast is calculated.
#'
#' @param df_nowcast a list of data frames containing columns `Distribution`
#' (containing the name of the observation model), `quantile_50`
#' (the point nowcasts), `quantile_2.5`, `quantile_25`, `quantile_75`,
#' `quantile_97.5` (bounds of the prediction intervals), `date` (x-axis
#' dates), `nowcast_date` (when the nowcast was issued) and `delay` (nowcast
#' horizon). The individual data frames in the list correspond to different
#' dates, when the nowcast is calculated.
#' @param df_total a list of data frames containing columns `date`, `counts` and
#' `data`. The last column `data` is an indicator, whether the values in the
#' `counts` column are the final sums of the counts, or the preliminary data
#' version. Needed to plot the observations alongside the nowcasts. The
#' individual data frames in the list correspond to different dates, when the
#' nowcast is calculated.
#' @param dates_to_show a selection of a few consecutive dates for which we want
#' to show the estimates.
#' @param data_origin a string indicating the data generating process of
#' simulated data, or whether the data correspond to the case study. Possible
#' values are "case_study", "NegBinX", "NegBin2D" and "NegBin1D"
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggplot2::ggsave()}
#'
#' @return a patchwork plot, or NULL if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2
#' @importFrom patchwork wrap_plots
#'
#' @export
plot_nowcast_example <- function(
  df_nowcast,
  df_total,
  dates_to_show,
  data_origin = c("case_study", "NegBinX", "NegBin2D", "NegBin1D"),
  save_plot = TRUE
) {
  data_origin <- match.arg(data_origin)

  # Bind rows of all the data frames that are in a list format
  df_total <- bind_rows(df_total)
  df_nowcast <- bind_rows(df_nowcast) |>
    # Create a faceting variable, which is the interaction of the observation
    # model and the fitting method
    mutate(
      model_method_interact = factor(
        interaction(.data$Distribution, .data$method),
        levels = names(get_interaction_names(nowcast_bands_ordering = TRUE)),
        labels = get_interaction_names(nowcast_bands_ordering = TRUE)
      )
    )
  # Find the axis limits. X-axis limits are defined as the beginning of the
  # first rolling window and the end of the last rolling window. Having fixed
  # x-axis limits helps with collecting the axis via `plot_layout()`, which
  # saves some vertical space, that would otherwise be occupied by x-axis
  # labels. Y-axis labels are set for the patchwork plot to have consistent
  # y-axis scale for all subplots.
  axis_limits <- list(
    x = range(df_total$date),
    y = c(
      min(df_total$counts),
      max(c(df_nowcast$quantile_97.5, df_total$counts))
    )
  )

  # Loop over the dates, on which the nowcasts are calculated
  patches <- vector("list", length(dates_to_show))
  names(patches) <- dates_to_show
  for (k in seq_along(patches)) {
    df_nowcast_filtered <- df_nowcast |>
      filter(.data$nowcast_date == dates_to_show[k])
    df_total_filtered <- df_total |>
      filter(.data$nowcast_date == dates_to_show[k])
    patches[[k]] <- plot_nowcast(
      df_nowcast_filtered,
      df_total_filtered,
      get_interaction_names(),
      dates_to_show[k],
      fitting_method = "both",
      data_origin = data_origin,
      sensitivity_sc = "",
      axis_limits = axis_limits,
      save_plot = FALSE
    )
  }
  # Arrange all the patches
  arranged <- patchwork::wrap_plots(
    patches,
    nrow = length(dates_to_show),
    ncol = 1,
    guides = "collect",
    axes = "collect"
  )
  # Save the plot if required, the width, height and path are hard-coded here
  if (save_plot) {
    save_figure(
      arranged,
      paste0(
        "inst/figure/nowcast_example_",
        data_origin
      ),
      width = 11.5,
      height = 7
    )
    ret <- NULL
  } else {
    ret <- arranged
  }
  ret
}

#' Plot the nowcasts, where the GLM method overshoots
#'
#' @description This function plots and possibly saves the plot of nowcasts and
#' estimate of the mean process \eqn{\lambda_t} for selected dates to highlight
#' the shortcomings of the GLM method leading to somewhat worse results compared
#' to the MCMC method.
#'
#' @param df_nowcast a data frame with columns `date` (date of the
#' nowcasting target), `delay` (reporting delay in weeks), `nowcast_date`,
#' `Distribution`, `quantile_2.5`,`quantile_25`, `quantile_50`, `quantile_75`,
#' `quantile_97.5` and `method`, that contains summarized nowcasting results
#' from both the MCMC and GLM method
#' @param df_lambda a data frame with columns `week`, `.value`,
#' `Distribution` and `nowcast_date`, that contains the distribution of the mean
#' process in a sample format
#' @param df_total a data frame containing columns `date`, `counts` and `data`.
#' The last column `data` is an indicator, whether the values in the `counts`
#' column are the final sums of the counts, or the preliminary data version.
#' Needed to plot the observations alongside the nowcasts.
#' @param dates_to_show a selection of a few consecutive dates for which we want
#' to show the estimates.
#' @param model_to_show a string indicating an observation model, from which we
#' want to show th estimates
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggplot2::ggsave()}
#'
#' @return a ggplot object, or NULL if \code{save_plot = TRUE}
#'
#' @import dplyr ggplot2 patchwork
#'
#' @export
plot_glm_overshoot <- function(
  df_nowcast,
  df_lambda,
  df_total,
  dates_to_show,
  model_to_show = "NegBinX",
  save_plot = TRUE
) {
  # Bind rows of all the data frames that are in a list format
  df_nowcast <- bind_rows(df_nowcast)
  df_total <- bind_rows(df_total) |>
    # Reverse the factor ordering to plot the colors in the correct ordering
    mutate(data = factor(data, levels = c("Preliminary", "Final")))
  df_lambda <- bind_rows(df_lambda) |>
    group_by(.data$week, .data$nowcast_date, .data$method) |>
    summarize(
      lambda_median = median(.data$.value),
      .groups = "drop_last"
    )
  # Derive the length of the rolling window
  window_len <- max(df_lambda$week)
  # Derive the date based on the week number and the nowcast date
  df_lambda <- df_lambda |>
    group_by(.data$nowcast_date, .data$method) |>
    mutate(
      date = as.Date(
        (.data$week - window_len) * 7,
        origin = .data$nowcast_date[1]
      )
    )

  # Set the colors and labels for the GLM and MCMC method in the plot of lambda
  # and the nowcasts
  method_colors <- c("glm" = "sienna3", "mcmc" = "turquoise3")
  method_labels <- c("glm" = "GAM", "mcmc" = "HMC")
  # Set the x-axis breaks
  x_axis_dates <- as.Date(sort(unique(df_total$date)))
  x_axis_breaks <- x_axis_dates[seq(1, length(x_axis_dates), by = 6)]
  # Format the facet titles
  facet_titles <- rlang::set_names(
    paste0("Nowcasts on ", format(as.Date(dates_to_show), "%d %b %Y")),
    dates_to_show
  )

  # Plot the nowcasts faceted by different rolling windows
  p_nowcast <- ggplot() +
    geom_line(
      df_total,
      mapping = aes(x = .data$date, y = .data$counts, color = .data$data)
    ) +
    geom_line(
      df_nowcast,
      mapping = aes(
        x = .data$date,
        y = .data$quantile_50,
        color = .data$method
      ),
      linetype = "dashed"
    ) +
    geom_ribbon(
      df_nowcast,
      mapping = aes(
        x = .data$date,
        ymin = .data$quantile_2.5,
        ymax = .data$quantile_97.5,
        fill = .data$method
      ),
      alpha = 0.2
    ) +
    scale_color_manual(
      values = c("Final" = "black", "Preliminary" = "gray", method_colors),
      breaks = c("Final", "Preliminary"),
      name = "Data"
    ) +
    scale_fill_manual(
      values = method_colors,
      name = "Method",
      labels = method_labels
    ) +
    scale_x_date(
      breaks = x_axis_breaks,
      date_labels = "%d %b",
      minor_breaks = x_axis_dates
    ) +
    labs(x = "Date", y = "Incidence", title = "Nowcast with 95% PIs") +
    facet_wrap(
      ~nowcast_date,
      labeller = as_labeller(facet_titles),
      nrow = 2
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  # Plot the mean process estimates faceted by different rolling windows
  p_lambda <- ggplot() +
    geom_line(
      df_total,
      mapping = aes(x = .data$date, y = .data$counts, color = .data$data)
    ) +
    geom_line(
      df_lambda,
      mapping = aes(
        x = .data$date,
        y = .data$lambda_median,
        color = .data$method
      )
    ) +
    scale_color_manual(
      values = c("Final" = "black", "Preliminary" = "gray", method_colors),
      breaks = c("Final", "Preliminary"),
      name = "Data"
    ) +
    labs(
      x = "Date",
      y = "Incidence",
      title = expression(Estimate~of~lambda[t]) # nolint
    ) +
    scale_x_date(
      breaks = x_axis_breaks,
      date_labels = "%d %b",
      minor_breaks = x_axis_dates
    ) +
    facet_wrap(
      ~nowcast_date,
      labeller = as_labeller(facet_titles),
      nrow = 2
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  # Compose the plots vertically
  p_combined <- (p_nowcast / p_lambda) +
    patchwork::plot_layout(nrow = 2, axes = "collect", guides = "collect")

  if (save_plot) {
    save_figure(
      p_combined,
      paste("inst/figure/glm_overshoot", model_to_show, sep = "_"),
      width = 7,
      height = 8
    )
    ret <- NULL
  } else {
    ret <- p_combined
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
