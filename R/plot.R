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
#' @param df_total a data frame
#' @param model_codes a vector of observation model names. Must be in the
#' correct order to label the models correctly. The order in case all models are
#' used is: "Poisson", "NegBinX", "NegBin2D", "NegBin1D", "NegBin2M",
#' "NegBin1M".
#' @param model_colors a named vector of the model colors corresponding to each
#' observation model
#' @param nowcast_date a date, when the nowcast is made
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggplot2::ggsave()}
#'
#' @return a ggplot object with one facet per observation model
#'
#' @import dplyr ggplot2
#'
#' @export
plot_nowcast <- function(
  df_summarized_nowcast,
  df_total,
  model_codes,
  model_colors,
  nowcast_date,
  save_plot = TRUE
) {
  # Convert and order factors so that they display with correct labels and in a
  # correct order
  df_summarized_nowcast <- df_summarized_nowcast |>
    mutate(Distribution = factor(.data$Distribution, labels = model_codes))
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
    scale_fill_manual(values = model_colors) +
    labs(x = "Date", y = "Incidence") +
    facet_wrap(~Distribution)
  # Save the plot if required, the width, height and path is hard-coded here
  if (save_plot) {
    save_figure(
      nowcasts_plot,
      path = paste0("inst/figure/nowcast_plots/nowcast_", nowcast_date),
      width = 9,
      height = 7
    )
  }
  nowcasts_plot
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
#' @param model_codes a vector of observation model names. Must be in the
#' correct order to label the models correctly. The order in case all models are
#' used is: "Poisson", "NegBinX", "NegBin2D", "NegBin1D", "NegBin2M",
#' "NegBin1M".
#' @param model_colors a named vector of the model colors corresponding to each
#' observation model
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggsave()}
#'
#' @return a ggplot object with one facet per nowcasting horizon
#'
#' @import dplyr ggplot2
#'
#' @export
plot_coverage <- function(
  df_summarized_nowcast,
  model_codes,
  model_colors,
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
    ) |>
    ungroup() |>
    # Pivot for easier definition of the alpha aesthetic
    pivot_longer(
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
    scale_fill_manual(values = model_colors) +
    labs(x = "Empirical coverage", title = "Empirical coverage by horizon") +
    facet_wrap(~delay)
  # Save the plot if required, the width, height and path is hard-coded here
  if (save_plot) {
    save_figure(
      coverage_plot,
      "inst/figure/coverage_plot",
      width = 9,
      height = 7
    )
  }
  coverage_plot
}

#' Plot and save the CRPS density
#'
#' @description This function plots and possibly saves the chart of CRPS
#' densities for all models.
#'
#' @param df_summarized_nowcast a data frame containing columns `Distribution`
#' (containing the name of the observation model), `CRPS` (the empirical
#' distribution of the CRPS) and `delay` (the nowcasting horizon)
#' @param model_colors a named vector of the model colors corresponding to each
#' observation model
#' @param save_plot logical indicator, whether to save the plot using
#' \code{ggsave()}
#'
#' @return a ggplot object with one facet per nowcasting horizon
#'
#' @import dplyr ggplot2
#'
#' @export
plot_crps <- function(
  df_summarized_nowcast,
  model_colors,
  save_plot = TRUE
) {
  crps_plot <- ggplot(
    df_summarized_nowcast,
    aes(x = .data$CRPS, color = .data$Distribution)
  ) +
    # Plot the density of the CRPS
    geom_line(stat = "density", alpha = 0.6) +
    scale_color_manual(values = model_colors) +
    labs(x = "CRPS", title = "CRPS distribution by horizon") +
    facet_wrap(~delay) +
    xlim(c(0, 2000))
  # Save the plot if required, the width, height and path is hard-coded here
  if (save_plot) {
    save_figure(
      crps_plot,
      "inst/figure/crps_plot",
      width = 7,
      height = 5.5
    )
  }
  crps_plot
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
#' @return a ggplot object
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
  # Save the plot if required, the width, height and path is hard-coded here
  if (save_plot) {
    save_figure(
      trajectory_plot,
      "inst/figure/SARI_trajectory",
      width = 9,
      height = 7
    )
  }
  trajectory_plot
}

#' Save a figure in the PDF and the PNG format
#'
#' @param figure a ggplot chart to be saved
#' @param path a file path indicating where to save the plot, typically starting
#'  with "inst/figure". No file extention should be included.
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
