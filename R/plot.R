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
    mutate(Distribution = factor(Distribution, labels = model_codes))
  df_total <- df_total |>
    mutate(data = factor(data, levels = c("Preliminary", "Final")))
  # Plot the nowcasts
  nowcasts_plot <- ggplot() +
    # Point prediction
    geom_line(
      data = df_summarized_nowcast,
      mapping = aes(
        x = date,
        y = quantile_50,
        color = Distribution,
        linetype = "Nowcast"
      )
    ) +
    # Different data versions - preliminary and final
    geom_line(
      data = df_total,
      mapping = aes(x = date, y = counts, color = data, linetype = data)
    ) +
    # 95% prediction intervals
    geom_ribbon(
      data = df_summarized_nowcast,
      mapping = aes(
        x = date,
        ymin = quantile_2.5,
        ymax = quantile_97.5,
        fill = Distribution,
        alpha = "PI_95"
      )
    ) +
    # 50% prediction intervals
    geom_ribbon(
      data = df_summarized_nowcast,
      mapping = aes(
        x = date,
        ymin = quantile_25,
        ymax = quantile_75,
        fill = Distribution,
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

plot_coverage <- function(
  df_summarized_nowcast,
  model_codes,
  model_colors,
  save_plot = TRUE
) {
  # Calculate the empirical coverage
  df_coverage <- df_summarized_nowcast |>
    group_by(delay, Distribution) |>
    summarize(
      coverage_50 = sum(true_val > quantile_25 & true_val < quantile_75) / n(),
      coverage_95 = sum(true_val > quantile_2.5 & true_val < quantile_97.5) /
        n(),
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
      x = empirical_coverage,
      y = Distribution,
      fill = Distribution,
      alpha = nominal_coverage
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

plot_crps <- function(
  df_summarized_nowcast,
  model_colors,
  save_plot = TRUE
) {
  crps_plot <- ggplot(
    df_summarized_nowcast,
    aes(x = CRPS, color = Distribution)
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
    pull(counts) |>
    # na.rm = TRUE is usually not needed, but it prevents the plot element to
    # disappear in the case of missing values
    max(na.rm = TRUE)
  last_window_max_cases <- totals |>
    filter(date >= last_window_beg) |>
    pull(counts) |>
    max(na.rm = TRUE)
  # 5% offset of the braces to avoid overplotting the trajectory
  bracket_offset <- first_window_max_cases * 0.05
  trajectory_plot <- ggplot(totals, aes(x = date, y = counts)) +
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
