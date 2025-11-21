plot_nowcast <- function(
  df_summarized_nowcast,
  df_total,
  model_codes,
  model_colors
) {
  df_summarized_nowcast <- df_summarized_nowcast |>
    mutate(Distribution = factor(Distribution, labels = model_codes))
  df_total <- df_total |>
    mutate(data = factor(data, levels = c("Preliminary", "Final")))
  browser()
  ggplot() +
    geom_line(
      data = df_summarized_nowcast,
      mapping = aes(
        x = date,
        y = quantile_50,
        color = Distribution,
        linetype = "Nowcast"
      )
    ) +
    geom_line(
      data = df_total,
      mapping = aes(x = date, y = counts, color = data, linetype = data)
    ) +
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
    scale_alpha_manual(
      values = c("PI_50" = 0.2, "PI_95" = 0.2),
      labels = c("50%", "95%"),
      name = "Prediction interval"
    ) +
    guides(
      linetype = guide_legend(
        override.aes = list(color = c("gray60", "black", "black"))
      ),
      alpha = guide_legend(override.aes = list(alpha = c(0.4, 0.2)))
    ) +
    scale_fill_manual(values = model_colors) +
    labs(x = "Date", y = "Incidence") +
    facet_wrap(~Distribution)
}
