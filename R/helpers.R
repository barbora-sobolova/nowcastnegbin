#' Get the names of the observation models
#'
#' @return a character vector containing the names of the 6 observational
#' models. They are ordered the same way as in the STAN implementation.
get_model_names <- function() {
  c("Poisson", "NegBinX", "NegBin2D", "NegBin1D", "NegBin2M", "NegBin1M")
}

#' Get the colors of the observation models for the plots
#'
#' @return a named character vector containing the colors of the 6 observational
#' models. They are ordered the same way as in the STAN implementation.
get_model_colors <- function() {
  c(
    "Poisson" = "#CC79A7",
    "NegBinX" = "#D55E00",
    "NegBin2D" = "#009E73",
    "NegBin1D" = "#56B4E9",
    "NegBin2M" = "#004282",
    "NegBin1M" = "#F0E442"
  )
}

#' Get the ggplot theme for figures
#'
#' @return a ggplot theme optimized for the coverage and CRPS plots
get_plot_theme <- function() {
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16),
    legend.spacing.y = unit(0.4, "cm"),
    legend.key.spacing.y = unit(0.2, "cm")
  )
}

#' Labeller function for ggplot facet titles
#'
#' This function takes the number of columns of the reporting triangle and
#' creates titles for the plots faceted by the nowcasting horizon.
#'
#' @param max_lag integer indicating the number of columns of the reporting
#' triangle.
#' @return a named vector of length \code{max_lag} with the factor levels and
#' their corresponding labels
label_horizon_facet <- function(max_lag) {
  horizons <- -rev(seq_len(max_lag) - 1)
  ret <- paste0("Horizon: ", horizons)
  names(ret) <- as.character(horizons)
  ret
}
