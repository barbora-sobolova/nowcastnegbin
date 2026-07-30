#' Get the names of the observation models
#'
#' @return a character vector containing the names of the 6 observation
#' models. They are ordered the same way as in the STAN implementation.
get_model_names <- function() {
  c("Poisson", "NegBinX", "NegBin2D", "NegBin1D", "NegBin2M", "NegBin1M")
}

#' Get the names of the interaction of an observation model and a fitting
#' method
#'
#' @description This function creates a named vector that controls the names
#' and ordering of the observation models as they appear in the aggregated
#' plots.
#'
#' @param nowcast_bands_ordering logical indicator switching between different
#' vector orderings. If \code{nowcast_bands_ordering = FALSE}, the ordering
#' corresponds to the coverage and CRPS plots. If TRUE, the models are ordered
#' as in the plot of the incidence trajectory alongside the nowcasts.
#'
#' @return a named character vector containing the names of the 9 models
#' (6 for the MCMC method and 3 for the GLM method). The elements of the vector
#' correspond to model names as they appear in the plot. The names of the vector
#' correspond to the model names as they are created by the
#' \code{interaction()} function
get_interaction_names <- function(nowcast_bands_ordering = FALSE) {
  ret <- c(
    "NegBinX.glm" = "NegBinX-GAM",
    "NegBin1D.glm" = "NegBin1D-GAM",
    "Poisson.glm" = "Poisson-GAM",
    "NegBinX.mcmc" = "NegBinX-HMC",
    "NegBin2D.mcmc" = "NegBin2D-HMC",
    "NegBin1D.mcmc" = "NegBin1D-HMC",
    "NegBin2M.mcmc" = "NegBin2M-HMC",
    "NegBin1M.mcmc" = "NegBin1M-HMC",
    "Poisson.mcmc" = "Poisson-HMC"
  )
  # For plotting the nowcasting bands, we need to change the ordering to have
  # the NegBin1M, NegBin2M and NegBin2D models, which don't have the GLM
  # alternative in one row, Poisson, NegBinX and NegBin1D in another row and
  # GAM models in the last row
  if (nowcast_bands_ordering) {
    ret <- ret[c(7, 8, 5, 9, 4, 6, 3, 1, 2)]
  }
  ret
}

#' Get the colors of the observation models for the plots
#'
#' @return a named character vector containing the colors of the 6 observation
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

#' Get the colors of the interaction of an observation model and a fitting
#' method for the plots
#'
#' @return a named character vector containing the colors of the 9 models
#' (6 for the MCMC method and 3 for the GLM method).
get_interaction_colors <- function() {
  c(
    "Poisson-HMC" = "#CC79A7",
    "NegBinX-HMC" = "#D55E00",
    "NegBin2D-HMC" = "#009E73",
    "NegBin1D-HMC" = "#56B4E9",
    "NegBin2M-HMC" = "#004282",
    "NegBin1M-HMC" = "#F0E442",
    "Poisson-GAM" = "#862D67",
    "NegBinX-GAM" = "#993700",
    "NegBin1D-GAM" = "#1D79B9"
  )
}



#' Shared ggplot theme for standardized text sizing and legend spacing
#'
#' @return a \code{ggplot2::theme()} object
#'
#' @import ggplot2
get_plot_theme <- function() {
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14),
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
