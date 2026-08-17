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
#' plots. We change the model labels for the plots. NegBin1 changes to NegBin-L
#' to indicate the linear mean-variance relationship. NegBin2 is replaced by
#' NegBin-Q to indicate the quadratic mean-variance relationship. "M" still
#' denotes the multinomial splitting distribution, while "D" denotes the
#' Dirichlet-Multinomial splitting.
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
    "NegBinX.glm" = "NegBin-X-GAM",
    "NegBin1D.glm" = "NegBin-LD-GAM",
    "Poisson.glm" = "Poisson-GAM",
    "NegBinX.mcmc" = "NegBin-X-HMM",
    "NegBin2D.mcmc" = "NegBin-QD-HMM",
    "NegBin1D.mcmc" = "NegBin-LD-HMM",
    "NegBin2M.mcmc" = "NegBin-QM-HMM",
    "NegBin1M.mcmc" = "NegBin-LM-HMM",
    "Poisson.mcmc" = "Poisson-HMM"
  )
  # For plotting the nowcasting bands, we need to change the ordering to have
  # the NegBin2M, NegBin1M and NegBin2D models, which don't have the GLM
  # alternative in one row, Poisson, NegBinX and NegBin1D in another row and
  # GAM models in the last row
  if (nowcast_bands_ordering) {
    ret <- ret[c(7, 8, 5, 9, 4, 6, 3, 1, 2)]
  }
  ret
}

#' Get the model y-axis labels as formatted expressions
#'
#' This function takes the internal name of the data generating process and
#' creates a list of labels for the ggplot y-axis, where the true data
#' generating process is highlighted in bold.
#'
#' @param data_origin string indicating the data generating process
#' @return a named list or vector of 9 elements containing the model labels to
#' show as ticks on the ggplot y-axis. If the \code{data_origin = "case_study"}
#' the function returns the same named vector the \code{get_y_axis_model_labels}
#' fuction would return. If the data generating process is known, the return
#' object is a list and we show the model labels aligned with it in bold. For
#' the NegBinX and NegBin1D models, we fit the models with both methods (GLM and
#' MCMC), so 2 labels will be highlighted in these cases.
get_y_axis_model_labels <- function(data_origin) {
  model_y_labels <- get_interaction_names()
  if (data_origin != "case_study") {
    which_to_highlight <- grepl(data_origin, names(model_y_labels))
    # Make the selected labels in the ggplot in bold
    model_y_labels <- lapply(
      seq_along(which_to_highlight),
      function(ind) {
        if (which_to_highlight[ind]) {
          bquote(bold(.(model_y_labels[ind])))
        } else {
          bquote(.(model_y_labels[ind]))
        }
      }
    )
  }
  model_y_labels
}

#' Get the plot title stating the data generating process
#'
#' This function takes the internal name of the data generating process and
#' converts it to a label consistent with the manuscript. This function is
#' needed to create the subplot titles of the CRPS and coverage plots for the
#' simulation study.
#'
#' @param data_origin string indicating the data generating process
#' @return a string to use as a ggplot title
get_dgp_title <- function(data_origin) {
  formatted_data_origin <- switch(
    data_origin,
    # We use only these three as a data generating process in the simulations
    NegBinX = "NegBin-X",
    NegBin1D = "NegBin-LD",
    NegBin2D = "NegBin-QD"
  )
  paste0(formatted_data_origin, " data generating process")
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
    "Poisson-HMM" = "#CC79A7",
    "NegBin-X-HMM" = "#D55E00",
    "NegBin-QD-HMM" = "#009E73",
    "NegBin-LD-HMM" = "#56B4E9",
    "NegBin-QM-HMM" = "#004282",
    "NegBin-LM-HMM" = "#F0E442",
    "Poisson-GAM" = "#862D67",
    "NegBin-X-GAM" = "#993700",
    "NegBin-LD-GAM" = "#1D79B9"
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
