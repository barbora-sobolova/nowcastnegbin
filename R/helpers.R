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
