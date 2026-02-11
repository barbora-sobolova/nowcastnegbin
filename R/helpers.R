#' Get the names of the observation models
#'
#' @return a character vector containing the names of the 6 observational
#' models. They are ordered the same way as in the STAN implementation.
get_model_names <- function() {
  c("Poisson", "NegBinX", "NegBin2D", "NegBin1D", "NegBin2M", "NegBin1M")
}
