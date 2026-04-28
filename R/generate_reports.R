#' Generate the reporting table of a nowcasting problem
#'
#' @description Generate the full reporting table, i.e. with no right censoring.
#' The mean process follows an exponential random walk, or can be specified by
#' a sequence of values.
#'
#' @param lgt Integer, the length of the simulated table.
#' @param max_lag Integer, the maximum reporting delay, the width of the table.
#' @param log_lambda0 Numeric, the initial value of the random walk.
#' @param rw_noise_sd Numeric, standard deviation of the geometric random walk.
#' @param probs Numeric, a numeric vector specifying the delay distribution.
#' @param nb_size Numeric, a positive real value specifying the size of the
#' negbin distribution. The lower, the more dispersed
#' @param model name of the observation model
#' @param fixed_lambda Numeric vector containing the values of the mean process.
#' @param seed An integer for seeding the simulation
#'
#' @return The reporting table as a matrix
#'
#' @importFrom stats rgamma rnbinom rnorm rpois
#'
#' @export
#'
#' @examples
#' # Example: Sample a reporting table from the Poisson model for 100 days,
#' # with the maximum delay 3.
#' generate_reports(100, 3, log_lambda0 = log(12), rw_noise_sd = 0.9,
#'                  probs = c(0.3, 0.4, 0.3), model = "Poisson")
#'
#' # Example: Sample a reporting table from the NegBin2D model for 10 days,
#' # with the maximum delay 3 and a pre-specified mean process
#' generate_reports(10, 3,
#'                  fixed_lambda = c(10.3, 23, 17, 9, 5.1, 2, 4, 7, 11.7, 13),
#'                  probs = c(0.3, 0.4, 0.3), nb_size = 0.8, model = "NegBin2D")
generate_reports <- function(
  lgt,
  max_lag,
  probs,
  log_lambda0 = NULL,
  rw_noise_sd = NULL,
  nb_size = NULL,
  model = c(
    "Poisson",
    "NegBinX",
    "NegBin2D",
    "NegBin1D",
    "NegBin2M",
    "NegBin1M"
  ),
  fixed_lambda = NULL,
  seed = 123456
) {
  model <- match.arg(
    model,
    c("Poisson", "NegBinX", "NegBin2D", "NegBin1D", "NegBin2M", "NegBin1M")
  )

  # If the mean process lambda is not specified already, we generate a random
  # walk
  set.seed(seed)
  if (is.null(fixed_lambda)) {
    if (is.null(log_lambda0) || is.null(rw_noise_sd)) {
      stop(
        "If no mean process is given, random walk parameters must be specified."
      )
    }
    wn <- rnorm(lgt - 1, 0, 1)
    log_lambda <- log_lambda0 + cumsum(c(0, wn)) * rw_noise_sd
    lambda <- exp(log_lambda)
  } else {
    if (length(fixed_lambda) != lgt) {
      stop("'fixed_lambda' must be of length 'lgt'")
    }
    lambda <- fixed_lambda
  }

  # The expected counts split by delay
  exp_obs <- lambda %*% t(probs)

  # Calculate the size parameter according to the model
  nb_sizes <- switch(
    model,
    NegBinX = rep(nb_size, lgt * max_lag),
    NegBin2D = rep(nb_size, lgt) %*% t(probs),
    NegBin1D = exp_obs * nb_size
  )

  # Generate reports
  if (model %in% c("NegBinX", "NegBin2D", "NegBin1D")) {
    samp <- rnbinom(lgt * max_lag, nb_sizes, mu = exp_obs) |>
      matrix(nrow = lgt, ncol = max_lag)
  } else {
    # Sample the random effect, when necessary
    if (model %in% c("NegBin2M", "NegBin1M")) {
      if (model == "NegBin1M") {
        random_effect_gamma_par <- lambda * nb_size
      } else {
        random_effect_gamma_par <- rep(nb_size, lgt)
      }
      random_effect <- rgamma(
        lgt,
        shape = random_effect_gamma_par,
        rate = random_effect_gamma_par
      )
      exp_obs <- sweep(exp_obs, 1, random_effect, "*")
    }
    samp <- rpois(lgt * max_lag, lambda = exp_obs) |>
      matrix(nrow = lgt, ncol = max_lag)
  }
  list(reports = samp, exp_obs_total = lambda)
}
