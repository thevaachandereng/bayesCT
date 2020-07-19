#' @title Simulating enrollment dates
#'
#' @description This function simulates enrollment dates using either Poisson
#'   distribution
#'
#' @param param vector. Lambda values for Poisson distribution.
#' @param time vector. Knots (of \code{length(param)} - 1) indicating end of
#'   time when a specific lambda is used.
#' @param N_total integer. Value of total sample size.
#'
#' @return A vector of enrollment times (from time of first patient enrollment)
#'   in days.
#'
#' @examples
#' enrollment(param = c(0.003, 0.7), 100, time = 10)
#'
#' enrollment(param = c(0.3, 0.5, 0.9, 1.2, 2.1), 200, c(20, 30, 40, 60))
#'
#' @importFrom stats rpois
#'
#' @export enrollment

enrollment <- function(param, N_total, time = NULL) {
  if (any(param <= 0)) {
    stop("The lambda(s) for poisson enrollment rate should be non-negative")
  }

  if (N_total <= 0) {
    stop("The sample size for enrollment needs to be greater than 0")
  }

  if (length(param) > 1  & (is.null(time) | (length(param) != (length(time) + 1)))) {
    stop("The cutoff time for lambda is not correct")
  }

  if (length(param) == 1 & !is.null(time)) {
    warning("The time input is not used!")
  }

  output <- NULL
  count <- 0
  # For constant lambda in Poisson distribution
  if (length(param) == 1) {
    while (length(output) < N_total) {
      count <- count + 1
      output <- c(output, rep(count, rpois(1, param)))
    }
  }
  # For different lambda values in Poisson distribution as a function of time
  else {
    while (length(output) < N_total) {
      count <- count + 1
      index <- min(c(which(time >= (count - 1)), length(param)))
      output <- c(output, rep(count, rpois(1, param[index])))
    }
  }

  # Adjusting trial and outputs
  output <- output[1:N_total]
  output <- output - output[1]
  return(output)
}
