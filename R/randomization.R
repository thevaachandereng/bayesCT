#' @title Randomization scheme
#'
#' @description Implements a randomization scheme for control and treatment arms with different randomization ratios and block sizes.
#'
#' @param N_total an integer value of total sample size for randomization scheme.
#' @param block an integer value of the block size for randomization. Note that it needs to be a multiple of the sum of \code{scheme}.
#' @param scheme a numeric vector of the randomization scheme in the order \code{c(control, treatment)}.
#' @description Implements randomization scheme for control and treatment with different randomization ratios and block.
#'
#' @param N_total a numeric value of total sample size for randomization scheme.
#' @param block a integer value of the block size for randomization, it needs to be a multiple of the sum of scheme.
#' @param scheme a numeric vector of the randomization scheme \code{c(control, treatment)}.
#'
#' @return the randomization scheme with 0, 1 for control and treatment
#'
#' @examples
#' # Implementing randomization scheme for control to treatment with 1:1.5 randomization ratio
#' randomization(100, block = 5, c(2, 3))
#'
#' # Randomization scheme with 2:1 for control to treatment
#' randomization(70, block = 9, c(2, 1))
#' @export randomization
randomization <- function(N_total, block = 2, scheme = c(1, 1)) {
  if (block %% 1 != 0 || block <= 0) {
    stop("'block' must be a non-negative integer")
  }

  if (block %% sum(scheme) != 0) {
    stop("The sum of 'scheme' must be a multiple of 'block'")
  }

  if (N_total < block) {
    stop("The number of subjects must be at least the size of 'block'")
  }

  if (any(scheme %% 1 != 0)) {
    stop("All values of 'scheme' must be integer values")
  }

  sampling <- NULL

  while (length(sampling) < (N_total - block)) {
    item <- rep(rep(0:1, times = scheme), each = block / sum(scheme))
    sampling <- c(sampling,
                  sample(x = item, size = block))
  }

  if (N_total > length(sampling)) {
    item <- rep(rep(0:1, times = scheme), each = block / sum(scheme))
    sampling <- c(sampling,
                  sample(item, size = (N_total - length(sampling))))
  }

  return(sampling)
}
