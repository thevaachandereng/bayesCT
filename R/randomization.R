#' @title Randomization allocation
#'
#' @description Implements a randomization allocation for control and treatment arms with different randomization ratios and block sizes.
#'
#' @param N_total an integer value of total sample size for randomization allocation.
#' @param block an integer value of the block size for randomization. Note that it needs to be a multiple of the sum of \code{allocation}.
#' @param allocation a numeric vector of the randomization allocation in the order \code{c(control, treatment)}.
#'
#' @return the randomization allocation with 0, 1 for control and treatment
#'
#' @examples
#' # Implementing randomization allocation for control to treatment with 1:1.5 randomization ratio
#' randomization(N_total = 100, block = 5, allocation = c(2, 3))
#'
#' # Randomization allocation with 2:1 for control to treatment
#' randomization(N_total = 70, block = 9, allocation = c(2, 1))
#'
#' # For complete randomization set the N_total to block size
#' randomization(N_total = 100, block = 100, allocation = c(1, 1))
#'
#' @export randomization
randomization <- function(N_total, block = 2, allocation = c(1, 1)) {
  if (block %% 1 != 0 || block <= 0) {
    stop("'block' must be a non-negative integer")
  }

  if (block %% sum(allocation) != 0) {
    stop("The sum of 'allocation' must be a multiple of 'block'")
  }

  if (N_total < block) {
    stop("The number of subjects must be at least the size of 'block'")
  }

  if (any(allocation %% 1 != 0)) {
    stop("All values of 'allocation' must be integer values")
  }

  sampling <- NULL

  while (length(sampling) < (N_total - block)) {
    item <- rep(rep(0:1, times = allocation), each = block / sum(allocation))
    sampling <- c(sampling,
                  sample(x = item, size = block))
  }

  if (N_total > length(sampling)) {
    item <- rep(rep(0:1, times = allocation), each = block / sum(allocation))
    sampling <- c(sampling,
                  sample(item, size = (N_total - length(sampling))))
  }

  return(sampling)
}
