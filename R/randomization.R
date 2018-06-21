#' @title Randomization scheme
#'
#' @description Implements randomization scheme for control and treatment with different randomization ratios and block
#'
#' @param N_total a numeric value of total sample size for randomization scheme
#' @param block a integer value of the block size for randomization, it needs to be a multiple of the sum of scheme
#' @param scheme a numeric vector of the randomization scheme \code{c(control, treatment)}
#'
#' @return the randomization scheme with 0, 1 for control and treatment
#'
#' @examples
#' #implementing randomization scheme for control to treatment with 1:1.5 randomization ratio
#' randomization(100, block = 5, c(2, 3))
#' #randomization scheme with 2:1 for control to treatment
#' randomization(70, block = 9, c(2, 1))
#'

randomization <- function(N_total, block = 2, scheme = c(1, 1)){
  stopifnot(block %% 1 == 0, block %% sum(scheme) == 0,
            N_total >= block, block > 0, all(scheme %% 1 == 0))
  sampling <- NULL
  while(length(sampling) < (N_total - block)){
    item <- rep(rep(0:1, scheme), each = block / sum(scheme))
    sampling <- c(sampling, sample(x = item, size = block))
  }
  if(N_total > length(sampling)){
  sampling <- c(sampling, sample(item, size = (N_total - length(sampling))))
  }
  return(sampling)
}
