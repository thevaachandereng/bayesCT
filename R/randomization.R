#' @title Randomization allocation
#'
#' @description Implements a randomization allocation for control and treatment arms with different randomization ratios and block sizes.
#'
#' @param N_total an integer value of total sample size for randomization allocation.
#' @param block a vector value of the block size for randomization. Note that it needs to be a multiple of the sum of \code{allocation}.
#' @param allocation a numeric vector of the randomization allocation in the order \code{c(control, treatment)}.
#'
#' @return the randomization allocation with 0, 1 for control and treatment
#'
#' @examples
#' # Implementing treatment allocation for control to treatment with 1:1.5 randomization ratio
#' randomization(N_total = 100, block = 5, allocation = c(2, 3))
#'
#' # Treatment allocation with 2:1 for control to treatment
#' randomization(N_total = 70, block = 9, allocation = c(2, 1))
#'
#' # Treatment allocation for control to treatment with 1:2 for control to treatment with
#' # multiple block sizes c(3, 9, 6)
#' randomization(N_total = 100, block = c(3, 9, 6), allocation = c(1, 2))
#'
#' # For complete randomization set the N_total to block size
#' randomization(N_total = 100, block = 100, allocation = c(1, 1))
#'
#' @export randomization

randomization <- function(N_total, block = 2, allocation = c(1, 1)) {
  if (any(block %% 1 != 0) | any(block <= 0)) {
    stop("'block' must be a non-negative integer")
  }

  if (any(block %% sum(allocation) != 0)) {
    stop("The sum of 'allocation' must be a multiple of 'block'")
  }

  if (N_total < sum(block)) {
    stop("The number of subjects must be at least the size of 'block'")
  }

  if (any(allocation %% 1 != 0)) {
    stop("All values of 'allocation' must be integer values")
  }

  sampling <- NULL
  next_block <- NULL

  # creating different block sizes for multiple blocks
  blocking <- rep(block, N_total %/% sum(block))

  for(k in 1:length(block)){
    if((sum(blocking) + block[k]) < N_total){
     blocking <- c(blocking, block[k])
     next_block <- block[k + 1]
    }
    else{
      break
    }
  }

  # making sure the next block is assigned
  if(is.null(next_block)){
    next_block <- block[1]
  }

  # within each block, randomize with the correct allocation
  for(m in 1:length(blocking)){
    item <- rep(rep(0:1, times = allocation), each = blocking[m] / sum(allocation))
    sampling <- c(sampling, sample(item))
  }

  # fill up the remainder of the allocation using next block
  if (N_total > length(sampling)) {
    item <- rep(rep(0:1, times = allocation), each = next_block / sum(allocation))
    sampling <- c(sampling,
                  sample(item, size = (N_total - length(sampling))))
  }

  # return sampling
  return(sampling)
}
