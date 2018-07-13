#' @title Simulating enrollment dates
#'
#' @description This function simulates enrollment dates using either poisson or binomial distribution
#'
#' @param param a vector of lambda in poisson
#' @param time a vector of the \code{length(param)} - 1 indicating end of time when a specific lambda is used
#' @param N_total a numeric value of total sample size
#' @return a vector of enrollment dates
#'
#' @examples
#' enrollment(param = c(0.003, 0.7), 100, time = 10)
#' enrollment(param = c(0.3, 0.5, 0.9, 1.2, 2.1), 200, c(20, 30, 40, 60))
#'
#' @importFrom stats rpois
#' @export enrollment
enrollment <- function(param, N_total, time = NULL){
  #checking each inputs
  stopifnot(all(param > 0), if(length(param) >= 1){!is.null(time)},
            N_total > 0, length(param) == (length(time) + 1))
  output <- NULL
  count <- 0
  # for constant lambda in poisson
  if(length(param) == 1){
    while(length(output) < N_total){
      count <- count + 1
      output <- c(output, rep(count, rpois(1, param)))
    }
  }
  # for different lambda in poisson a function of time
  else{
    while(length(output) < N_total){
      count <- count + 1
      index <- min(c(which(time >= (count - 1)), length(param)))
      output <- c(output, rep(count, rpois(1, param[index])))
    }
  }
  #adjusting trial and ouputs
  output <- output[1:N_total]
  output <- output - output[1]
  output
}
