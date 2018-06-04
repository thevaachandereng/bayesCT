#' Simulating enrollment dates
#'
#' This function simulates enrollment dates using either poisson or binomial distribution
#'
#' @param family a factor value of either "poisson" or "binomial"
#' @param param a vector of lambda in poisson
#' @param time a vector of the same size as param indicating end of time when a specific lambda is used
#' @param N_total a numeric value of
#' @return a vector of enrollment dates
#'
#' @examples enrollment(param = c(0.003, 0.7), 100, time = c(10, 400))
#'
#'
#'


enrollment <- function(param, N_total, time = NULL){
  #checking each inputs
  stopifnot(all(param > 0), if(length(param) >= 1){!is.null(time)})
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
      #if the time exceeds count, reuse the last lambda
      if(sum(time >= (count - 1)) == 0){
        index <- length(param)
      }
      #else pick the minimum time where time exceeds count
      else{
        index <- min(which(time >= (count - 1)))
      }
      output <- c(output, rep(count, rpois(1, param[index])))
    }
  }
  output <- output[1:N_total]
  output <- output - output[1]
  output
}
