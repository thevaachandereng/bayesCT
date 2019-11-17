#' @title Simulates time-to-event outcomes.
#'
#' @description Simulation of time-to-event outcomes using
#' the piecewise constant hazard exponential function.
#'
#' @param lambda vector. The constant hazard rates for exponential failures.
#' @param cutpoint vector. The change-point vector indicating time when the hazard rates change.
#' @param n scalar. The number of outcomes for simulation.
#' @param maxtime scalar. maximum time before end of study.
#'
#' @return a dataset with simulated follow-up time (time) and respective event indicator (1 = event,
#'          0 = censoring)
#'
#' @examples pw_exp_sim(c(0.02, 0.01, 0.005), 100, 100, c(10, 20))
#'           pw_exp_sim(0.015, 100, 100)
#'
#'
#' @importFrom stats rexp
#' @export pw_exp_sim


pw_exp_sim <- function(lambda, n, maxtime, cutpoint = NULL) {
  ## checking input of parameters
  #make sure lambda is positive
  if(any(lambda < 0)) {stop("At least one of the hazard rate(s) is less than 0!")}

  # make sure n is positive integer
  if(n <= 0 | n %% 1 != 0 ) {stop("The number of simulations need to be a positive integer!") }

  #make sure maxtime is positive
  if(maxtime <= 0 | length(maxtime) > 1) {stop("The maxtime needs to be greater than 0
                                               and the length of maxtime needs to be 1!")}

  #make sure lambda is positive
  if(length(lambda) > 1) {
    if(length(lambda) != (length(cutpoint) + 1)){
      stop("The length of cutpoint needs to be one less than the length of lambda!")
    }
  }

  # make sure cutpoint is null if length of lambda is 1.
  if(length(lambda) == 1){
    if(!is.null(cutpoint)){
      stop("cutpoint needs to be null if length of lambda is 1!")
    }
  }

  # generate exp(1)
  x <- rexp(n)

  if(is.null(cutpoint)) {
    time <- x / lambda
  }

  if(length(cutpoint) >= 1) {
    lamtau <- lambda[1] * cutpoint[1]
    if(length(cutpoint) > 1){
      lamtau <- c(lamtau, rep(NA, (length(lambda) - 1)))
      for(h in 2:length(cutpoint)){
        lamtau[h] <- lambda[h] * (cutpoint[h] - cutpoint[h - 1])
      }
    }
    lamtau  <- cumsum(lamtau)
    lamtau[length(lamtau) + 1] <- 99999
    cp <- function(p){
      time  <- rep(NA, length(p))
      for(m in 1:length(p)){
        index <- min(which(p[m] < lamtau))
      if(p[m] < lamtau[1]){
        time[m] <-  p[m] / lambda[1]
      }
      else if(p[m] >= lamtau[1]){
        time[m] <- ((p[m] - lamtau[(index - 1)]) / lambda[index] + cutpoint[index - 1])
        }
      }
      return(time)
    }
    time <- cp(x)
  }

  # if maxtime is lower than observed time, censor the data
  min_time  <- pmin(time, maxtime)
  event     <- as.numeric(time == min_time)
  dat       <- data.frame(time = min_time, event = event)

  # return dataset with time and event
  return(dat)
}


