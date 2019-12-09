#' @title Simulates time-to-event outcomes.
#'
#' @description Simulation of time-to-event outcomes using
#' the piecewise constant hazard exponential function.
#'
#' @param hazard vector. The constant hazard rates for exponential failures.
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


pw_exp_sim <- function(hazard, n, maxtime = NULL, cutpoint = NULL) {
  ## checking input of parameters
  #make sure hazard is positive
  if(any(hazard < 0)) {stop("At least one of the hazard rate(s) is less than 0!")}

  # make sure n is positive integer
  if(n <= 0 | n %% 1 != 0 ) {stop("The number of simulations need to be a positive integer!") }

  #make sure maxtime is positive or null
  if(!is.null(maxtime)){
    if(maxtime <= 0 | length(maxtime) > 1) {
      stop("The maxtime needs to be greater than 0 and
           the length of maxtime needs to be 1!")
      }
  }

  #make sure hazard is positive
  if(length(hazard) > 1) {
    if(length(hazard) != (length(cutpoint) + 1)){
      stop("The length of cutpoint needs to be one less than the length of hazard!")
    }
  }

  # make sure cutpoint is null if length of hazard is 1.
  if(length(hazard) == 1){
    if(!is.null(cutpoint)){
      stop("cutpoint needs to be null if length of hazard is 1!")
    }
  }

  # generate exp(1)
  x <- rexp(n)

  if(is.null(cutpoint)) {
    time <- x / hazard
  }

  if(length(cutpoint) >= 1) {
    lamtau <- hazard[1] * cutpoint[1]
    if(length(cutpoint) > 1){
      lamtau <- c(lamtau, rep(NA, (length(hazard) - 2)))
      for(h in 2:length(cutpoint)){
        lamtau[h] <- hazard[h] * (cutpoint[h] - cutpoint[h - 1])
      }
    }
    lamtau  <- cumsum(lamtau)
    lamtau[length(lamtau) + 1] <- 99999
    cp <- function(p){
      time  <- rep(NA, length(p))
      for(m in 1:length(p)){
        index <- min(which(p[m] < lamtau))
      if(p[m] < lamtau[1]){
        time[m] <-  p[m] / hazard[1]
      }
      else if(p[m] >= lamtau[1]){
        time[m] <- ((p[m] - lamtau[(index - 1)]) / hazard[index] + cutpoint[index - 1])
        }
      }
      return(time)
    }
    time <- cp(x)
  }

  # if maxtime is lower than observed time, censor the data
  if(!is.null(maxtime)){
    min_time  <- pmin(time, maxtime)
    event     <- as.numeric(time == min_time)
    dat       <- data.frame(time = min_time, event = event)
  }
  else{
    dat       <- data.frame(time = time, event = rep(1, length(time)))
  }

  # return dataset with time and event
  return(dat)
}

#' @title Imputes time-to-event outcomes.
#'
#' @description Imputation of time-to-event outcomes using
#' the piecewise constant hazard exponential function.
#'
#' @param time vector. The observed time for patient that have had no event or passed maxtime.
#' @param hazard vector. The constant hazard rates for exponential failures.
#' @param maxtime scalar. maximum time before end of study.
#' @param cutpoint vector. The change-point vector indicating time when the hazard rates change.
#'
#' @return a dataset with simulated follow-up time (time) and respective event indicator (1 = event,
#'          0 = censoring)
#'
#' @examples pw_exp_impute(time = c(120), c(0.005, 0.001), 110, 40)
#'           pw_exp_impute(time = c(10, 20, 30), c(0.005, 0.01, 0.02), 100, c(40, 80))
#'           pw_exp_impute(time = c(40, 30), c(0.005, 0.01), 120, c(50))
#'
#'
#' @importFrom stats rexp
#' @export pw_exp_impute


pw_exp_impute <- function(time, hazard, maxtime = NULL, cutpoint = NULL) {
  ## checking input of parameters
  #make sure hazard is positive
  if(any(hazard < 0)) {stop("At least one of the hazard rate(s) is less than 0!")}

  # make sure n is positive integer
  if(any(time < 0)) {stop("The time has to be positive!") }

  #make sure maxtime is positive or null
  if(!is.null(maxtime)){
    if(maxtime <= 0 | length(maxtime) > 1) {
      stop("The maxtime needs to be greater than 0 and
           the length of maxtime needs to be 1!")
      }
  }

  #make sure hazard is positive
  if(length(hazard) > 1) {
    if(length(hazard) != (length(cutpoint) + 1)){
      stop("The length of cutpoint needs to be one less than the length of hazard!")
    }
  }

  # make sure cutpoint is null if length of hazard is 1.
  if(length(hazard) == 1){
    if(!is.null(cutpoint)){
      stop("cutpoint needs to be null if length of hazard is 1!")
    }
  }

  # generate exp(1)
  x <- rexp(length(time))

  large <- max(time) + 1

  t <- rep(NA, length(time))

  for(k in 1:length(time)){

      if(is.null(cutpoint)) {
        t[k] <- x[k] / hazard + time[k]
      }

      if(length(cutpoint) >= 1) {
        index    <- min(which(time[k] < c(cutpoint, large)))
        interval <- c(time[k], cutpoint[index:length(cutpoint)])
        haz      <- hazard[index:length(hazard)]
        lamtau   <- NULL
        if(length(cutpoint) > 1){
          lamtau <- rep(NA, (length(haz) - 1))
          if(length(lamtau) > 0){
            for(h in 1:(length(haz) - 1)){
              lamtau[h] <- haz[h] * (interval[h + 1] - interval[h])
            }
          }
        }
        lamtau[length(lamtau) + 1] <- 99999
        lamtau  <- cumsum(lamtau)
        cp <- function(p){
          index <- min(which(p < lamtau))
          if(p < lamtau[1]){
            return(p / hazard[1])
          }
          else if(p >= lamtau[1]){
            return(((p - lamtau[(index - 1)]) / hazard[index] + cutpoint[index - 1]))
          }
        }
        t[k] <- cp(x[k]) + time[k]
      }
  }

  if(!is.null(maxtime)){
    # if maxtime is lower than observed time, censor the data
    min_time  <- pmin(t, maxtime)
    event     <- as.numeric(t == min_time)
    dat       <- data.frame(time = min_time, event = event)
  }
  else{
    dat       <- data.frame(time = t, event = rep(1, length(time)))
  }

  # return dataset with time and event
  return(dat)
}

