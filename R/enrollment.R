#' Simulating enrollment dates
#'
#' This function simulates enrollment dates using either poisson or binomial distribution
#'
#' @param family a factor value of either "poisson" or "binomial"
#' @param param a numeric value of lambda in poisson or value of a in uniform(0, a)
#' @param N_total a numeric value of
#' @return a vector of enrollment dates
#'
#' @examples enrollment("uniform", 3, 100)
#'
#'
#'


enrollment <- function(family = "poisson", param, N_total){
  stopifnot((family == "uniform" | family == "poisson"), is.numeric(param), param > 0)
  output <- NULL
  if(family == "uniform"){
    output <- round(runif(N_total, 0, param) + cumsum(c(0, rep(param, (N_total - 1)))))
  }
  else if(family == "poisson"){
    count <- 0
    while(length(output) < N_total){
      output <- c(output, rep(count, rpois(1, param)))
      count <- count + 1
    }
    output <- output[1:N_total]
    output <- output - output[1]
  }
  output
}
