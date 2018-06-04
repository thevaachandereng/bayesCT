library(magrittr)
input <- list(samplesize = NULL,
              n.looks = NULL)

looks <- function(n.looks = NULL,
                  samplesize = NULL,
                  data = NULL){
  stopifnot(n.looks >= 2, length(samplesize) == n.looks)
  data$n.looks <- n.looks
  data$samplesize <- samplesize
  data
}

enroll <- function(dist = "poisson",
                   param = NULL,
                   data = NULL){
  stopifnot(is.vector(param))
  data$dist <- dist
  data$param <- param
  data
}





input %>%
  looks(n.looks = 3, samplesize = c(20, 40, 60)) %>%
  enroll(dist = "poisson", param = 20)

