## ----setup, echo=FALSE, results="hide"-----------------------------------
knitr::opts_chunk$set(comment = "#>", collapse = TRUE)
library(magrittr)

## ---- eval = FALSE-------------------------------------------------------
#  install.packages("BACT")

## ---- eval = FALSE-------------------------------------------------------
#  devtools::install_github("thevaachandereng/BACT@vx.xx.x")
#  #OR
#  devtools::install_version("BACT", version = "x.x.x", repos = "http://cran.us.r-project.org")

## ---- results ='asis', eval = FALSE--------------------------------------
#  devtools::install_github("thevaachandereng/BACT")

## ----lib, results="asis", eval=TRUE--------------------------------------
library(BACT)

## ------------------------------------------------------------------------
set.seed(20000)
enrollment(param = c(0.3, 0.7, 0.9, 1.2), 50, time = c(5, 10, 15))

## ---- eval = FALSE-------------------------------------------------------
#  randomization(N_total = 140, block = 7, scheme = c(2, 1))

## ------------------------------------------------------------------------
randomization(N_total = 140, block = 6, scheme = c(2, 1))

## ------------------------------------------------------------------------
library(magrittr)
value <- 
  proportion(p_control = 0.12, p_treatment = 0.08) %>%
  enrollment_rate(lambda = c(0.3, 1), time = 25) %>%
  sample_size(sample.size = 300, end.of.study = 50) %>%
  looks(interim_look = c(210, 240, 270)) %>%
  impute(no_of_impute = 5) %>%
  randomize(block_size = 2, randomization_ratio = c(1, 1))
  


do.call(binomialBACT, value)

