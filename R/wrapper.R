

#' @title Sample size wrapper
#'
#' @description Wrapper function for sample size and length of study.
#'
#' @param sample.size integer. The number of sample size needed.
#' @param end.of.study integer. The length of the study.
#' @param data NULL. stores the sample size and length of study.
#'
#' @return a list with sample size and length of the study.
#'
#' @examples sample_size(sample.size = 300, end.of.study = 50)
#' @export sample_size
sample_size <- function(sample.size, end.of.study, data = NULL){
  data$N_total <- sample.size
  data$EndofStudy <- end.of.study
  data
}





#' @title Interim look wrapper
#'
#' @description Wrapper function for interim looks
#'
#' @param interim_look vector. Vector with interim looks.
#' @param data NULL. stores the impute function
#'
#' @return a list with interim look information
#'
#' @examples looks(interim_look = c(210, 240, 270))
#' @export looks
looks <- function(interim_look = NULL, data = NULL){
  data$interim_look <- interim_look
  data
}




#' @title Enrollment rate wrapper
#'
#' @description Wrapper function for enrollment rate
#'
#' @param lambda vector. Vector with different enrollment rate.
#' @param time vector. Vector with different cutoff time for lambda.
#' @param data NULL. stores the impute function
#'
#' @return a list with enrollment rate information
#'
#' @examples enrollment_rate(lambda = c(0.3, 1), time = 25)
#' @export enrollment_rate
enrollment_rate <- function(lambda = NULL, time = NULL, data = NULL){
  data$lambda <- lambda
  data$lambda_time <- time
  data
}

#' @title Imputation wrapper
#'
#' @description Wrapper function for no_of _impute
#'
#' @param no_of_impute integer. Number of Monte Carlo imputation for missing data
#' @param data NULL. stores the impute function
#'
#' @return a list with number of imputation
#'
#' @examples impute(no_of_impute = 100)
#' @export impute
impute <- function(no_of_impute, data = NULL){
  data$N_impute <- no_of_impute
  data
}


#' @title Randomization scheme wrapper
#'
#' @description Wrapper function for the randomization scheme in the trial
#'
#' @param block_size integer. Block size for the complete randomization in a block
#' @param randomization_ratio vector. The randomization ratio control to treatment.
#' @param data NULL. stores the randomization scheme function
#'
#' @return a list with randomization details (block size and ratio).
#'
#' @examples
#' randomize(block_size = 100, randomization_ratio = c(2, 3))
#' randomize(block_size = 10, randomization_ratio = c(1, 4))
#' @export randomize
randomize <- function(block_size, randomization_ratio, data = NULL){
  data$block <- block_size
  data$rand_ratio <- randomization_ratio
  data
}




#' @title Hypothesis wrapper
#'
#' @description Wrapper function for the hypothesis in the trial
#'
#' @param delta numeric. threshold set for margin in null hypothesis
#' @param futility_prob numeric. probability of futility
#' @param prob_ha numeric. posterior probability of accepting alternative.
#' @param expected_success_prob numeric. probability of expecting success.
#' @param data NULL. stores the randomization scheme function
#'
#' @return a list with information of hypothesis testing (threshold, futility probability,
#' probability of alternative and probability of expected success.)
#'
#' @examples
#' hypothesis(delta = 0, futility_prob = 0.05, prob_ha = 0.95, expected_success_prob = 0.90)
#' hypothesis(delta= 0.2, futility_prob = 0.1, prob_ha = 0.975, expected_success_prob = 0.80)
#' @export hypothesis
hypothesis <- function(delta, futility_prob, prob_ha, expected_success_prob, data = NULL){
  data$h0 <- delta
  data$futility_prob  <- futility_prob
  data$prob_ha <- prob_ha
  data$expected_success_prob <- expected_success_prob
  data
}

