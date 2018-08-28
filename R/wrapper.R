

#' @title Details of the clinical study
#'
#' @description Wrapper function for details of the clinical study simulation.
#'
#' @param total_sample_size integer. The number of sample size needed.
#' @param study_period integer. The length of the study.
#' @param interim_look vector. Vector with interim looks.
#' @param .data NULL. stores the sample size and length of study.
#'
#' @return a list with sample size and length of the study and interim looks.
#'
#' @examples study_details(total_sample_size = 300, study_period = 50, interim_look = c(210, 240, 270))
#' @export study_details
study_details <- function(total_sample_size, study_period, interim_look, .data = NULL){
  .data$N_total <- total_sample_size
  .data$EndofStudy <- study_period
  .data$interim_look <- interim_look
  .data
}



#' @title Enrollment rate wrapper
#'
#' @description Wrapper function for enrollment rate
#'
#' @param lambda vector. Vector with different enrollment rate.
#' @param time vector. Vector with different cutoff time for lambda.
#' @param .data NULL. stores the impute function
#'
#' @return a list with enrollment rate information
#'
#' @examples enrollment_rate(lambda = c(0.3, 1), time = 25)
#' @export enrollment_rate
enrollment_rate <- function(lambda = NULL, time = NULL, .data = NULL){
  .data$lambda <- lambda
  .data$lambda_time <- time
  .data
}


#' @title Imputation wrapper
#'
#' @description Wrapper function for no_of _impute
#'
#' @param no_of_impute integer. Number of Monte Carlo imputation for missing data
#' @param number_mcmc scalar. Number of Monte Carlo Markov Chain in sampling posterior.
#' @param .data NULL. stores the impute function
#'
#' @return a list with number of imputation
#'
#' @examples impute(no_of_impute = 100, number_mcmc = 1000)
#' @export impute
impute <- function(no_of_impute, number_mcmc, .data = NULL){
  .data$N_impute     <- no_of_impute
  .data$number_mcmc  <- number_mcmc
  .data
}


#' @title Randomization scheme wrapper
#'
#' @description Wrapper function for the randomization scheme in the trial
#'
#' @param block_size integer. Block size for the complete randomization in a block
#' @param randomization_ratio vector. The randomization ratio control to treatment.
#' @param .data NULL. stores the randomization scheme function
#'
#' @return a list with randomization details (block size and ratio).
#'
#' @examples
#' randomize(block_size = 100, randomization_ratio = c(2, 3))
#' randomize(block_size = 10, randomization_ratio = c(1, 4))
#' @export randomize
randomize <- function(block_size, randomization_ratio, .data = NULL){
  .data$block <- block_size
  .data$rand_ratio <- randomization_ratio
  .data
}




#' @title Hypothesis wrapper
#'
#' @description Wrapper function for the hypothesis in the trial
#'
#' @param delta numeric. threshold set for margin in null hypothesis
#' @param futility_prob numeric. probability of futility
#' @param prob_ha numeric. posterior probability of accepting alternative.
#' @param expected_success_prob numeric. probability of expecting success.
#' @param .data NULL. stores the randomization scheme function
#'
#' @return a list with information of hypothesis testing (threshold, futility probability,
#' probability of alternative and probability of expected success.)
#'
#' @examples
#' hypothesis(delta = 0, futility_prob = 0.05, prob_ha = 0.95, expected_success_prob = 0.90)
#' hypothesis(delta= 0.2, futility_prob = 0.1, prob_ha = 0.975, expected_success_prob = 0.80)
#' @export hypothesis
hypothesis <- function(delta, futility_prob, prob_ha, expected_success_prob, .data = NULL){
  .data$h0                     <- delta
  .data$futility_prob          <- futility_prob
  .data$prob_ha                <- prob_ha
  .data$expected_success_prob  <- expected_success_prob
  .data
}
