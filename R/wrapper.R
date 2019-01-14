

#' @title Details of the clinical study
#'
#' @description Wrapper function for details of the clinical trial simulation.
#'
#' @param total_sample_size integer. The number of sample size needed.
#' @param prop_loss_to_followup integer. The proportion of loss to follow up.
#' @param study_period integer. The length of the study.
#' @param interim_look vector. Vector with interim looks.
#' @param .data NULL. This should not be changed by the user.
#'
#' @return a list with sample size, length of the study, interim looks and proportion loss to follow up
#'
#' @examples study_details(total_sample_size = 300, study_period = 50, interim_look = c(210, 240, 270))
#' @export study_details
study_details <- function(total_sample_size, study_period, interim_look = NULL,
                          prop_loss_to_followup = 0.10, .data = NULL){
  .data$N_total <- total_sample_size
  .data$EndofStudy <- study_period
  .data$interim_look <- interim_look
  .data$prop_loss_to_followup <- prop_loss_to_followup
  .data
}


#' @title Enrollment rate wrapper
#'
#' @description Wrapper function for enrollment rate.
#'
#' @param lambda vector. Vector with different enrollment rate parameters.
#' @param time vector. Vector with different cut-off times for lambda.
#' @inheritParams study_details
#'
#' @return a list with enrollment rate information
#'
#' @examples enrollment_rate(lambda = c(0.3, 1), time = 25)
#' @export enrollment_rate
enrollment_rate <- function(lambda = 0.3, time = NULL, .data = NULL){
  .data$lambda <- lambda
  .data$lambda_time <- time
  .data
}


#' @title Imputation wrapper
#'
#' @description Wrapper function for no_of_impute.
#'
#' @param no_of_impute integer. Number of Monte Carlo imputation for missing
#'   data.
#' @param number_mcmc scalar. Number of Monte Carlo Markov Chain draws from
#'   posterior distribution.
#' @inheritParams study_details
#'
#' @return a list with number of imputation
#'
#' @examples impute(no_of_impute = 100, number_mcmc = 1000)
#' @export impute
impute <- function(no_of_impute = 10000,
                   number_mcmc  = 10000,
                   .data = NULL){
  .data$N_impute     <- no_of_impute
  .data$number_mcmc  <- number_mcmc
  .data
}


#' @title Randomization scheme wrapper
#'
#' @description Wrapper function for the randomization scheme in the trial.
#'
#' @param block_size integer. Block size for the complete randomization in a
#'   block.
#' @param randomization_ratio vector. The randomization allocation for control to
#'   treatment.
#' @inheritParams study_details
#'
#' @return a list with randomization details (block size and ratio).
#'
#' @examples
#' randomize(block_size = 100, randomization_ratio = c(2, 3))
#' randomize(block_size = 10, randomization_ratio = c(1, 4))
#' @export randomize
randomize <- function(block_size            = 2,
                      randomization_ratio   = c(1, 1),
                      .data = NULL){
  .data$block <- block_size
  .data$rand_ratio <- randomization_ratio
  .data
}


#' @title Hypothesis wrapper
#'
#' @description Wrapper function for the hypothesis in the trial.
#'
#' @param delta numeric. Threshold set for margin in null hypothesis. The default
#'   is set to 0.
#' @param futility_prob numeric. Probability of futility. The default is 0.05.
#' @param prob_accept_ha numeric. Posterior probability of accepting alternative
#'   hypothesis. The default is 0.95.
#' @param expected_success_prob numeric. Probability of expected success.
#' @param alternative character. The string specifying the alternative hypothesis,
#'        must be one of \code{"greater"} (default), \code{"less"} or
#'        \code{"two.sided"}.
#' @inheritParams study_details
#'
#' @return a list with information of hypothesis testing (threshold, futility
#'   probability, probability of accepting the alternative hypothesis, and probability of
#'   expected success).
#'
#' @examples
#' hypothesis(delta = 0, futility_prob = 0.05, prob_accept_ha = 0.95,
#'            expected_success_prob = 0.90, alternative = "greater")
#' hypothesis(delta= 0.2, futility_prob = 0.1, prob_accept_ha = 0.975,
#'            expected_success_prob = 0.80, alternative = "less")
#' @export hypothesis
hypothesis <- function(delta                 = 0,
                       futility_prob         = 0.05,
                       prob_accept_ha        = 0.95,
                       expected_success_prob = 0.90,
                       alternative           = "greater",
                       .data                 = NULL){
  .data$h0                     <- delta
  .data$futility_prob          <- futility_prob
  .data$prob_ha                <- prob_accept_ha
  .data$expected_success_prob  <- expected_success_prob
  .data$alternative            <- alternative
  .data
}

