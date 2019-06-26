#' @title Analyzing bayesian trial for time-to-event data
#'
#' @description Function to analyze bayesian trial for time-to-event data
#'  which allows early stopping and incorporation of historical data using
#'  the discount function approach
#'
#' @param time vector. exposure time for the subjects. It must be the same length
#'    as the treatment variable.
#' @param event vector. The status indicator, normally 0=alive, 1=dead.
#'   Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death).
#'   For censored data, the status indicator is 0=right censored, 1 = event at time.
#'   Although unusual, the event indicator can be omitted, in which case
#'   all subjects are assumed to have an event.
#' @param time0 vector. Historical exposure time for the subjects. It must be the same length
#'    as the treatment variable.
#' @param event vector. Hisc status indicator, normally 0=alive, 1=dead.
#'   Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death).
#'   For censored data, the status indicator is 0=right censored, 1 = event at time.
#'   Although unusual, the event indicator can be omitted, in which case
#'   all subjects are assumed to have an event.
#' @param prior scalar. Prior value for the gamma shape of the piecewise
#'   exponential hazards. Default is c(0.1, 0.1).
#' @inheritParams binomial_analysis
#'
#' @importFrom survival Surv
#' @importFrom dplyr mutate filter group_by bind_rows select n summarize
#' @importFrom bayesDP bdpsurvival
#'
#' @return a list of output for the bayesian trial for time-to-event.
#'
#' \describe{
#'   \item{\code{prob_of_accepting_alternative}}{
#'     scalar. The input parameter of probability of accepting the alternative.}
#'   \item{\code{margin}}{
#'     scalar. The margin input value of difference between mean estimate of treatment
#'      and mean estimate of the control.}
#'   \item{\code{alternative}}{
#'     character. The input parameter of alternative hypothesis. }
#'   \item{\code{N_treatment}}{
#'     scalar. The number of patients enrolled in the experimental group for
#'     each simulation.}
#'   \item{\code{N_control}}{
#'     scalar. The number of patients enrolled in the control group for
#'     each simulation.}
#'   \item{\code{N_enrolled}}{
#'     vector. The number of patients enrolled in the trial (sum of control
#'     and experimental group for each simulation. )}
#'   \item{\code{N_complete}}{
#'     scalar. The number of patients who completed the trial and had no
#'     loss to follow-up.}
#'   \item{\code{post_prob_accept_alternative}}{
#'     vector. The final probability of accepting the alternative
#'     hypothesis after the analysis is done.}
#'   \item{\code{est_final}}{
#'     scalar. The final estimate of the difference in posterior estimate of
#'     treatment and posterior estimate of the control group.}
#'   \item{\code{stop_futility}}{
#'     scalar. Did the trial stop for futility during imputation of patient
#'     who had loss to follow up? 1 for yes and 0 for no.}
#'   \item{\code{stop_expected_success}}{
#'     scalar. Did the trial stop for early success during imputation of patient
#'     who had loss to follow up? 1 for yes and 0 for no.}
#'
#' }
#'
#' @export survival_analysis

survival_analysis <- function(
  treatment,
  time,
  event                 = NULL,
  treatment_0           = NULL,
  time_0                = NULL,
  event_0               = NULL,
  surv_time             = NULL,
  alternative           = "greater",
  N_impute              = 10,
  h0                    = 0,
  number_mcmc           = 10000,
  prob_ha               = 0.95,
  futility_prob         = 0.10,
  expected_success_prob = 0.90,
  prior                 = c(1, 1),
  discount_function     = "identity",
  fix_alpha             = FALSE,
  alpha_max             = 1,
  weibull_scale         = 0.135,
  weibull_shape         = 3
){

  #if complete is NULL, assume the data is complete
  if(is.null(event)){
    event <- rep(1, length(time))
  }

  if(!is.null(time0)){
    data0 <- data.frame(cbind(time0, event0, treatment0))
  }

  #reading the data
  data_total <- data.frame(cbind(time, event, treatment))



  # analyze the data using bayesDp
  post <- bdpbinomial(y_t                    = sum(data$outcome[data$treatment == 1]),
                      N_t                    = length(data$outcome[data$treatment == 1]),
                      y_c                    = y_c,
                      N_c                    = N_c,
                      y0_t                   = y0_treatment,
                      N0_t                   = N0_treatment,
                      y0_c                   = y0_control,
                      N0_c                   = N0_control,
                      discount_function      = discount_function,
                      number_mcmc            = number_mcmc,
                      a0                     = prior[1],
                      b0                     = prior[2],
                      alpha_max              = alpha_max,
                      fix_alpha              = fix_alpha,
                      weibull_scale          = weibull_scale,
                      weibull_shape          = weibull_shape)


  # assigning stop_futility and expected success
  stop_futility         <- 0
  stop_expected_success <- 0
  expected_success_test <- 0

  for(i in 1:N_impute){
    data_control_success_impute <- data_interim %>%
      filter(treatment == 0) %>%
      mutate(outcome_impute = ifelse(futility,
                                     rbinom(n(), 1, prop$p_outcome[1]),
                                     outcome))
    # imputing success for treatment group
    data_treatment_success_impute  <- data_interim %>%
      filter(treatment == 1) %>%
      mutate(outcome_impute = ifelse(futility,
                                     rbinom(n(), 1, prop$p_outcome[2]),
                                     outcome))

    # combine the treatment and control imputed datasets
    data_success_impute <- bind_rows(data_control_success_impute,
                                     data_treatment_success_impute) %>%
      mutate(outcome = outcome_impute) %>%
      select(-outcome_impute)

    # Create enrolled subject data frame for discount function analysis
    data <- data_success_impute

    # assigning input for control arm given it is a single or double arm
    if(sum(data$treatment == 0) != 0){
      y_c <- sum(data$outcome[data$treatment == 0])
      N_c <- length(data$outcome[data$treatment == 0])
    }
    else{
      y_c <- NULL
      N_c <- NULL
    }

    # analyze complete+imputed data using discount funtion via binomial
    post_imp <- bdpbinomial(y_t                    = sum(data$outcome[data$treatment == 1]),
                            N_t                    = length(data$outcome[data$treatment == 1]),
                            y_c                    = y_c,
                            N_c                    = N_c,
                            y0_t                   = y0_treatment,
                            N0_t                   = N0_treatment,
                            y0_c                   = y0_control,
                            N0_c                   = N0_control,
                            discount_function      = discount_function,
                            number_mcmc            = number_mcmc,
                            a0                     = prior[1],
                            b0                     = prior[2],
                            alpha_max              = alpha_max,
                            fix_alpha              = fix_alpha,
                            weibull_scale          = weibull_scale,
                            weibull_shape          = weibull_shape)

    if(sum(data$treatment == 0) != 0){
      if(alternative == "two-sided"){
        effect_imp <- post_imp$posterior_treatment$posterior - post_imp$posterior_control$posterior
        success <- max(c(mean(effect_imp > h0), mean(-effect_imp > h0)))
      }
      else if(alternative == "greater"){
        effect_imp <- post_imp$posterior_treatment$posterior - post_imp$posterior_control$posterior
        success <- mean(effect_imp > h0)
      }
      else{
        effect_imp <- post_imp$posterior_treatment$posterior - post_imp$posterior_control$posterior
        success <- mean(-effect_imp > h0)
      }
    }

    else{
      effect_imp <- post_imp$final$posterior
      if(alternative == "two-sided"){
        success <- max(c(mean(effect_imp > h0), mean(effect_imp < h0)))
      }
      else if(alternative == "greater"){
        success <- mean(effect_imp > h0)
      }
      else{
        success <- mean(effect_imp < h0)
      }
    }

    if(success > prob_ha){
      expected_success_test <- expected_success_test + 1
    }

  }

  if(expected_success_test / N_impute < futility_prob){
    stop_futility       <- 1
  }

  # Test if expected success criteria met
  if(expected_success_test / N_impute > expected_success_prob ){
    stop_expected_success <- 1
  }

  data_final <- data_interim %>%
    filter(!futility)

  if(sum(data_final$treatment == 0) != 0){
    y_c <- sum(data_final$outcome[data_final$treatment == 0])
    N_c <- length(data_final$outcome[data_final$treatment == 0])
  }
  else{
    y_c <- NULL
    N_c <- NULL
  }

  # Analyze complete data using discount funtion via binomial
  post_final <- bdpbinomial(y_t                  = sum(data_final$outcome[data_final$treatment == 1]),
                            N_t                  = length(data_final$outcome[data_final$treatment == 1]),
                            y_c                  = y_c,
                            N_c                  = N_c,
                            y0_t                 = y0_treatment,
                            N0_t                 = N0_treatment,
                            y0_c                 = y0_control,
                            N0_c                 = N0_control,
                            number_mcmc          = number_mcmc,
                            discount_function    = discount_function,
                            a0                   = prior[1],
                            b0                   = prior[2],
                            alpha_max            = alpha_max,
                            fix_alpha            = fix_alpha,
                            weibull_scale        = weibull_scale,
                            weibull_shape        = weibull_shape)

  ### Format and output results
  # Posterior effect size: test vs control or treatment itself
  if(sum(data_final$treatment == 0) != 0){
    if(alternative == "two-sided"){
      effect <- post_final$posterior_treatment$posterior - post_final$posterior_control$posterior
      post_paa <- max(c(mean(effect > h0), mean(-effect > h0)))
    }
    else if(alternative == "greater"){
      effect <- post_final$posterior_treatment$posterior - post_final$posterior_control$posterior
      post_paa <- mean(effect > h0)
    }
    else{
      effect <- post_final$posterior_treatment$posterior - post_final$posterior_control$posterior
      post_paa <- mean(-effect > h0)
    }
  }

  else{
    effect <- post_final$final$posterior
    if(alternative == "two-sided"){
      post_paa <- max(c(mean(effect > h0), mean(effect < h0)))
    }
    else if(alternative == "greater"){
      post_paa <- mean(effect > h0)
    }
    else{
      post_paa <- mean(effect < h0)
    }
  }

  N_treatment  <- sum(data_final$treatment)         # Total sample size analyzed - test group
  N_control    <- sum(!data_final$treatment)        # Total sample size analyzed - control group
  N_enrolled   <- dim(data_total)[1]

  #estimating prop
  prop <- data %>%
    group_by(treatment) %>%
    summarize(p_outcome = mean(outcome))

  ## output
  results_list <- list(
    prob_of_accepting_alternative              = prob_ha,
    margin                                     = h0,                       # margin for error
    alternative                                = alternative,              # alternative hypothesis
    N_treatment                                = N_treatment,
    N_control                                  = N_control,
    N_complete                                 = N_treatment + N_control,
    N_enrolled                                 = N_enrolled,               # Total sample size enrolled when trial stopped
    post_prob_accept_alternative               = post_paa,                 # Posterior probability that alternative hypothesis is true
    est_final                                  = mean(effect),             # Posterior Mean of treatment effect
    stop_futility                              = stop_futility,            # Did the trial stop for futility
    stop_expected_success                      = stop_expected_success     # Did the trial stop for expected success
    #MLE_est                                   = MLE$coe[2],               # Treatment effect useing MLE
    #MLE_est_interim                           = MLE_int$coe[2]            # Treatment effect useing MLE at interim analysis
  )

  return(results_list)

}

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("complete", "outcome", "outcome_impute", "id",
                                                        "futility", "treatment",
                                                        "subject_impute_success",
                                                        "subject_impute_futility", "p_outcome"))
