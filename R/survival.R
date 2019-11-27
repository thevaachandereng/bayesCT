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
#' @param event0 vector. Historical status indicator, normally 0=alive, 1=dead.
#'   Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death).
#'   For censored data, the status indicator is 0=right censored, 1 = event at time.
#'   Although unusual, the event indicator can be omitted, in which case
#'   all subjects are assumed to have an event.
#' @param treatment0 vector. the historical treatment assignment for patients,
#'    1 for treatment group and 0 for control group.
#' @param prior scalar. Prior value for the gamma shape of the piecewise
#'   exponential hazards. Default is c(0.1, 0.1).
#' @param surv_time scalar. scalar. Survival time of interest for computing the probability
#'    of survival for a single arm (OPC) trial. Default is overall, i.e.,
#'    current+historical, median survival time.
#' @param breaks vector. Breaks (interval starts) used to compose the breaks of the piecewise
#'     exponential model. Do not include zero. Default breaks are the quantiles of the input
#'     times.
#' @inheritParams binomial_analysis
#' @inheritParams pw_exp_sim
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
#'   \item{\code{alpha_max}}{
#'     scalar. The alpha_max input. }
#'   \item{\code{N_treatment}}{
#'     scalar. The number of patients enrolled in the experimental group for
#'     each simulation.}
#'   \item{\code{event_treatment}}{
#'     scalar. The number of events in the experimental group for
#'     each simulation.}
#'   \item{\code{N_control}}{
#'     scalar. The number of patients enrolled in the control group for
#'     each simulation.}
#'   \item{\code{event_control}}{
#'     scalar. The number of events in the control group for
#'     each simulation.}
#'   \item{\code{N_enrolled}}{
#'     scalar. The number of patients enrolled in the trial (sum of control
#'     and experimental group for each simulation. )}
#'   \item{\code{N_complete}}{
#'     scalar. The number of patients whose time passes the surv_time.}
#'   \item{\code{alpha_discount}}{
#'     vector. The alpha discount funtion used for control and treatment.}
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
#' @importFrom stats quantile
#' @importFrom survival Surv
#' @importFrom bayesDP bdpsurvival
#' @importFrom dplyr mutate bind_cols
#'
#' @export survival_analysis

survival_analysis <- function(
  time,
  treatment,
  event                 = NULL,
  time0                 = NULL,
  treatment0            = NULL,
  event0                = NULL,
  surv_time             = NULL,
  h0                    = 0,
  breaks                = NULL,
  alternative           = "greater",
  N_impute              = 10,
  number_mcmc           = 10000,
  prob_ha               = 0.95,
  futility_prob         = 0.10,
  expected_success_prob = 0.90,
  prior                 = c(.1, .1),
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

  # combining the histrotical data
  if(!is.null(time0)){
    data0 <- data.frame(cbind(time = time0, event = event0, treatment = treatment0))
  }
  else{
    data0 <- NULL
  }

  if(is.null(surv_time)){
    surv_time <- median(c(time, time0))
  }


  #reading the data
  data_total <- data.frame(cbind(time, event, treatment))

  # added interim data for incomplete data
  data_interim <- data_total %>%
    mutate(futility = (time < surv_time & event == 0))

  # analyze the data using bayesDp
  post <- bdpsurvival(formula      = Surv(time, event) ~ treatment,
                      data         = data_total,
                      data0        = data0,
                      fix_alpha    = TRUE,
                      number_mcmc  = number_mcmc,
                      breaks       = breaks,
                      method       = "fixed")


  # assigning stop_futility and expected success
  stop_futility         <- 0
  stop_expected_success <- 0
  expected_success_test <- 0

  for(i in 1:N_impute){

    control_impute <-  data_interim %>%
      filter(treatment == 0 & futility)

    impute_control <- pw_exp_impute(time      = control_impute$time,
                                    hazard    = post$posterior_control$posterior_flat_hazard[i, ],
                                    maxtime   = surv_time,
                                    cutpoint  = post$args1$breaks)

    data_control_success_impute <- data_interim %>%
      filter(treatment == 0 & futility) %>%
      bind_cols(time_impute  = impute_control$time,
                event_impute = impute_control$event)


    treatment_impute <-  data_interim %>%
      filter(treatment == 1 & futility)

    impute_treatment <- pw_exp_impute(time      = treatment_impute$time,
                                      hazard    = post$posterior_treatment$posterior_flat_hazard[i, ],
                                      maxtime   = surv_time,
                                      cutpoint  = post$args1$breaks)

    data_treatment_success_impute <- data_interim %>%
      filter(treatment == 1 & futility) %>%
      bind_cols(time_impute  = impute_treatment$time,
                event_impute = impute_treatment$event)


    data_noimpute <- data_interim %>%
      filter(!futility) %>%
      mutate(time_impute = time, event_impute = event)


    # combine the treatment and control imputed datasets
    data_success_impute <- bind_rows(data_control_success_impute,
                                     data_treatment_success_impute,
                                     data_noimpute) %>%
      mutate(time = time_impute, event = event_impute) %>%
      select(-c(time_impute, event_impute))


    # analyze complete+imputed data using discount funtion via binomial
    post_imp <-bdpsurvival(formula      = Surv(time, event) ~ treatment,
                           data         = data_success_impute,
                           data0        = data0,
                           fix_alpha    = TRUE,
                           number_mcmc  = number_mcmc,
                           breaks       = breaks,
                           method       = "fixed")

    # Posterior effect size: test vs control or treatment itself
    if(sum(data_success_impute$treatment == 0) != 0){
      effect_imp <- post_imp$final$posterior_loghazard
      if(alternative == "two-sided"){
        success  <- max(c(mean(effect_imp > h0), mean(-effect_imp > h0)))
      }
      else if(alternative == "greater"){
        success  <- mean(effect_imp > h0)
      }
      else{
        success  <- mean(-effect_imp > h0)
      }
    }

    else{
      effect <- post_imp$posterior_treatment$posterior_survival
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

  data_final <- data_total


  # Analyze complete data using discount funtion via binomial
  post_final <- bdpsurvival(formula      = Surv(time, event) ~ treatment,
                            data         = data_final,
                            data0        = data0,
                            fix_alpha    = TRUE,
                            number_mcmc  = number_mcmc,
                            breaks       = breaks,
                            method       = "fixed")

  ### Format and output results
  # Posterior effect size: test vs control or treatment itself
  if(sum(data_final$treatment == 0) != 0){
    effect <- post_final$final$posterior_loghazard
    if(alternative == "two-sided"){
      post_paa <- max(c(mean(effect > h0), mean(-effect > h0)))
    }
    else if(alternative == "greater"){
      post_paa <- mean(effect > h0)
    }
    else{
      post_paa <- mean(-effect > h0)
    }
  }

  else{
    effect <- post_final$posterior_treatment$posterior_survival
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

  N_treatment           <- sum(data_final$treatment)         # Total sample size analyzed - test group
  N_control             <- sum(!data_final$treatment)        # Total sample size analyzed - control group
  N_enrolled            <- N_treatment + N_control
  N_complete            <- sum(!data_interim$futility)
  alpha_discount        <- c(post_final$posterior_control$alpha_discount, post_final$posterior_treatment$alpha_discount)
  event_treatment       <- sum(data_final$event[!!data_final$treatment])
  event_control         <- sum(data_final$event[!data_final$treatment])

  ## output
  results_list <- list(
    prob_of_accepting_alternative              = prob_ha,
    margin                                     = h0,                       # margin for error
    alternative                                = alternative,              # alternative hypothesis
    alpha_max                                  = alpha_max,
    N_treatment                                = N_treatment,
    event_treatment                            = event_treatment,
    N_control                                  = N_control,
    event_control                              = event_control,
    N_enrolled                                 = N_treatment + N_control,
    N_complete                                 = N_complete,
    alpha_discount                             = alpha_discount,
    post_prob_accept_alternative               = post_paa,                 # Posterior probability that alternative hypothesis is true
    est_final                                  = mean(effect),             # Posterior Mean of treatment effect
    stop_futility                              = stop_futility,            # Did the trial stop for futility
    stop_expected_success                      = stop_expected_success     # Did the trial stop for expected success
  )

  return(results_list)

}

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("time", "event", "event_impute",
                                                        "time_impute", "id", "futility", "treatment",
                                                        "subject_impute_success",
                                                        "subject_impute_futility", "p_outcome"))
