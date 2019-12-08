#' @title Time-to-event outcome for Bayesian Adaptive Trials
#'
#' @description Simulation for time-to-event outcome for Bayesian Adaptive trial with
#'   different inputs to control for power, sample size, type 1 error rate, etc.
#' @param hazard_treatment vector. Constant hazard rates under the treatment arm.
#' @param hazard_control vector. Constant hazard rates under the control arm.
#' @inheritParams pw_exp_sim
#' @inheritParams normalBACT
#' @inheritParams survival_analysis
#'
#' @return a list of output for a single trial simulation.
#' \describe{
#'   \item{\code{lambda_treatment}}{
#'     vector. The input parameter of constant hazard rates in the
#'     treatment group.}
#'   \item{\code{cutpoint_treatment}}{
#'     vector. The change-point vector when the constant hazard rate(s) changes for
#'     the treatment group.}
#'   \item{\code{lambda_control}}{
#'     vector. The input parameter of constant hazard rates in the
#'     control group.}
#'   \item{\code{cutpoint_control}}{
#'     vector. The change-point vector when the constant hazard rate(s) changes for
#'     the control group.}
#'   \item{\code{prob_of_accepting_alternative}}{
#'     scalar. The input parameter of probability threshold of accepting the
#'     alternative.}
#'   \item{\code{margin}}{
#'     scalar. The margin input value of difference between mean estimate of treatment
#'      and mean estimate of the control.}
#'   \item{\code{alternative}}{
#'     character. The input parameter of alternative hypothesis. }
#'   \item{\code{interim_look}}{
#'     vector. The sample size for each interim look.}
#'   \item{\code{N_treatment}}{
#'     scalar. The number of patients enrolled in the treatment group for
#'     each simulation.}
#'   \item{\code{event_treatment}}{
#'     scalar. The number of events in the treatment group for
#'     each simulation.}
#'   \item{\code{N_control}}{
#'     scalar. The number of patients enrolled in the control group for
#'     each simulation.}
#'   \item{\code{event_control}}{
#'     scalar. The number of events in the control group for
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
#'   \item{\code{est_interim}}{
#'     scalar. The interim estimate of the difference in posterior estimate of
#'     treatment and posterior estimate of the control group.}
#' }
#'
#' @examples survivalBACT(0.01, NULL, 0.014, 200, EndofStudy = 100)
#'
#'
#' @importFrom stats runif
#' @importFrom dplyr mutate filter group_by bind_cols bind_rows select ungroup
#' @importFrom bayesDP bdpsurvival
#' @export survivalBACT
#'
survivalBACT <- function(
  hazard_treatment,
  cutpoint              = NULL,
  hazard_control        = NULL,
  N_total,
  breaks                = NULL,
  time0                 = NULL,
  treatment0            = NULL,
  event0                = NULL,
  lambda                = 0.3,
  lambda_time           = NULL,
  interim_look          = NULL,
  EndofStudy,
  prior                 = c(.1, .1),
  block                 = 2,             # block size for randomization
  rand_ratio            = c(1, 1),       # randomization ratio in control to treatament (default 1:1)
  prop_loss_to_followup = 0.10,          # Proportion of loss in data
  alternative           = "greater",     # the alternative hypothesis (either two-sided, greater, less)
  h0                    = 0,             # Null hypothesis value
  futility_prob         = 0.05,          # Futility probability
  expected_success_prob = 0.9,           # Expected success probability
  prob_ha               = 0.95,          # Posterior probability of accepting alternative hypothesis
  N_impute              = 10,            # Number of imputation simulations for predictive distribution
  number_mcmc           = 10000,         # Number of posterior sampling
  discount_function     = "identity",
  alpha_max             = 1,             # max weight on incorporating historical data
  fix_alpha             = FALSE,         # fix alpha set weight of historical data to alpha_max
  weibull_scale         = 0.135,         # weibull parameter
  weibull_shape         = 3,             # weibull parameter
  method                = "fixed"
){


  # checking interim_look
  if(!is.null(interim_look)){
    stopifnot(all(N_total > interim_look))
  }

  # combining the histrotical data
  if(!is.null(time0)){
    data0 <- data.frame(cbind(time = time0, event = event0, treatment = treatment0))
  }
  else{
    data0 <- NULL
  }

  #checking if alternative is right
  if(alternative != "two-sided" & alternative != "greater" & alternative != "less"){
    stop("The input for alternative is wrong!")
  }

  # if cutpoint is not NULL, assign breaks to cutpoint
  if(!is.null(cutpoint)){
    breaks <- cutpoint
  }

  ## make sure breaks is not more than Endofstudy
  if(!is.null(breaks)){
    stopifnot(any(breaks < EndofStudy))
  }

  if(is.null(hazard_control) & h0 == 0){
    h0 <- 0.5
  }

  # assigining interim look and final look
  analysis_at_enrollnumber <- c(interim_look, N_total)

  # assignment of enrollment based on the enrollment function
  enrollment <- enrollment(param = lambda, N_total = N_total, time = lambda_time)

  # simulating group and treatment group assignment
  if(!is.null(hazard_control)){
    group <- randomization(N_total = N_total, block = block, allocation = rand_ratio)
  }
  else{
    group <- rep(1, N_total)
  }

  time  <- rep(NA, length = N_total)
  event <- rep(NA, length = N_total)
  # simulate binomial outcome
  if(!is.null(hazard_control)){
    sim_control          <- pw_exp_sim(hazard    = hazard_control,
                                       n         = sum(!group),
                                       maxtime   = EndofStudy,
                                       cutpoint  = cutpoint)
    time[which(!group)]  <- sim_control$time
    event[which(!group)] <- sim_control$event
  }

  sim_treatment        <- pw_exp_sim(hazard   = hazard_treatment,
                                     n        = sum(group),
                                     maxtime  = EndofStudy,
                                     cutpoint = cutpoint)
  time[which(!!group)] <- sim_treatment$time
  event[which(!!group)] <- sim_treatment$event

  if(is.null(breaks)){
    breaks <- mean(time)
  }

  # simulate loss to follow-up
  n_loss_to_fu <- ceiling(prop_loss_to_followup * N_total)
  loss_to_fu   <- rep(FALSE, N_total)
  loss_to_fu[sample(1:N_total, n_loss_to_fu)] <- TRUE

  # creating a new data.frame for all the variables
  data_total <- data.frame(
    time       = time,
    treatment  = group,
    event      = event,
    enrollment = enrollment,
    id         = 1:N_total,
    loss_to_fu = loss_to_fu)

  ## patient lost are uniformly distributed
  data_total$time[data_total$loss_to_fu]  <-
    runif(n_loss_to_fu, 0, data_total$time[data_total$loss_to_fu])
  data_total$event[data_total$loss_to_fu] <- rep(0, n_loss_to_fu)


  # assigning stop_futility and expected success
  stop_futility         <- 0
  stop_expected_success <- 0

  if(length(analysis_at_enrollnumber) > 1){
    for(i in 1:(length(analysis_at_enrollnumber) - 1)){

      # Analysis at the `analysis_at_enrollnumber` look
      # Indicators for subject type:
      # - subject_enrolled: subject has data present in the current look
      # - subject_impute_success: subject has data present in the current look but has not
      #                           reached end of study or subject is loss to follow-up: needs BP change imputed
      # - subject_impute_futility: subject has no data present in the current look;
      #                            needs baseline BP and BP change imputed
      data_interim <- data_total %>%
        mutate(subject_enrolled = id <= analysis_at_enrollnumber[i],
               subject_impute_futility = !subject_enrolled) %>%
        group_by(subject_enrolled) %>%
        mutate(subject_impute_success =
                 event * (enrollment[analysis_at_enrollnumber[i]] - enrollment <= EndofStudy & subject_enrolled) |
                 (1 - event) * (enrollment[analysis_at_enrollnumber[i]] - enrollment <= time & subject_enrolled) |
                 (subject_enrolled & loss_to_fu)) %>%
        mutate(real_time = ifelse(event * (enrollment[analysis_at_enrollnumber[i]] - enrollment <= EndofStudy & subject_enrolled) |
                                    (1 - event) * (enrollment[analysis_at_enrollnumber[i]] - enrollment <= time & subject_enrolled),
                                  enrollment[analysis_at_enrollnumber[i]] - enrollment + sd(time) / 10000, time)) %>%
        mutate(real_event = ifelse(event * (enrollment[analysis_at_enrollnumber[i]] - enrollment <= EndofStudy & subject_enrolled) |
                                     (1 - event) * (enrollment[analysis_at_enrollnumber[i]] - enrollment <= time & subject_enrolled),
                                   0, event))

      data_interim <- data_interim %>%
        mutate(time = real_time, event = real_event) %>%
        select(-c(real_time, real_event))


      # Carry out interim analysis on patients with complete data only
      # - Set-up `new data` data frame
      data <- data_interim %>%
        filter(subject_enrolled) %>%
        ungroup() %>%
        select(time, event, treatment)


      # Analyze data using discount function via binomial
      post <- suppressWarnings(
        bdpsurvival(formula            = Surv(time, event) ~ treatment,
                    data               = data,
                    data0              = data0,
                    breaks             = breaks,
                    a0                 = prior[1],
                    b0                 = prior[2],
                    surv_time          = EndofStudy,
                    discount_function  = discount_function,
                    alpha_max          = alpha_max,
                    fix_alpha          = fix_alpha,
                    number_mcmc        = number_mcmc,
                    weibull_scale      = weibull_scale,
                    weibull_shape      = weibull_shape,
                    method             = method))


      # Imputation phase futility and expected success - initialize counters
      # for the current imputation phase
      futility_test         <- 0
      expected_success_test <- 0

      for(j in 1:N_impute){
        #imputing the success for control group
        if(!is.null(hazard_control)){
          control_impute <-  data_interim %>%
            filter(treatment == 0 & subject_impute_success)

          impute_control <- pw_exp_impute(time      = control_impute$time,
                                          hazard    = post$posterior_control$posterior_flat_hazard[i, ],
                                          maxtime   = EndofStudy,
                                          cutpoint  = post$args1$breaks)

          data_control_success_impute <- data_interim %>%
            filter(treatment == 0 & subject_impute_success) %>%
            bind_cols(time_impute  = impute_control$time,
                      event_impute = impute_control$event)
        }
        else{
          data_control_success_impute <- NULL
        }

        # imputing success for treatment group

        treatment_impute <-  data_interim %>%
          filter(treatment == 1 & subject_impute_success)

        impute_treatment <- pw_exp_impute(time      = treatment_impute$time,
                                          hazard    = post$posterior_treatment$posterior_flat_hazard[i, ],
                                          maxtime   = EndofStudy,
                                          cutpoint  = post$args1$breaks)

        data_treatment_success_impute <- data_interim %>%
          filter(treatment == 1 & subject_impute_success) %>%
          bind_cols(time_impute  = impute_treatment$time,
                    event_impute = impute_treatment$event)


        # combine the treatment and control imputed datasets
        data_noimpute <- data_interim %>%
          filter(!subject_impute_success) %>%
          mutate(time_impute = time, event_impute = event)


        # combine the treatment and control imputed datasets
        data_success_impute <- bind_rows(data_control_success_impute,
                                         data_treatment_success_impute,
                                         data_noimpute) %>%
          mutate(time = time_impute, event = event_impute) %>%
          select(-c(time_impute, event_impute))

        # create enrolled subject data frame for discount function analysis
        data <- data_success_impute %>%
          filter(subject_enrolled)  %>%
          ungroup() %>%
          select(time, event, treatment)


        # analyze complete+imputed data using discount funtion via binomial
        # analyze complete+imputed data using discount funtion via binomial
        post_imp <- suppressWarnings(
          bdpsurvival(formula            = Surv(time, event) ~ treatment,
                      data               = data,
                      data0              = data0,
                      breaks             = breaks,
                      a0                 = prior[1],
                      b0                 = prior[2],
                      surv_time          = EndofStudy,
                      discount_function  = discount_function,
                      alpha_max          = alpha_max,
                      fix_alpha          = fix_alpha,
                      number_mcmc        = number_mcmc,
                      weibull_scale      = weibull_scale,
                      weibull_shape      = weibull_shape,
                      method             = method))

        # estimation of the posterior effect for difference between test and
        # control - If expected success, add 1 to the counter
        if(sum(data$treatment == 0) != 0){
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
          effect_imp <- post_imp$final$posterior_survival
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

        ##########################################################################
        ### Futility computations
        ##########################################################################

        # For patients not enrolled, impute the outcome
        # imputing the control group
        if(!is.null(hazard_control)){
          control_impute <-  data_success_impute %>%
            filter(treatment == 0 & subject_impute_futility)

          impute_control <- pw_exp_sim(hazard    = post$posterior_control$posterior_flat_hazard[i, ],
                                       n         = nrow(control_impute),
                                       maxtime   = EndofStudy,
                                       cutpoint  = post$args1$breaks)

          data_control_futility_impute <- data_success_impute %>%
            filter(treatment == 0 & subject_impute_futility) %>%
            bind_cols(time_impute  = impute_control$time,
                      event_impute = impute_control$event)
        }
        else{
          data_control_futility_impute <- NULL
        }

        # imputing the treatment group
        treatment_impute <-  data_success_impute %>%
          filter(treatment == 1 & subject_impute_futility)

        impute_treatment <- pw_exp_sim(hazard    = post$posterior_treatment$posterior_flat_hazard[i, ],
                                       n         = nrow(treatment_impute),
                                       maxtime   = EndofStudy,
                                       cutpoint  = post$args1$breaks)

        data_treatment_futility_impute <- data_success_impute %>%
          filter(treatment == 1 & subject_impute_futility) %>%
          bind_cols(time_impute  = impute_treatment$time,
                    event_impute = impute_treatment$event)

        # combine the treatment and control imputed datasets
        data_noimpute_futility <- data_success_impute %>%
          filter(!subject_impute_futility) %>%
          mutate(time_impute = time, event_impute = event)


        # combine the treatment and control imputed datasets
        data_futility_impute <- bind_rows(data_control_futility_impute,
                                          data_treatment_futility_impute,
                                          data_noimpute_futility) %>%
          mutate(time = time_impute, event = event_impute) %>%
          select(-c(time_impute, event_impute))


        # Create enrolled subject data frame for discount function analysis
        data <- data_futility_impute %>%
          ungroup() %>%
          select(time, event, treatment)

        # Analyze complete+imputed data using discount funtion via binomial
        post_imp <- suppressWarnings(
          bdpsurvival(formula            = Surv(time, event) ~ treatment,
                      data               = data,
                      data0              = data0,
                      breaks             = breaks,
                      a0                 = prior[1],
                      b0                 = prior[2],
                      surv_time          = EndofStudy,
                      discount_function  = discount_function,
                      alpha_max          = alpha_max,
                      fix_alpha          = fix_alpha,
                      number_mcmc        = number_mcmc,
                      weibull_scale      = weibull_scale,
                      weibull_shape      = weibull_shape,
                      method             = method))


        # Estimation of the posterior effect for difference between test and
        # control
        if(sum(data$treatment == 0) != 0){
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
          effect_imp <- post_imp$final$posterior_survival
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

        # Increase futility counter by 1 if P(effect_imp < h0) > ha
        if(success > prob_ha){
          futility_test <- futility_test + 1
        }

      }

      # print(analysis_at_enrollnumber[i])

      # Test if expected success criteria met
      if(expected_success_test / N_impute > expected_success_prob ){
        stop_expected_success <- 1
        stage_trial_stopped   <- analysis_at_enrollnumber[i]
        break
      }

      # Test if futility success criteria is met
      if(futility_test / N_impute < futility_prob){
        stop_futility       <- 1
        stage_trial_stopped <- analysis_at_enrollnumber[i]
        break
      }

      # Stop study if at last interim look
      if(analysis_at_enrollnumber[i + 1] == N_total){
        stage_trial_stopped <- analysis_at_enrollnumber[i + 1]
        break
      }

    }

    ##############################################################################
    ### Final analysis
    ##############################################################################

    # Estimation of the posterior of the difference
    if(sum(data$treatment == 0) != 0){
      if(alternative == "two-sided"){
        effect_int <- post$final$posterior_loghazard
      }
      else if(alternative == "greater"){
        effect_int <- post$final$posterior_loghazard
      }
      else{
        effect_int <- post$final$posterior_loghazard
      }
    }

    else{
      effect_int <- post_imp$posterior_treatment$posterior_survival
    }

    # Number of patients enrolled at trial stop
    N_enrolled <- nrow(data_interim[data_interim$id <= stage_trial_stopped, ])
  }

  # assigning stage trial stopped given no interim look
  else{
    N_enrolled            <- N_total
    stage_trial_stopped   <- N_total
    stop_futility         <- 0
    stop_expected_success <- 0
  }

  # All patients that have made it to the end of study
  # - Subset out patients loss to follow-up
  data_final <- data_total %>%
    filter(id <= stage_trial_stopped) %>%
    ungroup() %>%
    select(time, event, treatment)

  # Compute the final MLE for the complete data using GLM function
  # MLE <- glm(Y ~ treatment, data = data_final, family = "binomial")


  # Analyze complete data using discount funtion via binomial
  post_final <- suppressWarnings(
    bdpsurvival(formula            = Surv(time, event) ~ treatment,
                data               = data_final,
                data0              = data0,
                breaks             = breaks,
                a0                 = prior[1],
                b0                 = prior[2],
                surv_time          = EndofStudy,
                discount_function  = discount_function,
                alpha_max          = alpha_max,
                fix_alpha          = fix_alpha,
                number_mcmc        = number_mcmc,
                weibull_scale      = weibull_scale,
                weibull_shape      = weibull_shape,
                method             = method))

  data_final <- data_total %>%
    filter(id <= stage_trial_stopped)

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
    effect <- post_final$final$posterior_survival
    if(alternative == "two-sided"){
      post_paa <- max(c(mean(effect > h0), mean(effect > h0)))
    }
    else if(alternative == "greater"){
      post_paa <- mean(effect > h0)
    }
    else{
      post_paa <- mean(effect < h0)
    }
  }



  N_treatment  <- sum(!!data_final$treatment)         # Total sample size analyzed - test group
  N_control    <- sum(!data_final$treatment)          # Total sample size analyzed - control group

  ## output
  results_list <- list(
    hazard_treatment                           = hazard_treatment,              # probability of treatment in binomial
    hazard_control                             = hazard_control,                # probability of control in binomial
    cutpoint                                   = cutpoint,
    prob_of_accepting_alternative              = prob_ha,
    margin                                     = h0,                       # margin for error
    alternative                                = alternative,              # alternative hypothesis
    interim_look                               = interim_look,             # print interim looks
    N_treatment                                = N_treatment,
    N_control                                  = N_control,
    N_enrolled                                 = N_treatment + N_control,
    N_max                                      = N_total, 				         # Total potential sample size
    post_prob_accept_alternative               = post_paa,                 # Posterior probability that alternative hypothesis is true
    est_final                                  = mean(effect),             # Posterior Mean of treatment effect
    stop_futility                              = stop_futility,            # Did the trial stop for futility
    stop_expected_success                      = stop_expected_success     # Did the trial stop for expected success
  )

  if(length(analysis_at_enrollnumber) > 1){
    results_list                            <- c(results_list,
                                                 est_interim                             = mean(effect_int))           # Posterior Mean of treatment effect at interim analysis

  }

  #return results
  results_list

}

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("time", "event", "event_impute", "real_time",
                                                        "time_impute", "id", "futility", "treatment",
                                                        "subject_impute_success", "real_event",
                                                        "subject_impute_futility"))


#' @title Analyzing Bayesian trial for time-to-event data
#'
#' @description Function to analyze Bayesian trial for time-to-event data
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
#' @param surv_time scalar. scalar. Survival time of interest for computing the probability
#'    of survival for a single arm (OPC) trial. Default is overall, i.e.,
#'    current+historical, median survival time.
#' @param breaks vector. Breaks (interval starts) used to compose the breaks of the piecewise
#'     exponential model. Do not include zero. Default breaks are the quantiles of the input
#'     times.
#' @param prior vector. Prior values of the gamma rate, Gamma(a0, b0). The default is
#'   set to Gamma(.1, .1).
#' @inheritParams binomial_analysis
#' @inheritParams pw_exp_sim
#'
#'
#' @return a list of output for the Bayesian trial for time-to-event.
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
#' }
#'
#' @importFrom dplyr mutate bind_cols
#' @importFrom survival Surv
#' @importFrom dplyr mutate filter group_by bind_rows select n summarize
#' @importFrom bayesDP bdpsurvival
#' @importFrom stats median
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
  weibull_shape         = 3,
  method                = "fixed"
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
  post <- suppressWarnings(
    bdpsurvival(formula            = Surv(time, event) ~ treatment,
                data               = data_total,
                data0              = data0,
                breaks             = breaks,
                a0                 = prior[1],
                b0                 = prior[2],
                surv_time          = surv_time,
                discount_function  = discount_function,
                alpha_max          = alpha_max,
                fix_alpha          = fix_alpha,
                number_mcmc        = number_mcmc,
                weibull_scale      = weibull_scale,
                weibull_shape      = weibull_shape,
                method             = method))


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
    post_imp <- suppressWarnings(
      bdpsurvival(formula            = Surv(time, event) ~ treatment,
                  data               = data_success_impute,
                  data0              = data0,
                  breaks             = breaks,
                  a0                 = prior[1],
                  b0                 = prior[2],
                  surv_time          = surv_time,
                  discount_function  = discount_function,
                  alpha_max          = alpha_max,
                  fix_alpha          = fix_alpha,
                  number_mcmc        = number_mcmc,
                  weibull_scale      = weibull_scale,
                  weibull_shape      = weibull_shape,
                  method             = method))

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
  post_final <- suppressWarnings(
    bdpsurvival(formula            = Surv(time, event) ~ treatment,
                data               = data_final,
                data0              = data0,
                breaks             = breaks,
                a0                 = prior[1],
                b0                 = prior[2],
                surv_time          = surv_time,
                discount_function  = discount_function,
                alpha_max          = alpha_max,
                fix_alpha          = fix_alpha,
                number_mcmc        = number_mcmc,
                weibull_scale      = weibull_scale,
                weibull_shape      = weibull_shape,
                method             = method))

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
                                                        "subject_impute_futility"))


#' @title Historical data for survival analysis
#'
#' @description Wrapper function for historical data from time-to-event outcome.
#'
#' @inheritParams survival_analysis
#' @param .data NULL. stores the historical time, treatment and event , please do not fill it in.
#'
#' @return a list with historical data for time-to-event outcome with the discount function.
#'
#' @examples
#' historical_survival(time      = rexp(10, 0.01),
#'                     treatment = rep(10, 1),
#'                     event     = rep(10, 1))
#'
#' @export historical_survival
historical_survival <- function(time               = NULL,
                                treatment          = NULL,
                                event              = NULL,
                                discount_function  = "identity",
                                alpha_max          = 1,            # max weight on incorporating historical data
                                fix_alpha          = FALSE,        # fix alpha set weight of historical data to alpha_max
                                weibull_scale      = 0.135,        # weibull parameter
                                weibull_shape      = 3,            # weibull parameter
                                method             = "fixed",
                                .data              = NULL
){
  .data$time0              <- time
  .data$treatment0         <- treatment
  .data$event0             <- event
  .data$discount_function  <- discount_function
  .data$alpha_max          <- alpha_max
  .data$fix_alpha          <- fix_alpha
  .data$weibull_scale      <- weibull_scale
  .data$weibull_shape      <- weibull_shape
  .data$method             <- method
  .data
}


#' @title Gamma prior for for control and treatment group
#'
#' @description Wrapper function for gamma prior \code{Gamma(a0, b0)}.
#'
#' @param a0 numeric. The shape parameter in the gamma distribution
#'   (\code{beta(a0, b0)}).
#' @param b0 numeric. The scale parameter in the beta distribution
#'   (\code{beta(a0, b0)}).
#' @param .data NULL. stores the gamma prior rate, please do not fill it in.
#'
#' @return a list with vector of gamma rate for the gamma prior for treatment and control group.
#'
#' @examples gamma_prior(a0 = .1, b0 = .1)
#' @export gamma_prior
gamma_prior <- function(a0 = .1, b0 = .1, .data = NULL){
  .data$prior  <- c(a0, b0)
  .data
}


#' @title Piecewise constant hazard rates and the cutpoint for control and treatment group
#'
#' @description Wrapper function for the piecewise constant hazard rates and the cutpoint
#' for control and treatment group.
#'
#' @inheritParams survivalBACT
#' @param .data NULL. stores the hazard rates and cutpoint, please do not fill it in.
#'
#' @return a list with hazard rates and cutpoint for control and treatment group.
#'
#' @examples survival_outcome(hazard_treatment = 0.06,
#'                            hazard_control   = 0.08,
#'                            cutpoint         = NULL )
#' @export survival_outcome
survival_outcome <- function(hazard_treatment = NULL,
                             cutpoint         = NULL,
                             hazard_control   = NULL,
                             .data = NULL){
  .data$hazard_treatment  <- hazard_treatment
  .data$cutpoint          <- cutpoint
  .data$hazard_control    <- hazard_control
  .data
}


#' @title Data file for survival analysis
#'
#' @description Wrapper function for data file in survival analysis.
#'
#' @inheritParams survival_analysis
#' @param .data NULL. stores the survival data for analysis, please do not fill it in.
#'
#' @return a list with time, treatment, and event with time-to-event
#'   outcome.
#'
#' @examples data_survival(time       = c(6.2, 8.2, 8.0, 2.3),
#'                         treatment  = c(0, 1, 0, 1),
#'                         event      = c(1, 1, 1, 1))
#'
#' @export data_survival
data_survival <- function(time, treatment, event, .data = NULL){
  .data$time       <- time
  .data$treatment  <- treatment
  .data$event      <- event
  .data
}


