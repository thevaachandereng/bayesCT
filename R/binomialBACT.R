#' @title Binomial counts for Bayesian Adaptive Trial
#'
#'
#' @param p_treatment scalar. Proportion of an event under the treatment group.
#' @param p_control scalar. Proportion of an event under the contrl group.
#' @param N_total scalar. Total sample size
#' @param lambda vector. Lambda for different enrollment rates across times.
#' @param lambda_time vector. Same size as lambda, denote times at lambda changes.
#' @param analysis_at_enrollnumber vector. Sample size for interim looks,
#'                                 total size not included.
#' @param EndofStudy scalar. Length of the study.
#' @param weibull_scale vector. The weibull scale used in control and treatment group.
#' @param weibull_shape vector. The weibull shape used in control and treatment group.
#' @param discount_function character. The default is set to identity.
#'                          For more options look at bayesDP package.
#' @param bpd_method character. The default is set to fixed.
#'                   For more options look at bayesDP package.
#' @param block scalar. Block size for randomization to be implemented.
#' @param rand.ratio vector. Randomization ratio for control to treatment.
#'                    Integer values mapping the size of the block.
#' @param prop_loss_to_followup scalar. Proportion of loss in follow up.
#' @param h0 scalar. treshold for comparing two proportions. default at 0.
#' @param futility_prob scalar. Type 1 error rate.
#' @param expected_success_prob scalar. Power of the study.
#' @param prob_ha scalar. Probability of alternative hypothesis.
#' @param N_impute scalar. Number of imputation for monte carlo simulation.
#'
#' @return a list of output
#'
#'
#'
#'
#'


binomialBACT <- function(
  p_treatment,
  p_control,
  N_total,
  lambda,
  lambda_time,
  analysis_at_enrollnumber,
  EndofStudy,
  weibull_scale,
  weibull_shape,
  discount_function     = "identity",
  bpd_method            = "fixed",
  block                 = 2,              # block size for randomization
  rand.ratio            = c(1, 1),        # randomization ratio in control to treatament (default 1:1)
  prop_loss_to_followup = 0.15,           # Proportion of loss in data
  h0                    = 0,              # Null hypothesis value
  futility_prob         = 0.05,           # Futility probability
  expected_success_prob = 0.9,            # Expected success probability
  prob_ha               = 0.95,           # Posterior probability of accepting alternative hypothesis
  N_impute              = 100             # Number of imputation simulations for predictive distribution

  ){

  #assignment of enrollment based on the enrollment function
  enrollment <- enrollment(param = lambda, N_total = N_total, time = lambda_time)

  # simulating group and treatment group assignment
  group <- randomization(N_total = N, block = block, scheme = rand.ratio)

  #simulate binomial outcome
  sim <- rbinom(N_total, 1, prob = group * p_treatment + (1 - group) * p_control)

  #dividing treatment and control
  control <- sim[group == 0]
  treatment <- sim[group == 1]

  # Simulate loss to follow-up
  n_loss_to_fu <- ceiling(prop_loss_to_followup * N_total)
  loss_to_fu   <- rep(FALSE,N_total)
  loss_to_fu[sample(1:N_total, n_loss_to_fu)] <- TRUE


  # creating a new data.frame for all the variables
  data_total <- data.frame(
    Y          = sim,
    treatment  = group,
    enrollment = enrollment,
    id         = 1:N_total,
    loss_to_fu = loss_to_fu)

  #assigning stop_futility and expected success

  stop_futility         <- 0
  stop_expected_success <- 0



  for(i in 1:(length(analysis_at_enrollnumber)-1) ){

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
      mutate(subject_impute_sucess = (max(enrollment)-enrollment <= EndofStudy & subject_enrolled) |
               (subject_enrolled & loss_to_fu))

    # Carry out interim analysis on patients with complete data only
    # - Set-up `new data` data frame
    data <- data_interim %>%
      filter(subject_enrolled,
             !subject_impute_sucess)

    # MLE of data at interim analysis
    MLE_int        <- glm(Y ~ treatment, data = data, family = "binomial")

    # Analyze data using discount funtion via binomial
    post <- bdpbinomial(y_t                    = sum(data$Y[data$treatment == 1]),
                        N_t                    = length(data$Y[data$treatment == 1]),
                        y_c                    = sum(data$Y[data$treatment == 0]),
                        N_c                    = length(data$Y[data$treatment == 0]),
                        number_mcmc            = number_mcmc,
                        discount_function      = discount_function,
                        method                 = bdp_method,
                        weibull_scale          = weibull_scale,
                        weibull_shape          = weibull_shape)


    # Imputation phase futility and expected success - initialize counters
    # for the current imputation phase
    futility_test         <- 0
    expected_success_test <- 0

    for(j in 1:N_impute){
      #imputing the success for control group
      data_control_success_impute <- data_interim %>%
        filter(treatment == 0) %>%
        mutate(Y_impute = ifelse(subject_impute_sucess & subject_enrolled,
                                 rbinom(n(), 1, p_control),
                                 Y))
      # imputing success for treatment group
      data_treatment_success_impute  <- data_interim %>%
        filter(treatment == 1) %>%
        mutate(Y_impute = ifelse(subject_impute_sucess & subject_enrolled,
                                 rbinom(n(), 1, p_treatment),
                                 Y))

      # Combine the treatment and control imputed datasets
      data_success_impute <- bind_rows(data_control_success_impute,
                                       data_treatment_success_impute) %>%
        mutate(Y = Y_impute) %>%
        select(-Y_impute)

      #Create enrolled subject data frame for discount function analysis
      data <- data_success_impute %>%
        filter(subject_enrolled)


      # Analyze complete+imputed data using discount funtion via binomial
      post_imp <- bdpbinomial(y_t                    = sum(data$Y[data$treatment == 1]),
                              N_t                    = length(data$Y[data$treatment == 1]),
                              y_c                    = sum(data$Y[data$treatment == 0]),
                              N_c                    = length(data$Y[data$treatment == 0]),
                              number_mcmc            = number_mcmc,
                              discount_function      = discount_function,
                              method                 = bdp_method,
                              weibull_scale          = weibull_scale,
                              weibull_shape          = weibull_shape)


      # Estimation of the posterior effect for difference between test and control
      # - If expected success, add 1 to the counter
      post_control <- post_imp$posterior_control$posterior_flat
      post_treatment <- post_imp$posterior_treatment$posterior_flat
      if(mean(post_control - post_treatment > h0) > prob_ha){
        expected_success_test <- expected_success_test + 1
      }

      ##########################################################################
      ### Futility computations
      ##########################################################################
      ## For patients not enrolled, impute the outcome

      data_control_futility_impute <- data_success_impute %>%
        filter(treatment == 0) %>%
        mutate(Y_impute = ifelse(subject_impute_futility,
                                 rbinom(n(), 1, p_control),
                                 Y))

      data_treatment_futility_impute <- data_success_impute %>%
        filter(treatment == 1) %>%
        mutate(Y_impute = ifelse(subject_impute_futility,
                                 rbinom(n(), 1, p_treatment),
                                 Y))

      # Combine the treatment and control imputed datasets
      data_futility_impute <- bind_rows(data_control_futility_impute,
                                        data_treatment_futility_impute) %>%
        mutate(Y = Y_impute) %>%
        select(-Y_impute)


      #Create enrolled subject data frame for discount function analysis
      data <- data_futility_impute

      # Analyze complete+imputed data using discount funtion via binomial
      post_imp <- bdpbinomial(y_t                    = sum(data$Y[data$treatment == 1]),
                              N_t                    = length(data$Y[data$treatment == 1]),
                              y_c                    = sum(data$Y[data$treatment == 0]),
                              N_c                    = length(data$Y[data$treatment == 0]),
                              number_mcmc            = number_mcmc,
                              discount_function      = discount_function,
                              method                 = bdp_method,
                              weibull_scale          = weibull_scale,
                              weibull_shape          = weibull_shape)

      # Estimation of the posterior effect for difference between test and control
      post_control <- post_imp$posterior_control$posterior_flat
      post_treatment <- post_imp$posterior_treatment$posterior_flat

      # Increase futility counter by 1 if P(effect_imp < h0) > ha
      if(mean(post_control - post_treatment > h0) > prob_ha){
        futility_test <- futility_test + 1
      }

    }

    # Test if futility success criteria is met
    if(futility_test/N_imput < futility_prob){
      stop_futility       <- 1
      stage_trial_stopped <- analysis_at_enrollnumber[ii]
      break
    }

    # Test if expected success criteria met
    if(expected_success_test/N_imput > expected_success_prob ){
      stop_expected_success <- 1
      stage_trial_stopped   <- analysis_at_enrollnumber[ii]
      break
    }

    # Stop study if at last interim look
    if(analysis_at_enrollnumber[ii+1] == N_total){
      stage_trial_stopped <- analysis_at_enrollnumber[ii+1]
      break
    }

  }

  ##############################################################################
  ### 4) Final analysis
  ##############################################################################
  # Estimation of the posterior of the difference
  effect_int <- post$final$posterior

  # Number of patients enrolled at trial stop
  N_enrolled <- nrow(data_interim[data_interim$id <= stage_trial_stopped, ])

  # All patients that have made it to the end of study
  # - Subset out patients loss to follow-up
  data_final <- data_interim %>%
    filter(id <= stage_trial_stopped,
           !loss_to_fu)

  MLE <- glm(Y ~ treatment, data=data_final, family = "binomial")

  post <- bdpbinomial(y_t                 = sum(data$Y[data$treatment == 1]),
                      N_t                 = length(data$Y[data$treatment == 1]),
                      y_c                 = sum(data$Y[data$treatment == 0]),
                      N_c                 = length(data$Y[data$treatment == 0]),
                      number_mcmc         = number_mcmc,
                      discount_function   = discount_function,
                      method              = bdp_method,
                      weibull_scale       = weibull_scale,
                      weibull_shape       = weibull_shape)


  ### Format and output results
  effect       <- post$final$posterior       #Posterior effect size: test vs control
  N_treatment  <- sum(data_final$treatment)  #Total sample size analyzed - test group
  N_control    <- sum(!data_final$treatment) #Total sample size analyzed - control group


  ## output
  results_list <- list(
    prob_ha               = prob_ha,
    N_treatment           = N_treatment,
    N_control             = N_control,
    N_total               = N_enrolled,             # Total sample size enrolled when trial stopped
    N_max                 = N_total, 				        # Total potential sample size
    stop_futility         = stop_futility,          # Did the trial stop for futility
    stop_expected_success = stop_expected_success,  # Did the trial stop for expected success
    prob_ha_final         = mean(effect < h0),      # Posterior probability that alternative hypothesis is true
    est_final             = mean(effect),           # Posterior Mean of treatment effect
    est_int               = mean(effect_int),       # Posterior Mean of treatment effect at interim analysis
    MLE_est               = MLE$coe[2],             # Treatment effect useing MLE
    MLE_est_int           = MLE_int$coe[2],         # Treatment effect useing MLE at interim analysis
    discount_function     = discount_function,      # Discount function used
    weibull_scale         = weibull_scale,          # Discount Weibull function scale
    weibull_shape         = weibull_shape,          # Discount Weibull function schape
    p_treatment           = p_treatment,            # probability of treatment in binomial
    p_control             = p_control               # probability of control in binomial
  )

  #return results
  results_list

}

