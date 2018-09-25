#' @title Normal distribution for Bayesian Adaptive Trials
#'
#' @description Simulation for normal distribution for Bayesian Adaptive trial with different inputs to
#'              control for power, sample size, type I error rate, etc.
#'
#' @param mu_control scalar. Mean outcome in the control arm.
#' @param sd_control scalar. Standard deviation of outcome in the control arm.
#' @param mu_treatment scalar. Mean outcome in the treatment arm.
#' @param sd_treatment scalar. Standard deviation of outcome in the treatment arm.
#' @param N_total scalar. Total sample size.
#' @param lambda vector. Lambda for different enrollment rates across times.
#' @param lambda_time vector. Same size as lambda, denote times at lambda changes.
#' @param interim_look vector. Sample size for interim looks. Note: the overall trial size should not be included.
#' @param EndofStudy scalar. Length of the study.
#' @param prior vector. Prior value of beta rate, beta(a0, b0).
#' @param block scalar. Block size for randomization to be implemented.
#' @param rand_ratio vector. Randomization allocation for control to treatment.
#'                    Integer values mapping the size of the block.
#' @param alternative character. The string specifying the alternative hypothesis, must be one
#'        of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.
#' @param prop_loss_to_followup scalar. Proportion of subjects lost to follow-up.
#' @param h0 scalar. Treshold for comparing two proportions. Default is \code{h0=0}.
#' @param futility_prob scalar. Type I error rate.
#' @param expected_success_prob scalar. Power of the study.
#' @param prob_ha scalar. Probability of alternative hypothesis.
#' @param N_impute scalar. Number of imputations for Monte Carlo simulation of missing data.
#' @param number_mcmc scalar. Number of Monte Carlo Markov Chain draws in sampling posterior.
#'
#' @return a list of output
#'
#' @examples
#' normalBACT(mu_control = 10, mu_treatment = 8,
#'            sd_control = 0.8, sd_treatment = 1.2, N_total = 300,
#'            lambda = c(0.3, 1), lambda_time = c(25),
#'            interim_look = c(110, 140, 220, 270),
#'            EndofStudy = 50)
#'
#' @importFrom stats rnorm lm sd
#' @importFrom dplyr mutate filter group_by bind_rows select n
#' @importFrom bayesDP bdpnormal
#' @export normalBACT
normalBACT <- function(
  mu_control,
  sd_control,
  mu_treatment,
  sd_treatment,
  N_total,
  lambda,
  lambda_time,
  interim_look,
  EndofStudy,
  prior                 = c(1, 1),
  block                 = 2,            # block size for randomization
  rand_ratio            = c(1, 1),      # randomization ratio in control to treatament (default 1:1)
  alternative           = "two-sided",  # the alternative hypothesis (either two-sided, greater, less)
  prop_loss_to_followup = 0.15,         # Proportion of loss in data
  h0                    = 0,            # Null hypothesis value
  futility_prob         = 0.05,         # Futility probability
  expected_success_prob = 0.9,          # Expected success probability
  prob_ha               = 0.95,         # Posterior probability of accepting alternative hypothesis
  N_impute              = 100,          # Number of imputation simulations for predictive distribution
  number_mcmc           = 1000
){
  # checking inputs
  stopifnot((mu_control > 0 & sd_control > 0), (mu_treatment > 0 & sd_treatment > 0),
            all(N_total > interim_look), length(lambda) == (length(lambda_time) + 1),
            EndofStudy > 0, block %% sum(rand_ratio)  == 0,
            (prop_loss_to_followup >= 0 & prop_loss_to_followup < 0.75),
            (h0 >= 0 & h0 < 1), (futility_prob < 0.20 & futility_prob > 0),
            (expected_success_prob > 0.70 & expected_success_prob <= 1),
            (prob_ha > 0.70 & prob_ha < 1), N_impute > 0)

  #checking if alternative is right
  if(alternative != "two-sided" | alternative != "greater" | alternative != "less"){
    stop("The input for alternative is wrong!")
  }

  # assigining interim look and final look
  analysis_at_enrollnumber <- c(interim_look, N_total)

  # assignment of enrollment based on the enrollment function
  enrollment <- enrollment(param = lambda, N_total = N_total, time = lambda_time)

  # simulating group and treatment group assignment
  group <- randomization(N_total = N_total, block = block, allocation = rand_ratio)

  # simulate binomial outcome
  sim <- rnorm(N_total, mean = group * mu_treatment + (1 - group) * mu_control,
               sd = group * sd_treatment + (1 - group) * sd_control)

  # dividing treatment and control
  control <- sim[group == 0]
  treatment <- sim[group == 1]

  # simulate loss to follow-up
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

  # assigning stop_futility and expected success
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
      mutate(subject_impute_success = (max(enrollment) - enrollment <= EndofStudy & subject_enrolled) |
               (subject_enrolled & loss_to_fu))

    # Carry out interim analysis on patients with complete data only
    # - Set-up `new data` data frame
    data <- data_interim %>%
      filter(subject_enrolled,
             !subject_impute_success)

    # MLE of data at interim analysis
    MLE_int <- lm(Y ~ treatment, data = data)

    # analyze data using discount funtion via normal
    post <- bdpnormal(mu_t         = mean(data$Y[data$treatment == 1]),
                      sigma_t      = sd(data$Y[data$treatment == 1]),
                      N_t          = length(data$treatment == 1),
                      mu_c         = mean(data$Y[data$treatment == 0]),
                      sigma_c      = sd(data$Y[data$treatment == 0]),
                      N_c          = length(data$treatment == 0),
                      number_mcmc  = number_mcmc)

    # Imputation phase futility and expected success - initialize counters
    # for the current imputation phase
    futility_test         <- 0
    expected_success_test <- 0

    for(j in 1:N_impute){
      # imputing the success for control group
      data_control_success_impute <- data_interim %>%
        filter(treatment == 0) %>%
        mutate(Y_impute = ifelse(subject_impute_success & subject_enrolled,
                                 rnorm(n(), mu_control, sd_control),
                                 Y))
      # imputing success for treatment group
      data_treatment_success_impute  <- data_interim %>%
        filter(treatment == 1) %>%
        mutate(Y_impute = ifelse(subject_impute_success & subject_enrolled,
                                 rnorm(n(), mu_treatment, sd_treatment),
                                 Y))

      # combine the treatment and control imputed datasets
      data_success_impute <- bind_rows(data_control_success_impute,
                                       data_treatment_success_impute) %>%
        mutate(Y = Y_impute) %>%
        select(-Y_impute)

      # create enrolled subject data frame for discount function analysis
      data <- data_success_impute %>%
        filter(subject_enrolled)

      # analyze complete+imputed data using discount funtion via normal
      post_imp <- bdpnormal(mu_t         = mean(data$Y[data$treatment == 1]),
                            sigma_t      = sd(data$Y[data$treatment == 1]),
                            N_t          = length(data$treatment == 1),
                            mu_c         = mean(data$Y[data$treatment == 0]),
                            sigma_c      = sd(data$Y[data$treatment == 0]),
                            N_c          = length(data$treatment == 0),
                            number_mcmc  = number_mcmc)


      # Estimation of the posterior effect for difference between test and control
      # - If expected success, add 1 to the counter
      post_imp_final <- post_imp$final$posterior
      if(mean(post_imp_final < h0) > prob_ha){
        expected_success_test <- expected_success_test + 1
      }

      ##########################################################################
      ### Futility computations
      ##########################################################################

      # For patients not enrolled, impute the outcome
      data_control_futility_impute <- data_success_impute %>%
        filter(treatment == 0) %>%
        mutate(Y_impute = ifelse(subject_impute_futility,
                                 rnorm(n(), mu_control, sd_control),
                                 Y))

      data_treatment_futility_impute <- data_success_impute %>%
        filter(treatment == 1) %>%
        mutate(Y_impute = ifelse(subject_impute_futility,
                                 rnorm(n(), mu_treatment, sd_treatment),
                                 Y))

      # Combine the treatment and control imputed datasets
      data_futility_impute <- bind_rows(data_control_futility_impute,
                                        data_treatment_futility_impute) %>%
        mutate(Y = Y_impute) %>%
        select(-Y_impute)


      # Create enrolled subject data frame for discount function analysis
      data <- data_futility_impute

      # Analyze complete+imputed data using discount funtion via normal
      post_imp <- bdpnormal(mu_t         = mean(data$Y[data$treatment == 1]),
                            sigma_t      = sd(data$Y[data$treatment == 1]),
                            N_t          = length(data$treatment == 1),
                            mu_c         = mean(data$Y[data$treatment == 0]),
                            sigma_c      = sd(data$Y[data$treatment == 0]),
                            N_c          = length(data$treatment == 0),
                            number_mcmc  = number_mcmc)

      # Estimation of the posterior effect for difference between test and control
      post_imp_final <- post_imp$final$posterior

      # Increase futility counter by 1 if P(effect_imp < h0) > ha
      if(mean(post_imp_final < h0) > prob_ha){
        futility_test <- futility_test + 1
      }

    }

    print(analysis_at_enrollnumber[i])

    # Test if futility success criteria is met
    if(futility_test / N_impute < futility_prob){
      stop_futility       <- 1
      stage_trial_stopped <- analysis_at_enrollnumber[i]
      break
    }

    # Test if expected success criteria met
    if(expected_success_test / N_impute > expected_success_prob ){
      stop_expected_success <- 1
      stage_trial_stopped   <- analysis_at_enrollnumber[i]
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
  effect_int <- post$final$posterior


  # Number of patients enrolled at trial stop
  N_enrolled <- nrow(data_interim[data_interim$id <= stage_trial_stopped, ])

  print(N_enrolled)

  # All patients that have made it to the end of study
  # - Subset out patients loss to follow-up
  data_final <- data_interim %>%
    filter(id <= stage_trial_stopped,
           !loss_to_fu)


  # Compute the final MLE for the complete data using GLM function
  MLE <- lm(Y ~ treatment, data = data_final)

  # Analyze complete data using discount funtion via binomial
  post <- bdpnormal(mu_t         = mean(data_final$Y[data_final$treatment == 1]),
                    sigma_t      = sd(data_final$Y[data_final$treatment == 1]),
                    N_t          = length(data_final$treatment == 1),
                    mu_c         = mean(data_final$Y[data_final$treatment == 0]),
                    sigma_c      = sd(data_final$Y[data_final$treatment == 0]),
                    N_c          = length(data_final$treatment == 0),
                    number_mcmc  = number_mcmc)


  # Format and output results
  effect        <- post$final$posterior       #Posterior effect size: test vs control
  N_treatment   <- sum(data_final$treatment)  #Total sample size analyzed - test group
  N_control     <- sum(!data_final$treatment) #Total sample size analyzed - control group


  # output
  results_list <- list(
    mu_treatment_true                          = mu_treatment,            # mean of treatment in normal
    mu_control_true                            = mu_control,              # mean of control in normal
    sd_treatment_true                          = sd_treatment,            # sd of treatment in normal
    sd_control_true                            = sd_control,              # sd of control in normal
    prob_of_accepting_alternative              = prob_ha,
    N_treatment                                = N_treatment,
    N_control                                  = N_control,
    N_complete                                 = N_treatment + N_control,
    N_enrolled                                 = N_enrolled,              # Total sample size enrolled when trial stopped
    N_max                                      = N_total, 				        # Total potential sample size
    stop_futility                              = stop_futility,           # Did the trial stop for futility
    stop_expected_success                      = stop_expected_success,   # Did the trial stop for expected success
    post_prob_accept_alternative               = mean(effect < h0),       # Posterior probability that alternative hypothesis is true
    est_final                                  = mean(effect),            # Posterior Mean of treatment effect
    est_interim                                = mean(effect_int)         # Posterior Mean of treatment effect at interim analysis
    #MLE_est                                   = MLE$coe[2],              # Treatment effect useing MLE
    #MLE_est_interim                           = MLE_int$coe[2]           # Treatment effect useing MLE at interim analysis
  )

  # return results
  results_list

}






## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Y", "Y_impute", "id", "subject_enrolled",
                                                        "subject_impute_success", "subject_impute_futility"))



#' @title Parameters for treatment and control in normal case
#'
#' @description Wrapper function for mean and standard deviation with normal outcome.
#'
#' @param mu_control_true numeric. The mean for the control group.
#' @param sd_control_true numeric. The standard deviation for the control group.
#' @param mu_treatment_true numeric. The mean for the treatment group.
#' @param sd_treatment_true numeric. The standard deviation for the treatment group.
#' @param data NULL. stores the proportion of control and treatment, please do not fill it in.
#'
#' @return a list with proportion of control and treatment group.
#'
#' @examples normal_outcome(mu_control = 12, mu_treatment = 8, sd_treatment = 2.2, sd_control = 1.6)
#' @export normal_outcome
#'
normal_outcome <- function(mu_control_true = NULL, sd_control_true = NULL, mu_treatment_true = NULL, sd_treatment_true = NULL, data = NULL){
  data$mu_control <- mu_control_true
  data$sd_control <- sd_control_true
  data$mu_treatment <- mu_treatment_true
  data$sd_treatment <- sd_treatment_true
  data
}



#' @title Complete normal wrapper function
#'
#' @description Wrapper function for complete normal bayesCT function to compute power and type 1 error.
#'
#' @param input list. Input function for all normalBACT.
#' @param .data NULL. stores the proportion of control and treatment, please do not fill it in.
#'
#' @return a list with results of the clinical outcome.
#'
#' @importFrom stats rnorm lm sd
#' @importFrom dplyr mutate filter group_by bind_rows select n
#' @importFrom bayesDP bdpnormal
#'
#' @export BACTnormal
#'
BACTnormal <- function(input, .data = NULL){
  .data <- do.call(normalBACT, input)
  .data
}

#' @title Historical data for normal distribution
#'
#' @description Wrapper function for historical data from normal outcome.
#'
#' @param mu0_treatment numeric. Mean of the historical treatment group.
#' @param sd0_treatment numeric. The  Standard deviation of the historical treatment group.
#' @param N0_treatment numeric. scalar. Number of observations of the historical treatment group.
#' @param mu0_control numeric. Mean of the historical control group.
#' @param sd0_control numeric. The  Standard deviation of the historical control group.
#' @param N0_control numeric. umber of observations of the historical control group.
#' @param discount_function character. Specify the discount function to use. Currently supports weibull,
#'                          scaledweibull, and identity. The discount function scaledweibull scales the
#'                          output of the Weibull CDF to have a max value of 1. The identity discount function
#'                          uses the posterior probability directly as the discount weight. Default value is
#'                          "identity".
#' @param .data NULL. stores the proportion of control and treatment, please do not fill it in.
#'
#' @return a list with historical data for control and treatment group with the discount function.
#'
#' @examples historical_normal(mu0_treatment = 15, sd0_treatment = 2, N0_treatment = 10,
#'                             mu0_control = 17, sd0_control = 3, N0_control = 20)
#' @export historical_normal
historical_normal <- function(mu0_treatment       = NULL,
                              sd0_treatment       = NULL,
                              N0_treatment        = NULL,
                              mu0_control         = NULL,
                              sd0_control         = NULL,
                              N0_control          = NULL,
                              discount_function   = "identity",
                              .data               = NULL
){
  .data$mu0_treatment       <- mu0_treatment
  .data$sd0_treatment       <- sd0_treatment
  .data$N0_treatment        <- N0_treatment
  .data$mu0_control         <- mu0_control
  .data$sd0_control         <- sd0_control
  .data$N0_control          <- N0_control
  .data$discount_function   <- discount_function
  .data
}




