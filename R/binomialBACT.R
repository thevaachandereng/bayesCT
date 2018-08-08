#' @title Binomial counts for Bayesian Adaptive Trial
#'
#' @description Simulation for binomial counts for Bayesian Adaptive trial with different inputs to
#'              control for power, sample size, type 1 error rate, and etc.
#'
#' @param p_control scalar. Proportion of an event under the contrl group.
#' @param p_treatment scalar. Proportion of an event under the treatment group.
#' @param N_total scalar. Total sample size
#' @param lambda vector. Lambda for different enrollment rates across times.
#' @param lambda_time vector. Same size as lambda, denote times at lambda changes.
#' @param interim_look vector. Sample size for interim looks,
#'                                 total size not included.
#' @param EndofStudy scalar. Length of the study.
#' @param prior vector. Prior value of beta rate, beta(a0, b0)
#' @param block scalar. Block size for randomization to be implemented.
#' @param rand_ratio vector. Randomization ratio for control to treatment.
#'                    Integer values mapping the size of the block.
#' @param prop_loss_to_followup scalar. Proportion of loss in follow up.
#' @param h0 scalar. treshold for comparing two proportions. default at 0.
#' @param futility_prob scalar. Type 1 error rate.
#' @param expected_success_prob scalar. Power of the study.
#' @param prob_ha scalar. Probability of alternative hypothesis.
#' @param N_impute scalar. Number of imputation for monte carlo simulation for missing data.
#' @param number_mcmc scalar. Number of Monte Carlo Markov Chain in sampling posterior.
#'
#' @return a list of output
#'
#' @examples
#' binomialBACT(p_control = 0.12, p_treatment = 0.10, N_total = 300,
#'              lambda = c(0.3, 1), lambda_time = c(25),
#'              interim_look = c(110, 140, 220, 270),
#'              EndofStudy = 50)
#'
#' @importFrom magrittr %>%
#' @importFrom stats rbinom glm
#' @importFrom dplyr mutate filter group_by bind_rows select n
#' @importFrom bayesDP bdpbinomial
#' @export binomialBACT

binomialBACT <- function(
  p_treatment,
  p_control             = NULL,
  N_total,
  lambda,
  lambda_time,
  interim_look,
  EndofStudy,
  prior                 = c(1, 1),
  block                 = 2,            # block size for randomization
  rand_ratio            = c(1, 1),      # randomization ratio in control to treatament (default 1:1)
  prop_loss_to_followup = 0.10,         # Proportion of loss in data
  h0                    = 0,            # Null hypothesis value
  futility_prob         = 0.05,         # Futility probability
  expected_success_prob = 0.9,          # Expected success probability
  prob_ha               = 0.95,         # Posterior probability of accepting alternative hypothesis
  N_impute              = 1000,         # Number of imputation simulations for predictive distribution
  number_mcmc           = 1000
  ){
  #checking p_control
  if(!is.null(p_control)){stopifnot(p_control > 0 & p_control < 1)}
  #checking inputs
  stopifnot((p_treatment < 1 & p_treatment > 0), all(N_total > interim_look),
            length(lambda) == (length(lambda_time) + 1),
            EndofStudy > 0, block %% sum(rand_ratio)  == 0,
            (prop_loss_to_followup >= 0 & prop_loss_to_followup < 0.75),
            (h0 >= 0 & h0 < 1), (futility_prob < 0.20 & futility_prob > 0),
            (expected_success_prob > 0.70 & expected_success_prob <= 1),
            (prob_ha > 0.70 & prob_ha < 1), N_impute > 0)

  ##assigining interim look and final look
  analysis_at_enrollnumber <- c(interim_look, N_total)

  #assignment of enrollment based on the enrollment function
  enrollment <- enrollment(param = lambda, N_total = N_total, time = lambda_time)

  # simulating group and treatment group assignment
  if(!is.null(p_control)){
    group <- randomization(N_total = N_total, block = block, scheme = rand_ratio)
  }
  else{
    group <- rep(1, N_total)
  }

  #simulate binomial outcome
  if(!is.null(p_control)){
    sim <- rbinom(N_total, 1, prob = group * p_treatment + (1 - group) * p_control)
    #dividing treatment and control
    control <- sim[group == 0]
    treatment <- sim[group == 1]
  }
  else{
    sim <- rbinom(N_total, 1, p_treatment)
  }

  # Simulate loss to follow-up
  n_loss_to_fu <- ceiling(prop_loss_to_followup * N_total)
  loss_to_fu   <- rep(FALSE, N_total)
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
      mutate(subject_impute_success = (max(enrollment) - enrollment <= EndofStudy & subject_enrolled) |
               (subject_enrolled & loss_to_fu))


    # Carry out interim analysis on patients with complete data only
    # - Set-up `new data` data frame
    data <- data_interim %>%
      filter(subject_enrolled,
             !subject_impute_success)

    # MLE of data at interim analysis
    # MLE_int        <- glm(Y ~ treatment, data = data, family = "binomial")

    # assigning input for control arm given it is a single or double arm
    if(!is.null(p_control)){
      y_c <- sum(data$Y[data$treatment == 0])
      N_c <- length(data$Y[data$treatment == 0])
    }
    else{
      y_c <- NULL
      N_c <- NULL
    }

    # Analyze data using discount funtion via binomial
    post <- bdpbinomial(y_t                    = sum(data$Y[data$treatment == 1]),
                        N_t                    = length(data$Y[data$treatment == 1]),
                        y_c                    = y_c,
                        N_c                    = N_c,
                        number_mcmc            = number_mcmc,
                        a0                     = prior[1],
                        b0                     = prior[2])


    # Imputation phase futility and expected success - initialize counters
    # for the current imputation phase
    futility_test         <- 0
    expected_success_test <- 0

    for(j in 1:N_impute){
      #imputing the success for control group
      data_control_success_impute <- data_interim %>%
        filter(treatment == 0) %>%
        mutate(Y_impute = ifelse(subject_impute_success & subject_enrolled,
                                 rbinom(n(), 1, p_control),
                                 Y))
      # imputing success for treatment group
      data_treatment_success_impute  <- data_interim %>%
        filter(treatment == 1) %>%
        mutate(Y_impute = ifelse(subject_impute_success & subject_enrolled,
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


      # assigning input for control arm given it is a single or double arm
      if(!is.null(p_control)){
        y_c <- sum(data$Y[data$treatment == 0])
        N_c <- length(data$Y[data$treatment == 0])
      }
      else{
        y_c <- NULL
        N_c <- NULL
      }

      # Analyze complete+imputed data using discount funtion via binomial
      post_imp <- bdpbinomial(y_t                    = sum(data$Y[data$treatment == 1]),
                              N_t                    = length(data$Y[data$treatment == 1]),
                              y_c                    = y_c,
                              N_c                    = N_c,
                              number_mcmc            = number_mcmc,
                              a0                     = prior[1],
                              b0                     = prior[2])


      # Estimation of the posterior effect for difference between test and control
      # - If expected success, add 1 to the counter
      effect_imp <- post_imp$final$posterior


      if(mean(effect_imp < h0) > prob_ha){
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

      if(!is.null(p_control)){
        y_c <- sum(data$Y[data$treatment == 0])
        N_c <- length(data$Y[data$treatment == 0])
      }
      else{
        y_c <- NULL
        N_c <- NULL
      }

      # Analyze complete+imputed data using discount funtion via binomial
      post_imp <- bdpbinomial(y_t                    = sum(data$Y[data$treatment == 1]),
                              N_t                    = length(data$Y[data$treatment == 1]),
                              y_c                    = y_c,
                              N_c                    = N_c,
                              number_mcmc            = number_mcmc,
                              a0                     = prior[1],
                              b0                     = prior[2])

      # Estimation of the posterior effect for difference between test and control
      effect_imp <- post_imp$final$posterior

      # Increase futility counter by 1 if P(effect_imp < h0) > ha
      if(mean(effect_imp < h0) > prob_ha){
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
  ### 4) Final analysis
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
  # MLE <- glm(Y ~ treatment, data = data_final, family = "binomial")

  if(!is.null(p_control)){
    y_c <- sum(data_final$Y[data_final$treatment == 0])
    N_c <- length(data_final$Y[data_final$treatment == 0])
  }
  else{
    y_c <- NULL
    N_c <- NULL
  }

  # Analyze complete data using discount funtion via binomial
  post <- bdpbinomial(y_t                 = sum(data_final$Y[data_final$treatment == 1]),
                      N_t                 = length(data_final$Y[data_final$treatment == 1]),
                      y_c                 = y_c,
                      N_c                 = N_c,
                      number_mcmc         = number_mcmc,
                      a0                  = prior[1],
                      b0                  = prior[2])


  ### Format and output results
  effect       <- post$final$posterior        # Posterior effect size: test vs control or treatment itself
  N_treatment  <- sum(data_final$treatment)   # Total sample size analyzed - test group
  N_control    <- sum(!data_final$treatment)  # Total sample size analyzed - control group


  ## output
  results_list <- list(
    p_treatment                                = p_treatment,             # probability of treatment in binomial
    p_control                                  = p_control,               # probability of control in binomial
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
    est_interim                                = mean(effect_int)        # Posterior Mean of treatment effect at interim analysis
    #MLE_est                                    = MLE$coe[2],              # Treatment effect useing MLE
    #MLE_est_interim                            = MLE_int$coe[2]           # Treatment effect useing MLE at interim analysis
  )

  #return results
  results_list

}



## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Y", "Y_impute", "id", "subject_enrolled",
"subject_impute_success", "subject_impute_futility"))


#' @title Proportion of failure in control and treatment
#'
#' @description Wrapper function for proportion of failure in control and treatment group with binomial outcome.
#'
#' @param p_control_true numeric. The proportion of failure in the control group, 0 < $p_control$ < 1.
#' @param p_treatment_true numeric. The proportion of failure in the treatment group, 0 < $p_treatment$ < 1.
#' @param data NULL. stores the proportion of control and treatment, please do not fill it in.
#'
#' @return a list with proportion of control and treatment group.
#'
#' @examples binomial_outcome(p_control_true = 0.12, p_treatment_true = 0.08)
#' @export binomial_outcome
binomial_outcome <- function(p_control_true = NULL, p_treatment_true = NULL, data = NULL){
  data$p_control <- p_control_true
  data$p_treatment <- p_treatment_true
  data
}






