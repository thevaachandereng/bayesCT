#' @title Normal distribution for Bayesian Adaptive Trials
#'
#' @description Simulation of normally distributed data for Bayesian adaptive
#'   trials with various inputs to control for power, sample size, type I error
#'   rate, etc.
#'
#' @param mu_treatment scalar. Mean outcome in the treatment arm.
#' @param sd_treatment scalar. Standard deviation of outcome in the treatment
#' @param mu_control scalar. Mean outcome in the control arm.
#' @param sd_control scalar. Standard deviation of outcome in the control arm.
#'   arm.
#' @param mu0_treatment scalar. Mean of the historical treatment group.
#' @param sd0_treatment scalar. Standard deviation of the historical treatment group.
#' @param N0_treatment scalar. Number of observations of the historical treatment group.
#' @param mu0_control scalar. Mean of the historical control group.
#' @param sd0_control scalar. Standard deviation of the historical control group.
#' @param N0_control scalar. Number of observations of the historical control group.
#' @param N_total scalar. Total sample size.
#' @param lambda vector. Enrollment rates across simulated enrollment times. See
#'   \code{\link{enrollment}} for more details.
#' @param lambda_time vector. Enrollment time(s) at which the enrollment rates
#'   change. Must be same length as lambda. See \code{\link{enrollment}} for
#'   more details.
#' @param interim_look vector. Sample size for each interim look. Note: the
#'   maximum sample size should not be included.
#' @param EndofStudy scalar. Length of the study.
#' @param block scalar. Block size for generating the randomization schedule.
#' @param rand_ratio vector. Randomization allocation for the ratio of control
#'   to treatment. Integer values mapping the size of the block. See
#'   \code{\link{randomization}} for more details.
#' @param discount_function character. If incorporating historical data, specify
#'   the discount function. Currently supports the Weibull function
#'   (\code{discount_function="weibull"}), the scaled-Weibull function
#'   (\code{discount_function="scaledweibull"}), and the identity function
#'   (\code{discount_function="identity"}). The scaled-Weibull discount function
#'   scales the output of the Weibull CDF to have a max value of 1. The identity
#'   discount function uses the posterior probability directly as the discount
#'   weight. Default value is \code{"identity"}. See \code{\link{bdpnormal}} for
#'   more details.
#' @param alternative character. The string specifying the alternative
#'   hypothesis, must be one of \code{"greater"} (default), \code{"less"} or
#'   \code{"two.sided"}.
#' @param prop_loss_to_followup scalar. Overall oroportion of subjects lost to
#'   follow-up.
#' @param h0 scalar. Threshold for comparing two mean values. Default is
#'   \code{h0=0}.
#' @param futility_prob scalar. Probability of stopping early for futility.
#' @param expected_success_prob scalar. Probability of stopping early for success.
#' @param prob_ha scalar. Probability of alternative hypothesis.
#' @param N_impute scalar. Number of imputations for Monte Carlo simulation of
#'   missing data.
#' @param number_mcmc scalar. Number of Monte Carlo Markov Chain draws in
#'   sampling posterior.
#' @param alpha_max scalar. Maximum weight the discount function can apply.
#'   Default is 1. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to weight the historical treatment group and
#'   the second value is used to weight the historical control group.
#' @param fix_alpha logical. Fix alpha at alpha_max? Default value is FALSE.
#' @param weibull_scale scalar. Scale parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 0.135. For a two-arm trial, users may specify a vector of two
#'   values where the first value is used to estimate the weight of the
#'   historical treatment group and the second value is used to estimate the
#'   weight of the historical control group. Not used when discount_function =
#'   "identity".
#' @param weibull_shape scalar. Shape parameter of the Weibull discount function
#'   used to compute alpha, the weight parameter of the historical data. Default
#'   value is 3. For a two-arm trial, users may specify a vector of two values
#'   where the first value is used to estimate the weight of the historical
#'   treatment group and the second value is used to estimate the weight of the
#'   historical control group. Not used when discount_function = "identity".
#' @param method character. Analysis method with respect to estimation of the weight
#'   paramter alpha. Default method "\code{mc}" estimates alpha for each
#'   Monte Carlo iteration. Alternate value "\code{fixed}" estimates alpha once
#'   and holds it fixed throughout the analysis.  See the the
#'   \code{bdpsurvival} vignette \cr
#'   \code{vignette("bdpsurvival-vignette", package="bayesDP")} for more details.
#'
#'
#' @return a list of output for a single trial simulation.
#' \describe{
#'   \item{\code{mu_treatment}}{
#'     scalar. The input parameter of mean value of the outcome in the
#'     treatment group.}
#'   \item{\code{p_control}}{
#'     scalar. The input parameter of mean value of the outcome in the
#'     control group.}
#'   \item{\code{sd_treatment}}{
#'     scalar. The input parameter of standard deviation of the outcome
#'     in the control group.}
#'   \item{\code{sd_control}}{
#'     scalar. The input parameter of standard deviation of the outcome
#'     in the control group.}
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
#'   \item{\code{est_interim}}{
#'     scalar. The interim estimate of the difference in posterior estimate of
#'     treatment and posterior estimate of the control group.}
#' }
#'
#'
#' @importFrom stats rnorm lm sd
#' @importFrom dplyr mutate filter group_by bind_rows select n
#' @importFrom bayesDP bdpnormal
#' @export normalBACT
normalBACT <- function(
  mu_treatment,
  sd_treatment,
  mu_control            = NULL,
  sd_control            = NULL,
  mu0_treatment         = NULL,
  sd0_treatment         = NULL,
  N0_treatment          = NULL,
  mu0_control           = NULL,
  sd0_control           = NULL,
  N0_control            = NULL,
  N_total,
  lambda                = 0.3,
  lambda_time           = NULL,
  interim_look          = NULL,
  EndofStudy,
  block                 = 2,            # block size for randomization
  rand_ratio            = c(1, 1),      # randomization ratio in control to treatament (default 1:1)
  discount_function     = "identity",   # discount_function used in sampling
  alternative           = "greater",  # the alternative hypothesis (either two-sided, greater, less)
  prop_loss_to_followup = 0.15,         # Proportion of loss in data
  h0                    = 0,            # Null hypothesis value
  futility_prob         = 0.05,         # Futility probability
  expected_success_prob = 0.9,          # Expected success probability
  prob_ha               = 0.95,         # Posterior probability of accepting alternative hypothesis
  N_impute              = 10,          # Number of imputation simulations for predictive distribution
  number_mcmc           = 10000,         # Number of posterior sampling
  alpha_max             = 1,            # max weight on incorporating historical data
  fix_alpha             = FALSE,        # fix alpha set weight of historical data to alpha_max
  weibull_scale         = 0.135,        # weibull parameter
  weibull_shape         = 3,            # weibull parameter
  method                = "fixed"
 ){
  # checking inputs
  stopifnot((mu_treatment > 0 & sd_treatment > 0),
            all(N_total > interim_look), length(lambda) == (length(lambda_time) + 1),
            EndofStudy > 0, block %% sum(rand_ratio)  == 0,
            (prop_loss_to_followup >= 0 & prop_loss_to_followup < 0.75),
            (futility_prob < 0.20 & futility_prob >= 0),
            (expected_success_prob > 0.70 & expected_success_prob <= 1),
            (prob_ha > 0.70 & prob_ha < 1), N_impute > 0)

  #checking if alternative is right
  if(alternative != "two-sided" & alternative != "greater" & alternative != "less"){
    stop("The input for alternative is wrong!")
  }

  # assigining interim look and final look
  analysis_at_enrollnumber <- c(interim_look, N_total)

  # assignment of enrollment based on the enrollment function
  enrollment <- enrollment(param = lambda, N_total = N_total, time = lambda_time)


  # simulating group and treatment group assignment
  if(!is.null(mu_control)){
    group <- randomization(N_total = N_total, block = block, allocation = rand_ratio)
  }
  else{
    group <- rep(1, N_total)
  }

  # simulate normal outcome
  if(!is.null(mu_control)){
    sim <- rnorm(N_total, mean = group * mu_treatment + (1 - group) * mu_control,
                 sd = group * sd_treatment + (1 - group) * sd_control)
    # dividing treatment and control
    control <- sim[group == 0]
  }
  else{
    sim <- rnorm(N_total, mean = mu_treatment, sd = sd_treatment)
  }
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

  if(length(analysis_at_enrollnumber) > 1){
  for(i in 1:(length(analysis_at_enrollnumber)-1)){

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
      mutate(subject_impute_success = (enrollment[analysis_at_enrollnumber[i]] - enrollment <= EndofStudy & subject_enrolled) |
               (subject_enrolled & loss_to_fu))

    # Carry out interim analysis on patients with complete data only
    # - Set-up `new data` data frame
    data <- data_interim %>%
      filter(subject_enrolled,
             !subject_impute_success)

    # MLE of data at interim analysis
    #MLE_int <- lm(Y ~ treatment, data = data)

    # assigning input for control arm given it is a single or double arm
    if(!is.null(mu_control)){
      mu_c     <- mean(data$Y[data$treatment == 0])
      sd_c  <- sd(data$Y[data$treatment == 0])
      N_c      <- length(data$treatment == 0)
    }
    else{
      mu_c     <- NULL
      sd_c  <- NULL
      N_c      <- NULL
    }

    # analyze data using discount funtion via normal
    post <- bdpnormal(mu_t                = mean(data$Y[data$treatment == 1]),
                      sigma_t             = sd(data$Y[data$treatment == 1]),
                      N_t                 = length(data$treatment == 1),
                      mu_c                = mu_c,
                      sigma_c             = sd_c,
                      N_c                 = N_c,
                      mu0_t               = mu0_treatment,
                      sigma0_t            = sd0_treatment,
                      N0_t                = N0_treatment,
                      mu0_c               = mu0_control,
                      sigma0_c            = sd0_control,
                      N0_c                = N0_control,
                      number_mcmc         = number_mcmc,
                      discount_function   = discount_function,
                      alpha_max           = alpha_max,
                      fix_alpha           = fix_alpha,
                      weibull_scale       = weibull_scale,
                      weibull_shape       = weibull_shape,
                      method              = method)

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

      # assigning input for control arm given it is a single or double arm
      if(!is.null(mu_control)){
        mu_c     <- mean(data$Y[data$treatment == 0])
        sd_c     <- sd(data$Y[data$treatment == 0])
        N_c      <- length(data$treatment == 0)
      }
      else{
        mu_c     <- NULL
        sd_c     <- NULL
        N_c      <- NULL
      }

      # analyze complete+imputed data using discount funtion via normal
      post_imp <- bdpnormal(mu_t              = mean(data$Y[data$treatment == 1]),
                            sigma_t           = sd(data$Y[data$treatment == 1]),
                            N_t               = length(data$treatment == 1),
                            mu_c              = mu_c,
                            sigma_c           = sd_c,
                            N_c               = N_c,
                            mu0_t             = mu0_treatment,
                            sigma0_t          = sd0_treatment,
                            N0_t              = N0_treatment,
                            mu0_c             = mu0_control,
                            sigma0_c          = sd0_control,
                            N0_c              = N0_control,
                            number_mcmc       = number_mcmc,
                            discount_function = discount_function,
                            alpha_max         = alpha_max,
                            fix_alpha         = fix_alpha,
                            weibull_scale     = weibull_scale,
                            weibull_shape     = weibull_shape,
                            method            = method)


      # Estimation of the posterior effect for difference between test and control
      # - If expected success, add 1 to the counter

      if(!is.null(mu_control)){
        if(alternative == "two-sided"){
          effect_imp <- post_imp$posterior_treatment$posterior_mu - post_imp$posterior_control$posterior_mu
          success <- max(c(mean(effect_imp > h0), mean(-effect_imp > h0)))
        }
        else if(alternative == "greater"){
          effect_imp <- post_imp$posterior_treatment$posterior_mu - post_imp$posterior_control$posterior_mu
          success <- mean(effect_imp > h0)
        }
        else{
          effect_imp <- post_imp$posterior_treatment$posterior_mu - post_imp$posterior_control$posterior_mu
          success <- mean(-effect_imp > h0)
        }
      }

      else{
        effect_imp <- post_imp$posterior_treatment$posterior_mu
        if(alternative == "two-sided"){
          success <- max(c(mean(effect_imp - mu_treatment > h0), mean(mu_treatment - effect_imp > h0)))
        }
        else if(alternative == "greater"){
          success <- mean(effect_imp - mu_treatment > h0)
        }
        else{
          success <- mean(mu_treatment - effect_imp > h0)
        }
      }

      if(success > prob_ha){
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

      # assigning input for control arm given it is a single or double arm
      if(!is.null(mu_control)){
        mu_c     <- mean(data$Y[data$treatment == 0])
        sd_c     <- sd(data$Y[data$treatment == 0])
        N_c      <- length(data$treatment == 0)
      }
      else{
        mu_c     <- NULL
        sd_c     <- NULL
        N_c      <- NULL
      }

      # Analyze complete+imputed data using discount funtion via normal
      post_imp <- bdpnormal(mu_t                = mean(data$Y[data$treatment == 1]),
                            sigma_t             = sd(data$Y[data$treatment == 1]),
                            N_t                 = length(data$treatment == 1),
                            mu_c                = mu_c,
                            sigma_c             = sd_c,
                            N_c                 = N_c,
                            mu0_t               = mu0_treatment,
                            sigma0_t            = sd0_treatment,
                            N0_t                = N0_treatment,
                            mu0_c               = mu0_control,
                            sigma0_c            = sd0_control,
                            N0_c                = N0_control,
                            number_mcmc         = number_mcmc,
                            discount_function   = discount_function,
                            alpha_max           = alpha_max,
                            fix_alpha           = fix_alpha,
                            weibull_scale       = weibull_scale,
                            weibull_shape       = weibull_shape,
                            method              = method)


      # Estimation of the posterior effect for difference between test and
      # control

      if(!is.null(mu_control)){
        if(alternative == "two-sided"){
          effect_imp <- post_imp$posterior_treatment$posterior_mu - post_imp$posterior_control$posterior_mu
          success <- max(c(mean(effect_imp > h0), mean(-effect_imp > h0)))
        }
        else if(alternative == "greater"){
          effect_imp <- post_imp$posterior_treatment$posterior_mu - post_imp$posterior_control$posterior_mu
          success <- mean(effect_imp > h0)
        }
        else{
          effect_imp <- post_imp$posterior_treatment$posterior_mu - post_imp$posterior_control$posterior_mu
          success <- mean(-effect_imp > h0)
        }
      }

      else{
        effect_imp <- post_imp$posterior_treatment$posterior_mu
        if(alternative == "two-sided"){
          success <- max(c(mean(effect_imp - mu_treatment > h0), mean(mu_treatment - effect_imp > h0)))
        }
        else if(alternative == "greater"){
          success <- mean(effect_imp - mu_treatment > h0)
        }
        else{
          success <- mean(mu_treatment - effect_imp > h0)
        }
      }

      # Increase futility counter by 1 if P(effect_imp < h0) > ha
      if(success > prob_ha){
        futility_test <- futility_test + 1
      }

    }

    #print(analysis_at_enrollnumber[i])

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
    if(!is.null(mu_control)){
      if(alternative == "two-sided"){
        effect_int <- post$posterior_treatment$posterior_mu - post$posterior_control$posterior_mu
      }
      else if(alternative == "greater"){
        effect_int <- post$posterior_treatment$posterior_mu - post$posterior_control$posterior_mu
      }
      else{
        effect_int <- post$posterior_treatment$posterior_mu - post$posterior_control$posterior_mu
      }
    }

    else{
      effect_int <- post$posterior_treatment$posterior_mu
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

  #print(N_enrolled)

  # All patients that have made it to the end of study
  # - Subset out patients loss to follow-up
  data_final <- data_total %>%
    filter(id <= stage_trial_stopped,
           !loss_to_fu)


  # Compute the final MLE for the complete data using GLM function
  # MLE <- lm(Y ~ treatment, data = data_final)

  # assigning input for control arm given it is a single or double arm
  if(!is.null(mu_control)){
    mu_c     <- mean(data_final$Y[data_final$treatment == 0])
    sd_c     <- sd(data_final$Y[data_final$treatment == 0])
    N_c      <- length(data_final$treatment == 0)
  }
  else{
    mu_c     <- NULL
    sd_c     <- NULL
    N_c      <- NULL
  }

  # Analyze complete data using discount funtion via binomial
  post_final <- bdpnormal(mu_t                = mean(data_final$Y[data_final$treatment == 1]),
                          sigma_t             = sd(data_final$Y[data_final$treatment == 1]),
                          N_t                 = length(data_final$treatment == 1),
                          mu_c                = mu_c,
                          sigma_c             = sd_c,
                          N_c                 = N_c,
                          mu0_t               = mu0_treatment,
                          sigma0_t            = sd0_treatment,
                          N0_t                = N0_treatment,
                          mu0_c               = mu0_control,
                          sigma0_c            = sd0_control,
                          N0_c                = N0_control,
                          discount_function   = discount_function,
                          number_mcmc         = number_mcmc,
                          alpha_max           = alpha_max,
                          fix_alpha           = fix_alpha,
                          weibull_scale       = weibull_scale,
                          weibull_shape       = weibull_shape,
                          method              = method)


  ### Format and output results
  # Posterior effect size: test vs control or treatment itself
  if(!is.null(mu_control)){
    if(alternative == "two-sided"){
      effect <- post_final$posterior_treatment$posterior_mu - post_final$posterior_control$posterior_mu
      post_paa <- max(c(mean(effect > h0), mean(-effect > h0)))
    }
    else if(alternative == "greater"){
      effect <- post_final$posterior_treatment$posterior_mu - post_final$posterior_control$posterior_mu
      post_paa <- mean(effect > h0)
    }
    else{
      effect <- post_final$posterior_treatment$posterior_mu - post_final$posterior_control$posterior_mu
      post_paa <- mean(-effect > h0)
    }
  }

  else{
    effect <- post_final$posterior_treatment$posterior_mu
    if(alternative == "two-sided"){
      post_paa <- max(c(mean(effect - mu_treatment > h0), mean(mu_treatment - effect > h0)))
    }
    else if(alternative == "greater"){
      post_paa <- mean(effect - mu_treatment > h0)
    }
    else{
      post_paa <- mean(mu_treatment - effect > h0)
    }
  }

  N_treatment  <- sum(data_final$treatment)         # Total sample size analyzed - test group
  N_control    <- sum(!data_final$treatment)        # Total sample size analyzed - control group


  # output
  results_list <- list(
    mu_treatment                               = mu_treatment,             # mean of treatment in normal
    mu_control                                 = mu_control,               # mean of control in normal
    sd_treatment                               = sd_treatment,             # sd of treatment in normal
    sd_control                                 = sd_control,               # sd of control in normal
    prob_of_accepting_alternative              = prob_ha,
    margin                                     = h0,                       # margin for error
    alternative                                = alternative,              # alternative hypothesis
    interim_look                               = interim_look,             # print interim looks
    N_treatment                                = N_treatment,
    N_control                                  = N_control,
    N_complete                                 = N_treatment + N_control,
    N_enrolled                                 = N_enrolled,               # Total sample size enrolled when trial stopped
    N_max                                      = N_total, 				         # Total potential sample size
    post_prob_accept_alternative               = post_paa,                 # Posterior probability that alternative hypothesis is true
    est_final                                  = mean(effect),             # Posterior Mean of treatment effect
    stop_futility                              = stop_futility,            # Did the trial stop for futility
    stop_expected_success                      = stop_expected_success     # Did the trial stop for expected success
    #MLE_est                                   = MLE$coe[2],               # Treatment effect useing MLE
    #MLE_est_interim                           = MLE_int$coe[2]            # Treatment effect useing MLE at interim analysis
  )

  #if there is an interim look
  if(length(analysis_at_enrollnumber) > 1){
    results_list                            <- c(results_list,
    est_interim                             = mean(effect_int))           # Posterior Mean of treatment effect at interim analysis
    }

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
#' @param mu_control numeric. The mean for the control group.
#' @param sd_control numeric. The standard deviation for the control group.
#' @param mu_treatment numeric. The mean for the treatment group.
#' @param sd_treatment numeric. The standard deviation for the treatment group.
#' @param data NULL. stores the proportion of control and treatment, please do not fill it in.
#'
#' @return a list with proportion of control and treatment group.
#'
#' @examples normal_outcome(mu_control = 12, mu_treatment = 8, sd_treatment = 2.2, sd_control = 1.6)
#' @export normal_outcome
#'
normal_outcome <- function(mu_control = NULL, sd_control = NULL, mu_treatment = NULL,
                           sd_treatment = NULL, data = NULL){
  data$mu_control <- mu_control
  data$sd_control <- sd_control
  data$mu_treatment <- mu_treatment
  data$sd_treatment <- sd_treatment
  data
}


#' @title Historical data for normal distribution
#'
#' @description Wrapper function for historical data from normal outcome.
#'
#' @inheritParams normalBACT
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
                              alpha_max           = 1,            # max weight on incorporating historical data
                              fix_alpha           = FALSE,        # fix alpha set weight of historical data to alpha_max
                              weibull_scale       = 0.135,        # weibull parameter
                              weibull_shape       = 3,            # weibull parameter
                              method              = "fixed",
                              .data               = NULL){
  .data$mu0_treatment       <- mu0_treatment
  .data$sd0_treatment       <- sd0_treatment
  .data$N0_treatment        <- N0_treatment
  .data$mu0_control         <- mu0_control
  .data$sd0_control         <- sd0_control
  .data$N0_control          <- N0_control
  .data$discount_function   <- discount_function
  .data$alpha_max           <- alpha_max
  .data$fix_alpha           <- fix_alpha
  .data$weibull_scale       <- weibull_scale
  .data$weibull_shape       <- weibull_shape
  .data$method              <- method
  .data
}


#' @title Analyzing Bayesian trial for normal mean data
#'
#' @description Function to analyze Bayesian trial for normal mean data
#'  which allows early stopping and incorporation of historical data using
#'  the discount function approach
#'
#' @inheritParams normalBACT
#' @param treatment vector. treatment assignment for patients, 1 for treatment group and
#'    0 for control group
#' @param outcome vector. normal outcome of the trial.
#' @param complete vector. similar length as treatment and outcome variable,
#'    1 for complete outcome, 0 for loss to follow up. If complete is not provided,
#'    the dataset is assumed to be complete.
#'
#' @importFrom stats rnorm lm
#' @importFrom dplyr mutate filter group_by bind_rows select n summarize
#' @importFrom bayesDP bdpnormal
#'
#' @return a list of output for the analysis of Bayesian trial for normal mean.
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
#' @export normal_analysis

normal_analysis <- function(
  treatment,
  outcome,
  complete              = NULL,
  mu0_treatment         = NULL,
  sd0_treatment         = NULL,
  N0_treatment          = NULL,
  mu0_control           = NULL,
  sd0_control           = NULL,
  N0_control            = NULL,
  alternative           = "greater",
  N_impute              = 100,
  h0                    = 0,
  number_mcmc           = 10000,
  prob_ha               = 0.95,
  futility_prob         = 0.10,
  expected_success_prob = 0.90,
  discount_function     = "identity",
  fix_alpha             = FALSE,
  alpha_max             = 1,
  weibull_scale         = 0.135,
  weibull_shape         = 3,
  method                = "fixed"
){

  #if complete is NULL, assume the data is complete
  if(is.null(complete)){
    complete <- rep(1, length(outcome))
  }

  ##reading the data
  data_total <- data.frame(cbind(treatment, outcome, complete))

  data_interim <- data_total %>%
    mutate(futility = (complete == 0))

  data <- data_interim %>%
    filter(!futility)

  summary <- data %>%
    group_by(treatment) %>%
    summarize(mu_outcome = mean(outcome), sd_outcome = sd(outcome))


  if(sum(data$treatment == 0) != 0){
    mu_c <- mean(data$outcome[data$treatment == 0])
    sd_c <- sd(data$outcome[data$treatment == 0])
    N_c  <- length(data$outcome[data$treatment == 0])
  }
  else{
    mu_c <- NULL
    sd_c <- NULL
    N_c  <- NULL
  }

  # analyze the data using bayesDp
  post <-  bdpnormal(mu_t                = mean(data$outcome[data$treatment == 1]),
                     sigma_t             = sd(data$outcome[data$treatment == 1]),
                     N_t                 = length(data$treatment == 1),
                     mu_c                = mu_c,
                     sigma_c             = sd_c,
                     N_c                 = N_c,
                     mu0_t               = mu0_treatment,
                     sigma0_t            = sd0_treatment,
                     N0_t                = N0_treatment,
                     mu0_c               = mu0_control,
                     sigma0_c            = sd0_control,
                     N0_c                = N0_control,
                     number_mcmc         = number_mcmc,
                     discount_function   = discount_function,
                     alpha_max           = alpha_max,
                     fix_alpha           = fix_alpha,
                     weibull_scale       = weibull_scale,
                     weibull_shape       = weibull_shape,
                     method              = method)


  # assigning stop_futility and expected success
  stop_futility         <- 0
  stop_expected_success <- 0
  expected_success_test <- 0

  for(i in 1:N_impute){
    data_control_success_impute <- data_interim %>%
      filter(treatment == 0) %>%
      mutate(outcome_impute = ifelse(futility,
                                     rnorm(n(), summary$mu_outcome[1], summary$sd_outcome[1]),
                                     outcome))
    # imputing success for treatment group
    data_treatment_success_impute  <- data_interim %>%
      filter(treatment == 1) %>%
      mutate(outcome_impute = ifelse(futility,
                                     rnorm(n(), summary$mu_outcome[2], summary$sd_outcome[2]),
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
      mu_c <- mean(data$outcome[data$treatment == 0])
      sd_c <- sd(data$outcome[data$treatment == 0])
      N_c  <- length(data$outcome[data$treatment == 0])
    }
    else{
      mu_c <- NULL
      sd_c <- NULL
      N_c  <- NULL
    }

    # analyze complete+imputed data using discount funtion via binomial
    post_imp <- bdpnormal(mu_t                = mean(data$outcome[data$treatment == 1]),
                          sigma_t             = sd(data$outcome[data$treatment == 1]),
                          N_t                 = length(data$treatment == 1),
                          mu_c                = mu_c,
                          sigma_c             = sd_c,
                          N_c                 = N_c,
                          mu0_t               = mu0_treatment,
                          sigma0_t            = sd0_treatment,
                          N0_t                = N0_treatment,
                          mu0_c               = mu0_control,
                          sigma0_c            = sd0_control,
                          N0_c                = N0_control,
                          number_mcmc         = number_mcmc,
                          discount_function   = discount_function,
                          alpha_max           = alpha_max,
                          fix_alpha           = fix_alpha,
                          weibull_scale       = weibull_scale,
                          weibull_shape       = weibull_shape,
                          method              = method)

    if(sum(data$treatment == 0) != 0){
      if(alternative == "two-sided"){
        effect_imp <- post_imp$posterior_treatment$posterior_mu - post_imp$posterior_control$posterior_mu
        success <- max(c(mean(effect_imp > h0), mean(-effect_imp > h0)))
      }
      else if(alternative == "greater"){
        effect_imp <- post_imp$posterior_treatment$posterior_mu - post_imp$posterior_control$posterior_mu
        success <- mean(effect_imp > h0)
      }
      else{
        effect_imp <- post_imp$posterior_treatment$posterior_mu - post_imp$posterior_control$posterior_mu
        success <- mean(-effect_imp > h0)
      }
    }

    else{
      effect_imp <- post_imp$final$posterior_mu
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
    mu_c <- mean(data_final$outcome[data_final$treatment == 0])
    sd_c <- sd(data_final$outcome[data_final$treatment == 0])
    N_c  <- length(data_final$outcome[data_final$treatment == 0])
  }
  else{
    mu_c <- NULL
    sd_c <- NULL
    N_c  <- NULL
  }

  # Analyze complete data using discount funtion via binomial
  post_final <- bdpnormal(mu_t                = mean(data_final$outcome[data_final$treatment == 1]),
                          sigma_t             = sd(data_final$outcome[data_final$treatment == 1]),
                          N_t                 = length(data_final$treatment == 1),
                          mu_c                = mu_c,
                          sigma_c             = sd_c,
                          N_c                 = N_c,
                          mu0_t               = mu0_treatment,
                          sigma0_t            = sd0_treatment,
                          N0_t                = N0_treatment,
                          mu0_c               = mu0_control,
                          sigma0_c            = sd0_control,
                          N0_c                = N0_control,
                          discount_function   = discount_function,
                          number_mcmc         = number_mcmc,
                          alpha_max           = alpha_max,
                          fix_alpha           = fix_alpha,
                          weibull_scale       = weibull_scale,
                          weibull_shape       = weibull_shape,
                          method              = method)

  ### Format and output results
  # Posterior effect size: test vs control or treatment itself
  if(sum(data_final$treatment == 0) != 0){
    if(alternative == "two-sided"){
      effect <- post_final$posterior_treatment$posterior_mu - post_final$posterior_control$posterior_mu
      post_paa <- max(c(mean(effect > h0), mean(-effect > h0)))
    }
    else if(alternative == "greater"){
      effect <- post_final$posterior_treatment$posterior_mu - post_final$posterior_control$posterior_mu
      post_paa <- mean(effect > h0)
    }
    else{
      effect <- post_final$posterior_treatment$posterior_mu - post_final$posterior_control$posterior_mu
      post_paa <- mean(-effect > h0)
    }
  }

  else{
    effect <- post_final$final$posterior_mu
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
                                                        "subject_impute_success", "p_outcome"))




#' @title Data file for normal analysis
#'
#' @description Wrapper function for data file in normal analysis.
#'
#' @param treatment vector. treatment assignment for patients, 1 for treatment group and
#'    0 for control group
#' @param outcome vector. normal outcome of the trial.
#' @param complete vector. similar length as treatment and outcome variable,
#'    1 for complete outcome, 0 for loss to follow up. If complete is not provided,
#'    the dataset is assumed to be complete.
#' @param .data NULL. stores the normal data for analysis, please do not fill it in.
#'
#' @return a list with treatment, outcome and loss to follow up vector with normal
#'   outcome.
#'
#' @export data_normal
data_normal <- function(treatment, outcome, complete, .data = NULL){
  .data$treatment <- treatment
  .data$outcome   <- outcome
  .data$complete  <- complete
  .data
}





