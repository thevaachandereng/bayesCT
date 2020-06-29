#' @title Binomial counts for Bayesian Adaptive Trials
#'
#' @description Simulation for binomial counts for Bayesian Adaptive trial with
#'   different inputs to control for power, sample size, type 1 error rate, etc.
#' @param p_treatment scalar. Proportion of events under the treatment arm.
#' @param p_control scalar. Proportion of events under the control arm.
#' @param y0_treatment scalar. Number of events for the historical treatment
#'   arm.
#' @param N0_treatment scalar. Sample size of the historical treatment arm.
#' @param y0_control scalar. Number of events for the historical control arm.
#' @param N0_control scalar. Sample size of the historical control arm.
#' @param prior vector. Prior values of beta rate, Beta(a0, b0). The default is
#'   set to Beta(1, 1).
#' @inheritParams normalBACT
#'
#' @return a list of output for a single trial simulation.
#' \describe{
#'   \item{\code{p_treatment}}{
#'     scalar. The input parameter of proportion of events in the
#'     treatment group.}
#'   \item{\code{p_control}}{
#'     scalar. The input parameter of proportion of events in the
#'     control group.}
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
#' @importFrom stats rbinom glm
#' @importFrom dplyr mutate filter group_by bind_rows select n
#' @importFrom bayesDP bdpbinomial
#' @export binomialBACT
binomialBACT <- function(
  p_treatment,
  p_control             = NULL,
  y0_treatment          = NULL,
  N0_treatment          = NULL,
  y0_control            = NULL,
  N0_control            = NULL,
  N_total,
  lambda                = 0.3,
  lambda_time           = NULL,
  interim_look          = NULL,
  EndofStudy,
  prior                 = c(1, 1),
  block                 = 2,            # block size for randomization
  rand_ratio            = c(1, 1),      # randomization ratio in control to treatament (default 1:1)
  prop_loss_to_followup = 0.10,         # Proportion of loss in data
  alternative           = "greater",    # the alternative hypothesis (either two-sided, greater, less)
  h0                    = 0,            # Null hypothesis value
  futility_prob         = 0.05,         # Futility probability
  expected_success_prob = 0.9,          # Expected success probability
  prob_ha               = 0.95,         # Posterior probability of accepting alternative hypothesis
  N_impute              = 10,           # Number of imputation simulations for predictive distribution
  number_mcmc           = 10000,        # Number of posterior sampling
  discount_function     = "identity",
  alpha_max             = 1,            # max weight on incorporating historical data
  fix_alpha             = FALSE,        # fix alpha set weight of historical data to alpha_max
  weibull_scale         = 0.135,        # weibull parameter
  weibull_shape         = 3,            # weibull parameter
  method                = "fixed"
  ){

  # checking p_control
  if(!is.null(p_control)){stopifnot(p_control > 0 & p_control < 1)}

  # checking historical data for treatment group
  if(!is.null(y0_treatment) | !is.null(N0_treatment)){
    stopifnot((y0_treatment > 0 & N0_treatment > 0), y0_treatment <= N0_treatment,
              (discount_function == "identity" | discount_function == "weibull" |
               discount_function == "scaledweibull"))
  }

  # checking historical data for control group
  if(!is.null(y0_control) | !is.null(N0_control)){
    stopifnot((y0_control > 0 & N0_control > 0), y0_control <= N0_control)
  }

  # checking interim_look
  if(!is.null(interim_look)){
    stopifnot(all(N_total > interim_look))
  }

  # checking other inputs
  stopifnot((p_treatment < 1 & p_treatment > 0),
            length(lambda) == (length(lambda_time) + 1),
            EndofStudy > 0, block %% sum(rand_ratio)  == 0,
            (prop_loss_to_followup >= 0 & prop_loss_to_followup < 0.75),
            (h0 >= -1 & h0 < 1), (futility_prob < 0.20 & futility_prob >= 0),
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
  if(!is.null(p_control)){
    group <- randomization(N_total = N_total, block = block, allocation = rand_ratio)
  }
  else{
    group <- rep(1, N_total)
  }

  # simulate binomial outcome
  if(!is.null(p_control)){
    sim <- rbinom(N_total, 1, prob = group * p_treatment + (1 - group) * p_control)
    # dividing treatment and control
    control <- sim[group == 0]
  }
  else{
    sim <- rbinom(N_total, 1, p_treatment)
  }
  treatment <- sim[group == 1]

  # simulate loss to follow-up
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
      mutate(subject_impute_success = (enrollment[analysis_at_enrollnumber[i]] - enrollment <= EndofStudy & subject_enrolled) |
               (subject_enrolled & loss_to_fu))

    # Carry out interim analysis on patients with complete data only
    # - Set-up `new data` data frame
    data <- data_interim %>%
      filter(subject_enrolled,
             !subject_impute_success)

    # MLE of data at interim analysis
    # MLE_int <- glm(Y ~ treatment, data = data, family = "binomial")

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
                        weibull_shape          = weibull_shape,
                        method                 = method)

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

      # combine the treatment and control imputed datasets
      data_success_impute <- bind_rows(data_control_success_impute,
                                       data_treatment_success_impute) %>%
        mutate(Y = Y_impute) %>%
        select(-Y_impute)

      # create enrolled subject data frame for discount function analysis
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

      # analyze complete+imputed data using discount funtion via binomial
      post_imp <- bdpbinomial(y_t                    = sum(data$Y[data$treatment == 1]),
                              N_t                    = length(data$Y[data$treatment == 1]),
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
                              weibull_shape          = weibull_shape,
                              method                 = method)

      # estimation of the posterior effect for difference between test and
      # control - If expected success, add 1 to the counter

      if(!is.null(p_control)){
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
          success <- max(c(mean(effect_imp - p_treatment > h0), mean(p_treatment - effect_imp > h0)))
        }
        else if(alternative == "greater"){
          success <- mean(effect_imp - p_treatment > h0)
        }
        else{
          success <- mean(p_treatment - effect_imp > h0)
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


      # Create enrolled subject data frame for discount function analysis
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
                              weibull_shape          = weibull_shape,
                              method                 = method)

      # Estimation of the posterior effect for difference between test and
      # control
      if(!is.null(p_control)){
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
          success <- max(c(mean(effect_imp - p_treatment > h0), mean(p_treatment - effect_imp > h0)))
        }
        else if(alternative == "greater"){
          success <- mean(effect_imp - p_treatment > h0)
        }
        else{
          success <- mean(p_treatment - effect_imp > h0)
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
    if(!is.null(p_control)){
      if(alternative == "two-sided"){
        effect_int <- post$posterior_treatment$posterior - post$posterior_control$posterior
      }
      else if(alternative == "greater"){
        effect_int <- post$posterior_treatment$posterior - post$posterior_control$posterior
      }
      else{
        effect_int <- post$posterior_treatment$posterior - post$posterior_control$posterior
      }
    }

    else{
      effect_int <- post$final$posterior
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
  post_final <- bdpbinomial(y_t                  = sum(data_final$Y[data_final$treatment == 1]),
                            N_t                  = length(data_final$Y[data_final$treatment == 1]),
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
                            weibull_shape        = weibull_shape,
                            method               = method)

  ### Format and output results
  # Posterior effect size: test vs control or treatment itself
  if(!is.null(p_control)){
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
      post_paa <- max(c(mean(effect - p_treatment > h0), mean(p_treatment - effect > h0)))
    }
    else if(alternative == "greater"){
      post_paa <- mean(effect - p_treatment > h0)
    }
    else{
      post_paa <- mean(p_treatment - effect > h0)
    }
  }



  N_treatment  <- sum(data_final$treatment)         # Total sample size analyzed - test group
  N_control    <- sum(!data_final$treatment)        # Total sample size analyzed - control group

  ## output
  results_list <- list(
    p_treatment                                = p_treatment,              # probability of treatment in binomial
    p_control                                  = p_control,                # probability of control in binomial
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

  if(length(analysis_at_enrollnumber) > 1){
    results_list                            <- c(results_list,
    est_interim                             = mean(effect_int))           # Posterior Mean of treatment effect at interim analysis

  }

  #return results
  results_list

}


## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Y", "Y_impute", "id", "subject_enrolled",
"subject_impute_success", "subject_impute_futility"))


#' @title Proportion of an event in control and treatment
#'
#' @description Wrapper function for proportion of an event in control and treatment group with binomial outcome.
#'
#' @param p_treatment numeric. The proportion of an event in the treatment group, 0 < $p_treatment$ < 1.
#' @param p_control numeric. The proportion of an event in the control group, 0 < $p_control$ < 1.
#' @param .data NULL. stores the proportion of control and treatment, please do not fill it in.
#'
#' @return a list with proportion of control and treatment group.
#'
#' @examples binomial_outcome(p_control = 0.12, p_treatment = 0.08)
#' @export binomial_outcome
binomial_outcome <- function(p_treatment = NULL, p_control = NULL, .data = NULL){
  .data$p_treatment  <- p_treatment
  .data$p_control    <- p_control
  .data
}


#' @title Historical data for binomial distribution
#'
#' @description Wrapper function for historical data from binomial outcome.
#'
#' @inheritParams normalBACT
#' @inheritParams binomialBACT
#' @param .data NULL. stores the proportion of control and treatment, please do not fill it in.
#'
#' @return a list with historical data for control and treatment group with the discount function.
#'
#' @examples
#' historical_binomial(y0_treatment = 5, N0_treatment = 10, y0_control = 15, N0_control = 23)
#' historical_binomial(y0_treatment = 5, N0_treatment = 10, y0_control = 15, N0_control = 23,
#'                      discount_function = "weibull", alpha_max = 1, fix_alpha = FALSE,
#'                      weibull_scale = 0.135, weibull_shape = 3)
#' @export historical_binomial
historical_binomial <- function(y0_treatment       = NULL,
                                N0_treatment       = NULL,
                                discount_function  = "identity",
                                y0_control         = NULL,
                                N0_control         = NULL,
                                alpha_max          = 1,            # max weight on incorporating historical data
                                fix_alpha          = FALSE,        # fix alpha set weight of historical data to alpha_max
                                weibull_scale      = 0.135,        # weibull parameter
                                weibull_shape      = 3,            # weibull parameter
                                method             = "fixed",
                                .data              = NULL
                                ){
  .data$y0_treatment       <- y0_treatment
  .data$N0_treatment       <- N0_treatment
  .data$y0_control         <- y0_control
  .data$N0_control         <- N0_control
  .data$discount_function  <- discount_function
  .data$alpha_max          <- alpha_max
  .data$fix_alpha          <- fix_alpha
  .data$weibull_scale      <- weibull_scale
  .data$weibull_shape      <- weibull_shape
  .data$method             <- method
  .data
}


#' @title Beta prior for for control and treatment group
#'
#' @description Wrapper function for beta prior \code{beta(a0, b0)}.
#'
#' @param a0 numeric. The first shape parameter in the beta distribution
#'   (\code{beta(a0, b0)}).
#' @param b0 numeric. The second shape parameter in the beta distribution
#'   (\code{beta(a0, b0)}).
#' @param .data NULL. stores the beta prior rate, please do not fill it in.
#'
#' @return a list with vector of beta rate for the beta prior for treatment and control group.
#'
#' @examples beta_prior(a0 = 1, b0 = 1)
#' @export beta_prior
beta_prior <- function(a0 = 1, b0 = 1, .data = NULL){
  .data$prior  <- c(a0, b0)
  .data
}


#' @title Analyzing Bayesian trial for binomial counts
#'
#' @description Function to analyze Bayesian trial for binomial count data
#'  which allows early stopping and incorporation of historical data using
#'  the discount function approach
#'
#' @param treatment vector. treatment assignment for patients, 1 for treatment group and
#'    0 for control group
#' @param outcome vector. binomial outcome of the trial, 1 for response (success or
#'    failure), 0 for no response.
#' @param complete vector. similar length as treatment and outcome variable,
#'    1 for complete outcome, 0 for loss to follow up. If complete is not provided,
#'    the dataset is assumed to be complete.
#' @inheritParams normalBACT
#' @inheritParams binomialBACT
#'
#' @importFrom stats rbinom glm
#' @importFrom dplyr mutate filter group_by bind_rows select n summarize
#' @importFrom bayesDP bdpbinomial
#'
#' @return a list of output for the Bayesian trial for binomial count.
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
#' @export binomial_analysis

binomial_analysis <- function(
  treatment,
  outcome,
  complete              = NULL,
  y0_treatment          = NULL,
  N0_treatment          = NULL,
  y0_control            = NULL,
  N0_control            = NULL,
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
  weibull_shape         = 3,
  method                = "fixed"
){

  #if complete is NULL, assume the data is complete
  if(is.null(complete)){
   complete <- rep(1, length(outcome))
  }

  #reading the data
  data_total <- data.frame(cbind(treatment, outcome, complete))

  data_interim <- data_total %>%
    mutate(futility = (complete == 0))

  data <- data_interim %>%
    filter(!futility)

  prop <- data %>%
    group_by(treatment) %>%
    summarize(p_outcome = mean(outcome))


  if(sum(data$treatment == 0) != 0){
    y_c <- sum(data$outcome[data$treatment == 0])
    N_c <- length(data$outcome[data$treatment == 0])
  }
  else{
    y_c <- NULL
    N_c <- NULL
  }

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
                      weibull_shape          = weibull_shape,
                      method                 = method)


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
                            weibull_shape          = weibull_shape,
                            method                 = method)

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
                            weibull_shape        = weibull_shape,
                            method               = method)

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




#' @title Data file for binomial analysis
#'
#' @description Wrapper function for data file in binomial analysis.
#'
#' @inheritParams binomial_analysis
#' @param .data NULL. stores the binomial data for analysis, please do not fill it in.
#'
#' @return a list with treatment, outcome and loss to follow up vector with binomial
#'   outcome.
#'
#' @examples data_binomial(treatment = c(0, 1), outcome = c(1, 1), complete = c(1, 1))
#' @export data_binomial
data_binomial <- function(treatment, outcome, complete, .data = NULL){
  .data$treatment <- treatment
  .data$outcome   <- outcome
  .data$complete  <- complete
  .data
}


