#' @title Analysis wrapper function
#'
#' @description Wrapper function to analyze Bayesian trials.
#'
#' @param input list. Input function for all the analysis.
#' @param type character. Type of analysis to be ran (binomial (default),
#'   normal. etc.).
#' @param N_max_treatment integer. Maximum allowable sample size for the
#'   treatment arm (including the currently enrolled subjects). Default is NULL,
#'   meaning we are already at the final analysis.
#' @param N_max_control integer. Maximum allowable sample size for the control
#'   arm (including the currently enrolled subjects). Default is NULL, meaning
#'   we are already at the final analysis.
#' @param .data NULL. Stores the binomial data for analysis. Should not be
#'   edited by user.
#'
#' @return A list with results of the analysis of Bayesian trial.
#'
#' \describe{
#'   \item{\code{prob_of_accepting_alternative}}{
#'     scalar. The input parameter of probability of accepting the alternative.}
#'   \item{\code{margin}}{
#'     scalar. The margin input value of difference between mean estimate of treatment
#'      and mean estimate of the control.}
#'   \item{\code{alternative}}{
#'     character. The input parameter of alternative hypothesis.}
#'   \item{\code{N_treatment}}{
#'     scalar. The number of patients enrolled in the experimental group for
#'     each simulation.}
#'   \item{\code{N_control}}{
#'     scalar. The number of patients enrolled in the control group for
#'     each simulation.}
#'   \item{\code{N_enrolled}}{
#'     vector. The number of patients enrolled in the trial (sum of control
#'     and experimental group for each simulation.)}
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
#'     scalar. Did the trial stop for futility during imputation of patients
#'     who had loss to follow up? 1 for yes and 0 for no.}
#'   \item{\code{stop_expected_success}}{
#'     scalar. Did the trial stop for early success during imputation of patients
#'     who had loss to follow up? 1 for yes and 0 for no.}
#' }
#'
#' @importFrom stats rbinom glm rnorm lm runif
#' @importFrom dplyr mutate filter group_by bind_rows select n bind_cols
#' @importFrom bayesDP bdpbinomial bdpnormal bdpsurvival
#'
#' @export analysis

analysis <- function(input,
                     type            = "binomial",
                     N_max_treatment = NULL,
                     N_max_control   = NULL,
                     .data           = NULL
) {
  if (type == "binomial") {
    input <- c(input, "N_max_treatment" = N_max_treatment, "N_max_control" = N_max_control)
    do.call(binomial_analysis, input)
  } else if (type == "normal") {
    input <- c(input, "N_max_treatment" = N_max_treatment, "N_max_control" = N_max_control)
    do.call(normal_analysis, input)
  } else if (type == "survival") {
    do.call(survival_analysis, input)
  }
}
