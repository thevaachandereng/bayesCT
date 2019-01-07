#' @title Analysis wrapper function
#'
#' @description Wrapper function to analyze bayesian trials.
#'
#' @param input list. Input function for all the analysis.
#' @param type character. Type of analysis to be ran (binomial (default),
#'   normal. etc.).
#' @param .data NULL. stores the all the details, please do not fill it in.
#'
#' @return a list with results of the analysis of bayesian trial.
#'
#' @importFrom stats rbinom glm rnorm lm
#' @importFrom dplyr mutate filter group_by bind_rows select n
#' @importFrom bayesDP bdpbinomial bdpnormal
#'
#' @export analysis
#'

analysis <- function(input, type = "binomial", .data = NULL){
  if(type == "binomial"){
    do.call(binomial_analysis, input)
  }

  else if(type == "normal"){
    do.call(normal_analysis, input)
  }

}
