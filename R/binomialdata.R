#' Binomial dataset for analyzing adaptive Bayesian trials
#'
#' @description A dataset containing the results of 300 patients with binomial
#'   outcome, the dataset is filled with loss to follow up.
#'
#' @docType data
#'
#' @usage data(binomialdata)
#'
#' @format A data frame with 300 rows and 4 variables:
#'
#' \describe{
#'   \item{\code{id}}{
#'     Patient ID in the trial.}
#'   \item{\code{treatment}}{
#'     Treatment assignment for patients, 1 for treatment group 0 for control
#'     group.}
#'   \item{\code{outcome}}{
#'     Binomial outcome of the trial, 1 for response (success or failure), 0 for
#'     no response.}
#'   \item{\code{complete}}{
#'     1 for complete outcome, 0 for loss to follow-up.}
#' }
#'
#' @keywords datasets
#' @examples
#' data(binomialdata)

"binomialdata"
