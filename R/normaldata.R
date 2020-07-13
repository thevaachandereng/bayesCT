#' Gaussian dataset for analyzing adaptive Bayesian trials
#'
#' @description A dataset containing the results of 300 patients with continuous (normal) outcome,
#' the dataset is filled with loss to follow up.
#'
#' @docType data
#'
#' @keywords dataset
#'
#' @usage data(normaldata)
#'
#' @format A data frame with 300 rows and 4 variables:
#' \describe{
#'   \item{\code{id}}{
#'     Patient ID in the trial.}
#'   \item{\code{treatment}}{
#'     Treatment assignment for patients, 1 for treatment group 0 for control
#'     group.}
#'   \item{\code{outcome}}{
#'     Continuous outcome of the trial (Gaussian distributed).}
#'   \item{\code{complete}}{
#'     1 for complete outcome, 0 for loss to follow-up.}
#' }
#'
#' @keywords datasets
#' @examples
#' data(normaldata)

"normaldata"
