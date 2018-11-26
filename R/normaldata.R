#' Gaussian dataset for analyzing adaptive Bayesian trials
#'
#' A dataset containing the results of 300 patients with continuous (normal) outcome,
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
#'   \item{id}{Patient ID in the trial}
#'   \item{treatment}{treatment assignment for patients, 1 for treatment group
#'   0 for control group}
#'   \item{outcome}{continuous outcome of the trial (gaussian distribution)}
#'   \item{complete}{1 for complete outcome, 0 for loss to follow up}
#' }
#'
#' @examples
#' data(normaldata)
"normaldata"
