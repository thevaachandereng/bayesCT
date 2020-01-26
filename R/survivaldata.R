#' Time-to-event dataset for analyzing adaptive Bayesian trials
#'
#' A dataset containing the results of 100 patients with time-to-event outcome,
#' the dataset is filled with treatment assignment and status (0 = censored,
#' 1 = not censored).
#'
#' @docType data
#'
#' @keywords dataset
#'
#' @usage data(survivaldata)
#'
#' @format A data frame with 100 rows and 4 variables:
#' \describe{
#'   \item{id}{Patient ID in the trial}
#'   \item{treatment}{treatment assignment for patients, 1 for treatment group
#'   0 for control group}
#'   \item{time}{the follow up time for patients}
#'   \item{event}{The status indicator, normally 0=alive, 1=dead or
#'   0 = no event, 1 = event occurred}
#' }
#'
#' @examples
#' data(survivaldata)
"survivaldata"
