#' Time-to-event dataset for analyzing adaptive Bayesian trials
#'
#' @description A dataset containing the results of 100 patients with time-to-event outcome,
#' the dataset is filled with treatment assignment and status (0 = censored,
#' 1 = not censored).
#'
#' @docType data
#'
#' @usage data(survivaldata)
#'
#' @format A data frame with 100 rows and 4 variables:
#'
#' \describe{
#'   \item{\code{id}}{Patient ID in the trial.}
#'   \item{\code{treatment}}{Treatment assignment for patients, 1 for
#'   treatment group 0 for control group.}
#'   \item{\code{time}}{The follow up time for patients.}
#'   \item{\code{event}}{The status indicator, normally 0=alive, 1=dead or
#'   0 = no event, 1 = event occurred.}
#' }
#'
#' @keywords dataset
#' @examples
#' data(survivaldata)

"survivaldata"
