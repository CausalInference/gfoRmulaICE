#' Example Dataset for a Survival Outcome with Both Censoring and Competing Event
#'
#' A dataset with 26581 observations on 10000 individuals and 4 time points. The dataset is in long format with each row representing the record of one individual at one time point.
#'
#' @format A data frame with 26581 rows and 9 variables:
#' \describe{
#'   \item{t0}{Time index.}
#'   \item{id}{Unique identifier for each individual.}
#'   \item{L1}{Binary covariate.}
#'   \item{L2}{Continuous covariate.}
#'   \item{A1}{Categorical treatment variable with levels 1, 2, and 3.}
#'   \item{A2}{Binary treatment variable.}
#'   \item{C}{Censoring event indicator.}
#'   \item{D}{Competing event indicator.}
#'   \item{Y}{Outcome indicator.}
#' }
"compData"