#' Returns lagged and differenced data where input and output are same length.
#'
#' @param x A numeric, time series, or xts variable.
#' @param lag A number that represents the maximum lag order.
#' @param differences Number of times data should be differenced.
#' @param ... Optional input arguments.
#'
#' @return Vector of data that is the sample length as input.
#' @export
#'
#' @examples
#' dif(rnorm(100), 1, 1)

dif <- function(x, lag = 1, differences = 1, ...) {

  # number of observations
  num <- length(x)

  # calculations
  tmp0 <- diff(x, lag = lag, differences = differences)
  tmp1 <- c(rep(NA, num - length(tmp0)), tmp0)

  return(tmp1)
}
