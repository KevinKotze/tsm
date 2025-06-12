#' Perform Lagrange Multiplier Test for ARCH effect of a time series.
#'
#' @param data Time series variable.
#' @param m Selected AR order.
#' @return ARCH Lagrange Multiplier Statistics.
#' @export
#'
#' @examples
#' arch_test(rnorm(100))

arch_test <- function(data, m = 10) {

  # demean data and take the square
  y <- (data - mean(data))^2

  # set up lags for RHS variable
  T <- length(data)
  atsq <- y[(m + 1):T]
  x <- matrix(0, (T - m), m)
  for (i in 1:m) {
    x[, i] <- y[(m + 1 - i):(T - i)]
  }

  # estimate least squares regression
  result <- stats::lm(atsq ~ x)

  # return the result
  return(result)
}

