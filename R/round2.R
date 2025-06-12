#' Round numbers following normal convention.
#'
#' @param x A numeric scalar.
#' @param n The order for the rounding.
#'
#' @return The result that would be equivalent to Matlab.
#' @export
#'
#' @examples
#' round2(0.55, 1)

round2 <- function(x, n) {
  # where 0.5 rounds up to 1 and 0.4 to 0
  posneg <- sign(x)
  z <- abs(x) * 10^n
  z <- z + 0.5
  z <- trunc(z)
  z <- z / 10^n
  z * posneg
}
