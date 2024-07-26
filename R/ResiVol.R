#' Auxillary function for reporting of GARCH-M results.
#'
#' @param pars initial results from model.
#' @return GARCH results.
#' @examples

ResiVol <- function(pars) {
  rtn <- garchMdata
  mu <- pars[1]
  gamma <- pars[2]
  omega <- pars[3]
  alpha <- pars[4]
  beta <- pars[5]
  type <- Vtmp[1]
  T <- length(rtn)
  # use conditional variance
  if (type == 1) {
    ht <- Vtmp[2]
    et <- rtn[1] - mu - gamma * ht
    at <- c(et)
    for (i in 2:T) {
      sig2t <- omega + alpha * at[i - 1] ^ 2 + beta * ht[i - 1]
      ept <- rtn[i] - mu - gamma ** sig2t
      at <- c(at, ept)
      ht <- c(ht, sig2t)
    }
  }
  # use volatility
  if (type == 2) {
    ht <- Vtmp[2] ^ 2
    et <- rtn[1] - mu - gamma * Vtmp[2]
    at <- c(et)
    for (i in 2:T) {
      sig2t <- omega + alpha * at[i - 1] ^ 2 + beta * ht[i - 1]
      ept <- rtn[i] - mu - gamma * sqrt(sig2t)
      at <- c(at, ept)
      ht <- c(ht, sig2t)
    }
  }
  # use log(variance)
  if (type == 3) {
    ht = exp(Vtmp[2])
    et <- rtn[1] - mu - gamma * Vtmp[2]
    at <- c(et)
    for (i in 2:T) {
      sig2t <- omega + alpha * at[i - 1] + beta * ht[i - 1]
      ept <- rtn[i] - mu - gamma * log(abs(sig2t))
      at <- c(at, ept)
      ht <- c(ht, sig2t)
    }
  }
  #

  ResiVol <- list(residuals = at, sigma.t = sqrt(ht))
}

