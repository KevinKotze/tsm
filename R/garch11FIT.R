#' Auxillary function for fitting GARCH model.
#'
#' @param x time series vector.
#' @return GARCH parameter estimates.
#' @examples


garch11FIT = function(x) {
  # Step 1: Initialize Time Series Globally:
  tx <<- x
  # Step 2: Initialize Model Parameters and Bounds:
  Mean <- mean(tx)
  Var <- var(tx)
  S <- 1e-6
  params <- c(
    mu = Mean,
    omega = 0.1 * Var,
    alpha = 0.1,
    beta = 0.8
  )
  lowerBounds <- c(
    mu = -10 * abs(Mean),
    omega = S ^ 2,
    alpha = S,
    beta = S
  )
  upperBounds <- c(
    mu = 10 * abs(Mean),
    omega = 100 * Var,
    alpha = 1 - S,
    beta = 1 - S
  )
  # Step 3: Set Conditional Distribution Function:
  garch11Dist <- function(z, hh) {
    dnorm(x = z / hh) / hh
  }
  # Step 4: Compose log-Likelihood Function:
  garch11LLH <- function(parm) {
    mu = parm[1]
    omega = parm[2]
    alpha = parm[3]
    beta = parm[4]
    z = (tx - mu)
    Mean = mean(z ^ 2)
    # Use Filter Representation:
    e = omega + alpha * c(Mean, z[-length(tx)] ^ 2)
    h = stats::filter(e, beta, "r", init = Mean)
    hh = sqrt(abs(h))
    llh = -sum(log(garch11Dist(z, hh)))
    llh
  }
  # Step 5: Estimate Parameters and Compute Numerically Hessian:
  fit <- nlminb(
    start = params,
    objective = garch11LLH,
    lower = lowerBounds,
    upper = upperBounds
  )
  #
  epsilon <- 0.0001 * fit$par
  npar <- length(params)
  Hessian <- matrix(0, ncol = npar, nrow = npar)
  for (i in 1:npar) {
    for (j in 1:npar) {
      x1 = x2 = x3 = x4  = fit$par
      x1[i] = x1[i] + epsilon[i]
      x1[j] = x1[j] + epsilon[j]
      x2[i] = x2[i] + epsilon[i]
      x2[j] = x2[j] - epsilon[j]
      x3[i] = x3[i] - epsilon[i]
      x3[j] = x3[j] + epsilon[j]
      x4[i] = x4[i] - epsilon[i]
      x4[j] = x4[j] - epsilon[j]
      Hessian[i, j] = (garch11LLH(x1) - garch11LLH(x2) - garch11LLH(x3) +
                         garch11LLH(x4)) /
        (4 * epsilon[i] * epsilon[j])
    }
  }
  #
  se.coef <- sqrt(diag(solve(Hessian)))
  est <- fit$par
  # compute the sigma.t^2 series
  z <- tx - est[1]
  Mean <- mean(z ^ 2)
  e <- est[2] + est[3] * c(Mean, z[-length(tx)] ^ 2)
  h <- stats::filter(e, est[4], "r", init = Mean)

  garch11FIT <- list(par = est, separ = se.coef, ht = h)
}
