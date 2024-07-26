#' General-to-Specific application of Dickey-Fuller (1981) Test.
#'
#' @importFrom urca ur.df
#' @param series A vector of numeric values.
#' @return Summary of the results of the various tests.
#' @examples
#' gts_ur(rnorm(100))

gts_ur<-function(series)
{
  ur.trend <- urca::ur.df(series, type='trend', selectlags = c("AIC"))
  tstat.trend <- ur.trend@teststat
  cv.trend <- ur.trend@cval
  res.trend <- cbind(t(round(tstat.trend,2)),cv.trend)
  nam.trend <- rownames(res.trend)
  nam.trend[agrep("tau", nam.trend)] <- "pi"
  nam.trend[grep("phi2", nam.trend)] <- "varphi2"
  nam.trend[grep("phi3", nam.trend)] <- "varphi3"
  rownames(res.trend) <- nam.trend

  ur.drift <- urca::ur.df(series, type='drift', selectlags = c("AIC"))
  tstat.drift <- ur.drift@teststat
  cv.drift <- ur.drift@cval
  res.drift <- cbind(t(round(tstat.drift,2)),cv.drift)
  nam.drift <- rownames(res.drift)
  nam.drift[agrep("tau", nam.drift)] <- "pi"
  nam.drift[grep("phi1", nam.drift)] <- "varphi1"
  rownames(res.drift) <- nam.drift

  ur.none <- urca::ur.df(series, type='none', selectlags = c("AIC"))
  tstat.none <- ur.none@teststat
  cv.none <- ur.none@cval
  res.none <- cbind(t(round(tstat.none,2)),cv.none)
  nam.none <- rownames(res.none)
  nam.none[agrep("tau", nam.none)] <- "pi"
  rownames(res.none) <- nam.none

  # print summary
  cat(" ", "\n")
  cat("#################", "\n")
  cat("## ADF summary ##", "\n")
  cat("#################", "\n")
  cat(" ", "\n")

  if (res.trend[1,1] <= res.trend[1,3]) {
    cat("Able to reject null of unit root at 5% - with constant & trend", "\n")
  }else if (res.trend[1,1] > res.trend[1,3] && res.trend[2,1] >= res.trend[2,3]) {
    cat("Unable to reject null of unit root at 5% - with constant & trend", "\n")

  }else if (res.drift[1,1] <= res.drift[1,3]) {
    cat("Able to reject null of unit root at 5% - with constant", "\n")
  }else if (res.drift[1,1] > res.drift[1,3] && res.drift[2,1] >= res.drift[2,3]) {
    cat("Unable to reject null of unit root at 5% - with constant", "\n")

  }else if (res.none[1,1] <= res.none[1,3]) {
    cat("Able to reject null of unit root at 5% - no deterministic", "\n")
  }else {cat("Cannot reject null of unit root at 5% - no deterministic", "\n")}

  # print results
  cat(" ", "\n")
  cat("## ADF with constrant and time trend ##", "\n")
  print(res.trend)
  cat(" ", "\n")
  if (res.trend[1,1] > res.trend[1,3]) {
    cat("Cannot reject null of unit root at 5%", "\n")
  }else cat("Able to reject null of unit root at 5%", "\n")
    if (res.trend[2,1] < res.trend[2,3]) {
    cat("Cannot reject null of no constant and no trend at 5%", "\n")
  }else cat("Able to reject null of no constant and no trend at 5%", "\n")
  if (res.trend[3,1] < res.trend[3,3]) {
    cat("Cannot reject null of no trend at 5%", "\n")
  }else cat("Able to reject null of no trend at 5%", "\n")
  cat(" ", "\n")

  cat("## ADF with constrant ##", "\n")
  print(res.drift)
  cat(" ", "\n")
  if (res.drift[1,1] > res.drift[1,3]) {
    cat("Cannot reject null of unit root at 5%", "\n")
  }else cat("Able to reject null of unit root at 5%", "\n")
  if (res.drift[2,1] < res.drift[2,3]) {
    cat("Cannot reject null of no constant at 5%", "\n")
  }else cat("Able to reject null of no constant at 5%", "\n")
  cat(" ", "\n")

  cat("## ADF with no deterministic ##", "\n")
  print(res.none)
  cat(" ", "\n")
  if (res.none[1,1] > res.none[1,3]) {
    cat("Cannot reject null of unit root at 5%", "\n")
  }else cat("Able to reject null of unit root at 5%", "\n")
  cat(" ", "\n")

}
