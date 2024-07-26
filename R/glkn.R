#' Auxillary function for reporting of NGARCH results.
#'
#' @param par initial results from model.
#' @return GARCH results.
#' @examples

glkn <- function(par) {
  rtn = read.table("tmp.txt")[, 1]
  glkn = 0
  ht = var(rtn)
  T = length(rtn)
  if (T > 40)
    ht = var(rtn[1:40])
  at = rtn[1] - par[1]
  for (i in 2:T) {
    ept = rtn[i] - par[1]
    at = c(at, ept)
    sig2t = par[2] + par[3] * ht[i - 1] + par[4] * ht[i - 1] * (at[i - 1] /
                                                                  sqrt(ht[i - 1]) - par[5]) ^ 2
    ht = c(ht, sig2t)
    glkn = glkn + 0.5 * (log(sig2t) + ept ^ 2 / sig2t)
  }
  glkn
}
