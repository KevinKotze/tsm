#' Plot leading and lagging correlations
#' 
#' @param x1 A vector of numeric values
#' @param x2 A vector of numeric values
#' @param nlag A numberic scalar for the lag-order of the autoregressive (cyclical) part
#' @return Values for the stochastic trend and the stationary cycle
#' @examples

leadlag  <- function(x1,x2,nlag){

  # empty output matrix
  x <- rep(NaN,nlag*2+1)
  p <- x
  # make index to get the correct lead and lag structure
  index <- seq(nlag,-(nlag),-1)
  
  for (i in 1:(nlag*2+1)) {
    if (index[i]==0) {
      xi <- cor(x1, x2)
    } else if (index[i]>0) {
      xi <- cor(x1[1:(length(x1)-index[i])], x2[(index[i]+1):length(x2)])
    } else if (index[i]<0) {
      xi <- cor(x1[(1-index[i]):length(x1)], x2[1:(length(x2)+index[i])])
    }
    x[i] <- xi
  }
  return(x)
}
