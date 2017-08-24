#' Newey-West HAC covariance estimator.
#' 
#' @param y A T*k vector of numeric values.
#' @param qn A numberic scalar for the truncation lag.
#' @return Newey-West HAC covariance estimator (as derived in Hayashi).
#' @examples
#' 

nw <- function(y,qn){

T <- length(y)
ybar <- rep(1,T) * ((sum(y))/T)
dy <- y-ybar
G0 <- t(dy) %*% dy/T

for (j in 1:qn){
   gamma <- t(dy[(j+1):T]) %*% dy[1:T-j]/(T-1)
   G0 <- G0+(gamma+t(gamma))*(1-abs(j/qn))
  }

  return(as.numeric(G0))

}

