nw <- function(y,qn){
#input: y is a T*k vector and qn is the truncation lag
#output: the newey west HAC covariance estimator
#Formulas are from Hayashi
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

