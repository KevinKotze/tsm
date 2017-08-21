CWtest <- function(e.m1,e.m2,yf.m1,yf.m2){
  #input: e.m1 - errors model 1, e.m2 - errors model2, yf.m1 - forecast model1, yf.m2 - forecast model1
  #output: Clark-West (2007) Approximately normal tests for equal predictive accuracy in nested models
  
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
  
  P <- length(e.m1)
  froll.adj <- e.m1^2-(e.m2^2-(yf.m1-yf.m2)^2)
  varfroll.adj <- nw(froll.adj,1)
  CW <- sqrt(P)*(mean(froll.adj))/sqrt(varfroll.adj)
  pv <- 1-pnorm(CW,0,1)
  results=list(test=0,pvalue=0)
  results$test <- CW
  results$pvalue <- pv
  return(results)
}