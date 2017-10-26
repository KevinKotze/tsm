#' Estimation of a Gaussian GARCH-in-Mean with GARCH(1,1) model.
#'
#' @param rtn Time series variable.
#' @param type 1 for Variance-in-mean, 2 for volatility-in-mean, and 3 for log(variance)-in-mean.
#' @return GARCH-M results.
#' @examples

"garchM" <- function(rtn,type=1){

  if(is.matrix(rtn))rtn=c(rtn[,1])
  garchMdata <<- rtn 
  # obtain initial estimates
  m1=garch11FIT(garchMdata)
  est=as.numeric(m1$par); se.est=as.numeric(m1$separ);ht=m1$ht
  Mean=est[1]; cc=est[2]; ar=est[3]; ma=est[4]; S=1e-6
  seMean=se.est[1];secc=se.est[2];sear=se.est[3];sema=se.est[4]
  v1=ht
  if(type==2)v1=sqrt(v1)
  if(type==3)v1=log(v1)
  m2=lm(garchMdata~v1)
  Cnst=as.numeric(m2$coefficients[1])
  gam=as.numeric(m2$coefficients[2])
  params=c(mu=Cnst,gamma=gam, omega=cc, alpha=ar,beta=ma)
  lowBounds=c(mu=Mean-4*seMean,gamma=-10*abs(gam), omega=cc-2*secc, alpha=S, beta=ma-2*sema)
  uppBounds=c(mu=Mean+4*seMean,gamma=10*abs(gam), omega=cc+2*secc ,alpha=ar+2*sear,beta=1-S)
  Vtmp <<- c(type,v1[1])
  #
  fit=nlminb(start = params, objective= glkM, lower=lowBounds, upper=uppBounds,
             control=list(trace=3,rel.tol=1e-5))
  epsilon = 0.0001 * fit$par
  npar=length(params)
  Hessian = matrix(0, ncol = npar, nrow = npar)
  for (i in 1:npar) {
    for (j in 1:npar) {
      x1 = x2 = x3 = x4  = fit$par
      x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
      x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
      x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
      x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
      Hessian[i, j] = (glkM(x1)-glkM(x2)-glkM(x3)+glkM(x4))/
        (4*epsilon[i]*epsilon[j])
    }
  }
  cat("Maximized log-likehood: ",glkM(fit$par),"\n")
  # Step 6: Create and Print Summary Report:
  se.coef = sqrt(diag(solve(Hessian)))
  tval = fit$par/se.coef
  matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
  dimnames(matcoef) = list(names(tval), c(" Estimate",
                                          " Std. Error", " t value", "Pr(>|t|)"))
  cat("\nCoefficient(s):\n")
  printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
  
  m3=ResiVol(fit$par)
  
  garchM <- list(residuals=m3$residuals,sigma.t=m3$sigma.t)
}

glkM = function(pars){
  rtn <- garchMdata
  mu=pars[1]; gamma=pars[2]; omega=pars[3]; alpha=pars[4]; beta=pars[5]
  type=Vtmp[1]
  T=length(rtn)
  # use conditional variance
  if(type==1){
    ht=Vtmp[2]
    et=rtn[1]-mu-gamma*ht
    at=c(et)
    for (i in 2:T){
      sig2t=omega+alpha*at[i-1]^2+beta*ht[i-1]
      ept = rtn[i]-mu-gamma**sig2t
      at=c(at,ept)
      ht=c(ht,sig2t)
    }
  }
  # use volatility
  if(type==2){
    ht=Vtmp[2]^2
    et=rtn[1]-mu-gamma*Vtmp[2]
    at=c(et)
    for (i in 2:T){
      sig2t=omega+alpha*at[i-1]^2+beta*ht[i-1]
      ept=rtn[i]-mu-gamma*sqrt(sig2t)
      at=c(at,ept)
      ht=c(ht,sig2t)
    }
  }
  # use log(variance)
  if(type==3){
    ht=exp(Vtmp[2])
    et=rtn[1]-mu-gamma*Vtmp[2]
    at=c(et)
    for (i in 2:T){
      sig2t=omega+alpha*at[i-1]+beta*ht[i-1]
      ept=rtn[i]-mu-gamma*log(abs(sig2t))
      at=c(at,ept)
      ht=c(ht,sig2t)
    }
  }
  #
  hh=sqrt(abs(ht))
  glk=-sum(log(dnorm(x=at/hh)/hh))
  
  glk
}


ResiVol = function(pars){
  rtn <- garchMdata
  mu=pars[1]; gamma=pars[2]; omega=pars[3]; alpha=pars[4]; beta=pars[5]
  type=Vtmp[1]
  T=length(rtn)
  # use conditional variance
  if(type==1){
    ht=Vtmp[2]
    et=rtn[1]-mu-gamma*ht
    at=c(et)
    for (i in 2:T){
      sig2t=omega+alpha*at[i-1]^2+beta*ht[i-1]
      ept = rtn[i]-mu-gamma**sig2t
      at=c(at,ept)
      ht=c(ht,sig2t)
    }
  }
  # use volatility
  if(type==2){
    ht=Vtmp[2]^2
    et=rtn[1]-mu-gamma*Vtmp[2]
    at=c(et)
    for (i in 2:T){
      sig2t=omega+alpha*at[i-1]^2+beta*ht[i-1]
      ept=rtn[i]-mu-gamma*sqrt(sig2t)
      at=c(at,ept)
      ht=c(ht,sig2t)
    }
  }
  # use log(variance)
  if(type==3){
    ht=exp(Vtmp[2])
    et=rtn[1]-mu-gamma*Vtmp[2]
    at=c(et)
    for (i in 2:T){
      sig2t=omega+alpha*at[i-1]+beta*ht[i-1]
      ept=rtn[i]-mu-gamma*log(abs(sig2t))
      at=c(at,ept)
      ht=c(ht,sig2t)
    }
  }
  #
  
  ResiVol <- list(residuals=at,sigma.t=sqrt(ht))
}
#####################
garch11FIT = function(x){
  # Step 1: Initialize Time Series Globally:
  tx <<- x
  # Step 2: Initialize Model Parameters and Bounds:
  Mean = mean(tx); Var = var(tx); S = 1e-6
  params = c(mu = Mean, omega = 0.1*Var, alpha = 0.1, beta = 0.8)
  lowerBounds = c(mu = -10*abs(Mean), omega = S^2, alpha = S, beta = S)
  upperBounds = c(mu = 10*abs(Mean), omega = 100*Var, alpha = 1-S, beta = 1-S)
  # Step 3: Set Conditional Distribution Function:
  garch11Dist = function(z, hh) { dnorm(x = z/hh)/hh }
  # Step 4: Compose log-Likelihood Function:
  garch11LLH = function(parm) {
    mu = parm[1]; omega = parm[2]; alpha = parm[3]; beta = parm[4]
    z = (tx-mu); Mean = mean(z^2)
    # Use Filter Representation:
    e = omega + alpha * c(Mean, z[-length(tx)]^2)
    h = filter(e, beta, "r", init = Mean)
    hh = sqrt(abs(h))
    llh = -sum(log(garch11Dist(z, hh)))
    llh }
  # Step 5: Estimate Parameters and Compute Numerically Hessian:
  fit = nlminb(start = params, objective = garch11LLH,
               lower = lowerBounds, upper = upperBounds)
  #
  epsilon = 0.0001 * fit$par
  npar=length(params)
  Hessian = matrix(0, ncol = npar, nrow = npar)
  for (i in 1:npar) {
    for (j in 1:npar) {
      x1 = x2 = x3 = x4  = fit$par
      x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
      x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
      x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
      x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
      Hessian[i, j] = (garch11LLH(x1)-garch11LLH(x2)-garch11LLH(x3)+garch11LLH(x4))/
        (4*epsilon[i]*epsilon[j])
    }
  }
  # 
  se.coef = sqrt(diag(solve(Hessian)))
  est=fit$par
  # compute the sigma.t^2 series
  z=tx-est[1]; Mean=mean(z^2)
  e=est[2]+est[3]*c(Mean,z[-length(tx)]^2)
  h=filter(e,est[4],"r",init=Mean)
  
  garch11FIT <- list(par=est,separ=se.coef,ht=h)
}
