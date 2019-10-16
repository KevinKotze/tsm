#' Beveridge-Nelson decomposition
#' 
#' @param data A vector of first-order integrated numeric values
#' @param nlag A numberic scalar for the lag-order of the autoregressive (cyclical) part
#' @return Values for the stochastic trend and the stationary cycle
#' @examples
#' 

bnd  <- function(data,nlag){

y=matrix(data,ncol=1)
#nlag=8
  
yd=diff(y,lag=1)
yl=matrix(rep(0,length(y)*nlag),ncol=nlag)
yl[,1] = c(0,yd)
for (i in 2:nlag){
  yl[,i] = c(0,yl[1:(length(y)-1),i-1])
}
x=yl[(nlag+1):(length(y)-1),1:nlag]
yy=matrix(yd[(1+nlag):length(yd),],ncol=1)

# OLS
beta <- lm(yy ~ x)$coefficients

# Companion form of matrix
eye=diag(1,nlag)
coef.tmp=matrix(beta[2:length(beta)],nrow=1)
betac.tmp=rbind(coef.tmp,eye)
betac=betac.tmp[1:nlag,]


c1=betac %*% solve(diag(nlag)-betac)
ydd=c(rep(0,1+nlag),yd)-beta[1]
ydd.len=length(ydd)

# Construct matrix of historical lags
yD = matrix(rep(0,nlag*(ydd.len-nlag)),nrow=nlag)

for (i in 1:nlag){
  yD[i,]=matrix(ydd[(1+i):(ydd.len-(nlag)+i)],nrow=1)
}

yD.tmp=apply(yD,2,rev)
yD=yD.tmp[nrow(yD.tmp):1,]

# Selection vector
sel=rep(0,nlag)
sel[1]=1

# Compute trend and cycle
ytr=y+t(sel%*%c1%*%yD)
yc=t(sel%*%c1%*%yD)
out=cbind(ytr,yc)
colnames(out) = c("trend","cycle")
out=out
}
