library(Matrix)
library(numDeriv)
##########################
# Compute Satterthwaite degrees of freedom for an lmer object
satterthwaite.lmerMod = function(rslt,lvec){
  ####  checks
  if (class(rslt)!="lmerMod") stop("rslt must be of class lmerMod")
  p = length(fixef(rslt))
  if (!is.matrix(lvec)) {
    if (length(lvec) != p) {
      stop("length of lvec must be equal to number of fixed effects")
    } else {
      lvec = matrix(lvec,1,p) 
    }
  } else {
    if (ncol(lvec)!=p) stop("number of columns of lvec must be equal to number of fixed effects")
  }
  r = nrow(lvec)
  if (r != rankMatrix(lvec)) stop("contrasts in lvec are not independent")
  ####
  # get var-cov matrix of fixed effects
  C = vcov(rslt)
  # get var-cov matrix of random effects (theta parameterization)
  theta = getME(rslt,"theta")
  q = length(theta)
  if (any(theta==0)) {
    #    if (all(theta==0)){glm()$df.residual}
    warning("Unable to compute df with zero variance component, returning NA")
    return(NA)
  } 
  if (is.null(rslt@optinfo$derivs)){
    warning("Variance of theta parameters not available from lmer, returning NA")
    return(NA)
  }
  A = MASS::ginv(rslt@optinfo$derivs$Hessian/2)[1:q,1:q] # theta parameterization
  if (r==1){
    g = grad(func_lmer,theta,obj=rslt,lvec=lvec,side=rep(1,length(theta)),method="Richardson",method.args=list(r=6))
    g = matrix(g,q,1)
    df = (2*(lvec%*%C%*%t(lvec))^2) / (t(g)%*%A%*%g)
  } else{
    temp=eigen(lvec%*%C%*%t(lvec),only.values=TRUE)
    D = temp$values
    g = t(jacobian(func_lmer,theta,obj=rslt,lvec=lvec,side=rep(1,length(theta)),method="Richardson",method.args=list(r=6)))
    # g is q x r   
    E = 0
    for (j in 1:r){
      v = (2*D[j]^2) / (t(g[,j,drop=FALSE])%*%A%*%g[,j,drop=FALSE])
      if (v>2) E = E + v/(v-2)
    }
    df = ifelse(E>r,2*E/(E-r),0)
  }
  as.numeric(df)
}


func_lmer = function(newtheta,obj,lvec){
###
#identity link
   linkinv <- function(x) {x}
   gprime <- function(x) {rep(1,length(x))}
   vv <- function(x) {rep(1,length(x))}
   phi <- 1
   a <- 1
###
  r = nrow(lvec)
  #
  u = getME(obj,"u")
  X = getME(obj,"X")
  Z = getME(obj,"Z")
  L = getME(obj,"L")
  beta = fixef(obj)
  lambda = getME(obj,"Lambda")
  Lind = getME(obj,"Lind")
  lambda@x = newtheta[Lind]
  Lambdat = t(lambda)
#  Rx = getME(obj,"RX")
  sigma = getME(obj,"sigma")
  b = lambda%*%u  
  mu = as.vector(linkinv(X%*%beta + Z%*%b))
  W = diag(1/(phi*a*vv(mu)*gprime(mu)^2)) 
#
  ZtW = t(Z)%*%sqrt(W)
  ZtWX = t(Z)%*%W%*%X
  XtWX = t(X)%*%W%*%X
  L <- update(L, Lambdat %*% ZtW, mult = 1)
  RZX <- solve(L, solve(L, Lambdat %*% ZtWX, system = "P"),system = "L")
  RXtRX <- XtWX - t(RZX)%*%RZX
  C = sigma^2 * solve(RXtRX)
  if (r==1) {
    g = as.numeric(lvec%*%C%*%t(lvec))
  } else {
    g = rep(NA,r)
    temp=eigen(lvec%*%C%*%t(lvec))
    U = temp$vectors
    b = U%*%lvec
    for (j in 1:r){
      g[j] = b[j,,drop=FALSE]%*%C%*%t(b[j,,drop=FALSE])
    }
  }
  return(g)
}