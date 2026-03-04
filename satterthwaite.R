library(Matrix)
satterthwaite.glmerMod = function(rslt,lvec){
  #
  if (class(rslt)!="glmerMod") stop("rslt must be of class glmerMod")
  p = length(fixef(rslt))
  if (length(lvec)!=p) stop("length of lvec must be equal to number of fixed effects")
  lvec = matrix(lvec,p,1) 
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
  A = MASS::ginv(rslt@optinfo$derivs$Hessian/2)[1:q,1:q] # theta parameterization
  # get gradient of t(lvec)%*%C%*%lvec wrt theta
  g = grad(func,theta,obj=rslt,lvec=lvec)
  g = matrix(g,q,1)
  #
  df = (2*(t(lvec)%*%C%*%lvec)^2) / (t(g)%*%A%*%g)
  as.numeric(df)
}

# First attempt - based on Breslow and Clayton equations
func = function(newtheta,obj,lvec){
  #
  # this is not quite right because u also depends on theta
  expit <- function(x){exp(x)/(1 + exp(x))}
  gprime <- function(x){1/(x*(1-x))}
  vv <- function(x){x*(1-x)}
  #
  u = getME(obj,"u")
  X = getME(obj,"X")
  Z = getME(obj,"Z")
  beta = fixef(obj)
  lambda = getME(obj,"Lambda")
  Lind = getME(obj,"Lind")
  lambda@x = newtheta[Lind]
  #
  phi <- 1
  a <- 1
  b = lambda%*%u
  mu = as.vector(expit(X%*%beta + Z%*%b))
  w = 1/(phi*a*vv(mu)*gprime(mu)^2)
  V = diag(1/w) + Z%*%lambda%*%t(lambda)%*%t(Z)
  as.numeric(t(lvec)%*%solve(t(X)%*%solve(V)%*%X)%*%lvec)
}

# Second attempt 
func = function(newtheta,obj,lvec){
  #
  expit <- function(x){exp(x)/(1 + exp(x))}
  iter <- function(mu,a,phi,y,X,Z,beta,lambda) {
    eps = .0000001
    gprime <- function(x){1/(x*(1-x))}
    vv <- function(x){x*(1-x)}
    repeat{
      W = diag(1/(phi*a*vv(mu)*gprime(mu)^2))
      V = solve(W) + Z%*%lambda%*%t(lambda)%*%t(Z)
      b = lambda%*%t(lambda)%*%t(Z)%*%solve(V)%*%(y- X%*%beta)
      newmu = expit(X%*%beta + Z%*%b)
      if (sum(abs((newmu - mu)/mu)/length(mu) < eps) break()
          mu = newmu
    }
    V
  }
  #
  y = getME(obj,"y")
  mu = getME(obj,"mu")
  X = getME(obj,"X")
  Z = getME(obj,"Z")
  beta = fixef(obj)
  lambda = getME(obj,"Lambda")
  Lind = getME(obj,"Lind")
  lambda@x = newtheta[Lind]
  #
  phi <- 1
  a <- 1
  V = iter(mu,a,phi,y,X,Z,beta,lambda)
  as.numeric(t(lvec)%*%solve(t(X)%*%solve(V)%*%X)%*%lvec)
}

# Third attempt
func = function(newtheta,obj,lvec){
  #
  # this is not quite right because u also depends on theta
  expit <- function(x){exp(x)/(1 + exp(x))}
  gprime <- function(x){1/(x*(1-x))}
  vv <- function(x){x*(1-x)}
  #
  u = getME(obj,"u")
  X = getME(obj,"X")
  Z = getME(obj,"Z")
  beta = fixef(obj)
  lambda = getME(obj,"Lambda")
  Lind = getME(obj,"Lind")
  lambda@x = newtheta[Lind]
  #
  phi <- 1
  a <- 1
  b = lambda%*%u
  mu = as.vector(expit(X%*%beta + Z%*%b))
  w = 1/(phi*a*vv(mu)*gprime(mu)^2)
  V = diag(1/w) + Z%*%lambda%*%t(lambda)%*%t(Z)
  as.numeric(t(lvec)%*%solve(t(X)%*%solve(V)%*%X)%*%lvec)
}

# Fourth attempt
# I wish this worked, but it doesn't ... calling func with theta = getME(obj,"theta") 
# does not return as.numeric(t(lvec)%*%vcov(obj)%*%lvec)
ctrl <- glmerControl(optCtrl = list(maxfun = 1, checkConv = FALSE),
                     nAGQ0initStep = FALSE, optimizer = "bobyqa")   
func = function(newtheta,obj,lvec,ctrl){
  beta <- fixef(obj)
  newobj <- update(obj, control = ctrl,
                   start=list(theta=newtheta,fixef=beta))
  as.numeric(t(lvec)%*%vcov(newobj)%*%lvec)
}

##########################
# To get weights can do this
# W = weights(objb,type="working")

# For testing purposes
datac = swSim(swDsn(c(2,2,2,2)),family="gaussian",n=10,mu0=0,mu1=1,time.effect=0,
              sigma=1.2,tau=.2)
objc = lmer(response.var ~ as.factor(time.var) + tx.var + (1 | cluster.var),data=datac)

datab = swSim(swDsn(c(2,2,2,2)),family="binomial",n=10,mu0=0,mu1= -0.5,time.effect=0,
              tau=.4)
objb = glmer(response.var ~ as.factor(time.var) + tx.var + (1 | cluster.var),family=binomial,data=datab)
# Note: function below does not work correctly with glmer object but does with lmer

func = function(newtheta,obj,lvec){
  #
  # this is not quite right because u also depends on theta?
  if (family(obj)$link=="identity"){
  #identity
   linkinv <- function(x) {x}
   gprime <- function(x) {rep(1,length(x))}
   vv <- function(x) {rep(1,length(x))}
   phi <- 1
   a <- 1
  }
  if (family(obj)$link=="logit") {
  #logit
   linkinv <- function(x){exp(x)/(1 + exp(x))}
   gprime <- function(x){1/(x*(1-x))}
   vv <- function(x){x*(1-x)}
   phi <- 1
   a <- 1
  }
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
  V = sigma^2 * solve(RXtRX)
  as.numeric(t(lvec)%*%V%*%lvec)
}

satterthwaite.lmerMod = function(rslt,lvec){
  #
  if (class(rslt)!="lmerMod") stop("rslt must be of class lmerMod")
  p = length(fixef(rslt))
  if (length(lvec)!=p) stop("length of lvec must be equal to number of fixed effects")
  lvec = matrix(lvec,p,1) 
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
  A = MASS::ginv(rslt@optinfo$derivs$Hessian/2)[1:q,1:q] # theta parameterization
  # get gradient of t(lvec)%*%C%*%lvec wrt theta
  g = grad(func,theta,obj=rslt,lvec=lvec)
  g = matrix(g,q,1)
  #
  df = (2*(t(lvec)%*%C%*%lvec)^2) / (t(g)%*%A%*%g)
  as.numeric(df)
}