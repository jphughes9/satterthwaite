library(lme4)
library(Matrix)
library(numDeriv)
##########################
# Compute Satterthwaite degrees of freedom for an lmer object
# using code from Satterthwaite’s Method for Degrees of Freedom in Linear Mixed Models
# by Rune Haubo B Christensen
satterthwaite.lmerMod = function(model,L){
  #### Necessary functions
  devfun_varpar <- function(varpar, devfun, reml) {
    # Computes deviance as a function of 'varpar=c(theta, sigma)'
    # devfun: deviance function as a function of theta only.
    # reml: TRUE if REML; FALSE if ML
    nvarpar <- length(varpar)
    sigma2 <- varpar[nvarpar]^2
    theta <- varpar[-nvarpar]
    df_envir <- environment(devfun)
    devfun(theta) # Evaluate deviance function at varpar
    n <- nrow(df_envir$pp$V)
    # Compute deviance for ML:
    dev <- df_envir$pp$ldL2() + (df_envir$resp$wrss() + df_envir$pp$sqrL(1))/sigma2 +
      n * log(2 * pi * sigma2)
    if(!reml) return(dev)
    # Adjust of REML is used:
    RX <- df_envir$pp$RX() # X'V^{-1}X ~ crossprod(RX^{-1}) = cov(beta)^{-1} / sigma^2
    dev + 2*c(determinant(RX)$modulus) - ncol(RX) * log(2 * pi * sigma2)
  }
  get_covbeta <- function(varpar, devfun) {
    # Compute cov(beta) as a function of varpar
    # varpar: c(theta, sigma)
    # devfun: deviance function, ie. update(model, devFunOnly=TRUE) - a function of theta
    # return: cov(beta) given specified varpar
    #
    nvarpar <- length(varpar)
    sigma <- varpar[nvarpar] # residual std.dev.
    theta <- varpar[-nvarpar] # ranef var-par
    devfun(theta) # evaluate REML or ML deviance 'criterion'
    df_envir <- environment(devfun) # extract model environment
    covbeta <- sigma^2 * tcrossprod(df_envir$pp$RXi()) # vcov(beta)
#    covbeta <- sigma^2 * df_envir$pp$unsc() # gievs the same result
    return(covbeta)
  }
  check_boundary = function(theta){
    chk = !(abs(theta) < 1e-3) # Choose 0.001 because derivative calculations use +/- 0.0001
    if (all(!chk)) stop("All variance components are zero")
    return(chk)
  }
  ####  Checks
  if (class(model)!="lmerMod") stop("model must be of class lmerMod")
  if (!is.null(model@optinfo$conv$lme4$code)){
    if (grep("negative eigenvalues",rslt@optinfo$conv$lme4$messages)>0) stop("Model object has Hessian with negative eigenvalues; cannot compute df")
  }
  parlist <- list(beta=fixef(model),
                  theta=getME(model, "theta"),
                  sigma=sigma(model))
  p = length(parlist$beta)
  if (!is.matrix(L)) {
    if (length(L) != p) {
      stop("length of L must be equal to number of fixed effects")
    } else {
      L = matrix(L,1,p) 
    }
  } else {
    if (ncol(L)!=p) stop("number of columns of lvec must be equal to number of fixed effects")
  }
  q = nrow(L)
  if (q != rankMatrix(L)) stop("contrasts in L are not independent")
  if (any(!check_boundary(parlist$theta))) stop("One or more variance components are on the boundary on the theta scale; try fitting a simpler model")
  #### Set things up
  devfun <- update(model, devFunOnly=TRUE)
  is_reml <- getME(model, "is_REML")
  varpar_opt <- unname(c(parlist$theta, parlist$sigma))
  cov_beta <- as.matrix(vcov(model))
  #### Get Var-Cov matrix of variance parameters
  h <- hessian(func=devfun_varpar, x=varpar_opt, devfun=devfun, reml=is_reml)
  eig_h <- eigen(h, symmetric=TRUE)
  stopifnot(all(eig_h$values > 0))
  h_inv <- with(eig_h, vectors %*% diag(1/values) %*% t(vectors))
  cov_varpar <- 2 * h_inv
  #### Compute gradient of V(beta)
  Jac <- jacobian(func=get_covbeta, x=varpar_opt, devfun=devfun)
  Jac_list <- lapply(1:ncol(Jac), function(i) array(Jac[, i], dim=rep(length(parlist$beta), 2)))
  #### Compute df
  if (q==1){
  # 1 df case
  estimate <- sum(L * parlist$beta)
  var_Lbeta <- drop(L %*% cov_beta %*% t(L))
  se.estimate <- sqrt(var_Lbeta)
  grad_var_Lbeta <- vapply(Jac_list, function(x) {L %*% x %*% t(L)}, numeric(1L))
  satt_denom <- sum(grad_var_Lbeta * (cov_varpar %*% grad_var_Lbeta)) # g'Ag
  ddf <- drop(2 * var_Lbeta^2 / satt_denom) # denominator DF
  tstat <- estimate/se.estimate
  pvalue <- 2 * pt(abs(tstat), df = ddf, lower.tail = FALSE)
  data.frame(estimate, se=se.estimate, tstat, ddf, pvalue)
  } else {
  # multi-df case
  var_Lbeta <- L %*% cov_beta %*% t(L)
  eig_VLbeta <- eigen(var_Lbeta)
  P <- eig_VLbeta$vectors
  d <- eig_VLbeta$values
  PtL <- crossprod(P, L)
  t2 <- drop(PtL %*% parlist$beta)^2 / d
  Fvalue <- sum(t2) / q
  grad_PLcov <- lapply(1:q, function(m) {
    vapply(Jac_list, function(J) sum(PtL[m, ] * J %*% PtL[m, ]), numeric(1L))
  })
  nu_m <- vapply(1:q, function(m) {
    denom <- sum(grad_PLcov[[m]] * (cov_varpar %*% grad_PLcov[[m]])) # g'Ag
    2*(d[m])^2 / denom # 2d_m^2 / g'Ag
  }, numeric(1L))
  EQ <- sum(as.integer(nu_m>2)*nu_m / (nu_m - 2))
  ddf <- max(2, 2 * EQ / (EQ - q))
  pvalue <- pf(q=Fvalue, df1=q, df2=ddf, lower.tail=FALSE)
  data.frame('F value'=Fvalue, ndf=q, ddf=ddf, pvalue=pvalue,check.names = FALSE)
  }
  }

##########################
# Compute Satterthwaite degrees of freedom for an glmer object
# adapted from Satterthwaite’s Method for Degrees of Freedom in Linear Mixed Models
# by Rune Haubo B Christensen
satterthwaite.glmerMod = function(model,L,simplify=FALSE){
  #### Necessary functions
  get_covbeta <- function(theta, fixed, order, devfun) {
    # Compute cov(beta) as a function of theta
    # varpar: c(theta)
    # devfun: deviance function, ie. update(model, devFunOnly=TRUE) - a function of theta
    # return: cov(beta) given specified varpar
    #
    pars = c(theta,fixed)[order]
    devfun(pars) # evaluate REML or ML deviance 'criterion'
    df_envir <- environment(devfun) # extract model environment
    covbeta <- df_envir$pp$unsc()
    return(covbeta)
  }
  check_boundary = function(theta){
    chk = !(abs(theta) < 1e-3) # Choose 0.001 because derivative calculations use +/- 0.0001
    if (all(!chk)) stop("All variance components are zero")
    return(chk)
  }
  ####  Checks
  if (class(model)!="glmerMod") stop("model must be of class glmerMod")
  if (!is.null(model@optinfo$conv$lme4$code)){
    if (grep("negative eigenvalues",rslt@optinfo$conv$lme4$messages)>0) stop("Model object has Hessian with negative eigenvalues; cannot compute df")
  }
  parlist <- getME(model,c("theta","fixef"))
  p = length(parlist$fixef)
  r = length(parlist$theta)
  if (!is.matrix(L)) {
    if (length(L) != p) {
      stop("length of L must be equal to number of fixed effects")
    } else {
      L = matrix(L,1,p) 
    }
  } else {
    if (ncol(L)!=p) stop("number of columns of lvec must be equal to number of fixed effects")
  }
  q = nrow(L)
  if (q != rankMatrix(L)) stop("contrasts in L are not independent")
  #### Set things up
  devfun <- update(model, devFunOnly=TRUE)
  pars = unlist(parlist)
  #### Get Var-Cov matrix of variance parameters and fixed effects
  if (!is.null(model@optinfo$derivs$Hessian)){
    h = model@optinfo$derivs$Hessian
    dpars = pars[1:r]
    fpars = pars[-(1:r)]
    order = c(1:(r+p))
    cov_varpar = MASS::ginv(h/2)[1:r,1:r]
    cov_beta <- as.matrix(vcov(model))  
  } else {
    if (simplify){
     h = numDeriv::hessian(func=devfun, x=pars, method="Richardson",method.args=list(r=6))
     keep = check_boundary(parlist$theta)
     if (any(!keep)) warning("One or more variance components are zero on the theta scale; df calculation corresponds to a model with the zero variance components removed but the same beta - make sure this is what you want!")
     dpars = pars[1:r][keep]
     fpars = c(pars[1:r][!keep],pars[-(1:r)])
     order = c((1:r)[keep],(1:r)[!keep],(r+1):(r+p))
     r = sum(keep)
     keep = c(keep,rep(TRUE,p))
     cov_varpar = MASS::ginv(h[keep,keep]/2)[1:r,1:r]
     cov_beta = MASS::ginv(h[keep,keep]/2)[((r+1):(r+p)),((r+1):(r+p))]
    } else {
      stop("No hessian available; if this is due to a zero variance component, try fitting a simpler model")
    }
  }
  #### Compute gradient of V(beta) wrt non-boundary theta
  Jac <- jacobian(func=get_covbeta, x=dpars, fixed=fpars, order=order, devfun=devfun, method="Richardson")
  Jac_list <- lapply(1:ncol(Jac), function(i) array(Jac[, i], dim=rep(p, 2)))
  #### Compute df
  if (q==1){
    # 1 df case
    estimate <- sum(L * parlist$fixef)
    var_Lbeta <- drop(L %*% cov_beta %*% t(L))
    se.estimate <- sqrt(var_Lbeta)
    grad_var_Lbeta <- vapply(Jac_list, function(x) {L %*% x %*% t(L)}, numeric(1L))
    satt_denom <- sum(grad_var_Lbeta * (cov_varpar %*% grad_var_Lbeta)) # g'Ag
    ddf <- drop(2 * var_Lbeta^2 / satt_denom) # denominator DF
    tstat <- estimate/se.estimate
    pvalue <- 2 * pt(abs(tstat), df = ddf, lower.tail = FALSE)
    data.frame(estimate, se=se.estimate, tstat, ddf, pvalue)
  } else {
    # multi-df case
    var_Lbeta <- L %*% cov_beta %*% t(L)
    eig_VLbeta <- eigen(var_Lbeta)
    P <- eig_VLbeta$vectors
    d <- eig_VLbeta$values
    PtL <- crossprod(P, L)
    t2 <- drop(PtL %*% parlist$fixef)^2 / d
    Fvalue <- sum(t2) / q
    grad_PLcov <- lapply(1:q, function(m) {
      vapply(Jac_list, function(J) sum(PtL[m, ] * J %*% PtL[m, ]), numeric(1L))
    })
    nu_m <- vapply(1:q, function(m) {
      denom <- sum(grad_PLcov[[m]] * (cov_varpar %*% grad_PLcov[[m]])) # g'Ag
      2*(d[m])^2 / denom # 2d_m^2 / g'Ag
    }, numeric(1L))
    EQ <- sum(as.integer(nu_m>2)*nu_m / (nu_m - 2))
    ddf <- max(2, 2 * EQ / (EQ - q))
    pvalue <- pf(q=Fvalue, df1=q, df2=ddf, lower.tail=FALSE)
    data.frame('F value'=Fvalue, ndf=q, ddf=ddf, pvalue=pvalue,check.names = FALSE)
  }
}


