# These functions are based on code from the selectiveInference R package
# See https://github.com/selective-inference/R-software

#' Inference for the lasso, with a fixed lambda (selectiveInference package)
#'
#' Compute p-values and confidence intervals for a contrast of the lasso estimate,
#' at a fixed value of the tuning parameter lambda. This function was adapted
#' from the package selectiveInferece to perform the calculations on a contrast
#' of the estimates instead of a single coefficient
#'
#' @param x	Matrix of predictors (n by p);
#' @param y	Vector of outcomes (length n)
#' @param beta Estimated lasso coefficients (e.g., from glmnet). This is of length p (so the intercept is not included as the first component). Be careful! This function uses the "standard" lasso objective
#' 1/2 \|y - x β\|_2^2 + λ \|β\|_1.
#' In contrast, glmnet multiplies the first term by a factor of 1/n. So after running glmnet, to extract the beta corresponding to a value lambda, you need to use beta = coef(obj, s=lambda/n)[-1], where obj is the object returned by glmnet (and [-1] removes the intercept, which glmnet always puts in the first component)
#' @param lambda Value of lambda used to compute beta. See the above warning
#' @param family Response type: "gaussian" (default), "binomial", or "cox" (for censored survival data)
#' @param sigma	Estimate of error standard deviation. If NULL (default), this is estimated using the mean squared residual of the full least squares fit when n >= 2p, and using the standard deviation of y when n < 2p. In the latter case, the user should use estimateSigma function for a more accurate estimate. Not used for family= "binomial", or "cox"
#' @param alpha	Significance level for confidence intervals (target is miscoverage alpha/2 in each tail)
#' @param intercept	Was the lasso problem solved (e.g., by glmnet) with an intercept in the model? Default is TRUE. Must be TRUE for "binomial" family. Not used for 'cox" family, where no intercept is assumed.
#' @param status	Censoring status for Cox model; 1=failurem 0=censored
#' @param type	Contrast type for p-values and confidence intervals: default is "partial"—meaning that the contrasts tested are the partial population regression coefficients, within the active set of predictors; the alternative is "full"—meaning that the full population regression coefficients are tested. The latter does not make sense when p > n.
#' @param tol.beta	Tolerance for determining if a coefficient is zero
#' @param tol.kkt	Tolerance for determining if an entry of the subgradient is zero
#' @param gridrange	Grid range for constructing confidence intervals, on the standardized scale
#' @param bits Number of bits to be used for p-value and confidence interval calculations. Default is NULL, in which case standard floating point calculations are performed. When not NULL, multiple precision floating point calculations are performed with the specified number of bits, using the R package Rmpfr (if this package is not installed, then a warning is thrown, and standard floating point calculations are pursued). Note: standard double precision uses 53 bits so, e.g., a choice of 200 bits uses about 4 times double precision. The confidence interval computation is sometimes numerically challenging, and the extra precision can be helpful (though computationally more costly). In particular, extra precision might be tried if the values in the output columns of tailarea differ noticeably from alpha/2.
#' @param verbose	 Print out progress along the way? Default is FALSE
#' @param contrast A matrix with colnames as the Estimated lasso coefficients names to test and perform CIs
#'
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <-OneData(input,parameters)
#'
#' @export
fixedLassoInf <- function(x, y, beta, lambda, family=c("gaussian","binomial","cox"),
                          intercept=TRUE, status=NULL,
                          sigma=NULL, alpha=0.1,
                          type=c("partial","full"), tol.beta=1e-5, tol.kkt=0.1,
                          gridrange=c(-100,100), bits=NULL, verbose=FALSE) {

  family = match.arg(family)
  this.call = match.call()
  type = match.arg(type)

  if(family=="binomial")  {
    if(type!="partial") stop("Only type= partial allowed with binomial family")
    out=fixedLogitLassoInf(x,y,beta,lambda,alpha=alpha, type="partial", tol.beta=tol.beta, tol.kkt=tol.kkt,
                           gridrange=gridrange, bits=bits, verbose=verbose,this.call=this.call)
    return(out)
  }
  else if(family=="cox")  {
    if(type!="partial") stop("Only type= partial allowed with Cox family")
    out=fixedCoxLassoInf(x,y,status,beta,lambda,alpha=alpha, type="partial",tol.beta=tol.beta,
                         tol.kkt=tol.kkt, gridrange=gridrange, bits=bits, verbose=verbose,this.call=this.call)
    return(out)
  }

  else{
    checkargs.xy(x,y)
    if (missing(beta) || is.null(beta)) stop("Must supply the solution beta")
    if (missing(lambda) || is.null(lambda)) stop("Must supply the tuning parameter value lambda")
    checkargs.misc(beta=beta,lambda=lambda,sigma=sigma,alpha=alpha,
                   gridrange=gridrange,tol.beta=tol.beta,tol.kkt=tol.kkt)
    if (!is.null(bits) && !requireNamespace("Rmpfr",quietly=TRUE)) {
      warning("Package Rmpfr is not installed, reverting to standard precision")
      bits = NULL
    }

    n = nrow(x)
    p = ncol(x)
    beta = as.numeric(beta)
    if (length(beta) != p) stop("Since family='gaussian', beta must have length equal to ncol(x)")

    # If glmnet was run with an intercept term, center x and y
    if (intercept==TRUE) {
      obj = standardize(x,y,TRUE,FALSE)
      x = obj$x
      y = obj$y
    }

    # Check the KKT conditions
    g = t(x)%*%(y-x%*%beta) / lambda
    if (any(abs(g) > 1+tol.kkt * sqrt(sum(y^2))))
      warning(paste("Solution beta does not satisfy the KKT conditions",
                    "(to within specified tolerances)"))

    vars = which(abs(beta) > tol.beta / sqrt(colSums(x^2)))
    if(length(vars)==0){
      cat("Empty model",fill=T)
      return()
    }
    if (any(sign(g[vars]) != sign(beta[vars])))
      warning(paste("Solution beta does not satisfy the KKT conditions",
                    "(to within specified tolerances). You might try rerunning",
                    "glmnet with a lower setting of the",
                    "'thresh' parameter, for a more accurate convergence."))

    # Get lasso polyhedral region, of form Gy >= u
    out = fixedLasso.poly(x,y,beta,lambda,vars)
    G = out$G
    u = out$u

    # Check polyhedral region
    tol.poly = 0.01
    if (min(G %*% y - u) < -tol.poly * sqrt(sum(y^2)))
      warning(paste("Polyhedral constraints not satisfied; you must recompute beta",
                 "more accurately. With glmnet, make sure to use exact=TRUE in coef(),",
                 "and check whether the specified value of lambda is too small",
                 "(beyond the grid of values visited by glmnet).",
                 "You might also try rerunning glmnet with a lower setting of the",
                 "'thresh' parameter, for a more accurate convergence."))

    # Estimate sigma
    if (is.null(sigma)) {
      if (n >= 2*p) {
        oo = intercept
        sigma = sqrt(sum(lsfit(x,y,intercept=oo)$res^2)/(n-p-oo))
      }
      else {
        sigma = sd(y)
        warning(paste(sprintf("p > n/2, and sd(y) = %0.3f used as an estimate of sigma;",sigma),
                      "you may want to use the estimateSigma function"))
      }
    }

    k = length(vars)
    pv = vlo = vup = numeric(k)
    vmat = matrix(0,k,n)
    ci = tailarea = matrix(0,k,2)
    sign = numeric(k)

    if (type=="full" & p > n)
      warning(paste("type='full' does not make sense when p > n;",
                    "switching to type='partial'"))

    if (type=="partial" || p > n) {
      xa = x[,vars,drop=F]
      M = pinv(crossprod(xa)) %*% t(xa)
    }
    else {
      M = pinv(crossprod(x)) %*% t(x)
      M = M[vars,,drop=F]
    }

    for (j in 1:k) {
      if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))

      vj = M[j,]
      mj = sqrt(sum(vj^2))
      vj = vj / mj        # Standardize (divide by norm of vj)
      sign[j] = sign(sum(vj*y))
      vj = sign[j] * vj
      a = poly.pval(y,G,u,vj,sigma,bits)
      pv[j] = a$pv
      vlo[j] = a$vlo * mj # Unstandardize (mult by norm of vj)
      vup[j] = a$vup * mj # Unstandardize (mult by norm of vj)
      vmat[j,] = vj * mj * sign[j]  # Unstandardize (mult by norm of vj)

      a = poly.int(y,G,u,vj,sigma,alpha,gridrange=gridrange,
                   flip=(sign[j]==-1),bits=bits)
      ci[j,] = a$int * mj # Unstandardize (mult by norm of vj)
      tailarea[j,] = a$tailarea
    }

    out = list(type=type,lambda=lambda,pv=pv,ci=ci,
               tailarea=tailarea,vlo=vlo,vup=vup,vmat=vmat,y=y,
               vars=vars,sign=sign,sigma=sigma,alpha=alpha,
               sd=sigma*sqrt(rowSums(vmat^2)),
               coef0=vmat%*%y,
               call=this.call)
    class(out) = "fixedLassoInf"
    return(out)
  }
}

# Lasso inference function (for fixed lambda). Note: here we are providing inference
# for the solution of
# min 1/2 || y - \beta_0 - X \beta ||_2^2 + \lambda || \beta ||_1
#' Inference for the lasso, with a fixed lambda (selectiveInference package)
#'
#' Compute p-values and confidence intervals for a contrast of the lasso estimate,
#' at a fixed value of the tuning parameter lambda. This function was adapted
#' from the package selectiveInferece to perform the calculations on a contrast
#' of the estimates instead of a single coefficient
#'
#' @param x	Matrix of predictors (n by p);
#' @param y	Vector of outcomes (length n)
#' @param beta Estimated lasso coefficients (e.g., from glmnet). This is of length p (so the intercept is not included as the first component). Be careful! This function uses the "standard" lasso objective
#' 1/2 \|y - x β\|_2^2 + λ \|β\|_1.
#' In contrast, glmnet multiplies the first term by a factor of 1/n. So after running glmnet, to extract the beta corresponding to a value lambda, you need to use beta = coef(obj, s=lambda/n)[-1], where obj is the object returned by glmnet (and [-1] removes the intercept, which glmnet always puts in the first component)
#' @param lambda Value of lambda used to compute beta. See the above warning
#' @param family Response type: "gaussian" (default), "binomial", or "cox" (for censored survival data)
#' @param sigma	Estimate of error standard deviation. If NULL (default), this is estimated using the mean squared residual of the full least squares fit when n >= 2p, and using the standard deviation of y when n < 2p. In the latter case, the user should use estimateSigma function for a more accurate estimate. Not used for family= "binomial", or "cox"
#' @param alpha	Significance level for confidence intervals (target is miscoverage alpha/2 in each tail)
#' @param intercept	Was the lasso problem solved (e.g., by glmnet) with an intercept in the model? Default is TRUE. Must be TRUE for "binomial" family. Not used for 'cox" family, where no intercept is assumed.
#' @param status	Censoring status for Cox model; 1=failurem 0=censored
#' @param type	Contrast type for p-values and confidence intervals: default is "partial"—meaning that the contrasts tested are the partial population regression coefficients, within the active set of predictors; the alternative is "full"—meaning that the full population regression coefficients are tested. The latter does not make sense when p > n.
#' @param tol.beta	Tolerance for determining if a coefficient is zero
#' @param tol.kkt	Tolerance for determining if an entry of the subgradient is zero
#' @param gridrange	Grid range for constructing confidence intervals, on the standardized scale
#' @param bits Number of bits to be used for p-value and confidence interval calculations. Default is NULL, in which case standard floating point calculations are performed. When not NULL, multiple precision floating point calculations are performed with the specified number of bits, using the R package Rmpfr (if this package is not installed, then a warning is thrown, and standard floating point calculations are pursued). Note: standard double precision uses 53 bits so, e.g., a choice of 200 bits uses about 4 times double precision. The confidence interval computation is sometimes numerically challenging, and the extra precision can be helpful (though computationally more costly). In particular, extra precision might be tried if the values in the output columns of tailarea differ noticeably from alpha/2.
#' @param verbose	 Print out progress along the way? Default is FALSE
#' @param contrast A matrix with colnames as the Estimated lasso coefficients names to test and perform CIs
#'
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <-OneData(input,parameters)
#'
#' @export
fixedLassoInf_eta <- function(x, y, beta, lambda, family=c("gaussian","binomial","cox"),
                          intercept=TRUE, status=NULL,
                          sigma=NULL, alpha=0.1,
                          type=c("partial","full"), tol.beta=1e-5, tol.kkt=0.1,
                          gridrange=c(-100,100), bits=NULL, verbose=FALSE,
                          contrast) {

 family = match.arg(family)
  this.call = match.call()
  type = match.arg(type)

  if(family=="binomial")  {
      if(type!="partial") stop("Only type= partial allowed with binomial family")
       out=fixedLogitLassoInf(x,y,beta,lambda,alpha=alpha, type="partial", tol.beta=tol.beta, tol.kkt=tol.kkt,
                     gridrange=gridrange, bits=bits, verbose=verbose,this.call=this.call)
                      return(out)
                    }
else if(family=="cox")  {
    if(type!="partial") stop("Only type= partial allowed with Cox family")
     out=fixedCoxLassoInf(x,y,status,beta,lambda,alpha=alpha, type="partial",tol.beta=tol.beta,
          tol.kkt=tol.kkt, gridrange=gridrange, bits=bits, verbose=verbose,this.call=this.call)
                      return(out)
                    }

else{
  checkargs.xy(x,y)
  if (missing(beta) || is.null(beta)) stop("Must supply the solution beta")
  if (missing(lambda) || is.null(lambda)) stop("Must supply the tuning parameter value lambda")
  checkargs.misc(beta=beta,lambda=lambda,sigma=sigma,alpha=alpha,
                 gridrange=gridrange,tol.beta=tol.beta,tol.kkt=tol.kkt)
  if (!is.null(bits) && !requireNamespace("Rmpfr",quietly=TRUE)) {
    warning("Package Rmpfr is not installed, reverting to standard precision")
    bits = NULL
  }

  n = nrow(x)
  p = ncol(x)
  beta = as.numeric(beta)
  if (length(beta) != p) stop("Since family='gaussian', beta must have length equal to ncol(x)")

  # If glmnet was run with an intercept term, center x and y
  if (intercept==TRUE) {
    obj = standardize(x,y,TRUE,FALSE)
    x = obj$x
    y = obj$y
  }

  # Check the KKT conditions
  g = t(x)%*%(y-x%*%beta) / lambda
  if (any(abs(g) > 1+tol.kkt * sqrt(sum(y^2))))
    warning(paste("Solution beta does not satisfy the KKT conditions",
                  "(to within specified tolerances)"))

  vars = which(abs(beta) > tol.beta / sqrt(colSums(x^2)))
  if(length(vars)==0){
      cat("Empty model",fill=T)
      return()
  }
  if (any(sign(g[vars]) != sign(beta[vars])))
    warning(paste("Solution beta does not satisfy the KKT conditions",
                  "(to within specified tolerances). You might try rerunning",
                  "glmnet with a lower setting of the",
                  "'thresh' parameter, for a more accurate convergence."))

  # Get lasso polyhedral region, of form Gy >= u
  out = fixedLasso.poly(x,y,beta,lambda,vars)
  G = out$G
  u = out$u

  # Check polyhedral region
  tol.poly = 0.01
  if (min(G %*% y - u) < -tol.poly * sqrt(sum(y^2))){
    warning(paste("Polyhedral constraints not satisfied; you must recompute beta",
               "more accurately. With glmnet, make sure to use exact=TRUE in coef(),",
               "and check whether the specified value of lambda is too small",
               "(beyond the grid of values visited by glmnet).",
               "You might also try rerunning glmnet with a lower setting of the",
               "'thresh' parameter, for a more accurate convergence."))}

  # Estimate sigma
  if (is.null(sigma)) {
    if (n >= 2*p) {
      oo = intercept
      sigma = sqrt(sum(lsfit(x,y,intercept=oo)$res^2)/(n-p-oo))
    }
    else {
      sigma = sd(y)
      warning(paste(sprintf("p > n/2, and sd(y) = %0.3f used as an estimate of sigma;",sigma),
                    "you may want to use the estimateSigma function"))
    }
  }

  nc<-nrow(contrast) # Number of contrast specified in contrast=X.contrast
  k = length(vars)
  pv = vlo = vup = numeric(nc)
  exp. = var. = numeric(nc)
  vmat = matrix(0,nc,n)
  ci = tailarea = matrix(0,nc,2)
  ci_truncnorm = tailarea_truncnorm = matrix(0,nc,2)
  sign = numeric(nc)

  if (type=="full" & p > n)
      warning(paste("type='full' does not make sense when p > n;",
                    "switching to type='partial'"))

  if (type=="partial" || p > n) {
    xa = x[,vars,drop=F]
    M = pinv(crossprod(xa)) %*% t(xa)
  }  else {
    M = pinv(crossprod(x)) %*% t(x)
    M = M[vars,,drop=F]
  }

  ## The following lines implement contrast instead of inference for parameters.
  ## The contrast matrix is provided in the options of the function call
  ## However, this contrast matrix should be based on the full X matrix,
  ## i.e. taking all coefficients into account
  ## the Matrix c.i will adapt the contrast to the selected variables in the lasso.
  ## The matrix eta is then the collection of etas across individuals.
  ## Each row of eta correspond to the contrast provided in X1

  c.i <- t(contrast[, names(vars), drop = FALSE])
  eta <- t(c.i) %*% M
  for (j in 1:nc) {
    vj = eta[j,]
    if (all((vj)==0)) next
    mj = sqrt(sum(vj^2))
    vj = vj / mj        # Standardize (divide by norm of vj)
    sign[j] = sign(sum(vj*y))
    vj = sign[j] * vj
    a = poly.pval(y,G,u,vj,sigma,bits)
    pv[j] = a$pv
    vlo[j] = a$vlo * mj # Unstandardize (mult by norm of vj)
    vup[j] = a$vup * mj # Unstandardize (mult by norm of vj)
    vmat[j,] = vj * mj * sign[j]  # Unstandardize (mult by norm of vj)
    a = poly.int(y,G,u,vj,sigma,alpha,gridrange=gridrange,
                 flip=(sign[j]==-1),bits=bits)
    ci[j,] = a$int * mj # Unstandardize (mult by norm of vj)
    tailarea[j,] = a$tailarea
  }

  out = list(type=type,lambda=lambda,pv=pv,ci=ci,
             ci_truncnorm=ci_truncnorm,tailarea_truncnorm=tailarea_truncnorm,
             exp. = exp., var. = var.,
    tailarea=tailarea,vlo=vlo,vup=vup,vmat=vmat,y=y,
    vars=vars,sign=sign,sigma=sigma,alpha=alpha,
    sd=sigma*sqrt(rowSums(vmat^2)),
    coef0=vmat%*%y,
    call=this.call)
  class(out) = "fixedLassoInf_eta"
  return(out)
}
}

#############################-


fixedLasso.poly=
function(x, y, beta, lambda, a) {
  xa = x[,a,drop=F]
  xac = x[,-a,drop=F]
  # cat(xac)
  xai = pinv(crossprod(xa))
  xap = xai %*% t(xa)
  za = sign(beta[a])
  if (length(za)>1) dz = diag(za)
  if (length(za)==1) dz = matrix(za,1,1)
  #NOTE: inactive constraints not needed below!

  G = -rbind(
   #   1/lambda * t(xac) %*% P,
   # -1/lambda * t(xac) %*% P,
    -dz %*% xap
      )
     lambda2=lambda
     if(length(lambda)>1) lambda2=lambda[a]
  u = -c(
   #   1 - t(xac) %*% t(xap) %*% za,
   #   1 + t(xac) %*% t(xap) %*% za,
    -lambda2 * dz %*% xai %*% za)

  return(list(G=G,u=u))
}

##############################-

print.fixedLassoInf <- function(x, tailarea=TRUE, ...) {
  cat("\nCall:\n")
  dput(x$call)
  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))
  cat(sprintf("\nTesting results at lambda = %0.3f, with alpha = %0.3f\n",x$lambda,x$alpha))
  cat("",fill=T)
  tab = cbind(x$vars,
    round(x$coef0,3),
    round(x$coef0 / x$sd,3),
    round(x$pv,3),round(x$ci,3))
  colnames(tab) = c("Var", "Coef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
  if (tailarea) {
    tab = cbind(tab,round(x$tailarea,3))
    colnames(tab)[(ncol(tab)-1):ncol(tab)] = c("LowTailArea","UpTailArea")
  }
  rownames(tab) = rep("",nrow(tab))
  print(tab)

  cat(sprintf("\nNote: coefficients shown are %s regression coefficients\n",
              ifelse(x$type=="partial","partial","full")))
  invisible()
}


##############################-

print.fixedLassoInf_eta <- function(x, tailarea=TRUE, ...) {
  cat("\nCall:\n")
  dput(x$call)

  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))

  cat(sprintf("\nTesting results at lambda = %0.3f, with alpha = %0.3f\n",x$lambda,x$alpha))
  cat("",fill=T)
  tab = cbind(1:length(x$coef0),
              round(x$coef0,3),
              round(x$coef0 / x$sd,3),
              round(x$pv,3),round(x$ci,3))
  colnames(tab) = c("Subject Id", "Coef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
  if (tailarea) {
    tab = cbind(tab,round(x$tailarea,3))
    colnames(tab)[(ncol(tab)-1):ncol(tab)] = c("LowTailArea","UpTailArea")
  }
  rownames(tab) = rep("",nrow(tab))
  print(tab)

  cat(sprintf("\nNote: coefficients shown are %s regression coefficients\n",
              ifelse(x$type=="partial","partial","full")))
  invisible()
}
