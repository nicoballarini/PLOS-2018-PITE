# These functions are based on code from the selectiveInference R package
# See https://github.com/selective-inference/R-software

#' Perform selective inference for a specified contrast
#'
#' This function is a modified version of the selectiveInference::fixedCoxLassoInf.
#' This version allows an extra parameter 'contrast' which is a matrix with the
#' contrast to make inference for.
#'
#'
#' @export
fixedCoxLassoInf_eta=function(x, y, status, beta, lambda, alpha=.1, type=c("partial"), tol.beta=1e-5, tol.kkt=0.1,
                          gridrange=c(-100,100), bits=NULL, verbose=FALSE, this.call=NULL, contrast){


  checkargs.xy(x,y)
  if(is.null(status)) stop("Must supply `status' argument")
  if( sum(status==0)+sum(status==1)!=length(y)) stop("status vector must have values 0 or 1")
  if (missing(beta) || is.null(beta)) stop("Must supply the solution beta")
  if (missing(lambda) || is.null(lambda)) stop("Must supply the tuning parameter value lambda")
  checkargs.misc(beta=beta,lambda=lambda,alpha=alpha,
                 gridrange=gridrange,tol.beta=tol.beta,tol.kkt=tol.kkt)
  if (!is.null(bits) && !requireNamespace("Rmpfr",quietly=TRUE)) {
    warning("Package Rmpfr is not installed, reverting to standard precision")
    bits = NULL
  }

  n=nrow(x)
  p=ncol(x)
  nvar=sum(beta!=0)
  nc <- nrow(contrast) # Number of contrast specified in contrast=X.contrast
  pv=vlo=vup=sd=rep(NA, nc)
  ci=tailarea=matrix(NA,nc,2)
  m=beta!=0
  vars=which(m)
  l <- as.matrix(contrast[, vars])
  if(sum(m)>0){
    bhat=beta[beta!=0] #penalized coefs just for active variables
    s2=sign(bhat)
    lhat <- (l) %*% bhat
    sl=drop(sign(lhat))
    #check KKT

    aaa=coxph(Surv(y,status)~x[,m],init=bhat, iter.max=0) # this gives the Cox model at exactly bhat
    # so when we compute gradient and score
    # we are evaluating at the LASSO solution
    # naming of variables could be improved...
    res=residuals(aaa,type="score")
    if(!is.matrix(res)) res=matrix(res,ncol=1)
    scor=colSums(res)
    g=(scor+lambda*s2)/(2*lambda)
    #    cat(c(g,lambda,tol.kkt),fill=T)
    if (any(abs(g) > 1+tol.kkt) )
      warning(paste("Solution beta does not satisfy the KKT conditions",
                    "(to within specified tolerances)"))

    # Hessian of partial likelihood at the LASSO solution
    MM  = vcov(aaa)

    bbar = (bhat + lambda*MM%*%s2)

    A1= -(mydiag(s2))
    b1= -(mydiag(s2)%*%MM)%*%s2*lambda
    temp=max(A1%*%bbar-b1)
    # compute p-values

    for(jj in 1:nc){
      vj=l[jj,]
      junk=TG.pvalue(bbar, A1, b1, vj,MM)

      pv[jj] = junk$pv
      vlo[jj]=junk$vlo
      vup[jj]=junk$vup
      sd[jj]=junk$sd

      junk2=TG.interval(bbar, A1, b1, vj, MM, alpha, gridrange = gridrange)
      ci[jj,]=junk2$int
      tailarea[jj,] = junk2$tailarea

    }
    # JT: these don't seem to be the real one-step estimators
    fit0=coxph(Surv(y,status)~x[,m])
    coef0=l%*%fit0$coef
    se0=sqrt(diag(l%*%fit0$var%*%t(l)))
    zscore0=coef0/se0

    out = list(lambda=lambda,pv=pv,ci=ci,
               tailarea=tailarea,vlo=vlo,vup=vup,sd=sd,
               vars=vars,alpha=alpha,coef0=coef0,zscore0=zscore0,
               call=this.call)
    class(out) = "fixedCoxLassoInf_eta"
  }
  return(out)
}


#' @export
print.fixedCoxLassoInf_eta <- function(x, tailarea=TRUE, ...) {
  cat("\nCall:\n")
  dput(x$call)

  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))

  cat(sprintf("\nTesting results at lambda = %0.3f, with alpha = %0.3f\n",x$lambda,x$alpha))
  cat("",fill=T)
  tab = cbind(round(x$coef0,3),
              round(x$zscore0,3),
              round(x$pv,3),round(x$ci,3))
  colnames(tab) = c("Coef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
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






# Compute the truncation interval and SD of the corresponding Gaussian

TG.limits = function(Z, A, b, eta, Sigma=NULL) {

  target_estimate = sum(as.numeric(eta) * as.numeric(Z))

  if (max(A %*% as.numeric(Z) - b) > 0) {
    warning('Constraint not satisfied. A %*% Z should be elementwise less than or equal to b')
  }

  if (is.null(Sigma)) {
    Sigma = diag(rep(1, n))
  }

  # compute pvalues from poly lemma:  full version from Lee et al for full matrix Sigma

  n = length(Z)
  eta = matrix(eta, ncol=1, nrow=n)
  b = as.vector(b)
  var_estimate = sum(matrix(eta, nrow=1, ncol=n) %*% (Sigma %*% matrix(eta, ncol=1, nrow=n)))
  cross_cov = Sigma %*% matrix(eta, ncol=1, nrow=n)

  resid = (diag(n) - matrix(cross_cov / var_estimate, ncol=1, nrow=n) %*% matrix(eta, nrow=1, ncol=n)) %*% Z
  rho = A %*% cross_cov / var_estimate
  vec = (b - as.numeric(A %*% resid)) / rho

  vlo = suppressWarnings(max(vec[rho < 0]))
  vup = suppressWarnings(min(vec[rho > 0]))

  sd = sqrt(var_estimate)
  return(list(vlo=vlo, vup=vup, sd=sd, estimate=target_estimate))
}

TG.pvalue = function(Z, A, b, eta, Sigma=NULL, null_value=0, bits=NULL) {

  limits.info = TG.limits(Z, A, b, eta, Sigma)

  return(TG.pvalue.base(limits.info, null_value=null_value, bits=bits))
}

TG.interval = function(Z, A, b, eta, Sigma=NULL, alpha=0.1,
                       gridrange=c(-100,100),
                       gridpts=100,
                       griddepth=2,
                       flip=FALSE,
                       bits=NULL) {

  limits.info = TG.limits(Z, A, b, eta, Sigma)

  return(TG.interval.base(limits.info,
                          alpha=alpha,
                          gridrange=gridrange,
                          griddepth=griddepth,
                          flip=flip,
                          bits=bits))
}

TG.interval.base = function(limits.info, alpha=0.1,
                            gridrange=c(-100,100),
                            gridpts=100,
                            griddepth=2,
                            flip=FALSE,
                            bits=NULL) {

  # compute sel intervals from poly lemmma, full version from Lee et al for full matrix Sigma

  param_grid = seq(gridrange[1] * limits.info$sd, gridrange[2] * limits.info$sd, length=gridpts)

  pivot = function(param) {
    tnorm.surv(limits.info$estimate, param, limits.info$sd, limits.info$vlo, limits.info$vup, bits)
  }

  interval = grid.search(param_grid, pivot, alpha/2, 1-alpha/2, gridpts, griddepth)
  tailarea = c(pivot(interval[1]), 1- pivot(interval[2]))

  if (flip) {
    interval = -interval[2:1]
    tailarea = tailarea[2:1]
  }

  # int is not a good variable name, synonymous with integer...
  return(list(int=interval,
              tailarea=tailarea))
}

TG.pvalue.base = function(limits.info, null_value=0, bits=NULL) {
  pv = tnorm.surv(limits.info$estimate, null_value, limits.info$sd, limits.info$vlo, limits.info$vup, bits)
  return(list(pv=pv, vlo=limits.info$vlo, vup=limits.info$vup, sd=limits.info$sd))
}
