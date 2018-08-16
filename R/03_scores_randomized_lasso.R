#' Fit linear model to create a score using randomized Lasso
#'
#' Performs randomized lasso in a dataset with main effects and interactions.
#' calculates PITE and CI
#'
#' @param dataset A data.frame created with \code{OneData}
#' @param input A list of objects created with \code{MakeInput}
#' @param alpha	Significance level for confidence intervals (target is
#'  miscoverage alpha/2 in each tail)
#' @param parameters A \code{MakeParameters} object
#' @param verbose logical. whether to print a message with results.
#' @param alpha	Significance level for confidence intervals (target is
#'  miscoverage alpha/2 in each tail)
#' @param lam_frac multiplier for the lagrange parameter
#' @param lambda penalization parameter. one of "lagrange", "lambda.min", "lambda.1se" or a number
#' @param nfolds number of folds for crossvalidations. to be passed to glmnet
#' @param perturb_frac scale of Gaussian noise added to the response. this number is then multiplied to sd(y)
#' @param pite.ci logical. whether to calculate ci for pite
#' @param n.pite a number specifying the number of subjects in the dataset for which to calculate the ci for the pite if pite.ci is TRUE
#' @param ndraw  Number of samples of optimization variables to sample.
#' @param burnin How many samples of optimization variable to discard.
#'
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <- OneData(input, parameters)
#' score.lasso.added.noise(dataset, input, parameters)
#'
#' @export
score.lasso.added.noise <- function(dataset,
                                    input,
                                    parameters = NULL,
                                    verbose = TRUE,
                                    alpha = 0.05,
                                    lambda = "lagrange",
                                    lam_frac = 1,
                                    nfolds = 10,
                                    perturb_frac = 0.2,
                                    pite.ci = TRUE,
                                    n.pite = NULL,
                                    ndraw = 15000,
                                    burnin = 2000) {
  # Set Initial parameters -----------------------------------------------------
  start_python()
  selModel <- trueModel <- selModelTrue <- NA
  n_biom<-input$n_biom
  N <- n <- nrow(dataset)      # Total number of subjects
  if (is.null(n.pite)){
    n.pite <- N
  }
  scores = data.frame(Dx = numeric(N))
  K <- n_biom*2+2         # Number of coefficients in the model

  formula <- as.formula(paste0("y ~ treatment + ",
                               paste(paste0(" mkr",1:n_biom,
                                            collapse = " +"),
                                     paste0(" treatment*mkr", 1:n_biom,
                                            collapse = " +"), sep = " +")))
  X. <- model.matrix(formula, data = dataset)
  Y <- as.matrix(dataset[, "y"])
  # Standardize predictors
  x = scale(X.[, -1], TRUE, TRUE) # We dont want intercept in there.
  center.X <- attr(x, "scaled:center")
  scale.X  <- attr(x, "scaled:scale")
  # Standardize Response
  y = scale(Y, TRUE, TRUE)
  center.Y <- attr(y, "scaled:center")
  scale.Y  <- attr(y, "scaled:scale")

  if (!is.null(parameters)){
    n_biom_nopred <- input$n_biom - input$n_biom_pred
    prognosis <- c(parameters$prognosis, rep(0, n_biom_nopred))
    predictive <- c(parameters$predictive, rep(0, n_biom_nopred))
    coef.true <- c(input$a_, input$b_, prognosis, predictive)
    names(coef.true) <- colnames(X.)
    trueModel <- coef.true != 0
  }

  # Fit the Lasso --------------------------------------------------------------
  # X1 is the contrast matrix to create scores (Include intercept or not!)
  # The contrast matrix should have colnames as the names of the variables in the coef
  X1 <- cbind(0, 1, matrix(0, N, n_biom), X.[, 3:(n_biom + 2)])
  contrast = t(t(2 * X1) * scale.Y / c(1, scale.X))
  X = x
  out <- additive_noise(X = x,
                        y = y,
                        sigma = sd(y),
                        contrast = t(t(2 * X1) * scale.Y / c(1, scale.X)),
                        lambda = lambda,
                        nfolds = nfolds,
                        lam_frac = lam_frac,
                        perturb_frac = perturb_frac,
                        coverage = 1 - alpha,
                        ndraw = ndraw,
                        compute_intervals = TRUE,
                        burnin = burnin,
                        pite.ci = pite.ci,
                        n.pite = n.pite)
  bestlam <- out$bestlam
  coef <- out$coef
  selModel <- coef != 0
  coef.unsc <- coef * scale.Y / scale.X
  coef.dx <- coef.unsc[c(1, (n_biom + 1 + 1):(2 * n_biom + 1))]
  nvars <- sum(coef.dx != 0)

  # Null model Selected! -------------------------------------------------------
  if (sum(selModel) == 0 | nvars == 0) {
    warning("All estimated coefficients are 0. CIs are set to (0,0) since no model to condition for.")
    # Store score, variance and CIs in the dataset
    scores$Dx.ll <- 0
    # cat(dataset$DxLasso.ll.tibs)
    scores$Dx.ul <- 0
    # cat(dataset$DxLasso.ul.tibs)
    scores$Dx.cover <- 1 * (scores$Dx.ll <= dataset$TrueTrtEff &
                            scores$Dx.ul >= dataset$TrueTrtEff)

    scores$Dx.w  <- 0
    scores$PITE <- 0
    scores$Dx.cover <- 1
    return(list(scores = scores,
                null.model = TRUE,
                nvars = 0,
                Lasso.output = data.frame(term = names(coef),
                                          estimate = round(coef,6),
                                          statistic = NA,
                                          p.value = NA,
                                          LowConfPt = NA,
                                          UpConfPt = NA, row.names = NULL)))
  }
  selected     <- selModel[c(1, (n_biom+1+1):(2*n_biom+1))]
  selected.all <- selModel
  nvar <- sum(coef!=0)
  if (verbose == TRUE) cat("\n lambda:",round(bestlam, 5),
                           ", nvar:",nvar," - ",round(coef,2),
                           " - \n")
  # *Store results from regression ---------------------------------------------
  Lasso.output <- data.frame(term = names(coef.unsc),
                             estimate = round(coef.unsc, 6),
                             statistic = NA,
                             p.value = NA,
                             LowConfPt = NA,
                             UpConfPt = NA, row.names = NULL)
  Lasso.output$p.value[selected.all] <- out$Lasso.output2$pval
  Lasso.output$LowConfPt[selected.all] <- out$Lasso.output2$ll * scale.Y / scale.X[selected.all]
  Lasso.output$UpConfPt[selected.all] <- out$Lasso.output2$ul * scale.Y / scale.X[selected.all]
  ll <- out$Lasso.output2$ll
  ul <- out$Lasso.output2$ul

  # *Create estimated PITE ------------------------------------------------------
  # Retain only those coef in the score
  # Create matrix for calculating score
  X.dx <- as.matrix(cbind(1, dataset[, paste0("mkr", 1:n_biom)]))
  colnames(X.dx) <- c("treatment", paste0("treatment:mkr", 1:n_biom))
  Dx.m <- 2 * X.dx %*% coef.dx #Create score or PITE
  # Store score, variance and CIs in the dataset
  scores$Dx <- drop(Dx.m)
  out$out.contrast$observed

  if (!is.null(parameters)){
    coef.trueDx <- coef.true[c(2, (n_biom + 2 + 1):(2 * n_biom + 2))]  #Retain only those coef in the score
  }

  # X1 is the contrast matrix to create scores (Include intercept or not!)
  # The contrast matrix should have colnames as the names of the variables in the coef
  X1 <- cbind(1, matrix(0, N, n_biom), X[, 3:(n_biom + 2)])
  colnames(X1) <- names(coef.unsc)
  sigmahat = NULL
  scores$Dx    <- out$out.contrast$observed
  scores$Dx.ll <- out$out.contrast$ll
  scores$Dx.ul <- out$out.contrast$ul
  scores$Dx.cover <- 1*(scores$Dx.ll <= dataset$TrueTrtEff &
                        scores$Dx.ul >= dataset$TrueTrtEff)
  scores$Dx.w <- scores$Dx.ul - scores$Dx.ll

  ###--------------------------------------------------------------------------#
  ### Perform Expectation in Reduced Model ----
  ###--------------------------------------------------------------------------#
  if (!is.null(parameters)){
    # Flag for selected variables in the model
    selVars.Dx <- selected.all[grepl("treatment", names(selected.all))]
    selVars <- selVars.Dx[-1] # Excluding the treatment coefficient
    # First calculate the Term for selected biomarkers
    if (selVars[1] == FALSE & length(selVars) == 1) {
      SigmaS   <- 1
      SigmaScS <- 1
      muScS <- 0
      muScd <- 0
      XSd <- 0
    } else {
      S <- diag(1 * selVars)
      colnames(S) <- rownames(S) <- names(selVars)
      X.M.dx  <- X.dx[, -1, drop = FALSE][,  selVars, drop = FALSE] # matrix with selected variables
      X.Mc.dx <- X.dx[, -1, drop = FALSE][, !selVars, drop = FALSE] # matrix with no selected variables
      X_M <- X[, selected.all, drop = FALSE]
      X_MX <- solve(crossprod(X_M)) %*% t(X_M) %*% X.
      X.dx <- as.matrix(cbind(1, dataset[, paste0("mkr", 1:n_biom)]))
      X.dx.sel <- X.dx[, selected]
      # If they are not independent, we need to use the following.
      # contribution from the selected biomarkers
      XSd <- X.dx[, -1] %*% S %*% coef.trueDx[-1]

      # Obtain the mean for all variables from the parameters input
      # For the no predictive ones, it was set to 0
      n_biom_nopred <- input$n_biom - input$n_biom_pred
      mu <- c(parameters$means, rep(0, n_biom_nopred))
      names(mu) <- paste0("mkr", 1:n_biom)
      # Partition mean vector inot selected S and not selected Sc
      muSc <- mu[!selVars]
      muS  <- mu[selVars]
      if (sum(!selVars) > 0){ # If we dont select some variables
        if (sum(!selVars) == n_biom){ # But Only intercept was selected
          XSd <- 0
          muScd <- 0
        } else { # or at least one variable is selected
          covm <- BiomCorrMatrix(input, parameters)$covm # get the cov matrix
          colnames(covm) <- rownames(covm) <- names(selVars)
          # Get the partitions on Sigma needed for calculating conditional
          # expectation
          SigmaS   <- covm[selVars, selVars]
          SigmaScS <- covm[!selVars, selVars]
          # Calculate contribution of the non selected biomarkers
          # To sum the vector by colum, we use sweep with a -
          # First calculate the conditional mu
          muScS <- sweep(t(SigmaScS %*% solve(SigmaS) %*% t(sweep(X.M.dx, 2, muS))),
                         2, -muSc)
          # Then the Sc matrix
          Sc <- diag(1 * !selVars)[!selVars, ]
          # And then the full contribution of the non-selected
          muScd <- muScS %*% Sc %*% coef.trueDx[-1]
        }
      } else {# If all are selected, then no contribution from non-selected
        muScd <- 0
      }
    }
    # Now calculate the PITE
    PITE.M.LAN <- 2*(coef.trueDx[1] * selVars.Dx[1] + XSd + muScd)
    scores$PITE <- drop(PITE.M.LAN)
    scores$Dx.cover <- 1 * (scores$Dx.ll <= PITE.M.LAN &
                            scores$Dx.ul >= PITE.M.LAN)
    selModelTrue <- all.equal(selModel, trueModel)
  }

  # Output ---------------------------------------------------------------------
  if (verbose) cat("Lasso: Done - ")
  list(scores = scores,
       nvars = nvars,
       Lasso.output = Lasso.output,
       x = x,
       y = y,
       y_star = out$y_star,
       gamma = out$gamma,
       sigma = out$sigma,
       perturb_frac = perturb_frac,
       N = N,
       bestlam = bestlam,
       coef = coef,
       coef.unsc = coef.unsc,
       selected.all = selected.all,
       trueSel = selModelTrue,
       selModel = selModel,
       null.model = FALSE)
}

## Additive Noise ---------------------------------------------------------------------
additive_noise <- function(X,
                           y,
                           sigma,
                           contrast,
                           lambda = "lambda.min",
                           nfolds = 10,
                           lam_frac = 1.,
                           perturb_frac = 0.2,
                           coverage = 0.95,
                           ndraw = 8000,
                           compute_intervals = TRUE,
                           burnin = 2000,
                           pite.ci = TRUE,
                           n.pite = NULL) {
  # Additive noise LASSO.
  # This function implements randomized lasso as in the selectiveInference
  # python package
  # Parameters -
  # y : np.float Response vector
  # X : np.float Design matrix
  # sigma : np.float Noise variance
  # lam_frac : float (optional) Multiplier for choice of $\lambda$. Defaults to 2.
  # perturb_frac : float (optional) How much noise to add? Noise added has variance proportional to existing variance.
  # coverage : float Coverage for selective intervals. Defaults to 0.95.
  # ndraw : int (optional) How many draws to keep from Gibbs hit-and-run sampler. Defaults to 8000.
  # burnin : int (optional) Defaults to 2000.
  # compute_intervals : bool (optional) Compute selective intervals?
  #
  # Returns -
  #   results : [(variable, pvalue, interval)
  #              Indices of active variables,
  #              selected (twosided) pvalue and selective interval.
  #              If splitting, then each entry also includes
  #              a (split_pvalue, split_interval) using stage_two
  #              for inference.
  #              randomized_lasso : `lasso`
  #              Results of fitting LASSO to randomized data.
  n <- nrow(X)
  p <- ncol(X)
  alpha = 1 - coverage
  # Add some noise to y and fit the LASSO at a fixed lambda
  gamma = sqrt(perturb_frac) * sigma
  sigma_star = sqrt(sigma^2 + gamma^2)
  y_star = y + rnorm(n) * gamma
  y_star = y_star / sd(y_star)

  if (lambda == "lagrange"){
    standard_normal <- matrix(rnorm(n * 5000, 0, 1), nrow = n, ncol = 5000)
    XT.e <- abs(t(X) %*% standard_normal)
    lam = lam_frac * mean(apply(XT.e, 2, max)) * sigma_star ## Instead of doing this, we standardize the variable so that sigma_star=1
    lam
  }
  if (lambda == "lambda.min"){
    cv.out <- glmnet::cv.glmnet(x = X, y = y_star,
                                standardize = FALSE,
                                nfolds = nfolds)
    lam <- drop(cv.out[[lambda]] * n)
    cat(lam)
  }
  if (is.numeric(lambda)){
    lam <- lambda
  }
  randomized_lasso = glmnet::glmnet(x = X, y = y_star,
                                    standardize = FALSE, intercept = FALSE)
  coef <- coef(randomized_lasso, s = lam / n, exact = TRUE,
               x = X, y = y_star)[-1, 1]
  Lasso.output <- data.frame(term = names(coef),
                             estimate = coef,
                             statistic = NA,
                             p.value = NA,
                             LowConfPt = NA,
                             UpConfPt = NA, row.names = NULL)
  active <- coef != 0
  sum(active) == 0


  if (sum(active) == 0){
    Lasso.output2 <- data.frame(term = names(coef),
                                estimate = round(coef,6),
                                statistic = NA,
                                p.value = NA,
                                LowConfPt = NA,
                                UpConfPt = NA, row.names = NULL)
    return(list(coef = coef,
                Lasso.output2=Lasso.output2,
                bestlam = lam))
  }

  active_signs <- sign(coef[active])
  # Form the constraint matrix on (y, y^*)
  X_E = X[, active, drop = F]
  X_Ei = MASS::ginv(X_E)  ## This is equivalent to solve(t(X_E)%*%X_E)%*%t(X_E)
  Cov_E = X_Ei %*% t(X_Ei)    ## This is the same as solve(t(X_E)%*%X_E) which is (X'X)^(-1)
  W_E   = Cov_E %*% (active_signs)

  beta_E = X_Ei %*% (y)
  res <- list(0)
  # Coefficients --------------------------------------------------------------
  # compute each pvalue
  neta <- dim(X_E)[2]
  out.all <- data.frame(observed = numeric(neta),
                        pval = numeric(neta),
                        ll = numeric(neta),
                        ul = numeric(neta))
  for (j in 1:neta){
    s_obs = sum(active)
    keep = rep(1, s_obs)
    keep[j] = 0
    # form the 2s Gaussian vector we will condition on
    X_minus_j = X_E[, -j, drop = FALSE]
    if (dim(X_minus_j)[2] == 0){
      P_minus_j = matrix(0, n, n)
    } else {
      P_minus_j = X_minus_j %*% (MASS::ginv(X_minus_j))
    }

    theta_E = active_signs * (X_Ei %*% (P_minus_j %*% (y)) - lam * W_E)
    scale = sqrt(Cov_E[j, j])
    kappa = as.numeric(1. / scale**2)
    alpha_E = kappa * active_signs * Cov_E[, j]
    A = rbind(-t(alpha_E),
              diag(s_obs))
    con = list(linear_part = A, offset = theta_E, mean = rep(0, s_obs + 1))
    cov = matrix(0, s_obs + 1, s_obs + 1)
    cov[1,1] = scale**2 * sigma**2
    cov[2:(s_obs + 1), 2:(s_obs + 1)] =
      Cov_E * gamma**2 * outer(active_signs, active_signs)
    con$covariance = cov
    initial = rep(0, s_obs + 1)
    initial[1] = beta_E[j]
    initial[2:(s_obs + 1)] = -X_Ei %*% (y_star - y) * active_signs
    eta = rep(0, s_obs + 1)
    eta[1] = 1.
    observed = sum(initial * eta)

    ######################################################-
    empty.out <- c(observed = NA, pval = NA, ci = c(NA, NA))
    out = R.utils::evalWithTimeout(gibbs_test(con,
                                              initial,
                                              eta,
                                              ndraw = ndraw,
                                              burnin = burnin,
                                              alpha = alpha),
                                   timeout = 10, onTimeout = NULL)
    if (!is.null(out)){
      ci <- drop(eta %*% con$covariance %*% eta) * out$ci
      out.all[j, ] <- c(out$observed, out$pval, ci)
    } else {
      out.all[j, ] <- c(NA,NA,NA,NA)
    }
  }

  ncontrast <- (dim(contrast)[1])
  # Contrasts -----------------------------------------------------------------
  if (pite.ci){ # We only calculate ci for pite if pite.ci is TRUE
    out.contrast <- data.frame(observed = rep(NA, ncontrast),
                               pval = rep(NA, ncontrast),
                               ll = rep(NA, ncontrast),
                               ul = rep(NA, ncontrast))
    X1 <- contrast
    for (j in 1:n.pite){
      s_obs = sum(active)
      # form the 2s Gaussian vector we will condition on
      l <- X1[j, c(FALSE, active)]
      if (all(l == 0)) {
        out.contrast[j, ] <- c(0, 0, 0, 0)
        next()
      }
      ej = l
      A = ej
      At = t(A)
      beta <- solve(t(X_E) %*% X_E) %*% t(X_E) %*% y
      SIGMA = solve(t(X_E) %*% X_E)
      SIGMA_A = At %*% SIGMA %*% A
      C = SIGMA %*% A %*% solve(SIGMA_A)
      T_A = (diag(s_obs) - C %*% At) %*% beta
      theta_E = active_signs * (T_A - lam * W_E)
      scale = drop(sqrt((l) %*% Cov_E %*% (l)))
      kappa = 1. / scale**2
      alpha_E = drop(kappa) * active_signs *  drop(l %*% Cov_E)
      A = rbind(-t(alpha_E),
                diag(s_obs))
      con = list(linear_part = A, offset = theta_E, mean = rep(0, s_obs + 1))
      cov = matrix(0, s_obs + 1, s_obs + 1)
      cov[1, 1] = scale**2 * sigma**2
      cov[2:(s_obs + 1), 2:(s_obs + 1)] =
        Cov_E * gamma**2 * outer(active_signs, active_signs)
      con$covariance = cov
      initial = rep(0, s_obs + 1)
      initial[1] = l %*% beta_E
      initial[2:(s_obs + 1)] = -X_Ei %*% (y_star - y) * active_signs
      eta = rep(0, s_obs + 1)
      eta[1] = 1.
      observed = sum(initial * eta)
      ############################################################################
      empty.out <- c(observed = NA, pval = NA, ci = c(NA, NA))
      out = R.utils::evalWithTimeout(gibbs_test(con,
                                                initial,
                                                eta,
                                                alpha = alpha),
                                     timeout = 10, onTimeout = NULL)
      if (!is.null(out)){
        ci <- drop(eta %*% con$covariance %*% eta) * out$ci
        out.contrast[j, ] <- c(out$observed, out$pval, ci)
      } else {
        out.contrast[j, ] <- c(NA, NA, NA, NA)
      }
    }
  } else {
    out.contrast <- data.frame(observed = rep(NA, ncontrast),
                               pval = rep(NA, ncontrast),
                               ll = rep(NA, ncontrast),
                               ul = rep(NA, ncontrast))
  }
  list(out.all = out.all,
       y_star = y_star,
       gamma = gamma,
       sigma = sigma,
       out.contrast = out.contrast,
       coef = coef,
       Lasso.output = Lasso.output,
       Lasso.output2 = cbind(Lasso.output[which(Lasso.output$estimate!=0), ],
                             out.all),
       bestlam = lam)
}





gibbs_test <- function(con,
                       initial,
                       eta,
                       ndraw = 8000,
                       burnin = 2000,
                       alpha = 0.05){
  pySet("A", t(con$linear_part), useNumpy = T)
  pySet("theta_E", t(con$offset), useNumpy = T)
  pyExec("theta_E = theta_E.reshape(theta_E.shape[1])")
  pySet("cov", con$covariance, useNumpy = T)
  pyExec("con = constraints(A, theta_E)")
  pyExec("con.covariance[:] = cov")
  pySet("initial", t(initial), useNumpy = T)
  pyExec("initial = initial.reshape(initial.shape[1])")
  pySet("eta", t(eta), useNumpy = T)
  pyExec("eta = eta.reshape(eta.shape[1])")
  pyExec("observed = (initial * eta).sum()")
  pySet("ndraw", ndraw, useNumpy = T)
  pySet("burnin", burnin, useNumpy = T)
  pySet("alpha", alpha, useNumpy = T)
  pyExec("intervals = None
fail_count = 0

while intervals is None and fail_count < 5:
#print 'fail_count:', fail_count
  try:
    fail_count += 1
    _, _, _, family = gibbs_test(con,
    initial,
    eta,
    UMPU=False,
    sigma_known=True,
    ndraw=ndraw,
    burnin=burnin,
    how_often=5,
    tilt=con.covariance.dot(eta))
    pval = family.cdf(0, observed)
    pval = 2 * min(pval, 1 - pval)

    intervals=do_it(observed, alpha)[1]
    #intervals=family.equal_tailed_interval(observed, 1 - 0.95)

  except:
    pass")
  ci <- pyGet("np.array(intervals, dtype=float)", simplify = F)
  pyExec("A=(pval,observed)")
  pval <- pyGet("np.array(A,dtype=float)", simplify = F)
  list(observed = pval[2], pval = pval[1], ci = drop(ci))
}

start_python <- function(){
  suppressPackageStartupMessages(library(PythonInR))
  pyExec("import sys
import selection
import regreg
import numpy
import numpy as np
np.warnings.filterwarnings('ignore')
from functools import wraps
from selection.constraints.affine import gibbs_test
from selection.algorithms.lasso import standard_lasso
from selection.algorithms.lasso import constraints
from selection.algorithms.lasso import additive_noise
#from selection.tests.decorators import wait_for_return_value
import nose
import nose.tools
import time
import timeout_decorator")
pyExec("
def wait_for_return_value(max_tries=50, strict=True):
  def wait_for_decorator(test):

    @wraps(test)
    def _new_test(*args, **kwargs):
      count = 0
      while True:
        count += 1
        v = test(*args, **kwargs)
        if v is not None:
          return count, v
        if count >= max_tries:
          raise ValueError('test has not returned anything after %d tries' % max_tries)
    return nose.tools.make_decorator(test)(_new_test)

  return wait_for_decorator

@timeout_decorator.timeout(10, use_signals=False)
@wait_for_return_value(max_tries=10)
def do_it(observed, alpha = 0.05):
  intervals=family.equal_tailed_interval(observed, alpha)
  return intervals")
}
