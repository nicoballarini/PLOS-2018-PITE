#' Confidence Intervals for the full model
#'
#' Calculates the confidence intervals for the score given the full model
#'
#' @param dataset A data.frame created with \code{OneData}
#' @param ML.results The results obtained from score.lm
#' @param alpha	Significance level for confidence intervals (target is miscoverage alpha/2 in each tail)
#' @param input A list of objects created with \code{MakeInput}
#' @param parameters A \code{MakeParameters} object
#'
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <-OneData(input,parameters)
#' score.lm(dataset, input)
#'
#' @export
confidence.intervals.ml.test <- function(dataset, input, parameters,
                                         ML.results, alpha = 0.05){
  n_biom = input$n_biom
  K <- n_biom * 2 + 2         # Number of coefficients in the model

  scores <- data.frame(Dx=numeric(nrow(dataset)))
  if(ML.results$nvars > ML.results$N) { # Catch cases when ML is not possible because p>n
    warning("confidence.intervals.ml.test: No calculations were performed
            The number of parameters in the model is larger than
            the number of observations.")
    # Store score, variance and CIs in the dataset with NA
    scores$Dx    <- NA
    scores$Dx.v  <- NA
    scores$Dx.ll <- NA
    scores$Dx.ul <- NA
    scores$Dx.cover <- NA
    scores$Dx.w <- NA
    scores$PITE <- 0
    return(list(scores = scores))
  }

  # Calculate PITE in the test dataset
  X.dx <- as.matrix(cbind(1, dataset[, paste0("mkr", 1:n_biom)]))
  colnames(X.dx) <- c("treatment", paste0("treatment:mkr", 1:n_biom))
  # Retain only those coef in the score
  coef.dx <- ML.results$means
  # Calculate score or PITE for each individual
  Dx.m <- 2 * X.dx %*% coef.dx
  # Calculate Variance of the score
  ## Extract VAR and COV from coefficients in the score
  covb.dx <- ML.results$cov
  # Calculate variance of the scorev
  reg.0 <- ML.results$reg.0
  s2 <- reg.0$deviance/reg.0$df.residual
  vDx.m   <- 4 * (diag((X.dx) %*% covb.dx %*% t(X.dx)))
  vDx.m.pred <- s2 + vDx.m

  # Store score, variance and CIs in the dataset
  scores$Dx    <- drop(Dx.m)

  scores$Dx.v  <- drop(vDx.m)
  N<-nrow(dataset)      # Total number of subjects
  t.df <- reg.0$df.residual
  tsigma <- qt(1-alpha/2, df = t.df) * sqrt(vDx.m)
  scores$Dx.ll <- drop(Dx.m - tsigma)
  scores$Dx.ul <- drop(Dx.m + tsigma)

  # Calculate individual coverage ----------------------------------------------
  scores$Dx.cover <- 1*(scores$Dx.ll <= dataset$TrueTrtEff &
                          scores$Dx.ul >= dataset$TrueTrtEff)
  # Calculate width of the CI --------------------------------------------------
  scores$Dx.w <- scores$Dx.ul - scores$Dx.ll
  # if (verbose) cat("Linear model: Done - ")
  scores$PITE <- dataset$TrueTrtEff
  # Output ---------------------------------------------------------------------
  list(scores = scores)
}

#' Confidence Intervals for the model without interactions (ATE)
#'
#' Calculates the confidence intervals for the score given a maximum
#'  likelihood estimation
#'
#' @param dataset A data.frame created with \code{OneData}
#' @param ML.results The results obtained from score.lm
#' @param alpha	Significance level for confidence intervals (target is miscoverage alpha/2 in each tail)
#' @param input A list of objects created with \code{MakeInput}
#' @param parameters A \code{MakeParameters} object
#'
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <-OneData(input,parameters)
#' score.lm(dataset, input)
#'
#' @export
confidence.intervals.null.test <- function(dataset, input, parameters,
                                           ML.results, alpha=0.05) {
  n_biom = input$n_biom
  K <- n_biom + 2         # Number of coefficients in the model

  n_biom = input$n_biom
  N <- ML.results$N # Sample size in training data
  scores <- data.frame(Dx = numeric(nrow(dataset)))

  if(ML.results$nvars > ML.results$N) { # Catch cases when ML is not possible because p>n
    warning("confidence.intervals.null.test: No calculations were performed
            The number of parameters in the model is larger than
            the number of observations.")
    # Store score, variance and CIs in the dataset with NA
    scores$Dx    <- NA
    scores$Dx.v  <- NA
    scores$Dx.ll <- NA
    scores$Dx.ul <- NA
    scores$Dx.cover <- NA
    scores$Dx.w <- NA
    scores$PITE <- 0
    return(list(scores = scores))
  }

  # Calculate PITE in the test dataset
  # Retain only those coef in the score
  coef.dx <- ML.results$means
  # Calculate score or PITE for each individual
  Dx.m <- 2 * coef.dx
  # Calculate Variance of the score
  ## Extract VAR and COV from coefficients in the score
  covb.dx <- ML.results$cov
  # Calculate variance of the scorev
  reg.0 <- ML.results$reg.0
  s2 <- reg.0$deviance / reg.0$df.residual
  vDx.m   <- 4 * (covb.dx)
  vDx.m.pred <- s2 + vDx.m

  # Store score, variance and CIs in the dataset
  scores$Dx    <- drop(Dx.m)
  scores$Dx.v  <- drop(vDx.m)

  N<-nrow(dataset)      # Total number of subjects
  t.df <- reg.0$df.residual
  tsigma <- qt(1 - alpha / 2, df = t.df) * sqrt(vDx.m)
  scores$Dx.ll <- drop(Dx.m - tsigma)
  scores$Dx.ul <- drop(Dx.m + tsigma)

  # Calculate individual coverage ----------------------------------------------
  n_biom_nopred <- input$n_biom - input$n_biom_pred
  d <- c(parameters$predictive,rep(0, n_biom_nopred))
  mu <- c(parameters$means,rep(0, n_biom_nopred))
  # Now calculate the PITE
  PITE.M <- drop(2 * (input$b_ + sum(mu * d)))
  scores$Dx.cover <- 1 * (scores$Dx.ll <= PITE.M &
                          scores$Dx.ul >= PITE.M)
  scores$PITE <- drop(PITE.M)

  # Calculate width of the CI --------------------------------------------------
  scores$Dx.w <- scores$Dx.ul - scores$Dx.ll
  # Output ---------------------------------------------------------------------
  list(scores = scores)
}

#' Confidence Intervals for Lasso
#'
#' Calculates the confidence intervals for the score given a Lasso estimation
#'
#' @param dataset A data.frame created with \code{OneData}
#' @param lasso.results The results obtained from score.lasso
#' @param alpha	Significance level for confidence intervals (target is miscoverage alpha/2 in each tail)
#' @param input A list of objects created with \code{MakeInput}
#' @param parameters A \code{MakeParameters} object
#' @param gridrange_ Grid range for constructing confidence intervals, on the standardized scale. to be passed to the fixedLassoInf
#' @param tol.beta Tolerance for determining if a coefficient is zero. to be passed to the fixedLassoInf
#'
#'
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <-OneData(input,parameters)
#' score.lasso(dataset, input)
#'
#' @export
confidence.intervals.lasso.test <- function(dataset,
                                            input,
                                            parameters,
                                            lasso.results,
                                            alpha = 0.05,
                                            gridrange_ = 250,
                                            tol.beta = 0.01) {
  nbiom = input$n_biom
  scores = scoresML = scoresSch <- data.frame(Dx = numeric(nrow(dataset)))

  X <- lasso.results$x
  Y <- lasso.results$y
  # X1 is the contrast matrix to create scores (Include intercept or not!)
  # The contrast matrix should have colnames as the names of the variables in the coef
  N = nrow(dataset)
  X.test = dataset[paste0("mkr", 1:nbiom)]
  X1 <- cbind(0, 1, matrix(0, N, n_biom), X.test)

  colnames(X1) <- names(lasso.results$coef)
  sigmahat = NULL
  scale.X  <- attr(lasso.results$x, "scaled:scale")
  # Standardize Response
  scale.Y  <- attr(lasso.results$y, "scaled:scale")

  if (!is.null(parameters)){
    n_biom_nopred <- input$n_biom - input$n_biom_pred
    prognosis <- c(parameters$prognosis, rep(0, n_biom_nopred))
    predictive <- c(parameters$predictive, rep(0, n_biom_nopred))
    coef.true <- c(input$a_, input$b_, prognosis, predictive)
    names(coef.true) <- c("(Intercept)", colnames(X))
    trueModel <- coef.true != 0
    coef.trueDx <- coef.true[c(2, (n_biom + 2 + 1):(2 * n_biom + 2))]  #Retain only those coef in the score
  }

  # *Create estimated PITE ------------------------------------------------------
  # Retain only those coef in the score
  coef.dx <- lasso.results$coef.unsc[c(1, (n_biom+ 1 + 1):(2 * n_biom + 1))]
  # Create matrix for calculating score
  X.dx <- as.matrix(cbind(1, dataset[, paste0("mkr", 1:n_biom)]))
  colnames(X.dx) <- c("treatment", paste0("treatment:mkr", 1:n_biom))
  Dx.m <- 2 * X.dx %*% coef.dx #Create score or PITE
  # Store score, variance and CIs in the dataset

  scores$Dx <- drop(Dx.m)
  # **RETURN**  Handle cases when zero model is selected -------------------------------
  if(lasso.results$null.model == TRUE) {
    warning("All estimated coefficients are 0. CIs are set to (0,0) since no
            model to condition for.")
    # Store score, variance and CIs in the dataset
    scores$Dx    <- 0
    scores$Dx.v  <- 0
    scores$Dx.ll <- 0
    scores$Dx.ul <- 0
    scores$Dx.cover <- 1
    scores$Dx.w <- 0
    scores$PITE <- 0
    scoresML <- scoresSch <- scores
    return(list(scores = scores,
                scoresSch = scoresSch,
                scoresML = scoresML,
                tailarea = matrix(c(NA, NA, NA, NA), nrow = 1)))
  }

  # *Perform Selective Inference -----------------------------------------------
  out.ci <- fixedLassoInf_eta(x = lasso.results$x,
                              y = lasso.results$y,
                              beta = lasso.results$coef.noint,
                              lambda = lasso.results$bestlam * lasso.results$N,
                              sigma = NULL,
                              alpha = alpha,
                              tol.beta = tol.beta,
                              gridrange = c(-gridrange_, gridrange_),
                              contrast = t(t(2 * X1) * scale.Y / c(1, scale.X)))
  scores$Dx.ll <- out.ci$ci[, 1]
  scores$Dx.ul <- out.ci$ci[, 2]
  # Calculate individual coverage ----------------------------------------------
  # Calculate width of the CI --------------------------------------------------
  scores$Dx.w <- scores$Dx.ul - scores$Dx.ll

  ### Reduced Model -----------------------------------------------------------#
  ### Perform reduced model and estimate CI from there using PoSI Framework
  ###--------------------------------------------------------------------------#
  coef.M <- lasso.results$coef.M
  coef.M.dx <- coef.M[grepl("treatment", names(coef.M))]
  selVars.Dx <- names(coef.M[grepl("treatment", names(coef.M))])
  X.M.dx <- X.dx[, selVars.Dx, drop = FALSE]
  Dx.M.m <- 2 * as.matrix(X.M.dx) %*% coef.M.dx #Create score or PITE
  # Calculate Variance of the score
  # Extract VAR and COV from coefficients in the score
  covb.M.dx <- lasso.results$covb.M.dx
  # Calculate variance of the score
  vDx.M.m   <- (4 * diag(as.matrix(X.M.dx) %*% covb.M.dx %*% t(X.M.dx)))
  t.df <- lasso.results$reg.0$df.residual
  # Store score, variance and CIs in the dataset
  scoresML$Dx    <- drop(Dx.M.m)
  scoresML$Dx.v  <- drop(vDx.M.m)
  scoresML$Dx.ll <- drop(Dx.M.m - qt(1 - alpha / 2, df = t.df) * sqrt(vDx.M.m))
  scoresML$Dx.ul <- drop(Dx.M.m + qt(1 - alpha / 2, df = t.df) * sqrt(vDx.M.m))
  scoresML$Dx.w  <- scoresML$Dx.ul - scoresML$Dx.ll

  #############################################################################-
  #### *Scheffe confidence bounds ---------------------------------------------#
  # We need fist the full model matrix X
  # Get the coefficients for the scheffe conf interval
  selected.all <- lasso.results$selected.all
  X_M <- X[, selected.all[-1], drop = FALSE]
  XX.inv <- solve(t(X_M) %*% X_M)
  S2 <- lasso.results$S2
  if (!is.null(S2)){
    covb.s <- lasso.results$covb.s
    covb.s.dx <- covb.s[grepl("treatment", colnames(covb.s)),
                        grepl("treatment", colnames(covb.s))]
    # Calculate variance of the score, using S2 from full model
    # But X'X from reduced model
    vDx.M.s   <- (4 * diag(as.matrix(X.M.dx) %*% covb.s.dx %*% t(X.M.dx)))
    N. <- lasso.results$N
    p. <- input$n_biom # Number of biomarkers
    d. <- 2 * p. + 1     # Number of variables in design matrix
    d.M <- length(coef.M) - 1 # Degrees of freedom in reduced model
    r. <- N. - d. - 1
    alpha. <- alpha
    F. <- qf(p = 1 - alpha., df1 = d., df2 = r.)
    Ksch. <- sqrt(d. * F.)
    # Centered at Reducen model estimation via ML
    scoresSch$Dx    <- drop(Dx.M.m)
    scoresSch$Dx.ll <- drop(Dx.M.m - Ksch. * sqrt(vDx.M.s))
    scoresSch$Dx.ul <- drop(Dx.M.m + Ksch. * sqrt(vDx.M.s))
    scoresSch$Dx.w <- scoresSch$Dx.ul - scoresSch$Dx.ll
  } else {
    scoresSch$Dx    <- NA
    scoresSch$Dx.ll <- NA
    scoresSch$Dx.ul <- NA
    scoresSch$Dx.w  <- NA
  }

  ###--------------------------------------------------------------------------#
  ### Perform Expectation in Reduced Model -----
  ###--------------------------------------------------------------------------#
  if (!is.null(parameters)){
    # Flag for selected variables in the model
    selVars.Dx <- selected.all[grepl("treatment", names(selected.all))]
    selVars <- selVars.Dx[-1] # Excluding the treatment coefficient
    if (selVars[1] == FALSE & length(selVars) == 1) { # 1 biomarker case and not selected
      SigmaS   <- 1
      SigmaScS <- 1
      muScS <- 0
      muScd <- 0
      XSd <- 0
    } else {
      # First calculate the Term for selected biomarkers
      S <- diag(1 * selVars)
      colnames(S) <- rownames(S) <- names(selVars)
      X.M.dx  <- X.dx[, -1, drop = FALSE][,  selVars, drop = FALSE] # matrix with selected variables
      X.Mc.dx <- X.dx[, -1, drop = FALSE][, !selVars, drop = FALSE] # matrix with no selected variables
      X_M <- X[, selected.all[-1], drop = FALSE]
      X_MX <- solve(crossprod(X_M)) %*% t(X_M) %*% X
      # If they are not independent, we need to use the following.
      # contribution from the selected biomarkers
      XSd <- X.dx[, -1] %*% S %*% coef.trueDx[-1]

      # Obtain the mean for all variables from the parameters input
      # For the no predictive ones, it was set to 0
      n_biom_nopred <- input$n_biom - input$n_biom_pred
      mu <- c(parameters$means,rep(0, n_biom_nopred))
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
    PITE.M <- drop(2 * (coef.trueDx[1] * selVars.Dx[1] + XSd + muScd))
    scores$PITE <- scoresML$PITE <- scoresSch$PITE <- drop(PITE.M)
    scores$Dx.cover <- 1*(scores$Dx.ll <= PITE.M &
                          scores$Dx.ul >= PITE.M)
    scoresML$Dx.cover  <- 1*(scoresML$Dx.ll <= PITE.M &
                             scoresML$Dx.ul >= PITE.M)
    scoresSch$Dx.cover <- 1*(scoresSch$Dx.ll <= PITE.M &
                             scoresSch$Dx.ul >= PITE.M)
  }
  # Output --------------------------------------------------------------------#
  return(list(scores = scores,
              scoresML = scoresML,
              scoresSch = scoresSch,
              tailarea = cbind(out.ci$ci, out.ci$tailarea)))
}

#' Confidence Intervals for Lasso
#'
#' Calculates the confidence intervals for the score given a Lasso estimation
#'
#' @param dataset A data.frame created with \code{OneData}
#' @param input A list of objects created with \code{MakeInput}
#' @param alpha	Significance level for confidence intervals (target is
#'  miscoverage alpha/2 in each tail)
#' @param parameters A \code{MakeParameters} object
#' @param verbose logical. whether to print a message with results.
#' @param alpha	Significance level for confidence intervals (target is
#'  miscoverage alpha/2 in each tail)
#' @param lasso.results output from the score.lasso.added.noise
#' @param ndraw  Number of samples of optimization variables to sample.
#' @param burnin How many samples of optimization variable to discard.
#'
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <-OneData(input,parameters)
#' score.lasso(dataset, input)
#'
#' @export
confidence.intervals.lassoan.test <- function(dataset,
                                              input,
                                              parameters,
                                              lasso.results,
                                              ndraw = 8000,
                                              burnin = 2000,
                                              alpha = 0.05) {
  scores = data.frame(Dx = numeric(nrow(dataset)))
  if(lasso.results$null.model == TRUE) {
    warning("All estimated coefficients are 0. CIs are set to (0,0) since no
            model to condition for.")
    # Store score, variance and CIs in the dataset
    scores$Dx    <- 0
    scores$Dx.v  <- 0
    scores$Dx.ll <- 0
    scores$Dx.ul <- 0
    scores$Dx.cover <- 1
    scores$Dx.w <- 0
    scores$PITE <- 0
    scores
    return( list(scores=scores))
  }

  nbiom = input$n_biom
  X <- lasso.results$x
  y <- lasso.results$y
  y_star <- lasso.results$y_star

  # X1 is the contrast matrix to create scores (Include intercept or not!)
  # The contrast matrix should have colnames as the names of the variables in the coef
  N = nrow(dataset)
  X.test = dataset[paste0("mkr", 1:nbiom)]
  X1 <- cbind(0, 1, matrix(0, N, n_biom), X.test)
  colnames(X1) <- c("(Intercept)", names(lasso.results$coef))

  sigmahat = NULL
  scale.X  <- attr(lasso.results$x, "scaled:scale")
  # Standardize Response
  scale.Y  <- attr(lasso.results$y, "scaled:scale")

  if (!is.null(parameters)){
    n_biom_nopred <- input$n_biom - input$n_biom_pred
    prognosis <- c(parameters$prognosis, rep(0, n_biom_nopred))
    predictive <- c(parameters$predictive, rep(0, n_biom_nopred))
    coef.true <- c(input$a_, input$b_, prognosis, predictive)
    names(coef.true) <- colnames(X)
    trueModel <- coef.true != 0
    coef.trueDx <- coef.true[c(2, (n_biom + 2 + 1):(2 * n_biom + 2))]  #Retain only those coef in the score
  }

  lam  <- lasso.results$bestlam
  coef <- lasso.results$coef

  active <- coef != 0
  # **RETURN**  Handle cases when zero model is selected -----------------------
  if (sum(active) == 0){
    # Store score, variance and CIs in the dataset
    scores$Dx    <- 0
    scores$Dx.v  <- 0
    scores$Dx.ll <- 0
    scores$Dx.ul <- 0
    scores$Dx.cover <- 1
    scores$Dx.w <- 0
    scores$PITE <- 0
    scores
    return( list(scores = scores))
  }

  selModel <- coef != 0
  selected     <- selModel[c(1, (n_biom + 1 + 1):(2 * n_biom + 1))]
  selected.all <- selModel
  nvar <- sum(coef != 0)
  coef.unsc <- coef * scale.Y / scale.X
  # *Create estimated PITE ------------------------------------------------------
  # Retain only those coef in the score
  coef.dx <- coef.unsc[c(1, (n_biom + 1 + 1):(2 * n_biom + 1))]
  # Create matrix for calculating score
  X.dx <- as.matrix(cbind(1, dataset[, paste0("mkr",1:n_biom)]))
  colnames(X.dx) <- c("treatment", paste0("treatment:mkr", 1:n_biom))

  active_signs <- sign(coef[active])
  # Form the constraint matrix on (y,y^*)
  X_E = X[, active, drop = FALSE]
  X_Ei = MASS::ginv(X_E)  ## This is equivalent to solve(t(X_E)%*%X_E)%*%t(X_E)
  Cov_E = X_Ei %*% t(X_Ei)    ## This is the same as solve(t(X_E)%*%X_E) which is (X'X)^(-1)
  W_E = Cov_E %*% (active_signs)

  beta_E = X_Ei %*% (y)
  res <- list(0)

  contrast = t(t(2 * X1) * scale.Y / c(1, scale.X))
  X1 <- contrast
  ncontrast <- (dim(contrast)[1])

  gamma = lasso.results$gamma
  sigma = lasso.results$sigma

  # Contrasts -----------------------------------------------------------------
  out.contrast <- data.frame(observed = rep(NA, ncontrast),
                             pval = rep(NA, ncontrast),
                             ll = rep(NA, ncontrast),
                             ul = rep(NA, ncontrast))
  for (j in 1:nrow(contrast)){
    s_obs = sum(active)
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
    T_A = (diag(s_obs) - C%*%At) %*% beta
    theta_E = active_signs * (T_A - lam * W_E)
    scale = sqrt((l) %*% Cov_E %*% (l))
    kappa = as.numeric(1. / scale**2)
    alpha_E = kappa * active_signs *  drop(l %*% Cov_E)

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
    ###########################################################################-
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
      out.contrast[j,] <- c(out$observed, out$pval, ci)
    } else {
      out.contrast[j,] <- c(NA,NA,NA,NA)
    }
  }
  scores$Dx    <- out.contrast$observed
  scores$Dx.ll <- out.contrast$ll
  scores$Dx.ul <- out.contrast$ul
  scores$Dx.cover <- 1 * (scores$Dx.ll <= dataset$TrueTrtEff &
                          scores$Dx.ul >= dataset$TrueTrtEff)
  scores$Dx.w <- scores$Dx.ul - scores$Dx.ll

  ###--------------------------------------------------------------------------#
  ### Perform Expectation in Reduced Model ----
  ###--------------------------------------------------------------------------#
  if (!is.null(parameters)){
    # Flag for selected variables in the model
    selVars.Dx <- selected.all[grepl("treatment", names(selected.all))]
    selVars <- selVars.Dx[-1] # Excluding the treatment coefficient
    if (selVars[1] == FALSE & length(selVars) == 1) {
      SigmaS   <- 1
      SigmaScS <- 1
      muScS <- 0
      muScd <- 0
      XSd <- 0
    } else {# First calculate the Term for selected biomarkers
      S <- diag(1 * selVars)
      colnames(S) <- rownames(S) <- names(selVars)
      X.M.dx  <- X.dx[, -1, drop = FALSE][,  selVars, drop = FALSE] # matrix with selected variables
      X.Mc.dx <- X.dx[, -1, drop = FALSE][, !selVars, drop = FALSE]# matrix with no selected variables
      X_M <- X[, selected.all, drop = FALSE]
      X_MX <- solve(crossprod(X_M)) %*% t(X_M) %*% X
      # If they are not independent, we need to use the following.
      # contribution from the selected biomarkers
      XSd <- X.dx[, -1] %*% S %*% coef.trueDx[-1]

      # Obtain the mean for all variables from the parameters input
      # For the no predictive ones, it was set to 0
      n_biom_nopred <- input$n_biom - input$n_biom_pred
      mu <- c(parameters$means,rep(0, n_biom_nopred))
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
          SigmaS   <- covm[selVars,  selVars]
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
    PITE.M.LAN <- 2 * (coef.trueDx[1] * selVars.Dx[1] + XSd + muScd)
    scores$PITE <- drop(PITE.M.LAN)
    scores$Dx.cover  <- 1 * (scores$Dx.ll <= PITE.M.LAN &
                             scores$Dx.ul >= PITE.M.LAN)
  }
  # Output ---------------------------------------------------------------------
  return(list(scores = scores))
}
