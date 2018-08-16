#' Fit linear model to create a score using lm with no interactions
#'
#' Performs lasso in a dataset with main effects but no interactions. calculates
#' PITE and CI using the Average treatment effect (ATE)
#'
#' @param dataset A data.frame created with \code{OneData}
#' @param input A list of objects created with \code{MakeInput}
#' @param parameters A \code{MakeParameters} object
#' @param verbose logical. whether to print a message with results.
#' @param alpha	Significance level for confidence intervals (target is
#'  miscoverage alpha/2 in each tail)
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <- OneData(input,parameters)
#' score.null(dataset, input, parameters)
#'
#' @export
score.null <- function(dataset,
                       input,
                       parameters = NULL,
                       verbose = TRUE,
                       alpha = 0.05){
  # Set initial parameters -----------------------------------------------------
  n_biom<-input$n_biom  # n_biom
  N <- nrow(dataset)      # Total number of subjects
  K <- n_biom + 2         # Number of coefficients in the model
  scores = data.frame(Dx = numeric(N))
  ##---------------------------------------------------------------------------#
  ## Fit  ML -------------------------------------------------------------------
  # Creates formula to contain all main effects and interactions with treatments
  main <- paste0(" mkr", 1:n_biom, collapse = " +")
  formula <- as.formula(paste0("y ~ treatment + ", main))
  # Perform glm in the dataset
  reg.0 <- glm(formula, data = dataset, x = TRUE, y = TRUE)
  # Store estimated coefficients and sd
  ML.output <- broom::tidy(reg.0)
  ML.output$LowConfPt <- NA
  ML.output$UpConfPt  <- NA
  if(K > N) { # Catch cases when more parameters than sample size
    warning("No calculations were performed
            The number of parameters in the model is larger than
            the number of observations.")
    # Store score, variance and CIs in the dataset with NA
    scores$Dx    <- NA
    scores$Dx.v  <- NA
    scores$Dx.ll <- NA
    scores$Dx.ul <- NA
    scores$Dx.cover <- NA
    scores$Dx.w <- NA
    scores$PITE <- NA
    return(list(scores = scores,
                nvars = K,
                ML.output = ML.output,
                means = NULL,
                cov   = NULL,
                reg.0 = reg.0,
                N = N))
  }
  ML.output$LowConfPt <- ML.output$estimate - 1.96 * ML.output$std.error
  ML.output$UpConfPt  <- ML.output$estimate + 1.96 * ML.output$std.error

  ##---------------------------------------------------------------------------#
  ## Create score --------------------------------------------------------------
  s2 <- reg.0$deviance / reg.0$df.residual
  coef <- reg.0$coefficients # Get coefficients of the model
  covb <- vcov(reg.0)        # Get covariance matrix for the estimates
  # Create matrix X.dx for calculating score
  X.dx <- as.matrix(cbind(1, dataset[, paste0("mkr", 1:n_biom)]))
  colnames(X.dx) <- c("treatment", paste0("treatment:mkr", 1:n_biom))
  # Retain only the average treatment effect
  coef.dx <- coef[2]
  nvars   <- sum(coef.dx!=0)
  # Calculate score or PITE for each individual
  Dx.m <- 2 * coef.dx
  # Calculate Variance of the score
  ## Extract VAR and COV from coefficients in the score
  covb.dx <- covb[grepl("treatment", colnames(covb)),
                  grepl("treatment", colnames(covb))]
  # Calculate variance of the score
  vDx.m  <- 4 * covb.dx
  # Store score, variance and CIs in the dataset
  scores$Dx    <- drop(Dx.m)
  scores$Dx.v  <- drop(vDx.m)
  tsigma <- qt(1 - alpha / 2, df = N - K) * sqrt(vDx.m)
  scores$Dx.ll <- drop(Dx.m - tsigma)
  scores$Dx.ul <- drop(Dx.m + tsigma)

  # Calculate individual coverage ----------------------------------------------
  # Obtain the mean for all variables from the parameters input
  # For the no predictive ones, it was set to 0
  n_biom_nopred <- input$n_biom - input$n_biom_pred
  d  <- c(parameters$predictive, rep(0, n_biom_nopred))
  mu <- c(parameters$means, rep(0, n_biom_nopred))
  # Now calculate the PITE
  PITE.M <- drop(2 * (input$b_ + sum(mu * d)))
  scores$Dx.cover <- 1 * (scores$Dx.ll <= PITE.M &
                            scores$Dx.ul >= PITE.M)
  scores$PITE <- drop(PITE.M)

  # Calculate width of the CI --------------------------------------------------
  scores$Dx.w <- scores$Dx.ul - scores$Dx.ll

  if (verbose) cat("Null model: Done - ")
  # Output ---------------------------------------------------------------------
  list(scores = scores,
       nvars = nvars,
       ML.output = ML.output,
       means = coef.dx,
       cov   = covb.dx,
       reg.0 = reg.0,
       N = N)
}

#' Fit linear model to create a score using glm
#'
#' Performs glm in a dataset with main effects and interactions. Calculates
#' PITE and CI
#'
#' @param dataset A data.frame created with \code{OneData}
#' @param input A list of objects created with \code{MakeInput}
#' @param verbose logical. whether to print a message with results.
#' @param alpha	Significance level for confidence intervals (target is
#'  miscoverage alpha/2 in each tail)
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <- OneData(input,parameters)
#' score.lm(dataset, input)
#'
#' @export
score.lm <- function(dataset,
                     input,
                     verbose = TRUE,
                     alpha = 0.05) {
  # Set initial parameters -----------------------------------------------------
  n_biom <- input$n_biom  # n_biom
  N <- nrow(dataset)      # Total number of subjects
  K <- n_biom * 2 + 2         # Number of coefficients in the model
  scores <- data.frame(Dx = numeric(N))

  ## Fit glm -------------------------------------------------------------------
  # Creates formula to contain all main effects and interactions with treatments
  main <- paste0(" mkr",1:n_biom, collapse = " +")
  int  <- paste0(" treatment*mkr",1:n_biom, collapse = " +")
  formula <- as.formula(paste0("y ~ treatment + ",
                               paste(main, int, sep = " +")))
  # Perform glm in the dataset
  reg.0 <- glm(formula, data = dataset, x = TRUE, y = TRUE)
  # Store estimated coefficients and sd
  ML.output <- broom::tidy(reg.0)
  # And also store confidence intervals for coefficients
  ML.output$LowConfPt <- ML.output$estimate - qt(1-alpha/2, df=N-K) * ML.output$std.error
  ML.output$UpConfPt  <- ML.output$estimate + qt(1-alpha/2, df=N-K) * ML.output$std.error
  if(K > N) { # Catch cases when more terms than sample size
    warning("No calculations were performed
            The number of parameters in the model is larger than
            the number of observations.")
    # Store score, variance and CIs in the dataset with NA
    scores$Dx    <- NA
    scores$Dx.v  <- NA
    scores$Dx.ll <- NA
    scores$Dx.ul <- NA
    scores$Dx.cover <- NA
    scores$Dx.w  <- NA
    scores$PITE  <- 0
    return(list(scores = scores,
                nvars = K,
                ML.output = ML.output,
                means = NULL,
                cov   = NULL,
                reg.0 = reg.0,
                N = N))
  }
  ##---------------------------------------------------------------------------#
  ## Create score --------------------------------------------------------------
  s2 <- reg.0$deviance/reg.0$df.residual
  coef <- reg.0$coefficients # Get coefficients of the model
  covb <- vcov(reg.0)        # Get covariance matrix for the estimates
  # Create matrix X.dx for calculating score
  X.dx <- as.matrix(cbind(1, dataset[, paste0("mkr", 1:n_biom)]))

  colnames(X.dx) <- c("treatment", paste0("treatment:mkr", 1:n_biom))
  # Retain only those coef in the score
  coef.dx <- coef[c(2, (n_biom + 2 + 1):(2 * n_biom + 2))]
  nvars <- sum(coef.dx != 0)
  # Calculate score or PITE for each individual
  Dx.m <- 2 * X.dx %*% coef.dx
  # Calculate Variance of the score
  ## Extract VAR and COV from coefficients in the score
  covb.dx <- covb[grepl("treatment", colnames(covb)),
                  grepl("treatment", colnames(covb))]
  # Calculate variance of the score
  vDx.m   <- 4*(diag((X.dx)%*%covb.dx%*%t(X.dx)))
  # Store score, variance and CIs in the dataset
  scores$Dx    <- drop(Dx.m)
  scores$Dx.v  <- drop(vDx.m)

  tsigma <- qt(1-alpha/2, df=N-K)*sqrt(vDx.m)
  scores$Dx.ll <- drop(Dx.m - tsigma)
  scores$Dx.ul <- drop(Dx.m + tsigma)

  # Calculate individual coverage ----------------------------------------------
  scores$Dx.cover <- 1*(scores$Dx.ll <= dataset$TrueTrtEff &
                        scores$Dx.ul >= dataset$TrueTrtEff)
  # Calculate width of the CI --------------------------------------------------
  scores$Dx.w <- scores$Dx.ul - scores$Dx.ll
  scores$PITE <- dataset$TrueTrtEff

  if (verbose) cat("Linear model: Done - ")
  # Output ---------------------------------------------------------------------
  list(scores = scores,
       nvars = nvars,
       ML.output = ML.output,
       means = coef.dx,
       cov   = covb.dx,
       reg.0 = reg.0,
       N = N)
}

#' Fit linear model to create a score using LASSO
#'
#' Performs lasso in a dataset with main effects and interactions. calculates
#' PITE and CI
#'
#' @param dataset A data.frame created with \code{OneData}
#' @param input A list of objects created with \code{MakeInput}
#' @param parameters A \code{MakeParameters} object
#' @param verbose logical. whether to print a message with results.
#' @param alpha	Significance level for confidence intervals (target is
#'  miscoverage alpha/2 in each tail)
#' @param lam_frac multiplier for the lagrange parameter
#' @param lambda penalization parameter. one of "lagrange", "lambda.min", "lambda.1se" or a number
#' @param nfolds number of folds for crossvalidations. to be passed to glmnet
#' @param pite.ci logical. whether to calculate the ci for the pite.
#' @param gridrange_ Grid range for constructing confidence intervals, on the standardized scale. to be passed to the fixedLassoInf
#' @param tol.beta Tolerance for determining if a coefficient is zero. to be passed to the fixedLassoInf
#'
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <- OneData(input,parameters)
#' score.lasso(dataset, input)
#'
#' @export
score.lasso <- function(dataset,
                        input,
                        parameters = NULL,
                        verbose = TRUE,
                        alpha = 0.05,
                        lam_frac = 1,
                        lambda = "lambda.min",
                        nfolds = 10,
                        pite.ci = TRUE,
                        gridrange_ = 100,
                        tol.beta = 0.01){
  # Set Initial parameters -----------------------------------------------------
  selModel <- trueModel <- selModelTrue <- NA
  n_biom <- input$n_biom
  N <- n <- nrow(dataset)      # Total number of subjects
  K <- n_biom*2+2              # Number of terms in the model
  # Empty vectors for storage
  scores = scoresML = scoresSch = data.frame(Dx = numeric(N))
  S2 = NULL
  covb.s = NULL
  formula <- as.formula(paste0("y ~ treatment + ",
                               paste(paste0(" mkr",1:n_biom,
                                            collapse = " +"),
                                     paste0(" treatment*mkr", 1:n_biom,
                                            collapse = " +"), sep = " +")))
  X <- model.matrix(formula, data = dataset)
  Y <- as.matrix(dataset[, "y"])
  # Standardize predictors
  x <- scale(X[, -1], TRUE, TRUE) # We dont want intercept in there.
  center.X <- attr(x, "scaled:center")
  scale.X  <- attr(x, "scaled:scale")
  # Standardize Response
  y <- scale(Y, TRUE, TRUE)
  center.Y <- attr(y, "scaled:center")
  scale.Y  <- attr(y, "scaled:scale")

  if (!is.null(parameters)){
    n_biom_nopred <- input$n_biom - input$n_biom_pred
    prognosis  <- c(parameters$prognosis, rep(0,n_biom_nopred))
    predictive <- c(parameters$predictive, rep(0,n_biom_nopred))
    coef.true  <- c(input$a_, input$b_, prognosis, predictive)
    names(coef.true) <- colnames(X)
    trueModel  <- coef.true!=0
  }

  # Fit the Lasso --------------------------------------------------------------
  if (lambda == "lagrange"){
    standard_normal <- matrix(rnorm(n*5000, 0, 1),nrow = n, ncol = 5000)
    XT.e <- abs(t(x) %*% standard_normal)
    lam = lam_frac * mean(apply(XT.e, 2, max))
    bestlam <- lam / N
  }
  if (lambda == "universal"){
    bestlam = sqrt(2*log(ncol(x))/N)
    lam = bestlam * N
  }
  if (lambda %in% c("lambda.1se", "lambda.min")){
    cv.out <- glmnet::cv.glmnet(x = x, y = y, standardize = FALSE,
                                intercept = FALSE,
                                nfolds = nfolds)
    bestlam <- cv.out[[lambda]]
  }
  if (is.numeric(lambda)){
    bestlam <- lambda
  }
  cv.out <- glmnet::glmnet(x=x, y=y, standardize = FALSE, intercept = FALSE)

  # Store the estimated coefficients using the calculated lambda
  coef <- coef(cv.out, s = bestlam, exact = TRUE, x = x, y = y)[, 1]
  coef[-1][which(abs(coef[-1]) <= tol.beta/sqrt(colSums(x ^ 2)))] <- 0
  coef.noint <- coef[-1]
  coef.unsc <- coef.noint * scale.Y / scale.X

  coef.1se <- coef(cv.out, s = cv.out$lambda.1se, exact = TRUE)[,1]
  coef.min <- coef(cv.out, s = cv.out$lambda.min, exact = TRUE)[,1]
  selModel <- coef!=0
  selected     <- selModel[c(2, (n_biom+2+1):(2*n_biom+2))]
  selected.all <- selModel

  nvar.lam <- sum(coef!=0)
  nvar.min <- sum(coef.min!=0)
  nvar.1se <- sum(coef.1se!=0)

  if (verbose == TRUE) cat("\n lambda:", bestlam,
                           ", nvar:", sum(coef != 0),
                           " - ", round(coef, 2), " - \n")

  # *Store results from regression ---------------------------------------------
  Lasso.output <- data.frame(term = names(coef.unsc),
                             estimate = round(coef.unsc,6),
                             statistic = NA,
                             p.value = NA,
                             LowConfPt = NA,
                             UpConfPt = NA, row.names = NULL)
  ML.output.M <- data.frame(term = names(coef.unsc),
                             estimate = round(coef.unsc,6),
                              ste.err =NA,
                             statistic = NA,
                             p.value = NA,
                             LowConfPt = NA,
                             UpConfPt = NA, row.names = NULL)

  # **RETURN** Handle cases when more variables than observations --------------
  if(K > N/2){
    if(nvar.lam > N){
      if(nvar.1se < N){
        warning("Model is too big for N. Using lambda.1se instead")
        coef    <- coef.1se
        bestlam <- cv.out$lambda.1se
      } else {
        warning("Model is too big for N.
                More parameters than subjecs even with lambda.1se.
                Try increasing lambda")
        # Store score, variance and CIs in the dataset
        scores$Dx    <- NA
        scores$Dx.v  <- NA
        scores$Dx.ll <- NA
        scores$Dx.ul <- NA
        scores$Dx.cover <- NA
        scores$Dx.w <- NA
        scores$PITE <- NA
        scoresML <- scoresSch <- scores
        dataset$DxLasso.w <- 0

        return( list(Lasso.output = Lasso.output,
                     ML.output.M = ML.output.M,
                     nvars=0,
                     scores = scores,
                     scoresML = scoresML,
                     scoresSch = scoresSch,
                     x = x,
                     y = y,
                     bestlam = bestlam,
                     sigmahat = sigmahat,
                     coef = coef,
                     coef.noint = coef.noint,
                     coef.unsc = coef.unsc,
                     coef.M = coef.M,
                     covb.M.dx = covb.M.dx,
                     cv.out = cv.out,
                     S2 = S2,
                     selected.all = selected.all,
                     N = N,
                     trueSel = all.equal(selModel, trueModel),
                     selModel = selModel,
                     null.model = TRUE))
      }
    }
  }

  # *Create estimated PITE ------------------------------------------------------
  # Retain only those coef in the score
  coef.dx <- coef.unsc[c(1, (n_biom+1+1):(2*n_biom+1))]
  nvars <- sum(coef.dx!=0)
  # Create matrix for calculating score
  X.dx <- as.matrix(cbind(1, dataset[, paste0("mkr", 1:n_biom)]))
  colnames(X.dx) <- c("treatment", paste0("treatment:mkr", 1:n_biom))
  Dx.m <- 2*X.dx%*%coef.dx # Create score or PITE
  # Store score, variance and CIs in the dataset
  scores$Dx <- drop(Dx.m)
  if (!is.null(parameters)){
    coef.trueDx <- coef.true[c(2, (n_biom+2+1):(2*n_biom+2))]  #Retain only those coef in the score
  }

  # X1 is the contrast matrix to create scores
  # The contrast matrix should have colnames as the names of the variables in the coef
  X1 <- cbind(0, 1, matrix(0, N, n_biom), X[, 3:(n_biom + 2)])
  colnames(X1) <- names(coef)
  sigmahat = NULL
  vars = which(abs(coef.noint) > tol.beta / sqrt(colSums(x^2)))
  coef.noint_check <- coef.noint
  coef.noint_check[-vars] <- 0
  coef.dx_check <- coef.noint_check[c(1, (n_biom+1+1):(2*n_biom+1))]

  # **RETURN**  Handle cases when zero model is selected -----------------------
  if(all(abs(coef.dx_check) < tol.beta / sqrt(colSums(x[, c(1, (n_biom + 1 + 1):(2 * n_biom + 1))]^2))) | length(vars) == 0){
    warning("All estimated coefficients are 0.
            CIs are set to (0,0) since no model to condition for.")
    # Store score, variance and CIs in the dataset
    scores$Dx    <- 0
    scores$Dx.v  <- 0
    scores$Dx.ll <- 0
    scores$Dx.ul <- 0
    scores$Dx.cover <- 1
    scores$Dx.w <- 0
    scores$PITE <- 0
    scoresML <- scoresSch <- scores
    return( list(Lasso.output = Lasso.output,
                 ML.output.M = ML.output.M,
                 nvars = 0,
                 scores = scores,
                 scoresML = scoresML,
                 scoresSch = scoresSch,
                 x = x,
                 y = y,
                 bestlam = bestlam,
                 sigmahat = sigmahat,
                 coef = coef,
                 coef.noint = coef.noint,
                 coef.unsc = coef.unsc,
                 coef.M = 0,
                 covb.M.dx = 0,
                 cv.out = cv.out,
                 S2 = NULL,
                 covb.s = NULL,
                 selected.all = selected.all,
                 N = N,
                 trueSel = all.equal(selModel, trueModel),
                 selModel = selModel,
                 null.model = TRUE))
  }

  # *Perform Selective Inference -----------------------------------------------
  if(K > N/2){
    sigmahat = Inf
    while (!is.finite(sigmahat)){
      sigmahat = selectiveInference::estimateSigma(x, y,
                                                   standardize=FALSE)$sigmahat
    }
  }

  out.ci  <- fixedLassoInf_eta(x, y,
                               beta = coef.noint,
                               lambda = bestlam*N,
                               sigma=NULL,
                               alpha=alpha,
                               tol.beta = tol.beta,
                               gridrange = c(-gridrange_, gridrange_),
                               contrast = t(t(2*X1)*scale.Y/c(1, scale.X)))
  coef.ci <- selectiveInference::fixedLassoInf(x, y,
                                               beta = coef.noint,
                                               lambda = bestlam*N,
                                               sigma=NULL,
                                               tol.beta = tol.beta,
                                               gridrange = c(-gridrange_, gridrange_),
                                               alpha=alpha)
  Lasso.output[coef.ci$vars, "statistic"] <- coef.ci$coef0 / coef.ci$sd
  Lasso.output[coef.ci$vars, "p.value"]   <- coef.ci$pv
  ci <- coef.ci$ci
  Lasso.output[coef.ci$vars, "LowConfPt"] <- ci[, 1] * scale.Y / scale.X[coef.ci$vars]
  Lasso.output[coef.ci$vars, "UpConfPt"]  <- ci[, 2] * scale.Y / scale.X[coef.ci$vars]

  ci.c <- out.ci$ci
  scores$Dx.ll <- ci.c[, 1]
  scores$Dx.ul <- ci.c[, 2]
  scores$Dx.v <- out.ci$var.
  scores$Dx.w <- scores$Dx.ul - scores$Dx.ll


  ### Reduced Model ------------------------------------------------------------
  ### Perform reduced model and estimate CI from there using PoSI Framework
  ###--------------------------------------------------------------------------#

  #############################################################################-
  ## Calculate the coefficients in a reduced model with only
  ## those variables in the lasso model
  coef.M.l <- coef.noint[coef.noint != 0]
  formula.M <- as.formula(paste0("y ~ ", paste0(names(coef.M.l), collapse = " + ")))
  reg.0 <- glm(formula = formula.M,
               data = dataset, x = TRUE)

  ML.output.M <- broom::tidy(reg.0)
  K.red = length(coef.M.l)
  ML.output.M$LowConfPt <- ML.output.M$estimate - qt(1 - alpha / 2, df = reg.0$df.residual) * ML.output.M$std.error
  ML.output.M$UpConfPt  <- ML.output.M$estimate + qt(1 - alpha / 2, df = reg.0$df.residual) * ML.output.M$std.error
  coef.M <- reg.0$coefficients # Store coefficients of the model
  names(coef.M) <- c("(Intercept)",names(coef[coef!=0]))

  covb.M <- vcov(reg.0)        # Store covariance matrix for the estimates
  coef.M.dx <- coef.M[grepl("treatment", names(coef.M))]
  selVars.Dx <- names(coef.M[grepl("treatment", names(coef.M))])
  X.M.dx <- X.dx[, selVars.Dx]
  Dx.M.m <- 2 * as.matrix(X.M.dx) %*% coef.M.dx #Create score or PITE

  # Calculate Variance of the score
  # Extract VAR and COV from coefficients in the score
  covb.M.dx <- covb.M[grepl("treatment", colnames(covb.M)),
                      grepl("treatment", colnames(covb.M))]
  # Calculate variance of the score
  vDx.M.m   <- (4 * diag(as.matrix(X.M.dx) %*% covb.M.dx %*% t(X.M.dx)))
  t.df <- reg.0$df.residual
  # Store score, variance and CIs in the dataset
  scoresML$Dx    <- drop(Dx.M.m)
  scoresML$Dx.v  <- drop(vDx.M.m)
  scoresML$Dx.ll <- drop(Dx.M.m - qt(1 - alpha / 2, df = t.df) * sqrt(vDx.M.m))
  scoresML$Dx.ul <- drop(Dx.M.m + qt(1 - alpha / 2, df = t.df) * sqrt(vDx.M.m))
  scoresML$Dx.w  <- scoresML$Dx.ul - scoresML$Dx.ll
  #############################################################################-
  #### *Scheffe confidence bounds ----------------------------------------------
  # We need fist the full model matrix X
  # Get the coefficients for the scheffe conf interval

  if(K > N){
    scoresSch$Dx    <- drop(Dx.M.m)
    scoresSch$Dx.ll <- NA
    scoresSch$Dx.ul <- NA
    scoresSch$Dx.w    <- NA
  } else {
    X_M <- X[, selected.all, drop=FALSE]
    X_MX <- solve(crossprod(X_M))%*%t(X_M)%*%X
    XX.inv <- solve(t(X_M)%*%X_M)
    reg.1 <- glm(formula,
                 data = dataset)
    S2 <- (reg.1$deviance/reg.1$df.residual)
    covb.s <- XX.inv*S2       # Store covariance matrix for the estimates
    covb.s.dx <- covb.s[grepl( "treatment", colnames( covb.s ) ),
                        grepl( "treatment", colnames( covb.s ) )]
    # Calculate variance of the score, using S2 from full model
    # But X'X from reduced model
    vDx.M.s   <- (4*diag(as.matrix(X.M.dx)%*%covb.s.dx%*%t(X.M.dx)))
    N. <- N
    p. <- input$n_biom # Number of biomarkers
    d. <- 2*p. + 2     # Number of variables in design matrix
    r. <- N. - d.
    F. <- qf(p = 1-alpha, df1 = d., df2 = r.)
    Ksch. <- sqrt(d.*F.)
    scoresSch$Dx    <- drop(Dx.M.m)
    scoresSch$Dx.ll <- drop(Dx.M.m - Ksch. * sqrt(vDx.M.s))
    scoresSch$Dx.ul <- drop(Dx.M.m + Ksch. * sqrt(vDx.M.s))
    scoresSch$Dx.w    <- scoresSch$Dx.ul - scoresSch$Dx.ll
  }

  ###--------------------------------------------------------------------------#
  ### Perform Expectation in Reduced Model ----
  ###--------------------------------------------------------------------------#
  if (!is.null(parameters)){
    # Flag for selected variables in the model
    selVars.Dx <- selected.all[grepl("treatment",names(selected.all))]
    selVars <- selVars.Dx[-1] # Excluding the treatment coefficient
    if (selVars[1] == FALSE & length(selVars) == 1) {
      SigmaS   <- 1
      SigmaScS <- 1
      muScS <- 0
      muScd <- 0
      XSd <- 0
    } else {
      # First calculate the Term for selected biomarkers
      S <- diag(1 * selVars)
      colnames(S) <- rownames(S) <- names(selVars)
      X.M.dx  <- X.dx[, -1, drop = F][,  selVars, drop=F] # matrix with selected variables
      X.Mc.dx <- X.dx[, -1, drop = F][, !selVars, drop=F] # matrix with no selected variables

      X.dx <- as.matrix(cbind(1, dataset[, paste0("mkr", 1:n_biom)]))  # head(X.dx)
      X.dx.sel <- X.dx[, selected] # head(X.dx.sel)

      # If they are not independent, we need to use the following.
      # contribution from the selected biomarkers
      XSd <- X.dx[, -1] %*% S %*% coef.trueDx[-1]

      # Obtain the mean for all variables from the parameters input
      # For the no predictive ones, it was set to 0
      n_biom_nopred <- input$n_biom - input$n_biom_pred
      mu <- c(parameters$means, rep(0,n_biom_nopred))
      names(mu) <- paste0("mkr", 1:n_biom)
      # Partition mean vector inot selected S and not selected Sc
      muSc <- mu[!selVars]
      muS  <- mu[ selVars]
      if (sum(!selVars) > 0){ # If we dont select some variables
        if (sum(!selVars) == n_biom){ # But Only intercept was selected
          XSd   <- 0
          muScd <- 0
        } else { # or at least one variable is selected
          covm <- BiomCorrMatrix(input, parameters)$covm # get the cov matrix
          colnames(covm) <- rownames(covm) <- names(selVars)
          # Get the partitions on Sigma needed for calculating conditional
          # expectation
          SigmaS   <- covm[selVars,selVars]
          SigmaScS <- covm[!selVars,selVars]
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
    scores$Dx.cover    <- 1*(scores$Dx.ll <= PITE.M &
                             scores$Dx.ul >= PITE.M)
    scoresML$Dx.cover  <- 1*(scoresML$Dx.ll <= PITE.M &
                             scoresML$Dx.ul >= PITE.M)
    scoresSch$Dx.cover <- 1*(scoresSch$Dx.ll <= PITE.M &
                             scoresSch$Dx.ul >= PITE.M)
    selModelTrue <- all.equal(selModel, trueModel)
  }

  # Output ---------------------------------------------------------------------
  if (verbose) cat("Lasso: Done - ")
  list(Lasso.output = Lasso.output,
       ML.output.M = ML.output.M,
       reg.0 = reg.0,
       nvars = nvars,
       scores = scores,
       scoresML = scoresML,
       scoresSch = scoresSch,
       x = x,
       y = y,
       N = N,
       bestlam = bestlam,
       sigmahat = sigmahat,
       coef = coef,
       coef.noint = coef.noint,
       coef.unsc = coef.unsc,
       coef.M = coef.M,
       covb.M.dx = covb.M.dx,
       cv.out = cv.out,
       S2 = S2,
       covb.s = covb.s,
       selected.all = selected.all,
       trueSel = selModelTrue,
       selModel = selModel,
       null.model = FALSE)
}
