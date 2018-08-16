#' Creates one simulated dataset with survival outcome
#'
#' A function that will create a data.frame with simulated data
#'
#' @param alpha True value for intercept
#' @param beta True value for treatment effect
#' @param delta1 True value for predictive effect of biomarker 1
#' @param delta2 True value for predictive effect of biomarker 2
#' @param nsub number of subjects in dataset
#' @param sigma standard deviation for the biomarkers
#' @param maxtime maximum study time for censoring observations
#'
#'
#' @export
SurvOneData <- function(alpha,
                        beta,
                        delta1,
                        delta2,
                        nsub,
                        sigma = 0.2,
                        maxtime = 2){
  z  = c(rep(-1, nsub / 2), rep(1, nsub / 2))
  x1  = 1 * (rnorm(nsub, 0, sigma) > 0) * 2 - 1
  x2 = rnorm(nsub, 0, sigma)
  x3 = rnorm(nsub, 0, sigma)
  x4 = rnorm(nsub, 0, sigma)
  x5 = rnorm(nsub, 0, sigma)
  x6 = rnorm(nsub, 0, sigma)
  lambdat = exp(alpha +
                beta * z +
                0.17 * x1 + 0.18 * x2 + 0.11 * x3 + 0.14 * x4 + 0.21 * x5 + -0.12 * x6 +
                delta1 * z * x1 + delta2 * z * x2)
  Ex1 = 0
  Ex2 = 0
  pite.z    = 2 * (beta + delta1 * Ex1 + delta2 * Ex2)
  pite.x1   = 2 * (delta1 * x1  + delta2 * Ex2)
  pite.x2   = 2 * (delta1 * Ex1 + delta2 * x2)
  pite.x1x2 = 2 * (delta1 * x1  + delta2 * x2)
  pite.zx1  = 2 * (beta + delta1 * x1  + delta2 * Ex2)
  pite.zx2  = 2 * (beta + delta1 * Ex1 + delta2 * x2)
  pite      = 2 * (beta + delta1 * x1  + delta2 * x2)

  enter_time = runif(nsub, 0,1)
  t = rexp(nsub) / lambdat
  cens   = 1 * (t + enter_time > maxtime) # 1 if censored observation
  status = 1 - cens                       # 1 if dead
  time = apply(cbind(t, maxtime - enter_time), 1, min)
  data.frame(time, status, t, enter_time, lambdat, cens,
             z, x1, x2, x3, x4, x5, x6,
             pite.z, pite.x1, pite.x2,
             pite.x1x2, pite.zx1, pite.zx2, pite)
}


#' Creates one simulated dataset with survival outcome with only one subject
#' for test in the simulations
#'
#' A function that will create a data.frame with simulated data for one subject
#'
#' @param alpha True value for intercept
#' @param beta True value for treatment effect
#' @param delta1 True value for predictive effect of biomarker 1
#' @param delta2 True value for predictive effect of biomarker 2
#' @param sigma standard deviation for the biomarkers
#' @param maxtime maximum study time for censoring observations
#'
#' @export
SurvOneSubject <- function(alpha,
                           beta,
                           delta1,
                           delta2,
                           sigma = 0.2,
                           maxtime = 2){
  nsub = 1
  z    = 1
  x1  = 1 * (rnorm(nsub, 0, sigma) > 0) * 2 - 1
  x2 = rnorm(nsub, 0, sigma)
  x3 = rnorm(nsub, 0, sigma)
  x4 = rnorm(nsub, 0, sigma)
  x5 = rnorm(nsub, 0, sigma)
  x6 = rnorm(nsub, 0, sigma)
  lambdat = exp(alpha +
                beta * z +
                0.17 * x1 + 0.18 * x2 + 0.11 * x3 + 0.14 * x4 + 0.21 * x5 + -0.12 * x6 +
                delta1 * z * x1 + delta2 * z * x2)
  Ex1 = 0
  Ex2 = 0
  pite.z    = 2 * (beta + delta1 * Ex1 + delta2 * Ex2)
  pite.x1   = 2 * (delta1 * x1  + delta2 * Ex2)
  pite.x2   = 2 * (delta1 * Ex1 + delta2 * x2)
  pite.x1x2 = 2 * (delta1 * x1  + delta2 * x2)
  pite.zx1  = 2 * (beta + delta1 * x1  + delta2 * Ex2)
  pite.zx2  = 2 * (beta + delta1 * Ex1 + delta2 * x2)
  pite      = 2 * (beta + delta1 * x1  + delta2 * x2)

  enter_time = runif(nsub, 0, 1)
  t = rexp(nsub)/lambdat
  cens   = 1*(t + enter_time > maxtime) # 1 if censored observation
  status = 1 - cens                     # 1 if dead
  time = apply(cbind(t, maxtime - enter_time), 1, min)
  data.frame(time, status, t, enter_time, lambdat, cens,
             z, x1, x2, x3, x4, x5, x6,
             pite.z, pite.x1, pite.x2,
             pite.x1x2, pite.zx1, pite.zx2, pite)
}

#' Creates a list with nsim datasets using the SurvOneData function
#'
#' A function that will create list of data.frames with a simulated data
#'
#' @param alpha True value for intercept
#' @param beta True value for treatment effect
#' @param delta1 True value for predictive effect of biomarker 1
#' @param delta2 True value for predictive effect of biomarker 2
#' @param nsub number of subjects in dataset
#' @param sigma standard deviation for the biomarkers
#' @param maxtime maximum study time for censoring observations
#'
#' @export
SurvData <- function(alpha,
                     beta,
                     delta1,
                     delta2,
                     nsim,
                     nsub,
                     sigma = 1,
                     maxtime = 2){
  lapply(1:nsim,
         function(z) {
           SurvOneData(alpha, beta, delta1, delta2,
                                 nsub, sigma, maxtime)
           }
         )
}

#' Creates a list with nsim datasets using the SurvOneSubject function
#'
#' A function that will create list of data.frames with simulated data for one subject
#'
#' @param alpha True value for intercept
#' @param beta True value for treatment effect
#' @param delta1 True value for predictive effect of biomarker 1
#' @param delta2 True value for predictive effect of biomarker 2
#' @param nsub number of subjects in dataset
#' @param sigma standard deviation for the biomarkers
#' @param maxtime maximum study time for censoring observations
#'
#' @export
SurvSubjects <- function(alpha,
                         beta,
                         delta1,
                         delta2,
                         nsim,
                         sigma = 1,
                         maxtime = 2){
  lapply(1:nsim,
         function(z){
           SurvOneSubject(alpha, beta, delta1, delta2,
                          sigma, maxtime)})
}


#' Fits the 'null' model (ATE with no interactions)
#'
#' A function that fits the model with treatment effect and prognostic terms
#' but no interaction terms
#'
#' @param dataset a data.frame
#' @param alpha Significance level for confidence intervals (target is
#'  miscoverage alpha/2 in each tail)
#' @param calculate.ci logical. whether to calculate the confidence intervals for pite
#'
#' @export
score.cox.null = function(dataset,
                          alpha = 0.05,
                          calculate.ci = TRUE) {
  reg.results <- survival::coxph(Surv(time, status) ~
                                   z + x1 + x2 + x3 + x4 + x5 + x6,
                                 data = dataset, ties = "breslow")
  n_biom = 6
  reg.out    <- broom::tidy(reg.results, conf.level = 1 - alpha)
  coef       <- reg.results$coefficients
  cov.matrix <- reg.results$var
  all.vars   <- names(coef)
  coef.dx    <- coef[1]
  X.dx       <- 2
  X.dx.full  <- 2 * cbind(1, matrix(0, nrow(dataset), n_biom))
  alpha. <- alpha
  scores <- data.frame(Dx = numeric(nrow(dataset)))

  if (calculate.ci){
    scores$Dx   <- Dx.m  <- drop(as.matrix(X.dx) %*% coef.dx)
    scores$Dx.v <- vDx.m <- diag(as.matrix(X.dx.full) %*% cov.matrix %*% t(as.matrix(X.dx.full)))
    zsigma <- qnorm(1 - alpha. / 2) * sqrt(vDx.m)
    scores$Dx.ll <- drop(Dx.m - zsigma)
    scores$Dx.ul <- drop(Dx.m + zsigma)
    scores$pite <- dataset$pite.z
    # Store pite under full model
    scores$pite.full <- dataset$pite
    scores$cover <- (scores$Dx.ll < scores$pite) & (scores$pite < scores$Dx.ul)
    scores$width <- (scores$Dx.ul - scores$Dx.ll)
    scores$bias  <- (scores$Dx - scores$pite.full)
  } else {
    scores = NULL
  }
  list(reg.out = reg.out[, c(1, 2, 5, 6, 7)],
       scores = scores,
       vars = reg.out$term[which(reg.out$estimate != 0)],
       reg.results = reg.results)
}


#' Obtains the PITE and confidence interval for an additional subject
#'
#' A function that takes the results from score.cox.null to get the PITE
#' and the confidence intervals for one subject.
#'
#' @param dataset a data.frame
#' @param alpha Significance level for confidence intervals (target is
#'  miscoverage alpha/2 in each tail)
#' @param reg.results an object from the score.cox.null function
#'
#' @export
ci.cox.null <- function(dataset, reg.results, alpha = 0.05){
  n_biom = 6
  reg.out    <- broom::tidy(reg.results, conf.level = 1 - alpha)
  coef       <- reg.results$coefficients
  cov.matrix <- reg.results$var
  all.vars   <- names(coef)
  coef.dx    <- coef[1]
  X.dx       <- 2
  X.dx.full  <- 2 * cbind(1, matrix(0, nrow(dataset), n_biom))
  alpha. <- alpha
  scores <- data.frame(Dx = numeric(nrow(dataset)))
  scores$Dx   <- Dx.m  <- drop(as.matrix(X.dx) %*% coef.dx)
  scores$Dx.v <- vDx.m <- diag(as.matrix(X.dx.full) %*% cov.matrix %*% t(as.matrix(X.dx.full)))
  zsigma <- qnorm(1 - alpha. / 2) * sqrt(vDx.m)
  scores$Dx.ll <- drop(Dx.m - zsigma)
  scores$Dx.ul <- drop(Dx.m + zsigma)
  scores$pite <- dataset$pite.z
  # Store pite under full model
  scores$pite.full <- dataset$pite
  scores$cover <- (scores$Dx.ll < scores$pite) & (scores$pite < scores$Dx.ul)
  scores$width <- (scores$Dx.ul - scores$Dx.ll)
  scores$bias  <- (scores$Dx - scores$pite.full)
  # head(scores)
  list(reg.out = reg.out[, c(1, 2, 5, 6, 7)],
       scores = scores,
       vars = reg.out$term[which(reg.out$estimate != 0)])
}


#' Fits the full model (all treatment-covariate interactions)
#'
#' A function that fits the model with treatment effect, prognostic terms
#' and interaction terms
#' @param dataset a data.frame
#' @param alpha Significance level for confidence intervals (target is
#'  miscoverage alpha/2 in each tail)
#' @param calculate.ci logical. whether to calculate the confidence intervals for pite
#'
#' @export
score.cox <- function(dataset, alpha = 0.05, calculate.ci = TRUE){
  reg.results <- survival::coxph(Surv(time, status) ~
                                   z * (x1 + x2 + x3 + x4 + x5 + x6),
                                 data = dataset, ties = "breslow")
  n_biom = 6
  reg.out    <- broom::tidy(reg.results, conf.level = 1 - alpha)
  coef       <- reg.results$coefficients
  cov.matrix <- reg.results$var
  all.vars   <- names(coef)
  var.score  <- all.vars[2:(n_biom + 1)]
  coef.dx    <- coef[c(1, (n_biom + 2):(2 * n_biom + 1))]
  X.dx       <- 2 * cbind(1, dataset[var.score])
  X.dx.full  <- 2 * cbind(1, matrix(0, nrow(dataset), length(var.score)), dataset[var.score])
  alpha. <- alpha

  if (calculate.ci){
    scores <- data.frame(Dx = numeric(nrow(dataset)))
    scores$Dx <- Dx.m <- drop(as.matrix(X.dx) %*% coef.dx)
    scores$Dx.v <- vDx.m <- diag(as.matrix(X.dx.full) %*% cov.matrix %*% t(as.matrix(X.dx.full)))
    zsigma <- qnorm(1 - alpha. / 2) * sqrt(vDx.m)
    scores$Dx.ll <- drop(Dx.m - zsigma)
    scores$Dx.ul <- drop(Dx.m + zsigma)
    scores$pite <- dataset$pite
    # Store pite under full model
    scores$pite.full <- dataset$pite
    scores$cover <- (scores$Dx.ll < dataset$pite) & (dataset$pite < scores$Dx.ul)
    scores$width <- (scores$Dx.ul - scores$Dx.ll)
    scores$bias  <- (scores$Dx - scores$pite.full)
  } else {
    scores = NULL
  }
  list(reg.out = reg.out[, c(1, 2, 5, 6, 7)],
       scores  = scores,
       vars    = reg.out$term[which(reg.out$estimate != 0)],
       reg.results = reg.results)
}

#' Obtains the PITE and confidence interval for an additional subject
#'
#' A function that takes the results from score.cox to get the PITE
#' and the confidence intervals for one subject.
#'
#' @param dataset a data.frame
#' @param alpha Significance level for confidence intervals (target is
#'  miscoverage alpha/2 in each tail)
#' @param reg.results an object from the score.cox function
#'
#' @export
ci.cox <- function(dataset, reg.results, alpha = 0.05){
  n_biom = 6
  reg.out    <- broom::tidy(reg.results, conf.level = 1 - alpha)
  coef       <- reg.results$coefficients
  cov.matrix <- reg.results$var
  all.vars   <- names(coef)
  var.score  <- all.vars[2:(n_biom+1)]
  coef.dx    <- coef[c(1, (n_biom + 2):(2 * n_biom + 1))]
  X.dx       <- 2 * cbind(1, dataset[var.score])
  X.dx.full  <- 2 * cbind(1, matrix(0, nrow(dataset), length(var.score)), dataset[var.score])
  alpha. <- alpha

  scores <- data.frame(Dx = numeric(nrow(dataset)))
  scores$Dx   <- Dx.m  <- drop(as.matrix(X.dx) %*% coef.dx)
  scores$Dx.v <- vDx.m <- diag(as.matrix(X.dx.full) %*% cov.matrix %*% t(as.matrix(X.dx.full)))
  zsigma <- qnorm(1 - alpha. / 2) * sqrt(vDx.m)
  scores$Dx.ll <- drop(Dx.m - zsigma)
  scores$Dx.ul <- drop(Dx.m + zsigma)
  scores$pite <- dataset$pite
  # Store pite under full model
  scores$pite.full <- dataset$pite
  scores$cover <- (scores$Dx.ll < scores$pite) & (scores$pite < scores$Dx.ul)
  scores$width <- (scores$Dx.ul - scores$Dx.ll)
  scores$bias  <- (scores$Dx    - scores$pite.full)
  list(reg.out = reg.out[, c(1, 2, 5, 6, 7)],
       scores = scores,
       vars = reg.out$term[which(reg.out$estimate != 0)])
}


#' Printing function for fixedlassoInf
#'
#' @export
fixedLassoInfprint <- function(x, coef) {
  tab = data.frame(names(coef)[x$vars],
                   round(x$coef0, 3),
                   round(x$zscore0, 3),
                   round(x$pv, 3), round(x$ci, 3))
  colnames(tab) = c("Var", "Coef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
  tab = cbind(tab,round(x$tailarea, 3))
  colnames(tab)[(ncol(tab) - 1):ncol(tab)] = c("LowTailArea", "UpTailArea")
  rownames(tab) = x$vars
  tab
  tab %>%
    mutate(tmpLow = ifelse(Coef < 0, -UpConfPt,  LowConfPt),
           tmpUp  = ifelse(Coef < 0, -LowConfPt, UpConfPt),
           LowConfPt = tmpLow,
           UpConfPt  = tmpUp) %>%
    dplyr::select(-tmpLow, -tmpUp)
}

#' Fits the full model (all treatment-covariate interactions) with lasso penalty
#' for model selection
#'
#' A function that fits the model with treatment effect, prognostic terms
#' and interaction terms with glmnet
#' @param dataset a data.frame
#' @param alpha Significance level for confidence intervals (target is
#'  miscoverage alpha/2 in each tail)
#' @param calculate.ci logical. whether to calculate the confidence intervals for pite
#' @param lam_frac multiplier for the lagrange parameter
#' @param gridrange_ Grid range for constructing confidence intervals, on the standardized scale. to be passed to the fixedLassoInf
#' @param tol.beta Tolerance for determining if a coefficient is zero. to be passed to the fixedLassoInf
#'
#' @export
score.cox.lasso <- function(dataset,
                            alpha = 0.05,
                            calculate.ci = TRUE,
                            lam_frac = 0.5,
                            tol.beta = 0.1,
                            gridrange_ = 100) {
  y = Surv(dataset$time, dataset$status)
  X.m <- model.matrix( ~ -1 + z * (x1 + x2 + x3 + x4 + x5 + x6), data = dataset)
  n_biom = 6
  n = nrow(dataset)
  standard_normal <- matrix(rnorm(n * 5000, 0, 1),nrow = n, ncol = 5000)
  XT.e <- abs(t(X.m) %*% standard_normal)
  lam_frac = lam_frac
  lam = lam_frac * mean(apply(XT.e, 2, max))
  bestlam <- lam / n

  reg.results <- glmnet(X.m, y, family = "cox")
  coef.lasso  <- coef(reg.results, s = bestlam)[, 1]
  # compute fixed lambda p-values and selection intervals
  coef.lasso[which(abs(coef.lasso) <= tol.beta / sqrt(colSums(X.m^2)))] <- 0
  if (all(abs(coef.lasso) <= tol.beta/sqrt(colSums(X.m^2)))){
    scores <- data.frame(Dx = rep(0, n),
                         Dx.v = rep(0, n),
                         Dx.ll = rep(0, n),
                         Dx.ul = rep(0, n),
                         pite = rep(0, n),
                         pite.full = rep(0, n),
                         cover = rep(1, n),
                         width = rep(0, n),
                         bias = rep(0, n))
    lasso.results <- data.frame(term = names(coef.lasso),
                                estimate = coef.lasso,
                                p.value = rep(NA, length(coef.lasso)),
                                conf.low = rep(NA, length(coef.lasso)),
                                conf.high = rep(NA, length(coef.lasso)))
    return(list(reg.out = lasso.results,
                scores = scores,
                vars = NULL,
                reg.results = list(x = X.m,
                                   y = dataset$time,
                                   beta = coef.lasso,
                                   lambda = lam,
                                   status = dataset$status)))
  }
  out = selectiveInference::fixedLassoInf(x = X.m,
                                          y = dataset$time,
                                          beta = coef.lasso,
                                          lambda = lam,
                                          tol.beta = tol.beta,
                                          status = dataset$status,
                                          family = "cox",
                                          alpha = alpha)
  lasso.results1 <- data.frame(term = names(coef.lasso),
                               estimate = coef.lasso,
                               stringsAsFactors = FALSE)
  lasso.results2 <- fixedLassoInfprint(out, coef = coef.lasso)
  lasso.results2 <- lasso.results2[c("Var", "P-value", "LowConfPt", "UpConfPt")]
  names(lasso.results2) <- c("term", "p.value", "conf.low",   "conf.high")
  lasso.results2$term <- paste(lasso.results2$term)
  lasso.results <- dplyr::left_join(lasso.results1, lasso.results2, by = "term")

  vars = lasso.results$term[which(lasso.results$estimate != 0)]
  all.vars <- names(coef.lasso)
  var.score <- all.vars[2:(n_biom + 1)]
  coef.dx <- coef.lasso[c(1, (n_biom + 2):(2 * n_biom + 1))]
  if (all(abs(coef.lasso) <= tol.beta / sqrt(colSums(X.m^2))) |
      all(abs(coef.dx) <= tol.beta / sqrt(colSums(X.m[, c(1, (n_biom + 2):(2 * n_biom + 1))]^2)))){
    scores <- data.frame(Dx = rep(0, n),
                         Dx.v = rep(0, n),
                         Dx.ll = rep(0, n),
                         Dx.ul = rep(0, n),
                         pite = rep(0, n),
                         pite.full = rep(0, n),
                         cover = rep(1, n),
                         width = rep(0, n),
                         bias = rep(0, n))
    return(list(reg.out = lasso.results,
                scores = scores,
                vars = NULL,
                reg.results = list(x = X.m,
                                   y = dataset$time,
                                   beta = coef.lasso,
                                   lambda = lam,
                                   status = dataset$status)))
  }
  if (calculate.ci){
    X.dx <- cbind(1, dataset[var.score])
    X.dx.full <- 2 * cbind(1, matrix(0, nrow(dataset), length(var.score)), dataset[var.score])
    out.score <- PITE::fixedCoxLassoInf_eta(x = X.m,
                                                     y = dataset$time,
                                                     beta = coef.lasso,
                                                     lambda = lam,
                                                     status = dataset$status,
                                                     alpha = alpha,
                                                     tol.beta = tol.beta,
                                                     gridrange = c(-gridrange_, gridrange_),
                                                     contrast = X.dx.full)
    scores <- data.frame(Dx = as.matrix(X.dx.full) %*% coef.lasso)
    scores$Dx.v <- out.score$sd
    scores$Dx.ll <- out.score$ci[, 1]
    scores$Dx.ul <- out.score$ci[, 2]

    vars.score = vars[grep(pattern = "z", vars)]

    if (all(c("z:x1", "z:x2", "z") %in% vars.score)) {
      scores$pite <- dataset$pite
    } else   if (all(c("z:x1", "z:x2") %in% vars.score)) {
      scores$pite <- dataset$pite.x1x2
    } else   if (all(c("z:x1", "z") %in% vars.score)) {
      scores$pite <- dataset$pite.zx1
    } else   if (all(c("z:x2", "z") %in% vars.score)) {
      scores$pite <- dataset$pite.zx2
    } else   if (all(c("z") %in% vars.score)) {
      scores$pite <- dataset$pite.z
    } else   if (all(c("z:x2") %in% vars.score)) {
      scores$pite <- dataset$pite.x2
    } else   if (all(c("z:x1") %in% vars.score)) {
      scores$pite <- dataset$pite.x1
    } else  {
      scores$pite <- 0
    }
    # Store pite under full model
    scores$pite.full <- dataset$pite
    scores$cover <- (scores$Dx.ll < scores$pite) & (scores$pite < scores$Dx.ul)
    scores$width <- (scores$Dx.ul - scores$Dx.ll)
    scores$bias  <- (scores$Dx - scores$pite.full)
    head(scores)
  } else {
    scores = NULL
  }
  list(reg.out = lasso.results,
       scores = scores,
       vars = vars,
       reg.results = list(x = X.m,
                          y = dataset$time,
                          beta = coef.lasso,
                          lambda = lam,
                          status = dataset$status))
}

#' Obtains the PITE and confidence interval for an additional subject
#' using the Lasso output
#'
#' A function that takes the results from score.cox.lasso to get the PITE
#' and the confidence intervals for one subject.
#'
#' @param dataset a data.frame
#' @param alpha Significance level for confidence intervals (target is
#'  miscoverage alpha/2 in each tail)
#' @param reg.results an object from the score.cox.lasso function
#' @param gridrange_ Grid range for constructing confidence intervals, on the standardized scale. to be passed to the fixedLassoInf
#' @param tol.beta Tolerance for determining if a coefficient is zero. to be passed to the fixedLassoInf
#'
#' @export
ci.cox.lasso <- function(dataset,
                         alpha = 0.05,
                         reg.results,
                         reg.out,
                         tol.beta = 0.1,
                         gridrange_ = 100) {
  n_biom = 6
  n = 1
  n.train = nrow(reg.results$x)
  all.vars  <- names(reg.results$beta)
  var.score <- all.vars[2:(n_biom+1)]
  coef.lasso = reg.results$beta
  coef.dx <- coef.lasso[c(1,(n_biom+2):(2*n_biom+1))]
  if (all(abs(coef.lasso) <= tol.beta/sqrt(n.train)) | all(abs(coef.dx) <= tol.beta/sqrt(n.train))){
    scores <- data.frame(Dx = rep(0,n),
                         Dx.v = rep(0,n),
                         Dx.ll = rep(0,n),
                         Dx.ul = rep(0,n),
                         pite = rep(0,n),
                         pite.full = dataset$pite,
                         cover = rep(1,n),
                         width = rep(0,n),
                         bias = rep(0,n) - dataset$pite)
    return(list(reg.out = reg.out,
                scores = scores,
                vars = NULL))
  }

  X.dx <- cbind(1,dataset[var.score])
  X.dx.full <- 2*cbind(1,matrix(0, nrow(dataset), length(var.score)), dataset[var.score])
  head(X.dx.full)
  out.score <- PITE::fixedCoxLassoInf_eta(x = reg.results$x,
                                                   y = reg.results$y,
                                                   beta = reg.results$beta,
                                                   lambda = reg.results$lambda,
                                                   status = reg.results$status,
                                                   alpha  = alpha, tol.beta = tol.beta,
                                                   gridrange = c(-gridrange_, gridrange_),
                                                   contrast = X.dx.full)

  scores <- data.frame(Dx = as.matrix(X.dx.full)%*%coef.lasso)
  scores$Dx.v <- out.score$sd
  scores$Dx.ll <- out.score$ci[,1]
  scores$Dx.ul <- out.score$ci[,2]
  vars = reg.out$term[which(reg.out$estimate!=0)]
  vars.score = vars[grep(pattern = "z", vars)]
  vars.score
  if (all(c("z:x1", "z:x2", "z") %in% vars.score)) {
    scores$pite <- dataset$pite
  } else   if (all(c("z:x1", "z:x2") %in% vars.score)) {
    scores$pite <- dataset$pite.x1x2
  } else   if (all(c("z:x1", "z") %in% vars.score)) {
    scores$pite <- dataset$pite.zx1
  } else   if (all(c("z:x2", "z") %in% vars.score)) {
    scores$pite <- dataset$pite.zx2
  } else   if (all(c("z") %in% vars.score)) {
    scores$pite <- dataset$pite.z
  } else   if (all(c("z:x2") %in% vars.score)) {
    scores$pite <- dataset$pite.x2
  } else   if (all(c("z:x1") %in% vars.score)) {
    scores$pite <- dataset$pite.x1
  } else  {
    scores$pite <- 0
  }
  # Store pite under full model
  scores$pite.full <- dataset$pite
  scores$cover <- (scores$Dx.ll < scores$pite) & (scores$pite < scores$Dx.ul)
  scores$width <- (scores$Dx.ul - scores$Dx.ll)
  scores$bias  <- (scores$Dx - scores$pite.full)
  head(scores)

  # mean(scores$cover)
  list(reg.out = reg.out,
       scores = scores,
       vars = vars)
}

#' Fits the reduced model obtained with the lasso but with no penalization on
#' the paraemeters
#'
#' A function that fits the model with treatment effect, prognostic terms
#' and interaction terms with coxph using the model obtained by the lasso
#' @param dataset a data.frame
#' @param alpha Significance level for confidence intervals (target is
#'  miscoverage alpha/2 in each tail)
#' @param calculate.ci logical. whether to calculate the confidence intervals for pite
#' @param vars List of selected variables.
#' @param n_biom number of considered biomarkers
#'
#' @export
score.cox.reduced <- function(dataset,
                              vars,
                              n_biom = 6,
                              alpha = 0.05,
                              calculate.ci = TRUE) {
  X.m <- model.matrix( ~ -1 + z * (x1 + x2 + x3 + x4 + x5 + x6), data = dataset)
  n <- nrow(dataset)
  all.vars <- colnames(X.m)
  if (is.null(vars)){
    scores <- data.frame(Dx = rep(0, n),
                         Dx.v = rep(0, n),
                         Dx.ll = rep(0, n),
                         Dx.ul = rep(0, n),
                         pite = rep(0, n),
                         pite.full = rep(0, n),
                         cover = rep(1, n),
                         width = rep(0, n),
                         bias = rep(0, n))
    lasso.results <- data.frame(term = all.vars,
                                estimate = 0,
                                p.value = NA,
                                conf.low = NA,
                                conf.high = NA)
    return(list(reg.out = lasso.results,
                scores = scores,
                vars = NULL))
  }

  formula <- paste("Surv(time, status) ~", paste(vars, collapse = " + "))
  reg.results <- survival::coxph(eval(parse(text = formula)),
                                 data = dataset, ties = "breslow")
  ML.output.M <- broom::tidy(reg.results)[, c(1, 2, 5, 6, 7)]
  int <- ML.output.M$term[grepl(":z", ML.output.M$term)]
  if (length(int) != 0) {
    sub(":z", "", int)
    ML.output.M$term[grepl(":z", ML.output.M$term)] <- paste0("z:", sub(":z", "", int))
  }

  var.score <- all.vars[2:(n_biom + 1)]
  ml.results.m <- data.frame(term      = colnames(X.m),
                             estimate  = 0,
                             p.value   = NA,
                             conf.low  = NA,
                             conf.high = NA)
  ind <- match(ML.output.M$term, ml.results.m$term)
  ml.results.m[ind, 2:5] <- ML.output.M[2:5]

  if (calculate.ci){
    coef.m <- ml.results.m$estimate

    covb.M <- vcov(reg.results)        # Store covariance matrix for the estimates
    X.dx.full <- 2 * cbind(1, matrix(0, nrow(dataset), length(var.score)), dataset[var.score])
    names(X.dx.full) <- colnames(X.m)

    Dx.M.m <- as.matrix(X.dx.full) %*% coef.m #Create score or PITE
    covb.M.dx <- covb.M[grepl("z", colnames(covb.M)),
                        grepl("z", colnames(covb.M))]
    # Calculate variance of the score
    vDx.M.m   <- diag(as.matrix(X.dx.full[, coef.m!=0]) %*%
                        covb.M %*%
                        t(X.dx.full[, coef.m!=0]))

    scores <- data.frame(Dx = Dx.M.m)
    scores$Dx.v  <- vDx.M.m
    zsigma <- qnorm(1 - alpha / 2) * sqrt(vDx.M.m)
    scores$Dx.ll <- drop(Dx.M.m - zsigma)
    scores$Dx.ul <- drop(Dx.M.m + zsigma)

    vars.score = vars[grep(pattern = "z", vars)]
    vars.score
    if (all(c("z:x1", "z:x2", "z") %in% vars.score)) {
      scores$pite <- dataset$pite
    } else   if (all(c("z:x1", "z:x2") %in% vars.score)) {
      scores$pite <- dataset$pite.x1x2
    } else   if (all(c("z:x1", "z") %in% vars.score)) {
      scores$pite <- dataset$pite.zx1
    } else   if (all(c("z:x2", "z") %in% vars.score)) {
      scores$pite <- dataset$pite.zx1
    } else   if (all(c("z") %in% vars.score)) {
      scores$pite <- dataset$pite.z
    } else   if (all(c("z:x1") %in% vars.score)) {
      scores$pite <- dataset$pite.x1
    } else   if (all(c("z:x2") %in% vars.score)) {
      scores$pite <- dataset$pite.x2
    } else  {
      scores$pite <- 0
    }
    scores$pite.full <- dataset$pite
    scores$cover <- (scores$Dx.ll < scores$pite) & (scores$pite < scores$Dx.ul)
    scores$width <- (scores$Dx.ul - scores$Dx.ll)
    scores$bias  <- (scores$Dx - scores$pite)
  } else {
    scores = NULL
  }
  # Store pite under full model
  list(reg.out = ml.results.m,
       scores = scores,
       reg.results = reg.results)
}

#' Obtains the PITE and confidence interval for an additional subject
#' using the reduced model from the Lasso output
#'
#' A function that takes the results from score.cox.reduced to get the PITE
#' and the confidence intervals for one subject.
#'
#' @param dataset a data.frame
#' @param alpha Significance level for confidence intervals (target is
#'  miscoverage alpha/2 in each tail)
#' @param vars List of selected variables.
#' @param n_biom number of considered biomarkers
#' @param reg.results an object resulting from score.cox.reduced
#'
#' @export
ci.cox.reduced <- function(dataset,
                           vars,
                           reg.results,
                           reg.out,
                           n_biom = 6,
                           alpha = 0.05) {
  X.m <- model.matrix( ~ -1 + z * (x1 + x2 + x3 + x4 + x5 + x6), data = dataset)
  n <- nrow(dataset)
  all.vars <- colnames(X.m)
  if (is.null(vars)){
    scores <- data.frame(Dx = rep(0, n),
                         Dx.v = rep(0, n),
                         Dx.ll = rep(0, n),
                         Dx.ul = rep(0, n),
                         pite = rep(0, n),
                         pite.full = rep(0, n),
                         cover = rep(1, n),
                         width = rep(0, n),
                         bias = rep(0, n) - dataset$pite)
    lasso.results <- data.frame(term = all.vars,
                                estimate = 0,
                                p.value = NA,
                                conf.low = NA,
                                conf.high = NA)
    return(list(reg.out = lasso.results,
                scores = scores,
                vars = NULL))
  }


  var.score <- all.vars[2:(n_biom + 1)]
  coef.m <- reg.out$estimate
  covb.M <- vcov(reg.results)        # Store covariance matrix for the estimates
  X.dx.full <- 2 * cbind(1, matrix(0, nrow(dataset), length(var.score)), dataset[var.score])
  names(X.dx.full) <- colnames(X.m)

  Dx.M.m <- as.matrix(X.dx.full) %*% coef.m #Create score or PITE
  covb.M.dx <- covb.M[grepl("z", colnames(covb.M)),
                      grepl("z", colnames(covb.M))]

  # Calculate variance of the score
  vDx.M.m   <- diag(as.matrix(X.dx.full[, coef.m != 0]) %*%
                      covb.M %*%
                      t(X.dx.full[, coef.m != 0]))

  scores <- data.frame(Dx = Dx.M.m)
  scores$Dx.v  <- vDx.M.m
  zsigma <- qnorm(1 - alpha / 2) * sqrt(vDx.M.m)
  scores$Dx.ll <- drop(Dx.M.m - zsigma)
  scores$Dx.ul <- drop(Dx.M.m + zsigma)

  vars.score = vars[grep(pattern = "z", vars)]
  if (all(c("z:x1", "z:x2", "z") %in% vars.score)) {
    scores$pite <- dataset$pite
  } else   if (all(c("z:x1", "z:x2") %in% vars.score)) {
    scores$pite <- dataset$pite.x1x2
  } else   if (all(c("z:x1", "z") %in% vars.score)) {
    scores$pite <- dataset$pite.zx1
  } else   if (all(c("z:x2", "z") %in% vars.score)) {
    scores$pite <- dataset$pite.zx1
  } else   if (all(c("z") %in% vars.score)) {
    scores$pite <- dataset$pite.z
  } else   if (all(c("z:x1") %in% vars.score)) {
    scores$pite <- dataset$pite.x1
  } else   if (all(c("z:x2") %in% vars.score)) {
    scores$pite <- dataset$pite.x2
  } else  {
    scores$pite <- 0
  }
  scores$pite.full <- dataset$pite
  scores$cover <- (scores$Dx.ll < scores$pite) & (scores$pite < scores$Dx.ul)
  scores$width <- (scores$Dx.ul - scores$Dx.ll)
  scores$bias  <- (scores$Dx - scores$pite)
  # Store pite under full model
  list(reg.out = reg.out,
       scores = scores,
       vars = vars)
}
