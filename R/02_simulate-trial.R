#' Simulate datasets and Create score variables in each of them
#'
#' Creates \code{nsim} datasets with specified parameters,
#' performs lm, lasso and randomized lasso in this dataset and creates PITE and its CI.
#' Up to two biomarkers may have an effect on the response (mkr1 and mkr2).
#'
#' @param nsim number of simulated datasets to be generated
#' @param effects a 6-columns matrix of data.frame
#' with the parameters for a, b, g1, g2, d1 and d2 in the model.
#' The rows correspond to different cases to analyze
#' @param case Specify the row of 'effects' that needs to be evaluated
#' @param npergroup number of subjects per treatment arm
#' @param mc.cores number of cores for parallel computing. To be passed to mclapply
#' @param verbose logical. whether to display intermediate results.
#' @param n_biom number of biomarkers in the model
#' @param param Coding to be used. Only "EFFECT" is allowed now. Treatment and binary variables are coded -1/1
#' @param lambda penalization to be used for the lasso. Either "lambda.min", "lambda.1se", "lagrange", or a number
#' @param lam_frac numeric. a fraction to multiply the lagrange parameter in the fixed lambda case.
#' @param nfolds number of folds in crossvalidation for glmnet
#' @param typeCorr Correlation between biomarkers. See \link{input}. For "block" use simulateNtrial_block.
#' @param rho correlation between biomarkers when typeCorr is "AR1" or "CS"
#' @param alpha significance level to build the confidence intervals
#' @param test logical. whether to evaluate the score in a test dataset
#' @param n.test numeric. The size of the test dataset.
#' @param perturb_frac a numeric vector of length 2 that indicates the degree of noise to introduce in the randomized lasso problem.
#' @param par.method Only "mclapply" is allowed
#' @param pite.ci logical. whether to calculate the ci for pite in the train data.
#' @param n.pite number of subjects for which to analyze coverage in the train data if pite.ci = TRUE
#' @param tol.beta tolerance to determine a parameter estimate is 0. To be passed to the selectiveInference functions.
#'
#' @details
#' This function performs the simulations for the specified cases.
#' It assumes that up to two biomarkers are prognostic/predictive (one binary and one uniform if independent,
#' both normally distributed when considering correlated biomarkers).
#'
#'
#'
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <- OneData(input,parameters)
#' score.lm(dataset, input)
#'
#' @export
simulateNtrial <- function(nsim,
                           effects,
                           case,
                           npergroup,
                           mc.cores = 1,
                           verbose = TRUE,
                           n_biom = 2,
                           param = "EFFECT",
                           lambda = "lambda.min",
                           nfolds = 10,
                           typeCorr = "I",
                           rho = 0,
                           alpha = 0.05,
                           test = FALSE,
                           perturb_frac = c(0.2,0.8),
                           par.method = "mclapply",
                           pite.ci = TRUE,
                           n.pite = NULL,
                           n.test = NULL,
                           lam_frac = 1,
                           tol.beta = 0.1,
                           ... ){
  input <- MakeInput(nT = npergroup,
                     nC = npergroup,
                     n_biom = n_biom,
                     n_biom_pred = 2,
                     typeCorr = typeCorr,
                     rho = rho,
                     a_ = effects[case,"a"],
                     b_ = effects[case,"b"],
                     sigma = 1, ...)
  if (typeCorr!="I") { # If Correlated biomarkers, we use normal distribution
    types <- c(2, 2)
  } else {
    types <- c(1, 3)
  }
  parameters <- MakeParameters(prognosis = c(effects[case, "g1"], effects[case, "g2"]),
                               predictive = c(effects[case, "d1"], effects[case, "d2"]),
                               types = types,
                               prevalence = c(0.5, 0),
                               means = c(0, 0),
                               stddev = c(1, 1),
                               input = input)
  MakeInputTable(input)
  MakeParametersTable(input, parameters)

  if(par.method == "mclapply"){
    outN <- mclapply(1:nsim, function (x){
      seed. = round(runif(1, min = 1, max = 10000))
      set.seed(seed.)
      dataset <- OneData(input, parameters, standardize = FALSE, param = param)
      dataset.null <- score.null(dataset,
                                 input,
                                 verbose = verbose,
                                 alpha = alpha)
      dataset.lm <- score.lm(dataset,
                             input,
                             verbose = verbose,
                             alpha = alpha)
      dataset.lasso <- score.lasso(dataset,
                                   input,
                                   parameters = parameters,
                                   alpha = alpha,
                                   verbose = verbose,
                                   lambda = lambda,
                                   nfolds = nfolds,
                                   lam_frac = lam_frac,
                                   tol.beta = tol.beta)
      dataset.mlm <- list(scores = dataset.lasso$scoresML)
      dataset.sch <- list(scores = dataset.lasso$scoresSch)
      lam <- dataset.lasso$bestlam * dataset.lasso$N
      dataset.lasso.an1 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac = perturb_frac[1],
                                                   alpha = alpha,
                                                   pite.ci = pite.ci,
                                                   n.pite = n.pite)
      dataset.lasso.an2 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac = perturb_frac[2],
                                                   alpha = alpha,
                                                   pite.ci = pite.ci,
                                                   n.pite = n.pite)
      if(is.null(n.pite)){
        n <- 2 * npergroup
        n.pite <- 2 * npergroup
      } else {
        n <- n.pite
      }
      results <- rbind(
        cbind(method = "null",  summarize_scores(dataset.null$scores,  dataset,  n.pite = n.pite), nvars = dataset.null$nvars),
        cbind(method = "lm",    summarize_scores(dataset.lm$scores,    dataset,  n.pite = n.pite), nvars = dataset.lm$nvars),
        cbind(method = "lasso", summarize_scores(dataset.lasso$scores, dataset,  n.pite = n.pite), nvars = dataset.lasso$nvars),
        cbind(method = "mlm",   summarize_scores(dataset.mlm$scores,   dataset,  n.pite = n.pite), nvars = dataset.lasso$nvars),
        cbind(method = "sch",   summarize_scores(dataset.sch$scores,   dataset,  n.pite = n.pite), nvars = dataset.lasso$nvars),
        cbind(method = "an1",   summarize_scores(dataset.lasso.an1$scores, dataset, n.pite = n.pite), nvars = dataset.lasso.an1$nvars),
        cbind(method = "an2",   summarize_scores(dataset.lasso.an2$scores, dataset, n.pite = n.pite), nvars = dataset.lasso.an2$nvars))

      results.coef <- rbind(
        cbind(method = "null",  dataset.null$ML.output[, -3],    nvars = dataset.null$nvars),
        cbind(method = "lm",    dataset.lm$ML.output[, -3],      nvars = dataset.lm$nvars),
        cbind(method = "lasso", dataset.lasso$Lasso.output,      nvars = dataset.lasso$nvars),
        cbind(method = "mlm",   dataset.lasso$ML.output.M[, -3], nvars = dataset.lasso$nvars),
        cbind(method = "an1",   dataset.lasso.an1$Lasso.output,  nvars = dataset.lasso.an1$nvars),
        cbind(method = "an2",   dataset.lasso.an2$Lasso.output,  nvars = dataset.lasso.an2$nvars))

      if (test){
        input.test <- input
        if(n.test == 1){
          dataset.test <- OneSubject(input.test, parameters, standardize = FALSE, param = param)
        } else {
          if(is.null(n.test)){
            n.test <- 2 * npergroup
          }
          input.test$nT <- n.test / 2 + 1
          input.test$nC <- n.test / 2 + 1
          dataset.test <- OneData(input.test, parameters, standardize = FALSE, param = param)
        }

        dataset.test.null  <- confidence.intervals.null.test(dataset = dataset.test, input = input.test,
                                                             parameters = parameters, ML.results = dataset.null,
                                                             alpha = alpha)
        dataset.test.ml    <- confidence.intervals.ml.test(dataset = dataset.test, input = input.test,
                                                           parameters = parameters, ML.results = dataset.lm,
                                                           alpha = alpha)
        dataset.test.lasso <- confidence.intervals.lasso.test(dataset = dataset.test, input = input.test,
                                                              parameters = parameters,
                                                              lasso.results = dataset.lasso,
                                                              alpha = alpha, tol.beta = tol.beta)
        dataset.test.mlm <- list(scores = dataset.test.lasso$scoresML)
        dataset.test.sch <- list(scores = dataset.test.lasso$scoresSch)
        dataset.test.lasso.an1 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input.test,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an1,
                                                                    alpha = alpha)
        dataset.test.lasso.an2 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input.test,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an2,
                                                                    alpha = alpha)
        results.test <- rbind(
          cbind(method = "null", summarize_scores(dataset.test.null$scores, dataset.test, n.pite = n.test)),
          cbind(method = "lm",   summarize_scores(dataset.test.ml$scores,   dataset.test, n.pite = n.test)),
          cbind(method = "lasso",summarize_scores(dataset.test.lasso$scores,dataset.test, n.pite = n.test)),
          cbind(method = "mlm",  summarize_scores(dataset.test.mlm$scores,  dataset.test, n.pite = n.test)),
          cbind(method = "sch",  summarize_scores(dataset.test.sch$scores,  dataset.test, n.pite = n.test)),
          cbind(method = "an1",  summarize_scores(dataset.test.lasso.an1$scores, dataset.test, n.pite = n.test)),
          cbind(method = "an2",  summarize_scores(dataset.test.lasso.an2$scores, dataset.test, n.pite = n.test)))

        out1<- list(results = results,
                    results.test = cbind(results.test, sim = x),
                    results.coef = results.coef,
                    tailarea = cbind(dataset.test.lasso$tailarea, sim = x),
                    seed = seed.)
      } else {
        out1<- list(results = results,
                    results.coef=results.coef,
                    tailarea = cbind(dataset.test.lasso$tailarea, sim = x),
                    seed = seed.)
      }
      class(out1) <- "simulate1trial"
      out1
    }, mc.cores = mc.cores, mc.preschedule = FALSE)
  }
  out <- list(dataset = outN)
  class(out) <- "simulateNtrial"
  out
}

#' Simulate datasets and Create score variables in each of them
#'
#' Creates \code{nsim} datasets with specified parameters,
#' performs lm, lasso and randomized lasso in this dataset and creates PITE and its CI. This function
#' assumes up to one prognostic/predictive biomarkers.
#' Up to one biomarkers may have an effect on the response (mkr1).
#'
#' @param nsim number of simulated datasets to be generated
#' @param effects a 6-columns matrix of data.frame
#' with the parameters for a, b, g1, g2, d1 and d2 in the model.
#' The rows correspond to different cases to analyze
#' @param case Specify the row of 'effects' that needs to be evaluated
#' @param npergroup number of subjects per treatment arm
#' @param mc.cores number of cores for parallel computing. To be passed to mclapply
#' @param verbose logical. whether to display intermediate results.
#' @param n_biom number of biomarkers in the model
#' @param param Coding to be used. Only "EFFECT" is allowed now. Treatment and binary variables are coded -1/1
#' @param lambda penalization to be used for the lasso. Either "lambda.min", "lambda.1se", "lagrange", or a number
#' @param lam_frac numeric. a fraction to multiply the lagrange parameter in the fixed lambda case.
#' @param nfolds number of folds in crossvalidation for glmnet
#' @param typeCorr Correlation between biomarkers. See \link{input}. For "block" use simulateNtrial_block.
#' @param rho correlation between biomarkers when typeCorr is "AR1" or "CS"
#' @param alpha significance level to build the confidence intervals
#' @param test logical. whether to evaluate the score in a test dataset
#' @param n.test numeric. The size of the test dataset.
#' @param perturb_frac a numeric vector of length 2 that indicates the degree of noise to introduce in the randomized lasso problem.
#' @param par.method Only "mclapply" is allowed
#' @param pite.ci logical. whether to calculate the ci for pite in the train data.
#' @param n.pite number of subjects for which to analyze coverage in the train data if pite.ci = TRUE
#' @param tol.beta tolerance to determine a parameter estimate is 0. To be passed to the selectiveInference functions.
#'
#' @details
#' This function performs the simulations for the specified cases.
#' It assumes that up to one biomarker is prognostic/predictive (binary if independent of other biomarkers,
#' normally distributed when considering correlated biomarkers).
#'
#'
#'
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <- OneData(input,parameters)
#' score.lm(dataset, input)
#'
#' @export
simulateNtrial_1pred <- function(nsim,
                                 effects,
                                 case,
                                 npergroup,
                                 mc.cores = 1,
                                 verbose = TRUE,
                                 n_biom = 2,
                                 param = "EFFECT",
                                 lambda = "lambda.min",
                                 nfolds = 10,
                                 typeCorr = "I",
                                 rho = 0,
                                 alpha = 0.05,
                                 test=FALSE,
                                 perturb_frac = c(0.2, 0.8),
                                 par.method = "mclapply",
                                 pite.ci = FALSE,
                                 n.pite = NULL,
                                 lam_frac = 1){
  input <- MakeInput(nT = npergroup,
                     nC = npergroup,
                     n_biom = n_biom,
                     n_biom_pred = 1,
                     typeCorr = typeCorr,
                     rho = rho,
                     a_ = effects[case,"a"],
                     b_ = effects[case,"b"],
                     sigma = 1)
  if (typeCorr!="I") { # If the Correlated biomarkers, we use normal distribution
    types <- c(2)
  } else {
    types <- c(1)
  }
  parameters <- MakeParameters(prognosis = c(effects[case, "g1"]),
                               predictive = c(effects[case, "d1"]),
                               types = types,
                               prevalence = c(0.5),
                               means = c(0),
                               stddev = c(1),
                               input = input)
  MakeInputTable(input)
  MakeParametersTable(input, parameters)
  if(par.method=="mclapply"){
    outN <- mclapply(1:nsim, function (x){
      seed. = round(runif(1, min = 1, max = 10000))
      set.seed(seed.)
      dataset<-OneData(input, parameters, standardize = FALSE, param = param)
      dataset.lm <- score.lm(dataset,
                             input,
                             verbose = verbose,
                             alpha = alpha)
      dataset.lasso <- score.lasso(dataset,
                                   input,
                                   parameters = parameters,
                                   alpha = alpha,
                                   verbose = verbose,
                                   lambda = lambda,
                                   nfolds = nfolds,
                                   lam_frac = lam_frac)

      dataset.mlm <- list(scores = dataset.lasso$scoresML)
      dataset.sch <- list(scores = dataset.lasso$scoresSch)
      lam <- dataset.lasso$bestlam*dataset.lasso$N
      dataset.lasso.an1 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac = perturb_frac[1],
                                                   pite.ci = pite.ci,
                                                   n.pite = n.pite)
      dataset.lasso.an2 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac = perturb_frac[2],
                                                   pite.ci = pite.ci,
                                                   n.pite = n.pite)
      if(is.null(n.pite)){
        n <- 2 * npergroup
        n.pite <- 2 * npergroup
      } else {
        n <- n.pite
      }
      results <- rbind(
        cbind(method = "lm",    summarize_scores(dataset.lm$scores,    dataset, n.pite = n.pite), nvars = dataset.lm$nvars),
        cbind(method = "lasso", summarize_scores(dataset.lasso$scores, dataset, n.pite = n.pite), nvars = dataset.lasso$nvars),
        cbind(method = "mlm",   summarize_scores(dataset.mlm$scores,   dataset, n.pite = n.pite), nvars = dataset.lasso$nvars),
        cbind(method = "sch",   summarize_scores(dataset.sch$scores,   dataset, n.pite = n.pite), nvars = dataset.lasso$nvars),
        cbind(method = "an1",   summarize_scores(dataset.lasso.an1$scores, dataset, n.pite = n.pite), nvars = dataset.lasso.an1$nvars),
        cbind(method = "an2",   summarize_scores(dataset.lasso.an2$scores, dataset, n.pite = n.pite), nvars = dataset.lasso.an2$nvars))

      results.coef <- rbind(
        cbind(method = "lm",    dataset.lm$ML.output[,-3],  nvars = dataset.lm$nvars),
        cbind(method = "lasso", dataset.lasso$Lasso.output, nvars = dataset.lasso$nvars),
        cbind(method = "mlm",   dataset.lasso$ML.output.M[,-3], nvars = dataset.lasso$nvars),
        cbind(method = "an1",   dataset.lasso.an1$Lasso.output, nvars = dataset.lasso.an1$nvars),
        cbind(method = "an2",   dataset.lasso.an2$Lasso.output, nvars = dataset.lasso.an2$nvars))

      if (test){
        input.test <- input
        input.test$nT <- n.pite / 2 + 1
        input.test$nC <- n.pite / 2 + 1
        dataset.test <- OneData(input.test, parameters, standardize = FALSE, param = param)

        dataset.test.ml <- confidence.intervals.ml.test(dataset = dataset.test, input = input,
                                                        parameters = parameters, ML.results = dataset.lm,
                                                        alpha = alpha)
        dataset.test.lasso <- confidence.intervals.lasso.test(dataset = dataset.test,
                                                              input = input,
                                                              parameters = parameters,
                                                              lasso.results = dataset.lasso,
                                                              alpha = alpha)
        dataset.test.mlm <- list(scores = dataset.test.lasso$scoresML)
        dataset.test.sch <- list(scores = dataset.test.lasso$scoresSch)
        dataset.test.lasso.an1 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an1,
                                                                    alpha = alpha)
        dataset.test.lasso.an2 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an2,
                                                                    alpha = alpha)

        results.test <- rbind(
          cbind(method = "lm",   summarize_scores(dataset.test.ml$scores,    dataset.test, n.pite = n.pite)),
          cbind(method = "lasso",summarize_scores(dataset.test.lasso$scores, dataset.test, n.pite = n.pite)),
          cbind(method = "mlm",  summarize_scores(dataset.test.mlm$scores,   dataset.test, n.pite = n.pite)),
          cbind(method = "sch",  summarize_scores(dataset.test.sch$scores,   dataset.test, n.pite = n.pite)),
          cbind(method = "an1",  summarize_scores(dataset.test.lasso.an1$scores, dataset.test, n.pite = n.pite)),
          cbind(method = "an2",  summarize_scores(dataset.test.lasso.an2$scores, dataset.test, n.pite = n.pite)))

        out1<- list(results = results,
                    results.test = results.test,
                    results.coef = results.coef,
                    seed = seed.)
      } else {
        out1<- list(results = results,
                    results.coef = results.coef,
                    seed = seed.)
      }
      class(out1) <- "simulate1trial"
      out1
    }, mc.cores = mc.cores, mc.preschedule = FALSE)
  }
  out <- list(dataset = outN)
  class(out) <- "simulateNtrial"
  out
}





#' Simulate datasets and Create score variables in each of them
#'
#' Creates \code{nsim} datasets with specified parameters,
#' performs lm, lasso and randomized lasso in this dataset and creates PITE and its CI.
#' This function allows comparing randomized lasso with different degrees of noise.
#' It takes perturb_frac a vector of length 7 and performs seven model fits with the randomized
#' lasso.
#' Up to two biomarkers may have an effect on the response (mkr1 and mkr2).
#'
#' @param nsim number of simulated datasets to be generated
#' @param effects a 6-columns matrix of data.frame
#' with the parameters for a, b, g1, g2, d1 and d2 in the model.
#' The rows correspond to different cases to analyze
#' @param case Specify the row of 'effects' that needs to be evaluated
#' @param npergroup number of subjects per treatment arm
#' @param mc.cores number of cores for parallel computing. To be passed to mclapply
#' @param verbose logical. whether to display intermediate results.
#' @param n_biom number of biomarkers in the model
#' @param param Coding to be used. Only "EFFECT" is allowed now. Treatment and binary variables are coded -1/1
#' @param lambda penalization to be used for the lasso. Either "lambda.min", "lambda.1se", "lagrange", or a number
#' @param lam_frac numeric. a fraction to multiply the lagrange parameter in the fixed lambda case.
#' @param nfolds number of folds in crossvalidation for glmnet
#' @param typeCorr Correlation between biomarkers. See \link{input}. For "block" use simulateNtrial_block.
#' @param rho correlation between biomarkers when typeCorr is "AR1" or "CS"
#' @param alpha significance level to build the confidence intervals
#' @param test logical. whether to evaluate the score in a test dataset
#' @param n.test numeric. The size of the test dataset.
#' @param perturb_frac a numeric vector of length 7 that indicates the degree of noise to introduce in the randomized lasso problem.
#' @param par.method Only "mclapply" is allowed
#' @param pite.ci logical. whether to calculate the ci for pite in the train data.
#' @param n.pite number of subjects for which to analyze coverage in the train data if pite.ci = TRUE
#' @param tol.beta tolerance to determine a parameter estimate is 0. To be passed to the selectiveInference functions.
#'
#' @details
#' This function performs the simulations for the specified cases.
#'
#'
#'
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <- OneData(input,parameters)
#' score.lm(dataset, input)
#'
#' @export
simulateNtrial_randomized <- function(nsim,
                                      effects,
                                      case,
                                      npergroup,
                                      mc.cores = 1,
                                      verbose = TRUE,
                                      n_biom = 2,
                                      param = "EFFECT",
                                      lambda = "lambda.min",
                                      nfolds = 10,
                                      typeCorr = "I",
                                      rho = 0,
                                      alpha = 0.05,
                                      test = FALSE,
                                      perturb_frac = c(0, 0.1, 0.2, 0.8, 1, 1.5, 2),
                                      par.method = "mclapply",
                                      pite.ci = TRUE,
                                      n.pite = NULL,
                                      n.test = NULL,
                                      lam_frac = 1){
  input <- MakeInput(nT = npergroup,
                     nC = npergroup,
                     n_biom = n_biom,
                     n_biom_pred = 2,
                     typeCorr = typeCorr,
                     rho = rho,
                     a_ = effects[case, "a"],
                     b_ = effects[case, "b"],
                     sigma = 1)
  if (typeCorr!="I") { # If the Correlated biomarkers, we use normal distribution
    types <- c(2, 2)
  } else {
    types <- c(1, 3)
  }

  parameters <- MakeParameters(prognosis = c(effects[case, "g1"], effects[case, "g2"]),
                               predictive = c(effects[case, "d1"], effects[case, "d2"]),
                               types = types,
                               prevalence = c(0.5, 0),
                               means = c(0, 0),
                               stddev = c(1, 1),
                               input = input)
  MakeInputTable(input)
  MakeParametersTable(input, parameters)
  if(par.method=="mclapply"){
    outN <- mclapply(1:nsim, function (x){
      seed. = round(runif(1, min = 1, max = 10000))
      set.seed(seed.)
      dataset <-OneData(input, parameters, standardize = FALSE, param = param)
      dataset.null <- score.null(dataset,
                                 input,
                                 verbose = verbose,
                                 alpha = alpha)
      dataset.lasso <- score.lasso(dataset,
                                   input,
                                   parameters = parameters,
                                   alpha = alpha,
                                   verbose = verbose,
                                   lambda = lambda,
                                   nfolds = nfolds,
                                   lam_frac = lam_frac)
      lam<-dataset.lasso$bestlam*dataset.lasso$N
      dataset.lasso.an1 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac = perturb_frac[1],
                                                   pite.ci = pite.ci, alpha = alpha,
                                                   n.pite = n.pite)
      dataset.lasso.an2 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac = perturb_frac[2],
                                                   pite.ci = pite.ci, alpha = alpha,
                                                   n.pite = n.pite)
      dataset.lasso.an3 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac = perturb_frac[3],
                                                   pite.ci = pite.ci, alpha = alpha,
                                                   n.pite = n.pite)
      dataset.lasso.an4 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac = perturb_frac[4],
                                                   pite.ci = pite.ci, alpha = alpha,
                                                   n.pite = n.pite)
      dataset.lasso.an5 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac = perturb_frac[5],
                                                   pite.ci = pite.ci, alpha = alpha,
                                                   n.pite = n.pite)
      dataset.lasso.an6 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac = perturb_frac[6],
                                                   pite.ci = pite.ci, alpha = alpha,
                                                   n.pite = n.pite)
      dataset.lasso.an7 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac = perturb_frac[7],
                                                   pite.ci = pite.ci, alpha = alpha,
                                                   n.pite = n.pite)
      if(is.null(n.pite)){
        n<-2 * npergroup
        n.pite<-2 * npergroup
      } else {
        n<-n.pite
      }
      results <- rbind(
        cbind(method = "null", summarize_scores(dataset.null$scores, dataset,  n.pite = n.pite), nvars = dataset.null$nvars),
        cbind(method = "lasso",summarize_scores(dataset.lasso$scores,dataset,  n.pite = n.pite), nvars = dataset.lasso$nvars),
        cbind(method = "an1",  summarize_scores(dataset.lasso.an1$scores, dataset, n.pite = n.pite), nvars = dataset.lasso.an1$nvars),
        cbind(method = "an2",  summarize_scores(dataset.lasso.an2$scores, dataset, n.pite = n.pite), nvars = dataset.lasso.an2$nvars),
        cbind(method = "an3",  summarize_scores(dataset.lasso.an3$scores, dataset, n.pite = n.pite), nvars = dataset.lasso.an3$nvars),
        cbind(method = "an4",  summarize_scores(dataset.lasso.an4$scores, dataset, n.pite = n.pite), nvars = dataset.lasso.an4$nvars),
        cbind(method = "an5",  summarize_scores(dataset.lasso.an5$scores, dataset, n.pite = n.pite), nvars = dataset.lasso.an5$nvars),
        cbind(method = "an6",  summarize_scores(dataset.lasso.an6$scores, dataset, n.pite = n.pite), nvars = dataset.lasso.an6$nvars),
        cbind(method = "an7",  summarize_scores(dataset.lasso.an7$scores, dataset, n.pite = n.pite), nvars = dataset.lasso.an7$nvars))

      results.coef <- rbind(
        cbind(method = "null",   dataset.null$ML.output[,-3],  nvars = dataset.null$nvars),
        cbind(method = "Lasso",  dataset.lasso$Lasso.output,   nvars = dataset.lasso$nvars),
        cbind(method = "an1",  dataset.lasso.an1$Lasso.output, nvars = dataset.lasso.an1$nvars),
        cbind(method = "an2",  dataset.lasso.an2$Lasso.output, nvars = dataset.lasso.an2$nvars),
        cbind(method = "an3",  dataset.lasso.an3$Lasso.output, nvars = dataset.lasso.an3$nvars),
        cbind(method = "an4",  dataset.lasso.an4$Lasso.output, nvars = dataset.lasso.an4$nvars),
        cbind(method = "an5",  dataset.lasso.an5$Lasso.output, nvars = dataset.lasso.an5$nvars),
        cbind(method = "an6",  dataset.lasso.an6$Lasso.output, nvars = dataset.lasso.an6$nvars),
        cbind(method = "an7",  dataset.lasso.an7$Lasso.output, nvars = dataset.lasso.an7$nvars))

      if (test){
        input.test <- input
        if(is.null(n.test)){
          n.test <- 2 * npergroup
        }
        input.test$nT <- n.test / 2 + 1
        input.test$nC <- n.test / 2 + 1
        dataset.test <- OneData(input.test, parameters, standardize = FALSE, param = param)

        dataset.test.null <- confidence.intervals.null.test(dataset = dataset.test, input = input.test,
                                                            parameters = parameters, ML.results = dataset.null,
                                                            alpha = alpha)
        dataset.test.lasso <- confidence.intervals.lasso.test(dataset = dataset.test,
                                                              input = input.test,
                                                              parameters = parameters,
                                                              lasso.results = dataset.lasso,
                                                              alpha = alpha)
        dataset.test.lasso.an1 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input.test,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an1,
                                                                    alpha = alpha)
        dataset.test.lasso.an2 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input.test,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an2,
                                                                    alpha = alpha)
        dataset.test.lasso.an3 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input.test,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an3,
                                                                    alpha = alpha)
        dataset.test.lasso.an4 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input.test,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an4,
                                                                    alpha = alpha)
        dataset.test.lasso.an5 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input.test,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an5,
                                                                    alpha = alpha)
        dataset.test.lasso.an6 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input.test,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an6,
                                                                    alpha = alpha)
        dataset.test.lasso.an7 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input.test,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an7,
                                                                    alpha = alpha)

        results.test <- rbind(
          cbind(method = "null", summarize_scores(dataset.test.null$scores,  dataset.test, n.pite = n.test)),
          cbind(method = "lasso",summarize_scores(dataset.test.lasso$scores, dataset.test, n.pite = n.test)),
          cbind(method = "an1",  summarize_scores(dataset.test.lasso.an1$scores, dataset.test, n.pite = n.test)),
          cbind(method = "an2",  summarize_scores(dataset.test.lasso.an2$scores, dataset.test, n.pite = n.test)),
          cbind(method = "an3",  summarize_scores(dataset.test.lasso.an3$scores, dataset.test, n.pite = n.test)),
          cbind(method = "an4",  summarize_scores(dataset.test.lasso.an4$scores, dataset.test, n.pite = n.test)),
          cbind(method = "an5",  summarize_scores(dataset.test.lasso.an5$scores, dataset.test, n.pite = n.test)),
          cbind(method = "an6",  summarize_scores(dataset.test.lasso.an6$scores, dataset.test, n.pite = n.test)),
          cbind(method = "an7",  summarize_scores(dataset.test.lasso.an7$scores, dataset.test, n.pite = n.test)))

        out1<- list(results = results,
                    results.test = results.test,
                    results.coef = results.coef,
                    seed = seed.)
      } else {
        out1<- list(results = results,
                    results.coef = results.coef,
                    seed = seed.)
      }
      class(out1) <- "simulate1trial"
      out1
    }, mc.cores = mc.cores, mc.preschedule = FALSE)
  }
  out <- list(dataset = outN)
  class(out) <- "simulateNtrial"
  out
}


#' Simulate datasets and Create score variables in each of them
#'
#' Creates \code{nsim} datasets with specified parameters,
#' performs lm, lasso and randomized lasso in this dataset and creates PITE and its CI.
#' This function allows comparing the penalization parameters for the lasso
#'  with cross-validation and with fixed lambda.
#' Up to two biomarkers may have an effect on the response (mkr1 and mkr2).
#'
#'
#' @param nsim number of simulated datasets to be generated
#' @param effects a 6-columns matrix of data.frame
#' with the parameters for a, b, g1, g2, d1 and d2 in the model.
#' The rows correspond to different cases to analyze
#' @param case Specify the row of 'effects' that needs to be evaluated
#' @param npergroup number of subjects per treatment arm
#' @param mc.cores number of cores for parallel computing. To be passed to mclapply
#' @param verbose logical. whether to display intermediate results.
#' @param n_biom number of biomarkers in the model
#' @param param Coding to be used. Only "EFFECT" is allowed now. Treatment and binary variables are coded -1/1
#' @param lambda penalization to be used for the lasso. Either "lambda.min", "lambda.1se", "lagrange", or a number
#' @param lam_frac numeric. a fraction to multiply the lagrange parameter in the fixed lambda case.
#' @param nfolds number of folds in crossvalidation for glmnet
#' @param typeCorr Correlation between biomarkers. See \link{input}. For "block" use simulateNtrial_block.
#' @param rho correlation between biomarkers when typeCorr is "AR1" or "CS"
#' @param alpha significance level to build the confidence intervals
#' @param test logical. whether to evaluate the score in a test dataset
#' @param n.test numeric. The size of the test dataset.
#' @param perturb_frac a numeric vector of length 2 that indicates the degree of noise to introduce in the randomized lasso problem.
#' @param par.method Only "mclapply" is allowed
#' @param pite.ci logical. whether to calculate the ci for pite in the train data.
#' @param n.pite number of subjects for which to analyze coverage in the train data if pite.ci = TRUE
#' @param tol.beta tolerance to determine a parameter estimate is 0. To be passed to the selectiveInference functions.
#'
#' @details
#' This function performs the simulations for the specified cases.
#'
#'
#'
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <- OneData(input,parameters)
#' score.lm(dataset, input)
#'
#' @export
simulateNtrial_lasso <- function(nsim,
                                 effects,
                                 case,
                                 npergroup,
                                 mc.cores = 1,
                                 verbose = TRUE,
                                 n_biom = 2,
                                 param = "EFFECT",
                                 lambda = "lambda.min",
                                 nfolds = 10,
                                 typeCorr = "I",
                                 rho = 0,
                                 alpha = 0.05,
                                 test = FALSE,
                                 perturb_frac = c(0.2, 0.8),
                                 par.method = "mclapply",
                                 pite.ci = TRUE,
                                 n.pite = NULL,
                                 n.test = NULL,
                                 lam_frac = 1,
                                 tol.beta = 0.1){
  input <- MakeInput(nT = npergroup,
                     nC = npergroup,
                     n_biom = n_biom,
                     n_biom_pred=2,
                     typeCorr = typeCorr,
                     rho = rho,
                     a_ = effects[case,"a"],
                     b_ = effects[case,"b"],
                     sigma = 1)
  if (typeCorr!="I") { # If the Correlated biomarkers, we use normal distribution
    types <- c(2, 2)
  } else {types <- c(1, 3)}
  parameters <- MakeParameters(prognosis = c(effects[case,"g1"],effects[case,"g2"]),
                               predictive = c(effects[case,"d1"],effects[case,"d2"]),
                               types = types,
                               prevalence = c(0.5, 0),
                               means = c(0, 0),
                               stddev = c(1, 1),
                               input=input)
  MakeInputTable(input)
  MakeParametersTable(input, parameters)
  if(par.method=="mclapply"){
    outN <- mclapply(1:nsim, function (x){
      seed. = round(runif(1, min = 1, max = 10000))
      set.seed(seed.)
      dataset <-OneData(input, parameters, standardize = FALSE, param = param)
      dataset.lasso_min <- score.lasso(dataset,
                                       input,
                                       parameters = parameters,
                                       alpha = alpha,
                                       verbose = verbose,
                                       lambda = "lambda.min",
                                       nfolds = nfolds,
                                       pite.ci = pite.ci,
                                       lam_frac = lam_frac,
                                       tol.beta = tol.beta)
      dataset.lasso_1se <- score.lasso(dataset,
                                       input,
                                       parameters = parameters,
                                       alpha = alpha,
                                       verbose = verbose,
                                       lambda = "lambda.1se",
                                       nfolds = nfolds,
                                       pite.ci = pite.ci,
                                       lam_frac = lam_frac,
                                       tol.beta = tol.beta)
      dataset.lasso_lag <- score.lasso(dataset,
                                       input,
                                       parameters = parameters,
                                       alpha = alpha,
                                       verbose = verbose,
                                       lambda = "lagrange",
                                       nfolds = nfolds,
                                       pite.ci = pite.ci,
                                       lam_frac = lam_frac,
                                       tol.beta = tol.beta)
      out = c(min = dataset.lasso_min$bestlam,
              se  = dataset.lasso_1se$bestlam,
              lag = dataset.lasso_lag$bestlam,
              seed = seed.)
      out
    }, mc.cores = mc.cores, mc.preschedule = FALSE)
  }
  outN
}


#' Simulate datasets and Create score variables in each of them
#'
#' Creates \code{nsim} datasets with specified parameters,
#' performs lm, lasso and randomized lasso in this dataset and creates PITE and its CI.
#' This function is for the case where p>>n. In this case we do not perform the full model.
#' Only the lasso and randomized lasso are fitted.
#' Up to two biomarkers may have an effect on the response (mkr1 and mkr2).
#'
#' @param nsim number of simulated datasets to be generated
#' @param effects a 6-columns matrix of data.frame
#' with the parameters for a, b, g1, g2, d1 and d2 in the model.
#' The rows correspond to different cases to analyze
#' @param case Specify the row of 'effects' that needs to be evaluated
#' @param npergroup number of subjects per treatment arm
#' @param mc.cores number of cores for parallel computing. To be passed to mclapply
#' @param verbose logical. whether to display intermediate results.
#' @param n_biom number of biomarkers in the model
#' @param param Coding to be used. Only "EFFECT" is allowed now. Treatment and binary variables are coded -1/1
#' @param lambda penalization to be used for the lasso. Either "lambda.min", "lambda.1se", "lagrange", or a number
#' @param lam_frac numeric. a fraction to multiply the lagrange parameter in the fixed lambda case.
#' @param nfolds number of folds in crossvalidation for glmnet
#' @param typeCorr Correlation between biomarkers. See \link{input}. For "block" use simulateNtrial_block.
#' @param rho correlation between biomarkers when typeCorr is "AR1" or "CS"
#' @param alpha significance level to build the confidence intervals
#' @param test logical. whether to evaluate the score in a test dataset
#' @param n.test numeric. The size of the test dataset.
#' @param perturb_frac a numeric vector of length 2 that indicates the degree of noise to introduce in the randomized lasso problem.
#' @param par.method Only "mclapply" is allowed
#' @param pite.ci logical. whether to calculate the ci for pite in the train data.
#' @param n.pite number of subjects for which to analyze coverage in the train data if pite.ci = TRUE
#' @param tol.beta tolerance to determine a parameter estimate is 0. To be passed to the selectiveInference functions.
#'
#' @details
#' This function performs the simulations for the specified cases.
#'
#'
#'
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <- OneData(input,parameters)
#' score.lm(dataset, input)
#'
#' @export
simulateNtrial_many <- function(nsim,
                                effects,
                                case,
                                npergroup,
                                mc.cores = 1,
                                verbose = TRUE,
                                n_biom = 2,
                                param = "EFFECT",
                                lambda = "lambda.min",
                                nfolds = 10,
                                typeCorr = "I",
                                rho = 0,
                                alpha = 0.05,
                                test = FALSE,
                                perturb_frac = c(0.2, 0.8),
                                par.method = "mclapply",
                                pite.ci = TRUE,
                                n.pite = NULL,
                                n.test = NULL,
                                lam_frac = 1,
                                tol.beta = 0.1){
  input <- MakeInput(nT = npergroup,
                     nC = npergroup,
                     n_biom = n_biom,
                     n_biom_pred=2,
                     typeCorr = typeCorr,
                     rho = rho,
                     a_ = effects[case,"a"],
                     b_ = effects[case,"b"],
                     sigma = 1)
  if (typeCorr!="I") { # If the Correlated biomarkers, we use normal distribution
    types <- c(2, 2)
  } else {
    types <- c(1, 3)
  }
  parameters <- MakeParameters(prognosis  = c(effects[case, "g1"], effects[case, "g2"]),
                               predictive = c(effects[case, "d1"], effects[case, "d2"]),
                               types = types,
                               prevalence = c(0.5, 0),
                               means = c(0, 0),
                               stddev = c(1, 1),
                               input=input)
  MakeInputTable(input)
  MakeParametersTable(input, parameters)
  if(par.method=="mclapply"){
    outN <- mclapply(1:nsim, function (x){
      seed. = round(runif(1, min = 1, max = 10000))
      set.seed(seed.)
      dataset <- OneData(input, parameters, standardize = FALSE, param = param)
      dataset.lasso <- score.lasso(dataset,
                                   input,
                                   parameters = parameters,
                                   alpha = alpha,
                                   verbose = verbose,
                                   lambda = lambda,
                                   nfolds = nfolds,
                                   gridrange_ = 300,
                                   lam_frac=lam_frac, tol.beta = tol.beta)
      lam<-dataset.lasso$bestlam*dataset.lasso$N
      dataset.lasso.an1 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac = perturb_frac[1],
                                                   pite.ci = pite.ci, alpha = alpha,
                                                   n.pite = n.pite)
      dataset.lasso.an2 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac = perturb_frac[2],
                                                   pite.ci = pite.ci, alpha = alpha,
                                                   n.pite = n.pite)
      if(is.null(n.pite)){
        n <- 2 * npergroup
        n.pite <- 2 * npergroup
      } else {
        n <- n.pite
      }
      results <- rbind(
        cbind(method = "lasso",summarize_scores(dataset.lasso$scores,     dataset, n.pite = n.pite), nvars = dataset.lasso$nvars),
        cbind(method = "an1",  summarize_scores(dataset.lasso.an1$scores, dataset, n.pite = n.pite), nvars = dataset.lasso.an1$nvars),
        cbind(method = "an2",  summarize_scores(dataset.lasso.an2$scores, dataset, n.pite = n.pite), nvars = dataset.lasso.an2$nvars))

      results.coef <- rbind(
        cbind(method = "lasso",dataset.lasso$Lasso.output,     nvars = dataset.lasso$nvars),
        cbind(method = "an1",  dataset.lasso.an1$Lasso.output, nvars = dataset.lasso.an1$nvars),
        cbind(method = "an2",  dataset.lasso.an2$Lasso.output, nvars = dataset.lasso.an2$nvars))

      if (test){
        input.test <- input
        if(n.test == 1){
          dataset.test <- OneSubject(input.test, parameters, standardize = FALSE, param = param)
        } else {
          if(is.null(n.test)){
            n.test <- 2 * npergroup
          }
          input.test$nT <- n.test / 2 + 1
          input.test$nC <- n.test / 2 + 1
          dataset.test <- OneData(input.test, parameters, standardize = FALSE, param = param)
        }

        dataset.test.lasso <- confidence.intervals.lasso.test(dataset = dataset.test,
                                                              input = input.test,
                                                              parameters = parameters,
                                                              lasso.results = dataset.lasso,
                                                              gridrange_ = 300,
                                                              alpha = alpha, tol.beta = tol.beta)
        dataset.test.lasso.an1 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input.test,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an1,
                                                                    ndraw = 10000,
                                                                    burnin = 5000,
                                                                    alpha = alpha)
        dataset.test.lasso.an2 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input.test,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an2,
                                                                    ndraw = 10000,
                                                                    burnin = 5000,
                                                                    alpha = alpha)

        results.test <- rbind(
          cbind(method = "lasso",summarize_scores(dataset.test.lasso$scores,     dataset.test, n.pite = n.test)),
          cbind(method = "an1",  summarize_scores(dataset.test.lasso.an1$scores, dataset.test, n.pite = n.test)),
          cbind(method = "an2",  summarize_scores(dataset.test.lasso.an2$scores, dataset.test, n.pite = n.test)))

        out1<- list(results = cbind(results, sim = x),
                    results.test = cbind(results.test, sim = x),
                    results.coef = cbind(results.coef, sim = x),
                    tailarea = cbind(dataset.test.lasso$tailarea, sim = x),
                    seed = seed.)
      } else {
        out1<- list(results = results,
                    results.coef = results.coef)
      }
      class(out1) <- "simulate1trial"
      out1
    }, mc.cores = mc.cores, mc.preschedule = FALSE)
  }
  out <- list(dataset = outN)
  class(out) <- "simulateNtrial"
  out
}


#' Simulate datasets and Create score variables in each of them
#'
#' Creates \code{nsim} datasets with specified parameters,
#' performs lm, lasso and randomized lasso in this dataset and creates PITE and its CI.
#' This function allows simulating biomarkers with block correlation structure.
#' Up to two biomarkers may have an effect on the response (mkr1 and mkr6).
#'
#' @param nsim number of simulated datasets to be generated
#' @param effects a 6-columns matrix of data.frame
#' with the parameters for a, b, g1, g2, d1 and d2 in the model.
#' The rows correspond to different cases to analyze
#' @param case Specify the row of 'effects' that needs to be evaluated
#' @param npergroup number of subjects per treatment arm
#' @param mc.cores number of cores for parallel computing. To be passed to mclapply
#' @param verbose logical. whether to display intermediate results.
#' @param n_biom number of biomarkers in the model
#' @param param Coding to be used. Only "EFFECT" is allowed now. Treatment and binary variables are coded -1/1
#' @param lambda penalization to be used for the lasso. Either "lambda.min", "lambda.1se", "lagrange", or a number
#' @param lam_frac numeric. a fraction to multiply the lagrange parameter in the fixed lambda case.
#' @param nfolds number of folds in crossvalidation for glmnet
#' @param typeCorr Correlation between biomarkers. See \link{input}. For simulateNtrial_block this should be "block.
#' @param blocks Number of blocks. Set to 2
#' @param nvarblocks1_ Number of variables in block 1
#' @param nvarblocks2_ Number of variables in block 2
#' @param rhoblocks1_ Correlation in block 1
#' @param rhoblocks2_ Correlation in block 2
#' @param rho correlation between biomarkers when typeCorr is "AR1" or "CS"
#' @param alpha significance level to build the confidence intervals
#' @param test logical. whether to evaluate the score in a test dataset
#' @param n.test numeric. The size of the test dataset.
#' @param perturb_frac a numeric vector of length 2 that indicates the degree of noise to introduce in the randomized lasso problem.
#' @param par.method Only "mclapply" is allowed
#' @param pite.ci logical. whether to calculate the ci for pite in the train data.
#' @param n.pite number of subjects for which to analyze coverage in the train data if pite.ci = TRUE
#' @param tol.beta tolerance to determine a parameter estimate is 0. To be passed to the selectiveInference functions.
#'
#' @details
#' This function performs the simulations for the specified cases.
#' It assumes that up to two biomarker is prognostic/predictive (binary if independent of other biomarkers,
#' normally distributed when considering correlated biomarkers).
#'
#'
#'
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <- OneData(input,parameters)
#' score.lm(dataset, input)
#'
#' @export
simulateNtrial_block <- function(nsim,
                                 effects,
                                 case,
                                 npergroup,
                                 mc.cores = 1,
                                 verbose = TRUE,
                                 n_biom = 2,
                                 param = "EFFECT",
                                 lambda = "lambda.min",
                                 nfolds = 10,
                                 typeCorr = "block",
                                 rho = 0,
                                 alpha = 0.05,
                                 test = FALSE,
                                 perturb_frac = c(0.2, 0.8),
                                 par.method = "mclapply",
                                 blocks = 2,
                                 nvarblocks1_ = 5,
                                 nvarblocks2_ = 5,
                                 rhoblocks1_ = 0.5,
                                 rhoblocks2_ = 0.5,
                                 pite.ci = TRUE,
                                 n.pite = NULL,
                                 n.test = NULL,
                                 lam_frac = 1,
                                 tol.beta = 0.1, ...){
  input <- MakeInput(nT = npergroup,
                     nC = npergroup,
                     n_biom = n_biom,
                     n_biom_pred = n_biom,
                     typeCorr = typeCorr,
                     rho = rho,
                     a_ = effects[case,"a"],
                     b_ = effects[case,"b"],
                     blocks = blocks,
                     nvarblocks1_ = nvarblocks1_,
                     nvarblocks2_ = nvarblocks2_,
                     rhoblocks1_ = rhoblocks1_,
                     rhoblocks2_ = rhoblocks2_,
                     sigma = 1, ...)
  if (typeCorr != "block") warning("This function is for typeCorr = 'block'")
  types <- rep(2, n_biom)
  parameters <- MakeParameters(prognosis  = c(effects[case, "g1"], rep(0, nvarblocks1_ - 1),
                                              effects[case, "g2"], rep(0, nvarblocks2_ - 1)),
                               predictive = c(effects[case, "d1"], rep(0, nvarblocks1_ - 1),
                                              effects[case, "d2"], rep(0, nvarblocks2_ - 1)),
                               prevalence = rep(0.5, n_biom),
                               means  = rep(0, n_biom),
                               stddev = rep(1, n_biom),
                               types  = types,
                               input  = input)
  MakeInputTable(input)
  MakeParametersTable(input, parameters)
  if(par.method=="mclapply"){
    outN <- mclapply(1:nsim, function (x){
      seed. = round(runif(1, min = 1, max = 10000))
      set.seed(seed.)
      dataset <- OneData(input, parameters, standardize = FALSE, param = param)
      dataset.null <- score.null(dataset,
                                 input,
                                 verbose = verbose,
                                 alpha = alpha)
      dataset.lm <- score.lm(dataset,
                             input,
                             verbose = verbose,
                             alpha = alpha)
      dataset.lasso <- score.lasso(dataset,
                                   input,
                                   parameters = parameters,
                                   alpha = alpha,
                                   verbose = verbose,
                                   lambda = lambda,
                                   nfolds = nfolds,
                                   lam_frac = lam_frac, tol.beta = tol.beta)
      dataset.mlm <- list(scores = dataset.lasso$scoresML)
      dataset.sch <- list(scores = dataset.lasso$scoresSch)
      lam <- dataset.lasso$bestlam*dataset.lasso$N
      dataset.lasso.an1 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac = perturb_frac[1],
                                                   pite.ci = pite.ci, alpha = alpha,
                                                   n.pite = n.pite)
      dataset.lasso.an2 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac = perturb_frac[2],
                                                   pite.ci = pite.ci, alpha = alpha,
                                                   n.pite = n.pite)
      if(is.null(n.pite)){
        n <- 2 * npergroup
        n.pite <- 2 * npergroup
      } else {
        n <- n.pite
      }
      results <- rbind(
        cbind(method = "null", summarize_scores(dataset.null$scores, dataset,  n.pite = n.pite), nvars = dataset.null$nvars),
        cbind(method = "lm",   summarize_scores(dataset.lm$scores,   dataset,  n.pite = n.pite), nvars = dataset.lm$nvars),
        cbind(method = "lasso",summarize_scores(dataset.lasso$scores,dataset,  n.pite = n.pite), nvars = dataset.lasso$nvars),
        cbind(method = "mlm",  summarize_scores(dataset.mlm$scores,  dataset,  n.pite = n.pite), nvars = dataset.lasso$nvars),
        cbind(method = "sch",  summarize_scores(dataset.sch$scores,  dataset,  n.pite = n.pite), nvars = dataset.lasso$nvars),
        cbind(method = "an1",  summarize_scores(dataset.lasso.an1$scores, dataset, n.pite = n.pite), nvars = dataset.lasso.an1$nvars),
        cbind(method = "an2",  summarize_scores(dataset.lasso.an2$scores, dataset, n.pite = n.pite), nvars = dataset.lasso.an2$nvars))

      results.coef <- rbind(
        cbind(method = "null", dataset.null$ML.output[,-3],nvars = dataset.null$nvars),
        cbind(method = "lm",   dataset.lm$ML.output[,-3],  nvars = dataset.lm$nvars),
        cbind(method = "lasso",dataset.lasso$Lasso.output, nvars = dataset.lasso$nvars),
        cbind(method = "mlm",  dataset.lasso$ML.output.M[,-3], nvars = dataset.lasso$nvars),
        cbind(method = "an1",  dataset.lasso.an1$Lasso.output, nvars = dataset.lasso.an1$nvars),
        cbind(method = "an2",  dataset.lasso.an2$Lasso.output, nvars = dataset.lasso.an2$nvars))

      if (test){
        input.test <- input
        if(n.test == 1){
          dataset.test <- OneSubject(input.test, parameters, standardize = FALSE, param = param)
        } else {
          if(is.null(n.test)){
            n.test <- 2 * npergroup
          }
          input.test$nT <- n.test / 2 + 1
          input.test$nC <- n.test / 2 + 1
          dataset.test <- OneData(input.test, parameters, standardize = FALSE, param = param)
        }

        dataset.test.null <- confidence.intervals.null.test(dataset = dataset.test, input = input.test,
                                                            parameters = parameters, ML.results = dataset.null,
                                                            alpha = alpha)
        dataset.test.ml   <- confidence.intervals.ml.test(dataset = dataset.test, input = input.test,
                                                          parameters = parameters, ML.results = dataset.lm,
                                                          alpha = alpha)
        dataset.test.lasso <- confidence.intervals.lasso.test(dataset = dataset.test,
                                                              input = input.test,
                                                              parameters = parameters,
                                                              lasso.results = dataset.lasso,
                                                              alpha = alpha, tol.beta = tol.beta)
        dataset.test.mlm <- list(scores = dataset.test.lasso$scoresML)
        dataset.test.sch <- list(scores = dataset.test.lasso$scoresSch)
        dataset.test.lasso.an1 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input.test,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an1,
                                                                    alpha = alpha)
        dataset.test.lasso.an2 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input.test,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an2,
                                                                    alpha = alpha)

        results.test <- rbind(
          cbind(method = "null", summarize_scores(dataset.test.null$scores, dataset.test, n.pite = n.test)),
          cbind(method = "lm",   summarize_scores(dataset.test.ml$scores,   dataset.test, n.pite = n.test)),
          cbind(method = "lasso",summarize_scores(dataset.test.lasso$scores,dataset.test, n.pite = n.test)),
          cbind(method = "mlm",  summarize_scores(dataset.test.mlm$scores,  dataset.test, n.pite = n.test)),
          cbind(method = "sch",  summarize_scores(dataset.test.sch$scores,  dataset.test, n.pite = n.test)),
          cbind(method = "an1",  summarize_scores(dataset.test.lasso.an1$scores, dataset.test, n.pite = n.test)),
          cbind(method = "an2",  summarize_scores(dataset.test.lasso.an2$scores, dataset.test, n.pite = n.test)))

        out1<- list(results = results,
                    results.test = cbind(results.test, sim = x),
                    results.coef = results.coef,
                    tailarea = cbind(dataset.test.lasso$tailarea, sim = x),
                    seed = seed.)
      } else {
        out1<- list(results = results,
                    results.coef=results.coef,
                    tailarea = cbind(dataset.test.lasso$tailarea, sim = x),
                    seed = seed.)
      }
      class(out1) <- "simulate1trial"
      out1
    }, mc.cores = mc.cores, mc.preschedule = FALSE)
  }
  out <- list(dataset = outN)
  class(out) <- "simulateNtrial"
  out
}
