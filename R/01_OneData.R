#' Input parameters for simulations
#'
#' A function that will create a list of input values to be used in other functions
#' of this package. Mainly to specify the parameters in the simulation.
#'
#' @param y_type type of variable for the outcome variable: 1=Binary; 2=Normal
#' @param nT Number of subjects in the treatment group
#' @param nC Number of subjects in the control group
#' @param nsim Number of simulations to be performed
#' @param n_biom Number of biomarkers in the dataset
#' @param n_biom_pred Number of prognostic/predictive biomarkers
#' @param typeCorr Type of correlation matrix for the biomarkers: one of "I"
#' for independent biomarkers, "CS" for compound symetry,
#'   "AR1" for autocorrelated, "block" for a block matrix, or "manual"
#'   which allows specifying each correlation (to be used in shiny)
#' @param rho  Correlation coefficient if the matrix is specified to be AR1 or CS.
#'   must be in the range (-1, 1)
#' @param blocks Number of blocks in the matrix in case type.corr="block".
#' @param nvarblocks1_ If typeCorr = "block"; number of variables in block 1
#' @param nvarblocks2_ If typeCorr = "block"; number of variables in block 2
#' @param nvarblocks3_ If typeCorr = "block"; number of variables in block 3
#' @param rhoblocks1_  If typeCorr = "block"; correlation coefficient for
#'   variables in block 1
#' @param rhoblocks2_  If typeCorr = "block"; correlation coefficient for
#'   variables in block 2
#' @param rhoblocks3_  If typeCorr = "block"; correlation coefficient for
#'   variables in block 3
#' @param correlations TRUE/FALSE. Only used in shiny to display an input for
#'   manual correlation matrix
#' @param a_ A number indicating the baseline response in the patients
#' @param b_ A number indicating the treatment effect in the patients
#' @param sigma Standard deviation for the response variable y
#'
#' @examples
#' input <- MakeInput(y_type = 2, nT = 500,
#'                  nC = 500, n_biom = 30, n_biom_pred = 2)
#' input
#'
#' @export
MakeInput <- function(y_type = 2,   # 1=Binary; 2=Normal
                      nT = 50,
                      nC = 50,
                      nsim = 100,
                      n_biom = 1,
                      n_biom_pred = 1,
                      typeCorr = c("I","CS", "AR1", "block", "manual"),
                      rho = 0,           # Only if typeCorr in c("CS","AR1")
                      blocks = 1,        # If typeCorr = "block"
                      nvarblocks1_ = 2,  # If typeCorr = "block"
                      nvarblocks2_ = NA, # If typeCorr = "block"
                      nvarblocks3_ = NA, # If typeCorr = "block"
                      rhoblocks1_ = 0,   # If typeCorr = "block"
                      rhoblocks2_ = 0,   # If typeCorr = "block"
                      rhoblocks3_ = 0,   # If typeCorr = "block"
                      correlations = FALSE, # to display the corr in the app
                      a_ = 0,
                      b_ = 0,
                      sigma = 1) {
  typeCorr = match.arg(typeCorr)
  list(y_type = y_type,  # 1 = Binary; 2 = Normal
       nT = nT,
       nC = nC,
       nsim = nsim,
       n_biom = n_biom,
       n_biom_pred = n_biom_pred,
       typeCorr = typeCorr,           # "I","CS", "AR1", "block", "manual"
       rho = rho,                     # Only if typeCorr in c("CS","AR1")
       blocks = blocks,               # If typeCorr = "block"
       nvarblocks1_ = nvarblocks1_,   # If typeCorr = "block"
       nvarblocks2_ = nvarblocks2_,   # If typeCorr = "block"
       nvarblocks3_ = nvarblocks3_,   # If typeCorr = "block"
       rhoblocks1_  = rhoblocks1_,    # If typeCorr = "block"
       rhoblocks2_  = rhoblocks2_,    # If typeCorr = "block"
       rhoblocks3_  = rhoblocks3_,    # If typeCorr = "block"
       correlations = correlations,
       a_ = a_,
       b_ = b_,
       sigma = sigma)

}


#' Creates a nice table with the information of all the parameters.
#'
#' \code{MakeInputTable} returns a table with the information of all the
#'  parameters to create the dataset
#'
#' @param input A list of objects created with \code{MakeInput}
#'
#' @examples
#' input <- MakeInput()
#' MakeInputTable(input)
#'
#' @export
MakeInputTable <- function (input) {
  m1 <- as.data.frame(unlist(input))
  m1$p <- rownames(m1)
  m2 <- m1[c(1, 2, 3, 4, 5, 17, 18, 19), ]
  rownames(m2) <- NULL
  m2$Parameter <- c("Outcome Type",
                    "n for treatment group",
                    "n for control group",
                    "number of simulated datasets",
                    "number of biomarkers",
                    "Overall mean response (a)",
                    "Overall treatment effect (b)",
                    "Std dev of outcome")

  m3 <- m2[, c(3, 1)]
  names(m3) <- c("Parameter", "Value")
  knitr::kable(m3, align = "lr",
               row.names = FALSE,
               caption = "Input parameters for the dataset")
}



#' Parameters for biomarkers
#'
#' A function that will create a data.frame with the parameters to be used in
#' the simulations. It is used to specify distributions and effects of biomarkers
#'
#' @param prognosis A vector of length \code{n_biom_pred} indicating the
#'   prognostic coefficients in the model for the generation of the data
#' @param predictive A vector of length \code{n_biom_pred} indicating the
#'   predictive coefficients in the model for the generation of the data
#' @param types A vector of length \code{n_biom_pred} indicating the type of
#'   biomarker to be generated 1:Binary; 2:Normal, 3:Uniform. For uniform biomarkers
#'   we use U(-sqrt(3), sqrt(3)) to have std.dev=1
#' @param prevalence A vector of length \code{n_biom_pred} indicating the
#'   prevalence for the binary biomarkers. If the biomarker in the position i is
#'   not binary, the ith element of this vector is ignored
#' @param means A vector of length \code{n_biom_pred} indicating the mean for
#'   the NORMAL biomarkers. If the biomarker in the position i is not normal,
#'   the ith element of this vector is ignored
#' @param stddev A vector of length \code{n_biom_pred} indicating the standard
#'   deviation for the NORMAL biomarkers. If the biomarker in the position i is
#'   not normal, the ith element of this vector is ignored
#' @param input an object created with \code{MakeInput} function
#'
#' @examples
#' input <- MakeInput(y_type = 2, nT = 500, nC = 500,
#'                    n_biom = 30, n_biom_pred = 2)
#' MakeParameters(prognosis = c(1, 0), predictive = c(0, 1), types = c(1, 1),
#'   prevalence = c(0.5, 0), means = c(0, 0), stddev = c(1, 1), input = input)
#'
#' @export
MakeParameters <- function(prognosis = numeric(),
                           predictive = numeric(),
                           types = numeric(),
                           prevalence = numeric(),
                           means = numeric(),
                           stddev = numeric(),
                           input) {
  countbiom <- length(unique(lengths(list(prognosis,predictive,prevalence,means,
                                          stddev))))

  if (countbiom != 1) {
    stop("Check length of vectors")
  }

  if (length(prognosis) == 0)   prognosis <- rep(0, input$n_biom_pred)
  if (length(predictive) == 0)  predictive <- rep(0, input$n_biom_pred)
  if (length(types) == 0)       types <- rep(1,input$n_biom_pred)
  if (length(prevalence) == 0)  prevalence <- rep(0.5, input$n_biom_pred)
  if (length(means) == 0)       means <- rep(0, input$n_biom_pred)
  if (length(stddev) == 0)      stddev <- rep(0, input$n_biom_pred)

  if (length(prognosis) != input$n_biom_pred) {
    stop("The number of prognostic/predictive biomarkers specified in the MakeInput object is different than the number of parameters specified in this function")
  }

  means[types == 1] <- 0
  stddev[types == 1] <- 1
  prevalence[types == 2] <- 0.5
  m <- data.frame(prognosis = prognosis,
                  predictive = predictive,
                  types = types,
                  prevalence = prevalence,
                  means = means,
                  stddev = stddev)
  m
}

#' Creates a nice table with the information of all the parameters
#'
#' \code{MakeParametersTable} returns a table with the information of all the
#'  biomarkers created in the dataset
#'
#' @param input A list of objects created with \code{MakeInput}
#' @param parameters A list of objects created with \code{MakeParameters}
#'
#' @examples
#' input<-MakeInput()
#' parameters<-MakeParameters(input = input)
#' MakeParametersTable(input,parameters)
#'
#' @export
MakeParametersTable <- function (input, parameters) {
  biom_type<-rep(".", input$n_biom_pred)
  for (i in 1:input$n_biom_pred){
    if (parameters$types[i] == 2) {
      biom_type[i] = "Normal"
    }
    if (parameters$types[i] == 1) {
      biom_type[i] = "Binary"
    }
    if (parameters$types[i] == 3) {
      biom_type[i] = "Uniform"
    }
  }
  m3 <- data.frame(Biomarker = (1:input$n_biom_pred),
                   c = parameters$prognosis,
                   d = parameters$predictive,
                   Type = biom_type,
                   Prevalence = parameters$prevalence,
                   Mean = parameters$means,
                   Stddev = parameters$stddev)
  m3[which(m3[, "Type"] == "Binary"), "Mean"] <- NA
  m3[which(m3[, "Type"] == "Binary"), "Stddev"] <- NA
  m3[which(m3[, "Type"] == "Normal"), "Prevalence"] <- NA
  m3[which(m3[, "Type"] == "Uniform"), "Mean"] <- 0
  m3[which(m3[, "Type"] == "Uniform"), "Stddev"] <- sqrt((sqrt(3)-(-sqrt(3)))^2/12)

  n_biom_nopred <- input$n_biom-input$n_biom_pred
  if (n_biom_nopred > 0){
    m4 <- data.frame(Biomarker = ((input$n_biom_pred + 1):input$n_biom),
                     c = rep(0, n_biom_nopred),
                     d = rep(0, n_biom_nopred),
                     Type = rep("Normal", n_biom_nopred),
                     Prevalence = rep(NA, n_biom_nopred),
                     Mean = rep(0, n_biom_nopred),
                     Stddev = rep(1, n_biom_nopred))
    m3 <- rbind(m3, m4)
  }
  knitr::kable(m3,
               row.names = FALSE,
               caption = "Input parameters for prognostic and predictive effects of the biomarkers")
}


#' Generate Correlation Matrix
#'
#' A generic function that will create a correlation matrix depending on the
#' specified structure and the parameters
#'
#' @param nvar Number of variables that will determine the matrix dimension
#'   (\code{nvar} x \code{nvar})
#' @param type.corr Type of correlation matrix to be generated. Options:
#'   \code{"AR1"} create a first-order autoregressive correlation matrix;
#'   \code{"I"} creates an identity matrix; \code{"CS"} creates a Compound
#'   simmetry matrix; \code{"block"} creates a matrix with blocks of equally
#'   correlated variables; and \code{"step"} is used if the 3 cuartiles are used
#'   to dichotomize a variable and form 3 different variables.
#' @param rho Correlation coefficient if the matrix is specified to be AR1 or
#'  CS. must be in the range (-1, 1)
#' @param blocks Number of blocks in the matrix in case
#'   type.corr=\code{"block"}
#' @param nvar.blocks Vector of length = \code{blocks} with the number of
#'   variables in each block
#' @param rho.blocks Vector of length = \code{blocks} with the correlation
#'   coefficients for each block
#'
#' @examples
#' GenerateCorrMatrix(nvar = 3, type.corr = "I")
#' GenerateCorrMatrix(nvar = 5, type.corr = "AR1", rho = 0.3)
#' GenerateCorrMatrix(nvar = 5, type.corr = "CS", rho = 0.3)
#' GenerateCorrMatrix(nvar = 5, type.corr = "block", blocks = 2,
#'                    nvar.blocks = c(3, 2), rho.blocks = c(0.5, 0.3))
#'
#' @export
GenerateCorrMatrix <- function(nvar, type.corr = "I",
                               rho, blocks, nvar.blocks, rho.blocks) {
  #############################################################-
  ## AR(1)
  #############################################################-
  if (type.corr == "AR1"){
    times <- 1:nvar
    H <- abs(outer(times, times, "-"))
    V <- rho^H
  }

  #############################################################-
  ## Compound Symmetry (Equally correlated)
  #############################################################-
  if (type.corr == "CS"){
    V <- matrix(rep(rho, nvar*nvar), ncol=nvar)
    V[cbind(1:nvar, 1:nvar)] <- 1
  }

  #############################################################-
  ## independent
  #############################################################-
  if (type.corr == "I"){
    V <- matrix(rep(0, nvar*nvar), ncol=nvar)
    V[cbind(1:nvar, 1:nvar)] <- 1
  }

  #############################################################-
  ## Blocks of correlated biomarkers
  #############################################################-
  if (type.corr == "block"){
    if (nvar < sum(nvar.blocks)) {
      V <- NA
      stop(paste0("The total number of biomarkers cannot be smaller than",
                   "the sum of the number of biomarkers in each block"))
    } else {
      pos.blocks <- rep(1, blocks + 1)
      for (i in 1:blocks){
        pos.blocks[i + 1] <- pos.blocks[i] + nvar.blocks[i]
      }
      V <- matrix(rep(0, nvar*nvar), ncol=nvar)
      for (i in 1:blocks) {
        Vsmall <- matrix(rep(rho.blocks[i], nvar.blocks[i] * nvar.blocks[i]),
                         ncol=nvar.blocks[i])
        V[pos.blocks[i]:(pos.blocks[i + 1] - 1),
          pos.blocks[i]:(pos.blocks[i + 1] - 1)] <- Vsmall
      }
      V[cbind(1:nvar, 1:nvar)] <- 1
    }
  }

  #############################################################-
  ## Step corr matrix.
  #############################################################-
  if (type.corr == "step") {
    p1 <- 0.25
    p2 <- 0.5
    p3 <- 0.75
    corr12 <- ((1 - p2) - (1 - p1) * (1 - p2)) /
      sqrt(p1 * (1 - p1) * p2 * (1 - p2))
    corr13 <- ((1 - p3) - (1 - p1) * (1 - p3)) /
      sqrt(p1 * (1 - p1) * p3 * (1 - p3))
    pos.blocks <- rep(1, (nvar / 3) + 1)
    for (i in 1:(nvar / 3)) {
      pos.blocks[i + 1] <- pos.blocks[i] + 3
    }
    V <- matrix(rep(0, nvar * nvar), ncol = nvar)
    stepm <- matrix(c(1, corr12, corr13, corr12,
                      1, corr12, corr13, corr12, 1), ncol=3)
    for (i in 1:blocks) {
      V[pos.blocks[i]:(pos.blocks[i + 1] - 1),
        pos.blocks[i]:(pos.blocks[i + 1] - 1)] <- stepm
    }
    V[cbind(1:nvar, 1:nvar)] <- 1
  }
  V
}

#' Generate the Correlation Matrix for the biomarkers in a dataset
#'
#' A function that will create a correlation matrix using GenerateCorrM to be
#' used in the generation of a dataset
#'
#' @param input A list of objects created with \code{MakeInput}
#' @param parameters A data.frame of objects created with \code{MakeParameters}
#'
#' @examples
#' input <- MakeInput(n_biom = 5, typeCorr = "AR1", rho = 0.4)
#' parameters <- MakeParameters(input=input)
#' BiomCorrMatrix(input, parameters)
#'
#' @export
BiomCorrMatrix <- function(input, parameters) {
  corm <- diag(input$n_biom)
  if (input[["typeCorr"]] == 'CS') {
    corm <- GenerateCorrMatrix(nvar = input[["n_biom"]],
                               rho = input[["rho"]],
                               type.corr = "CS")
  }
  if (input[["typeCorr"]] == 'AR1') {
    corm <- GenerateCorrMatrix(nvar = input[["n_biom"]],
                               rho = input[["rho"]],
                               type.corr = "AR1")
  }
  if (input[["typeCorr"]] == 'block') {
    lista <- list("nvarblocks1_", "nvarblocks2_", "nvarblocks3_")
    a <- sapply(1:input$blocks, function(x) lista[[x]])
    b <- lapply(a, function(x) input[[x]])
    c <- unlist(b)
    nvarblocks <- c
    lista <- list("rhoblocks1_", "rhoblocks2_", "rhoblocks3_")
    a <- sapply(1:input$blocks,function(x) lista[[x]])
    b <- lapply(a,function(x) input[[x]])
    c <- unlist(b)
    rhoblocks <- c
    corm <- GenerateCorrMatrix(nvar = input[["n_biom"]], type.corr = 'block',
                               blocks = input$blocks, nvar.blocks = nvarblocks,
                               rho.blocks = rhoblocks)
  }
  if (input[["typeCorr"]] == 'step') {
    corm <- GenerateCorrMatrix(nvar = input[["n_biom"]],
                               type.corr = 'step',
                               blocks = input[["n_biom"]] / 3)
  }
  # The option typeCorr="Manual" is only used when executed in shiny
  if (input[["typeCorr"]] == 'manual' & input$correlations) {
    for (i in 1:input$n_biom) {
      for (j in 1:input$n_biom) {
        if (j > i) {
          corm[i, j]<-input[[paste0("corr", i, j)]]
        }
        if (j < i) {
          corm[i, j] <- input[[paste0("corr", j, i)]]
        }
      }
    }
  }

  # To generate now the Cov matrix we need the stddev parameters
  if (input$n_biom > input$n_biom_pred) {
    means  <- c(parameters$means,  rep(0, input$n_biom - input$n_biom_pred))
    stddev <- c(parameters$stddev, rep(1, input$n_biom - input$n_biom_pred))
  } else {
    means  <- parameters$means
    stddev <- parameters$stddev
  }
  covm <- diag(stddev) %*% corm %*% diag(stddev)
  out <- list(corm = corm,covm = covm)
  out
}





#' Creates one simulated dataset
#'
#' A function that creates a data.frame with a simulated dataset using the
#' specified input
#'
#' @param input A list of objects created with \code{MakeInput}
#' @param parameters A data.frame of objects created with \code{MakeParameters}
#' @param standardize logical. If TRUE all covariates are scaled to have 0 mean
#' and sd 1
#' @param param parametrization for the binary variables. Either "REF" or "EFFECT"
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' OneData(input,parameters)
#'
#' @export
OneData <- function(input, parameters, standardize = FALSE, param = "REF") {
  covm <- BiomCorrMatrix(input, parameters)$covm # covariance matrix
  # names of the variables to be used later
  mkrs <- paste0("mkr", 1:input[["n_biom"]])
  mydata <- data.frame(ID = 1:(input$nT + input$nC),
                       treatment = as.integer(c(rep(0, input$nT),
                                                rep(1, input$nC))))
  if (param == "EFFECT") {
    mydata$treatment <- mydata$treatment * 2 - 1
  }

  threshold <- qnorm(p = 1 - parameters$prevalence,
                     mean = 0,
                     sd = 1)
  if (input$n_biom > input$n_biom_pred) {
    means  <- c(parameters$means,
                rep(0, input$n_biom - input$n_biom_pred))
    stddev <- c(parameters$stddev,
                rep(1, input$n_biom - input$n_biom_pred))
  } else {
    means  <- parameters$means
    stddev <- parameters$stddev
  }
  if (input$n_biom > 1) {
    allMkrCols <- MASS::mvrnorm(n = (input$nT + input$nC),
                                mu = means,
                                Sigma = covm)
    for (z in  1:input$n_biom_pred) {
      if (parameters$types[z] == 1) {
        threshold <- qnorm(p    = 1 - parameters$prevalence[z],
                           mean = 0,
                           sd   = 1)
        allMkrCols[, z] <- 1 * (allMkrCols[, z] > threshold)
        if (param == "EFFECT") {
          allMkrCols[, z] <- allMkrCols[, z] * 2 - 1
        }
      }
      if (parameters$types[z] == 3) {
        # sqrt(3) gives mean=0 and V=1
        allMkrCols[, z] <- runif(n = input$nT + input$nC,
                                 min = -sqrt(3), max = sqrt(3))

      }
    }
  }

  if (input$n_biom == 1) {
    if (parameters$types[1] == 1) {
      allMkrCols <- data.frame("mkr1" = rbinom(n    = input$nT + input$nC,
                                               size = 1,
                                               prob = parameters$prevalence[1]))
      if (param == "EFFECT") {
        allMkrCols[, 1] <- allMkrCols[, 1] * 2 - 1
      }
    }
    if (parameters$types[1] == 2) {
      allMkrCols <- data.frame("mkr1" = rnorm(n = input$nT + input$nC,
                                              mean = parameters$means[1],
                                              sd   = parameters$stddev[1]))
    }
    if (parameters$types[1] == 3) {
      allMkrCols <- data.frame("mkr1" = runif(n = input$nT + input$nC,
                                              min = -sqrt(3), max = sqrt(3)))
    }
  }

  if (standardize == TRUE) {
    allMkrCols <- data.frame(scale(allMkrCols, center = TRUE, scale = TRUE))
  }

  if (input$n_biom_pred > 1) {
    wt1 <- diag(parameters$prognosis)
    wt2 <- diag(parameters$predictive)
    mkrCols <- as.matrix(allMkrCols[, 1:(input$n_biom_pred)])
    prog_term <- rowSums(mkrCols %*% wt1)
    pred_term <- rowSums(mkrCols %*% wt2) * mydata$treatment
    pred_termEff <- rowSums(mkrCols %*% wt2)
  } else {
    prog_term <- allMkrCols[, 1] * parameters$prognosis
    pred_term <- allMkrCols[, 1] * parameters$predictive * mydata$treatment
    pred_termEff <- allMkrCols[, 1] * parameters$predictive
  }
  mydata$mean <- input$a_ + input$b_ * mydata$treatment + prog_term + pred_term
  mydata$TrueTrtEff <- input$b_ +  pred_termEff
  if (param == "EFFECT") {
    mydata$TrueTrtEff <- mydata$TrueTrtEff * 2
  }
  if (input$y_type == 2) {
    mydata$y <- rnorm(n    = input$nT + input$nC,
                      mean = mydata$mean,
                      sd   = input$sigma)
  }
  if (input$y_type == 1) {
    # p.logis = plogis(drop(mydata$mean))
    # plogis.manual = exp(drop(mydata$mean))/(1 + exp(drop(mydata$mean)))
    mydata$y <- rbinom(n    = input$nT + input$nC,
                       size = 1,
                       prob = plogis(drop(mydata$mean)))
  }
  colnames(allMkrCols) <- mkrs
  allMkrCols <- data.frame(ID = mydata[, "ID"], allMkrCols)
  mydata <- merge(mydata, allMkrCols)
  mydata
}

#' Creates one simulated subject
#'
#' A function that will create a data.frame with one subject with simulated
#' covariates
#'
#' @param input A list of objects created with \code{MakeInput}
#' @param parameters A data.frame of objects created with \code{MakeParameters}
#' @param standardize logical. If TRUE all covariates are scaled to have 0 mean
#' and sd 1
#' @param param parametrization for the binary variables. Either "REF" or "EFFECT"
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#'
#' @export
OneSubject <- function(input, parameters, standardize = FALSE, param = "REF") {
  covm <- BiomCorrMatrix(input, parameters)$covm # cov matrix
  # names of the variables to be used later
  mkrs <- paste0("mkr", 1:input[["n_biom"]])
  mydata <- data.frame(ID = 1,
                       treatment = as.integer(1))
  if (param == "EFFECT") {
    mydata$treatment <- mydata$treatment * 2 - 1
  }

  threshold <- qnorm(p    = 1 - parameters$prevalence,
                     mean = 0,
                     sd   = 1)
  if (input$n_biom > input$n_biom_pred) {
    means  <- c(parameters$means,
                rep(0, input$n_biom - input$n_biom_pred))
    stddev <- c(parameters$stddev,
                rep(1, input$n_biom - input$n_biom_pred))
  } else {
    means  <- parameters$means
    stddev <- parameters$stddev
  }
  if (input$n_biom > 1) {
    allMkrCols <- MASS::mvrnorm(n     = 1,
                                mu    = means,
                                Sigma = covm)
    for (z in  1:input$n_biom_pred) {
      if (parameters$types[z] == 1) {
        threshold <- qnorm(p    = 1 - parameters$prevalence[z],
                           mean = 0,
                           sd   = 1)
        allMkrCols[z] <- 1 * (allMkrCols[z] > threshold)
        if (param == "EFFECT") {
          allMkrCols[z] <- allMkrCols[z] * 2 - 1
        }
      }
      if (parameters$types[z] == 3) {
        #sqrt(3) gives mean=0 and V=1
        allMkrCols[z] <- runif(n = 1,
                               min = -sqrt(3), max = sqrt(3))
      }
    }
  }

  if (input$n_biom == 1) {
    if (parameters$types[1] == 1) {
      allMkrCols <- data.frame("mkr1" = rbinom(n = 1,
                                               size = 1,
                                               prob = parameters$prevalence[1]))
      if (param == "EFFECT") {
        allMkrCols[1] <- allMkrCols[1] * 2 - 1
      }
    }
    if (parameters$types[1] == 2) {
      allMkrCols <- data.frame("mkr1" = rnorm(n = 1,
                                              mean = parameters$means[1],
                                              sd   = parameters$stddev[1]))
    }
    if (parameters$types[1] == 3) {
      allMkrCols <- data.frame("mkr1" = runif(n = 1,
                                              min = -sqrt(3), max = sqrt(3)))
    }
  }

  if (standardize == TRUE) {
    allMkrCols <- data.frame(scale(allMkrCols, center = TRUE, scale = TRUE))
  }

  if (input$n_biom_pred > 1) {
    wt1 <- diag(parameters$prognosis)
    wt2 <- diag(parameters$predictive)
    mkrCols <- matrix(allMkrCols[1:(input$n_biom_pred)], nrow = 1)
    prog_term <- rowSums(mkrCols %*% wt1)
    pred_term <- rowSums(mkrCols %*% wt2) * mydata$treatment
    pred_termEff <- rowSums(mkrCols %*% wt2)
  } else {
    prog_term <-    allMkrCols[1] * parameters$prognosis
    pred_term <-    allMkrCols[1] * parameters$predictive * mydata$treatment
    pred_termEff <- allMkrCols[1] * parameters$predictive
  }
  mydata$mean <- input$a_ + input$b_ * mydata$treatment + prog_term + pred_term
  mydata$TrueTrtEff <- input$b_ +  pred_termEff
  if (param == "EFFECT") {
    mydata$TrueTrtEff <- mydata$TrueTrtEff * 2
  }
  if (input$y_type == 2) {
    mydata$y <- rnorm(n    = 1,
                      mean = mydata$mean,
                      sd   = input$sigma)
  }
  if (input$y_type == 1) {
    # p.logis = plogis(drop(mydata$mean))
    # plogis.manual = exp(drop(mydata$mean))/(1 + exp(drop(mydata$mean)))
    mydata$y <- rbinom(n = 1,
                       size = 1,
                       prob = plogis(drop(mydata$mean)))
  }
  names(allMkrCols) <- mkrs
  allMkrCols <- data.frame(t(c(ID = mydata[, "ID"], allMkrCols)))
  mydata <- merge(mydata, allMkrCols)
  mydata
}
