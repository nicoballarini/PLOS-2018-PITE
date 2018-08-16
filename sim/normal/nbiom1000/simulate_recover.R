###############################################################################-
## This program is used to recover one simulation run
##
## Running in interactive mode
###############################################################################-
# Load needed libraries
library(PITE)
library(ggplot2)
library(monomvn)
library(parallel)
library(dplyr)
library(ggplot2)
library(selectiveInference)
RNGkind("L'Ecuyer-CMRG")
sim. = 11
j = 2 # For sample size
i = 1 # For Scenario (case)

# Number of biomarkers in this simulation
n_biom   =   1000
onlyPlots = F # Set to TRUE if simulations already run, and only plots want to be calculated
# Set the environment to store results by creating relative paths --------------
rproj <- rprojroot::find_package_root_file()
path  <- paste0(rproj,sprintf("/sim/normal/nbiom%s/",n_biom))
if (!dir.exists(paste0(path, "results/")))       dir.create(paste0(path, "results/"))
if (!dir.exists(paste0(path, "results/data/")))  dir.create(paste0(path, "results/data/"))
if (!dir.exists(paste0(path, "results/plots/"))) dir.create(paste0(path, "results/plots/"))
path.data <- paste0(path, "results/data/")
path.plot <- paste0(path, "results/plots/")

# Set initial parameters -------------------------------------------------------
alpha.   =    0.05
seed.    = 2005
nsim     = 1000
mc.cores =   35 # Number of cores for simulations

sample.size = c(40, 100, 220, 350) # Sample sizes considered in simulations
npergroup.  = sample.size/2

# Calculate effecs to be used in simulations
# First we suppose effects
effects <- expand.grid(b= c(0,1),   d1=c(0,1), d2=c(0,1))[c(1,5,3,7,2,8),]
# effects <- expand.grid(b= c(0,1/4), d1=c(0,1/2), d2=c(0,1/3))[c(1,5,3,7,2,8),]
all.ref <- cbind(case=1:6, a=0,effects["b"], g1=0, g2=0, effects[c('d1','d2')])
# The following table show the parameters in the simulation if 0-1 coding was used
knitr::kable(all.ref, row.names = F)
# We now translate these values for the case when -1,1 coding is used
all.eff <- data.frame(case=1:6,
                      a = all.ref["a"]+all.ref["b"]/2+all.ref["g1"]/2+all.ref["d1"]/4,
                      b = all.ref["b"]/2+all.ref["d1"]/4,
                      g1 = all.ref["g1"]/2+all.ref["d1"]/4,
                      g2 = sqrt(3)*(all.ref["g2"]+all.ref["d2"]/2),
                      d1 = all.ref["d1"]/4,
                      d2 = sqrt(3)*all.ref["d2"]/2)
knitr::kable(all.eff, row.names = F)

effects = all.eff
param = "EFFECT"
verbose = T

set.seed(seed.)


npergroup = npergroup.[j]
case = i
nsim = nsim
effects = all.eff
case = i
npergroup = npergroup
mc.cores = mc.cores
verbose = FALSE
n_biom = n_biom
param = "EFFECT"
lambda = "lagrange"
alpha = alpha.
test = T
perturb_frac = c(0.2, 0.8)
pite.ci = FALSE
n.pite = NULL
n.test = 1
tol.beta = 1e-5
lam_frac=.75

nfolds = 10
typeCorr = "I"
rho = 0

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
    types <- c(2,2)
} else {
    types <- c(1,3)
}
parameters <- MakeParameters(prognosis = c(effects[case,"g1"],effects[case,"g2"]),
                               predictive = c(effects[case,"d1"],effects[case,"d2"]),
                               types = types,
                               prevalence = c(0.5,0),
                               means = c(0,0),
                               stddev = c(1,1),
                               input=input)
MakeInputTable(input)
MakeParametersTable(input, parameters)


seed. = readRDS(
  sprintf("sim/normal/nbiom1000/results/data/nsim100-case%s-npergroup%s-seed.rds",
          case, npergroup))[[sim.]]







set.seed(seed.)
dataset <- OneData(input, parameters, standardize = FALSE, param = param)
dataset[1:10,1:10]
# tol.beta = 0.1
# tol.beta = 1e-5
dataset.lasso <- score.lasso(dataset,
                             input, infinity=T,
                             parameters = parameters,
                             alpha = alpha,
                             verbose = verbose,
                             lambda = lambda,
                             nfolds = nfolds, gridrange_ = 2000,
                             lam_frac=lam_frac, tol.beta = tol.beta)

output. = readRDS(
  sprintf("sim/normal/nbiom1000/results/data/nsim100-case%s-npergroup%s-coef.rds",
          case, npergroup))
output. %>%
  filter(method=="lasso") %>%
  filter(sim == sim.) %>%
  filter(estimate != 0)
dataset.lasso$Lasso.output %>%
  filter(estimate != 0)
dataset.lasso$ML.output.M %>%
  filter(estimate != 0)
lam<-dataset.lasso$bestlam*dataset.lasso$N
      dataset.lasso.an1 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac=perturb_frac[1],
                                                   pite.ci = pite.ci, alpha = alpha,
                                                   n.pite=n.pite)
      dataset.lasso.an2 <- score.lasso.added.noise(dataset,
                                                   input,
                                                   parameters = parameters,
                                                   verbose = verbose,
                                                   lambda = lam,
                                                   perturb_frac=perturb_frac[2],
                                                   pite.ci = pite.ci, alpha = alpha,
                                                   n.pite=n.pite)

  dataset.lasso.an1$Lasso.output %>%
    filter(estimate != 0)
  dataset.lasso.an2$Lasso.output %>%
    filter(estimate != 0)

      if(is.null(n.pite)){
        n<-2*npergroup
        n.pite<-2*npergroup
      } else {
        n<-n.pite
      }
      results <- rbind(cbind(method = "lasso",summarize_scores(dataset.lasso$scores,dataset,  n.pite=n.pite), nvars=dataset.lasso$nvars),
        cbind(method = "an1",  summarize_scores(dataset.lasso.an1$scores, dataset, n.pite=n.pite), nvars=dataset.lasso.an1$nvars),
        cbind(method = "an2",  summarize_scores(dataset.lasso.an2$scores, dataset, n.pite=n.pite), nvars=dataset.lasso.an2$nvars))

      results.coef <- rbind(
        cbind(method = "lasso",dataset.lasso$Lasso.output, nvars=dataset.lasso$nvars),
        cbind(method = "an1",  dataset.lasso.an1$Lasso.output, nvars=dataset.lasso.an1$nvars),
        cbind(method = "an2",  dataset.lasso.an2$Lasso.output, nvars=dataset.lasso.an2$nvars))

      input.test <- input
        if(n.test == 1){
          dataset.test <- OneSubject(input.test, parameters, standardize = FALSE, param = param)
        } else {
          if(is.null(n.test)){
            n.test<-2*npergroup
          }
          input.test$nT <- n.test/2+1
          input.test$nC <- n.test/2+1
          dataset.test <- OneData(input.test, parameters, standardize = FALSE, param = param)
        }
#   dataset.test. = readRDS(
#     sprintf("sim/normal/nbiom1000/results/data/nsim100-case%s-npergroup%s-dataset-test.rds",
#           case, npergroup))
dataset.test[1, 1:10]
# dataset.test.[sim., 1:10]
        dataset.test.lasso <- confidence.intervals.lasso.test(dataset = dataset.test,
                                                              input = input.test,
                                                              parameters = parameters,
                                                              lasso.results = dataset.lasso,
                                                              gridrange_ = 250,
                                                              alpha = alpha, tol.beta = tol.beta)
        dataset.test.lasso.an1 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input.test,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an1,
                                                                    ndraw=10000,
                                                                    burnin=5000,
                                                                    alpha = alpha)
        dataset.test.lasso.an2 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                                    input = input.test,
                                                                    parameters = parameters,
                                                                    lasso.results = dataset.lasso.an2,
                                                                    ndraw=10000,
                                                                    burnin=5000,
                                                                    alpha = alpha)

results.test <- rbind(cbind(method = "lasso",summarize_scores_many(dataset.test.lasso$scores,dataset.test, n.pite=n.test)),
          cbind(method = "an1",  summarize_scores_many(dataset.test.lasso.an1$scores, dataset.test, n.pite=n.test)),
          cbind(method = "an2",  summarize_scores_many(dataset.test.lasso.an2$scores, dataset.test, n.pite=n.test)))
results.test
results.test. = readRDS(
    sprintf("sim/normal/nbiom1000/results/data/nsim100-case%s-npergroup%s-test.rds",
          case, npergroup))
results.test. %>%
  filter(sim == sim.)

###############################################################################-
###############################################################################-

sessionInfo(package = NULL)
