# Make sure we are working in a clean workspace for reproducibility
rm(list=ls())
cat("\014") # Clear console
# Load libraries
library(parallel)
library(PITE)
library(survival)
library(glmnet)
library(dplyr)
library(ggplot2)

# Define simulation parameters
seed     = 285641
lam_frac =      0.50
alpha    =      0.05
delta1   =     -0.08
delta2   =      0.18
beta     =     -0.10
nsim     =   1000
mc.cores =     35 # Number of cores for parallel computing

# This file should be run in BATCH mode using the following line:
# nohup R CMD BATCH "--args 2000" --vanilla 8_SurvData.R
# The args argument correspond to the number of subjects for the simulation
# If it is not specified, we set it to 100
# Specify number of subjects in the dataset
if (length(commandArgs(TRUE))==0){
  # If running in interactive mode, then we should use
  args <- list(nsub=100)
} else{
# reading arguments
  args <- commandArgs(TRUE)
}

# Print args to check values
nsub  <- as.numeric(args[1])
print(args)
print(nsub)

# Create directories for storing results if they dont exist
if (!dir.exists("results")){
  dir.create("results")
}
if (!dir.exists("results/data")){
  dir.create("results/data")
}
if (!dir.exists("results/plots")){
  dir.create("results/plots")
}



### Test the function SurvOneData ---------------------------------------------
# alpha=0; beta=0; gamma=0; delta=0.10; nsub=10; sigma = 1
dataset <- SurvOneData(alpha=0, beta=0,  delta1=-0.2, delta2=0.1,
                       nsub=1000, sigma = 1, maxtime = 2)

res.null  <- score.cox.null(dataset)
# res.null$reg.out
res       <- score.cox(dataset)
res.lasso <- score.cox.lasso(dataset)
res.reduced <- score.cox.reduced(dataset, vars = res.lasso$vars)

dataset <- SurvOneSubject(alpha=0, beta=0,  delta1=-0.2, delta2=0.1,
                          sigma = 1, maxtime = 2)
ci.cox.null(dataset, reg.results = res.null$reg.results)$scores
ci.cox(dataset, reg.results = res$reg.results)$scores
ci.cox.lasso(dataset, reg.results = res.lasso$reg.results, reg.out = res.lasso$reg.out)$scores
ci.cox.reduced(dataset, reg.results = res.reduced$reg.results, reg.out = res.reduced$reg.out,
               vars = res.lasso$vars)$scores


# results = result.null
simplifyResults <- function(results){
  mean(unlist(lapply(results, function(x) {
    scores <- x$scores
    mean(scores$cover, na.rm=T)
  }))) -> coverage
  mean(unlist(lapply(results, function(x) {
    scores <- x$scores
    mean(scores$bias, na.rm=T)
  }))) -> bias
  mean(unlist(lapply(results, function(x) {
    scores <- x$scores
    scores$width <- ifelse(!is.finite(scores$width), NA, scores$width)
    mean(scores$width, na.rm=T)
  }))) -> width
  mean(unlist(lapply(results, function(x) {
    scores <- x$scores
    mean(scores$bias^2, na.rm=T)
  }))) -> mse
  c(coverage = coverage, width = width, bias = bias, mse = mse)
}


## Simulate data ------
set.seed(1234)
#
# delta <- -0.2
# nsub  <- 1000
data <- SurvOneData(alpha=1,
                    beta=beta,  delta1=delta1, delta2=delta2,
                    nsub=1000, sigma = 1)
alpha = 0.05
s.n <- score.cox.null(data, alpha = alpha, calculate.ci = FALSE)
s.c <- score.cox(data, alpha = alpha, calculate.ci = FALSE)
s.l <- score.cox.lasso(data, lam_frac = lam_frac, alpha = alpha, calculate.ci = FALSE)
s.r <- score.cox.reduced(data, vars = s.l$reg.out$term[which(s.l$reg.out$estimate!=0)], alpha = alpha, calculate.ci = FALSE)

s.n$reg.out
s.c$reg.out
s.l$reg.out
s.r$reg.out

subject <- SurvOneSubject(alpha=1,
                          beta=beta,  delta1=delta1, delta2=delta2,
                          sigma = 1)
c.n <- ci.cox.null(subject,
            reg.results = s.n$reg.results)$scores
c.c <- ci.cox(subject,
            reg.results = s.c$reg.results)$scores
c.l <- ci.cox.lasso(subject,
                    reg.results = s.l$reg.results,
                    reg.out = s.l$reg.out)$scores
c.r <- ci.cox.reduced(subject,
            reg.results = s.r$reg.results, reg.out = s.r$reg.out,
            vars = s.l$vars)$scores
c.n
c.c
c.l
c.r




set.seed(seed)
data.sim    <- SurvData(alpha=1,
                     beta=beta,
                     delta1=delta1, delta2=delta2,
                     nsim=nsim,
                     nsub=nsub, sigma = 1)
subject.sim <- SurvSubjects(alpha=1,
                     beta=beta,
                     delta1=delta1, delta2=delta2,
                     nsim=nsim,
                     sigma = 1)

result.null     <- mclapply(1:nsim, function(x){
    s.n = score.cox.null(data.sim[[x]], alpha = alpha, calculate.ci = FALSE)
    ci.cox.null(subject.sim[[x]], alpha = alpha,
                reg.results = s.n$reg.results)
  },
  mc.preschedule = FALSE,
  mc.cores = mc.cores)

result.cox     <- mclapply(1:nsim, function(x){
    s.c = score.cox(data.sim[[x]], alpha = alpha, calculate.ci = FALSE)
    ci.cox(subject.sim[[x]], alpha = alpha,
                reg.results = s.c$reg.results)
  },
  mc.preschedule = FALSE,
  mc.cores = mc.cores)


result.lasso   <- mclapply(1:nsim, function(x){
    s.l = score.cox.lasso(data.sim[[x]], lam_frac = lam_frac, alpha = alpha,  calculate.ci = FALSE)
    ci.cox.lasso(subject.sim[[x]], alpha = alpha,
                 reg.results = s.l$reg.results, reg.out = s.l$reg.out)
  },
  mc.preschedule = FALSE,
  mc.cores = mc.cores)

result.reduced <- mclapply(1:nsim, function(x){
    s.r = score.cox.reduced(data.sim[[x]], vars = result.lasso[[x]]$vars, alpha = alpha,  calculate.ci = FALSE)
    ci.cox.reduced(subject.sim[[x]], alpha = alpha,
                          reg.results = s.r$reg.results, reg.out = s.r$reg.out,
                          vars = result.lasso[[x]]$vars)
  },
  mc.preschedule = FALSE,
  mc.cores = mc.cores)

save(result.null,    file = paste0("results/data/result.null",nsub,".Rda"))
save(result.cox,     file = paste0("results/data/result.cox",nsub,".Rda"))
save(result.lasso,   file = paste0("results/data/result.lasso",nsub,".Rda"))
save(result.reduced, file = paste0("results/data/result.reduced",nsub,".Rda"))

simplifyResults(result.null)
simplifyResults(result.cox)
simplifyResults(result.lasso)
simplifyResults(result.reduced)

#------------------------------------------------------------------------------#
## Plot coefficients -----------------------------------------------------------
pdf(file = paste0("results/plots/n", nsub , "alpha", alpha, ".pdf"))
do.call(rbind,
        lapply(1:length(result.null), function(x){
          data.frame(result.null[[x]]$reg.out, sim = x, method="cox-null")})
) -> n.o
do.call(rbind,
        lapply(1:length(result.lasso), function(x){
          data.frame(result.lasso[[x]]$reg.out, sim = x, method="cox-lasso")})
) -> l.o
do.call(rbind,lapply(1:length(result.cox), function(x){
  data.frame(result.cox[[x]]$reg.out, sim = x, method="cox")
})) -> c.o
do.call(rbind,lapply(1:length(result.reduced), function(x){
  data.frame(result.reduced[[x]]$reg.out, sim = x, method="cox-reduced")
})) -> r.o
a.o <- rbind(l.o, c.o, r.o)
a.o$method <- factor(a.o$method, levels = c("cox", "cox-reduced", "cox-lasso"))
aa.o <- rbind(n.o, l.o, c.o, r.o)
aa.o$method <- factor(aa.o$method, levels = c("cox-null", "cox", "cox-reduced", "cox-lasso"))

head(a.o)
a.o[1:13, "term"]

true <- c(beta,
          0.17, 0.18, 0.11, 0.14, 0.21, -0.12,
          delta1, delta2,
          rep(0,4))

a.o %>%
  mutate(cover = 1*((conf.low < true)&(true < conf.high))) %>%
  group_by(term, method) %>%
  summarize(cover = mean(cover, na.rm=T),
            mean.est = mean(estimate, na.rm=T)) %>%
  ggplot() +
  ylim(c(0,1)) +
  geom_hline(yintercept= 1 - alpha) +
  geom_bar(aes(x= term, y=cover, fill=method), stat="identity", position="dodge") +
  ggtitle("Coverage of confidence intervals for coefficients")

aa.o %>%
  group_by(term, method) %>%
  ggplot() +
  geom_hline(yintercept=0) +
  geom_boxplot(aes(x= term, y=estimate, fill=method),position="dodge") +
  ggtitle("Estimates for coefficients") +
  coord_cartesian(ylim = c(-1,1))


aa.o %>%
  mutate(included = 1*(estimate != 0 )) %>%
  group_by(term, method) %>%
  summarize(included = mean(included, na.rm=T)) %>%
  ggplot() +
  # ylim(c(0,1)) +
  geom_hline(yintercept=0) +
  geom_bar(aes(x= term, y=included, fill=method),
           stat="identity", position="dodge") +
  ggtitle("Inclusion by terms")

c.o %>%
  filter(term=="z") %>%
  mutate(reject.pos = conf.low >0,
         reject.neg = conf.high<0) %>%
  summarize(reject.pos = mean(reject.pos),
            reject.neg = mean(reject.neg))


#------------------------------------------------------------------------------#
## Plot PITE ------------------------------------------------------------------
do.call(rbind,
        lapply(1:length(result.null), function(x){
          data.frame(result.null[[x]]$scores, sim = x, method="cox-null")})
) -> n.o
do.call(rbind,
        lapply(1:length(result.lasso), function(x){
          data.frame(result.lasso[[x]]$scores, sim = x, method="cox-lasso")})
) -> l.o
do.call(rbind,lapply(1:length(result.cox), function(x){
  data.frame(result.cox[[x]]$scores, sim = x, method="cox")
})) -> c.o
do.call(rbind,lapply(1:length(result.reduced), function(x){
  data.frame(result.reduced[[x]]$scores, sim = x, method="cox-reduced")
})) -> r.o
a.o <- rbind(n.o, l.o, c.o, r.o)
a.o$method <- factor(a.o$method, levels = c("cox-null", "cox", "cox-reduced", "cox-lasso"))
head(a.o)
a.o %>%
  group_by(method) %>%
  summarize(cover = mean(cover, na.rm=T)) %>%
  ggplot() +
  ylim(c(0,1)) +
  geom_hline(yintercept= 1 - alpha) +
  geom_bar(aes(x= method, y=cover, fill=method), stat="identity", position="dodge") +
  ggtitle("Coverage for PITE")
a.o %>%
  group_by(method) %>%
  ggplot() +
  geom_hline(yintercept=0) +
  geom_boxplot(aes(x= 1, y=Dx, fill=method),position="dodge") +
  ggtitle("Estimate for PITE")
a.o %>%
  group_by(method) %>%
  ggplot() +
  geom_hline(yintercept=0) +
  geom_boxplot(aes(x= 1, y=bias, fill=method),position="dodge") +
  ggtitle("BIAS for PITE")
a.o %>%
  mutate(bias2 = bias^2) %>%
  group_by(method, sim) %>%
  summarise(sqrt.mse = sqrt(mean(bias2))) %>%
  ggplot() +
  geom_hline(yintercept=0) +
  geom_boxplot(aes(x= 1, y = sqrt.mse, fill=method),position="dodge") +
  ggtitle("sqrt(MSE) for PITE")

dev.off() # Close pdf device

# Clean everything
rm(list=ls())
