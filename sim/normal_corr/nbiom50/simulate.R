###############################################################################-
## This program is used to run simulations and generate a summary
##
## Ideally, this is run in batch mode in a clean session via
## nohup R CMD BATCH --vanilla simulate.R
##
## If running in interactive mode, Control+Shift+F10 to start a new session
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

# Number of biomarkers in this simulation
n_biom    =   50
# Set onlyPlots to TRUE if simulations already run, and only plots want to be created
onlyPlots = FALSE

# Set the environment to store results by creating relative paths --------------
rproj <- rprojroot::find_package_root_file()
path  <- paste0(rproj,sprintf("/sim/normal_corr/nbiom%s/",n_biom))
if (!dir.exists(paste0(path, "results/")))       dir.create(paste0(path, "results/"))
if (!dir.exists(paste0(path, "results/data/")))  dir.create(paste0(path, "results/data/"))
if (!dir.exists(paste0(path, "results/plots/"))) dir.create(paste0(path, "results/plots/"))
path.data <- paste0(path, "results/data/")
path.plot <- paste0(path, "results/plots/")

# Set initial parameters -------------------------------------------------------
alpha.   =    0.05
seed.    = 35613
nsim     = 1000
mc.cores =   30 # Number of cores for simulations

sample.size = c(40, 100, 220, 350) # Sample sizes considered in simulations
npergroup.  = sample.size / 2

# Calculate effecs to be used in simulations
# First we suppose effects under the 0,1 coding
effects <- expand.grid(b = c(0, 1 / 4), d1 = c(0, 1 / 2),
                       d2 = c(0 , 1 / 3))[c(1, 5, 3, 7, 2, 8),]
all.ref <- cbind(case = 1:6, a = 0, effects["b"], g1 = 0, g2 = 0, effects[c('d1', 'd2')])
# The following table show the parameters in the simulation if 0-1 coding was used
knitr::kable(all.ref, row.names = F)
# We now translate these values for the case when -1,1 coding is used
all.eff <- data.frame(case = 1:6,
                      a  = all.ref["a"] + all.ref["b"]/2 + all.ref["g1"] / 2 + all.ref["d1"] / 4,
                      b  = all.ref["b"] / 2 + all.ref["d1"] / 4,
                      g1 = all.ref["g1"] / 2 + all.ref["d1"] / 4,
                      g2 = sqrt(3) * (all.ref["g2"] + all.ref["d2"] / 2),
                      d1 = all.ref["d1"] / 4,
                      d2 = sqrt(3) * all.ref["d2"] / 2)
knitr::kable(all.eff, row.names = F)

effects = all.eff
param   = "EFFECT"
verbose = T

## Now Simulate! -------------------------------------------------------------------
# We now run the scenarios by looping across sample sizes and parameter cases
# To try in one scenario, just run j=1;i=1 and the code inside the loop
set.seed(seed.)
if (!onlyPlots){
  for (j in c(4, 3, 2, 1)) {
    for (i in 1:6) {
      npergroup = npergroup.[j]
      case = i
      listdata = simulateNtrial(nsim = nsim,
                                effects = all.eff,
                                case = i,
                                npergroup = npergroup,
                                mc.cores = mc.cores,
                                verbose = FALSE,
                                n_biom = n_biom,
                                typeCorr = "AR1", rho = 0.9,
                                param = "EFFECT",
                                lambda = "lagrange",
                                alpha = alpha.,
                                test = T,
                                perturb_frac = c(0.2, 0.8),
                                pite.ci=FALSE,
                                n.pite=NULL,
                                n.test = 1,
                                lam_frac=0.5)

      filename <- paste0("results/data/nsim", nsim, "-case", case, "-npergroup", npergroup)
      cat(print(pryr::mem_used()))
      # For a better use of memory, we save the results and delete them from the
      # workspace. Then they are loaded again to create the plots
      saveRDS(listdata, file = paste0(filename, ".rds"))
      listdata1 <- do.call(rbind, lapply(listdata$dataset, function(x) (x[[1]])))
      listdata2 <- do.call(rbind, lapply(listdata$dataset, function(x) (x[[2]])))
      listdata3 <- do.call(rbind, lapply(listdata$dataset, function(x) (x[[3]])))
      listdata4 <- do.call(rbind, lapply(listdata$dataset, function(x) (x[[4]])))
      listdata2 %>%
        dplyr::group_by(method) %>%
        dplyr::summarize_all(.funs = function(x) mean(x, na.rm=TRUE)) %>%
        print()
      saveRDS(listdata1, file = paste0(filename,"-train.rds"))
      saveRDS(listdata2, file = paste0(filename,"-test.rds"))
      saveRDS(listdata3, file = paste0(filename,"-coef.rds"))
      saveRDS(listdata4, file = paste0(filename,"-tailarea.rds"))
      rm(listdata)
      rm(listdata1)
      rm(listdata2)
      rm(listdata3)
      rm(listdata4)
      rm(filename)
      cat("\n\n\nCOMPLETED: case (i)",case,"npergroup (j)",npergroup,"- \n\n\n")
      warnings()
      cat(print(pryr::mem_used()))
    }
  }
}

###############################################################################-
## Produce Plots and summary ---------------------------------------------------
files.2 <- list.files(path=path.data) # All results files
methods <- data.frame(method = c("null","lm","lasso","mlm","sch","an1", 'an2'),
                      Method = factor(c("ATE", "full", "Lasso","reduced", "reduced-Scheffe", "rLasso-1", "rLasso-2"),
                                      levels=c("ATE", "full", "Lasso","reduced", "reduced-Scheffe", "rLasso-1", "rLasso-2")))
myPal <- c("black", "#EE4035", "#F3A530", "#56B949", "#82b6cc", "#0f5e7f", "#024460")
myShape <- c(NA, 16, 17, 15, 3, 12, 13)
names(myPal) <- methods$Method
names(myShape) <- methods$Method
colorize.fill <- scale_fill_manual(name = "Method",  values = myPal)
colorize      <- scale_color_manual(name = "Method", values = myPal)
outlier.size = -1

##### Results at the Coefficients level ---------------------------------------
## The following results analyze selection of biomarkers, coverage of confidence
## intervals for the coefficients in the model, and their width

files.2.coef <- files.2[grep(files.2, pattern = "coef")]
results.coeficients <- data.frame()
for(i in files.2.coef){
  listdata3 = readRDS(paste0(path.data,i))
  listdata3 %>% filter(term %in% c("treatment", "mkr1", "mkr2",
                                   "treatment:mkr2",
                                   "treatment:mkr1",
                                   "mkr1:treatment",
                                   "mkr2:treatment")) -> listdata3
  results.coeficients <- rbind(results.coeficients,
                               data.frame(listdata3, scenario = i))
}

results.coeficients %>%
  left_join(methods, by = "method") -> results.coeficients

pos = do.call(rbind, gregexpr('-', results.coeficients$scenario))
results.coeficients$case <- as.numeric(substr(results.coeficients$scenario,
                                              start = pos[, 2] - 1, stop = pos[, 2] - 1))

results.coeficients$n    <- 2 * as.numeric(substr(results.coeficients$scenario,
                                                  start = pos[, 2] + 10, stop = pos[, 3] - 1))

results.coeficients <- results.coeficients %>%
  mutate(term = ifelse(term=="mkr1:treatment", "treatment:mkr1",
                       ifelse(term=="mkr2:treatment", "treatment:mkr2", term)))

##### *Plots by terms -----------------------------------------------------------
results.coeficients %>%
  group_by(n, case, Method, term) %>%
  summarize(estimate = mean(estimate)) %>%
  arrange(case, n) -> results.coeficients2
unique(results.coeficients2$term)

# Extract the true values
tmp.dat <- data.frame(variable=c("a", "b", "g1", "g2", "d1", "d2"),
                      term=c("(Intercept)", "treatment",
                             "mkr1", "mkr2",
                             "treatment:mkr1", "treatment:mkr2"))
reshape2::melt(effects, "case") %>%
  left_join(tmp.dat) -> tmp.dat2

# Coverage
results.coeficients %>%
  left_join(tmp.dat2[,-2]) %>%
  mutate(cover = 1*(LowConfPt <= value & value <= UpConfPt),
         width = -(LowConfPt - UpConfPt),
         bias = estimate - value,
         bias2 = (estimate - value)^2) -> res.coef
res.coef %>%
  group_by(n,case, Method, term) %>%
  mutate(sel=estimate!=0) %>%
  summarize(estimate = mean(estimate, na.rm=T),
            cover = mean(cover, na.rm=T),
            width = mean(width, na.rm=T),
            bias = mean(bias, na.rm=T),
            mse = mean(bias2, na.rm=T),
            value = mean(value, na.rm=T),
            sel = mean(sel) * 100) %>%
  arrange(case,n) -> res.coef.summary
res.coef.summary




## *Selection probabilities -----------------------------------------------------
head(res.coef.summary)
res.coef.summary %>%
  filter(Method!="full" & Method!="reduced")%>%
  ggplot() +
  geom_bar(aes(x=term, y=sel, fill=Method),
           stat="identity", position="dodge", color="black", size=0.1) +
  colorize.fill + colorize +
  ylab("Percentage of selection") +
  xlab("Term") +
  facet_grid(case ~ n, labeller = label_both) +
  theme_bw() +
  theme(axis.text.x =  element_text(angle=45, hjust = 1)) +
  ylim(0,100) -> plot.sel
# plot.sel

ggsave(filename = paste0(path.plot, "/percent_selected.pdf"),
       plot = plot.sel,
       width = 16, height = 20, units = "cm")


## *coverage of Confidence intervals --------------------------------------------
p = 0.05
q = 1 - p
e = 1.96 * sqrt(p * q / nsim)
head(res.coef.summary)
res.coef.summary %>%
  filter(term!="(Intercept)")%>%
  ggplot() +
  geom_bar(aes(x=term, y=100-cover*100, fill=Method),
           stat="identity", position="dodge", color="black", size=0.1) +
  colorize.fill + colorize +
  xlab("Term") +
  ylab("Miscoverage") +
  facet_grid(case ~ n, labeller = label_both) +
  theme_bw() +
  theme(axis.text.x =  element_text(angle=45, hjust = 1)) +
  geom_hline(yintercept=5+e*100, linetype=2) +
  geom_hline(yintercept=5) +
  geom_hline(yintercept=5-e*100, linetype=2) +
  ylim(0,30) -> plot.cover
# plot.cover

ggsave(filename = paste0(path.plot,"/cover.pdf"),
       plot = plot.cover,
       width = 16, height = 20, units = "cm")

## *Width of Confidence intervals -----------------------------------------------
head(res.coef)
res.coef %>%
  ggplot() +
  geom_boxplot(aes(x=term, y=width, fill=Method),
               position="dodge", size=0.2, outlier.size = outlier.size) +
  colorize.fill + colorize +
  ylab("Width of the confidence interval") +
  xlab("Term") +
  coord_cartesian(ylim=c(0,7.5))+
  facet_grid(case ~ n, labeller = label_both) +
  theme_bw() +
  theme(axis.text.x =  element_text(angle=45, hjust = 1))   -> plot.width
# plot.width


ggsave(filename = paste0(path.plot,"/width.pdf"),
       plot = plot.width,
       width = 16, height = 20, units = "cm")





###############################################################################-
###############################################################################-
##  Train data results ---------------------------------------------------------
files.2 <- list.files(path = path.data)
files.2.train <- files.2[grep(files.2, pattern = "-train")]
results.train <- data.frame()
for(i in files.2.train){
  listdata1 = readRDS(paste0(path.data,i))
  results.train <- rbind(results.train, data.frame(listdata1, scenario=i))
}

pos = do.call(rbind,gregexpr('-', results.train$scenario))
results.train$case <- as.numeric(substr(results.train$scenario,
                                        start = pos[,2] - 1, stop = pos[,2] - 1))

results.train$n <- 2*as.numeric(substr(results.train$scenario,
                                       start = pos[,2] + 10, stop = pos[,3] - 1))

results.train %>%
  group_by(n,case, method) %>%
  summarize(cover = mean(cover)) %>%
  arrange(case,n) -> results.train2
results.train2

results.train %>%
  left_join(methods, by = "method") %>%
  mutate(method = Method)-> results.train
results.train$Method


## *Number of variables in the score. Nvars -------------------------------------
head(results.train)
results.train %>%
  ggplot() +
  geom_boxplot(aes(x=as.factor(n), y=nvars, fill=method),
               position="dodge", size=0.2, outlier.size = outlier.size) +
  colorize.fill + colorize +
  ylab(expression("Number of variables in score")) +
  xlab("N") +
  facet_grid( ~ case, labeller = label_both) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(c(1.5), "lines"),
        strip.text.x = element_text(face = "bold"),
        strip.background = element_blank(),
        axis.text.x =  element_text(angle=45, hjust = 1)) -> plot.nvars.pite
# plot.nvars.pite

ggsave(filename = paste0(path.plot,"/train-pite-nvars.pdf"),
       plot = plot.nvars.pite,
       width = 16, height = 10, units = "cm")


## *Number of variables in the score = 0. Nvars = 0 -------
head(results.train)
results.train %>%
  group_by(case,n, method) %>%
  summarise(null = sum(nvars == 0)/1000*100) %>%
  ggplot() +
  geom_bar(aes(x=as.factor(n), y=null, fill=method),
           position="dodge", size=0.2,
           stat="identity") +
  colorize.fill + colorize +
  ylab(expression("Null model")) +
  xlab("N") +
  facet_grid( ~ case, labeller = label_both) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(c(1.5), "lines"),
        strip.text.x = element_text(face = "bold"),
        strip.background = element_blank(),
        axis.text.x =  element_text(angle=45, hjust = 1)) -> plot.null.pite
# plot.null.pite

ggsave(filename = paste0(path.plot,"/train-pite-null.pdf"),
       plot = plot.null.pite,
       width = 16, height = 10, units = "cm")



## *Function for plotting the Sensitivity and specificity plot-------
plotsens <- function(results, case.=4){
  sample.size <- c(unique(results$n))
  results %>%
    group_by(case,n, method) %>%
    filter(case==case.) %>%
    summarise(sensitivity.ll = mean(sensitivity.ll, na.rm = T),
              sensitivity.dx = mean(sensitivity.dx, na.rm = T),
              sensitivity.ul = mean(sensitivity.ul, na.rm = T),
              specificity.ll = mean(specificity.ll, na.rm = T),
              specificity.dx = mean(specificity.dx, na.rm = T),
              specificity.ul = mean(specificity.ul, na.rm = T)) %>%
    reshape2::melt(c("case","n","method")) %>%
    mutate(which = factor(substr(variable,start = 13,stop = 14),
                          levels = c("ll","dx","ul")),
           measure = factor(substr(variable,start = 1,stop = 11),
                            levels = c("sensitivity", "specificity"))) %>%
    rename(sensitivity = value) -> datplot
  levels(datplot$which) = c("bold(hat(B)[l])", "bold(hat(B))", "bold(hat(B)[u])")
  datplot %>%
    filter(method %in% c("ATE", "full", "Lasso", "reduced-Scheffe", "rLasso-1", "rLasso-2")) %>%
    ggplot() +
    geom_line(aes(y=sensitivity,  color = method, x=n), linetype  = 2) +
    geom_point(aes(y=sensitivity,  color = method, shape=method, x=n)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8)) +
    colorize.fill + colorize +
    scale_shape_manual(name = "Method", values = myShape) +
    scale_y_continuous(name = "", limits = c(0,1)) +
    scale_x_continuous(name = "n",
                       limits = c(0,350), breaks = c(0,sample.size)) -> p1
  p1 + facet_grid(measure ~ which, switch = "y", labeller = label_parsed) +
    theme(strip.background = element_blank(),
          panel.grid = element_blank(),
          strip.text   = element_text(size = 7),
          strip.text.y = element_text(size = 7),
          strip.text.x = element_text(size = 9),
          legend.text  = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.title.align = 0.5,
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 7),
          axis.text    = element_text(size = 7, color = "black"),
          axis.text.x  = element_text(size = 7, color = "black"),
          strip.placement = "outside")
}

##  Test data results ---------------------------------------------------------
files.2.test <- files.2[grep(files.2, pattern = "-test")]
results.test <- data.frame()
for(i in files.2.test){
  listdata2 = readRDS(paste0(path.data,i))
  results.test <- rbind(results.test, data.frame(listdata2, scenario=i))
}

# results.test
pos = do.call(rbind,gregexpr('-', results.test$scenario))
results.test$case <- as.numeric(substr(results.test$scenario,
                                       start = pos[,2] - 1, stop = pos[,2] - 1))

results.test$n <- 2*as.numeric(substr(results.test$scenario,
                                      start = pos[,2] + 10, stop = pos[,3] - 1))

p = n_biom*2+2
results.test %>%
  filter(!((method == "sch"|method == "lm") & n < p)) -> results.test

results.test %>%
  left_join(methods, by = "method") %>%
  mutate(method = Method)-> results.test

results.test %>% head()

results.test %>%
  group_by(n,case, method) %>%
  summarize(cover = mean(cover)) %>%
  arrange(case,n) -> results.test2
results.test2

## *PITE width ------------------------------------------------------------------
head(results.test2)
results.test %>%
  ggplot() +
  geom_boxplot(aes(x=as.factor(n), y=w, fill=method), position="dodge",
               size=0.2, outlier.size = outlier.size) +
  colorize.fill + colorize +
  ylab("Width") +
  xlab("n") +
  facet_wrap( ~ case, labeller = label_both, nrow = 2) +
  theme_bw() +
  coord_cartesian(ylim=c(0,20)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(c(1.5), "lines"),
        strip.text.x = element_text(face = "bold"),
        strip.background = element_blank(),
        axis.text.x =  element_text(angle=45, hjust = 1)) -> plot.width.pite

# plot.width.pite
ggsave(filename = paste0(path.plot,"/test-pite-width.pdf"),
       plot = plot.width.pite,
       width = 16, height = 20, units = "cm")

## *PITE coverage ---------------------------------------------------------------
results.test %>%
  group_by(n,case, method) %>%
  summarize(coverage = mean(cover, na.rm=T)) %>%
  arrange(case,n) -> data.coverage.pite
head(data.coverage.pite)
nsim=1000
p=0.05
q=1-p
e=1.96*sqrt(p*q/nsim)
data.coverage.pite %>%
  ggplot() +
  geom_bar(aes(x=as.factor(n), y=100-coverage*100, fill=method),
           position="dodge", stat = "identity", color="black", size=0.2) +
  colorize.fill + colorize +
  ylab("Miscoverage") +
  xlab("N") +
  facet_wrap( ~ case, labeller = label_both, nrow = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(c(1.5), "lines"),
        strip.text.x = element_text(face = "bold"),
        strip.background = element_blank(),
        axis.text.x =  element_text(angle=45, hjust = 1))+
  geom_hline(yintercept=5+e*100, linetype=2) +
  geom_hline(yintercept=5) +
  geom_hline(yintercept=5-e*100, linetype=2) +
  ylim(0,40)  -> plot.coverage
# plot.coverage

ggsave(filename = paste0(path.plot,"/test-pite-coverage-bar.pdf"),
       plot = plot.coverage,
       width = 16, height = 20, units = "cm")


## *PITE BIAS TRUE---------------------------------------------------------------
results.test %>%
  group_by(n,case, method) %>%
  arrange(case,n) -> data.bias.pite

results.test %>%
  group_by(n,case, method) %>%
  summarize(bias = mean(bias.true, na.rm=T)) %>%
  arrange(case,n) -> data.bias.pite.mean
outlier.size=-1
head(data.bias.pite)
data.bias.pite %>%
  ggplot() +
  geom_boxplot(aes(x=as.factor(n), y=bias.true, fill=method),
               position="dodge", size=0.2, outlier.size = outlier.size) +
  geom_hline(yintercept = 0, color="red") +
  geom_point(data = data.bias.pite.mean,
             aes(x=as.factor(n), y=bias, fill=method),
             position=position_dodge(width=0.75), size= 0.5)+
  colorize.fill + colorize +
  ylab("Bias") +
  xlab("N") +
  facet_wrap( ~ case, labeller = label_both, nrow = 2) +
  theme_bw() +
  coord_cartesian(ylim=c(-2.5,2.5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(c(1.5), "lines"),
        strip.text.x = element_text(face = "bold"),
        strip.background = element_blank(),
        axis.text.x =  element_text(angle=45, hjust = 1)) -> plot.bias.pite.true
# plot.bias.pite.true

ggsave(filename = paste0(path.plot,"/test-pite-bias-true.pdf"),
       plot = plot.bias.pite.true,
       width = 16, height = 20, units = "cm")

## *PITE MSE TRUE----------------------------------------------------------------
head(results.test)
results.test %>%
  ggplot() +
  geom_boxplot(aes(x=as.factor(n), y=sqrt(mse.true), fill=method),
               position="dodge", size=0.2, outlier.size = outlier.size) +
  colorize.fill + colorize +
  ylab(expression(sqrt(MSE))) +
  xlab("N") +
  facet_wrap( ~ case, labeller = label_both, nrow = 2) +
  theme_bw() +
  coord_cartesian(ylim=c(0,7.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(c(1.5), "lines"),
        strip.text.x = element_text(face = "bold"),
        strip.background = element_blank(),
        axis.text.x =  element_text(angle=45, hjust = 1)) -> plot.mse.pite.true
# plot.mse.pite.true

ggsave(filename = paste0(path.plot,"/test-pite-mse-true.pdf"),
       plot = plot.mse.pite.true,
       width = 16, height = 20, units = "cm")

# Sensitivity and specificity --------------------------------------------------
plotsens(results.test, case. = 1) -> ss1
ggsave(paste0(path.plot,"/test-ss-1.pdf"), plot = ss1, width = 17, height = 8, units = "cm")
plotsens(results.test, case. = 2) -> ss2
ggsave(paste0(path.plot,"/test-ss-2.pdf"), plot = ss2, width = 17, height = 8, units = "cm")
plotsens(results.test, case. = 3) -> ss3
ggsave(paste0(path.plot,"/test-ss-3.pdf"), plot = ss3, width = 17, height = 8, units = "cm")
plotsens(results.test, case. = 4) -> ss4
ggsave(paste0(path.plot,"/test-ss-4.pdf"), plot = ss4, width = 17, height = 8, units = "cm")
plotsens(results.test, case. = 5) -> ss5
ggsave(paste0(path.plot,"/test-ss-5.pdf"), plot = ss5, width = 17, height = 8, units = "cm")
plotsens(results.test, case. = 6) -> ss6
ggsave(paste0(path.plot,"/test-ss-6.pdf"), plot = ss6, width = 17, height = 8, units = "cm")


## Save all plots in only one document -----------------------------------------
pdf(file= paste0(path.plot,"/all-plots.pdf"), width = 8, height = 5)
ss1 + ggtitle("Sensitivity and Specificity. Case 1")
ss2 + ggtitle("Sensitivity and Specificity. Case 2")
ss3 + ggtitle("Sensitivity and Specificity. Case 3")
ss4 + ggtitle("Sensitivity and Specificity. Case 4")
ss5 + ggtitle("Sensitivity and Specificity. Case 5")
ss6 + ggtitle("Sensitivity and Specificity. Case 6")
plot.coverage       + ggtitle("Coverage of Confidence Intervals")
plot.width.pite     + ggtitle("Width of Confidence Intervals for PITE")
plot.bias.pite.true + ggtitle("BIAS for PITE")
plot.mse.pite.true  + ggtitle("sqrt(MSE) for PITE")
plot.sel            + ggtitle("Percent of inclusion by variable in the Lasso")
dev.off()

# Print session info
sessionInfo(package = NULL)
