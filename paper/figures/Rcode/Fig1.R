###############################################################################-
# pdf(NULL)

# Erase everything in environment for reproducible results and set seed
cat("\014")
print(pryr::mem_used())
rm(list = ls())
print(pryr::mem_used())
###############################################################################-
# Load required packages -------------------------------------------------------
library(dplyr)
library(ggplot2)
library(MASS)
library(Matrix)
library(matrixcalc)
library(matrixStats)
library(monomvn)
library(mvtnorm)
library(scales)
library(ggpubr)
library(PITE)


set.seed(123444)
###############################################################################-
# We use the dataset from Schnell (2016) https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.12522
# to illustrate the methods proposed in this work.
# The data comes from a clinical trial of an Alzheimer's disease treatment
# developed by AbbVie.
# Data is loaded from the article's supplementary material
# Load Dataset and data transformation -----------------------------------------
data_url = "https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fbiom.12522&attachmentId=164892786"
temp <- tempfile()
download.file(data_url,temp)
data <- read.csv(unz(temp, "alzheimers.csv"))
# Recategorize predictors to the 'PARAM' parametrization
# And rename them to use in our function
data$TREATMENT <- relevel(data$TREATMENT, ref="placebo")
data$CARRIER   <- relevel(data$CARRIER,   ref="NON-CARRIER")
data$treatment <- 2*(data$TREATMENT=="low")-1
data$mkr1      <- scale(data$SEVERITY)
data$mkr2      <- scale(data$AGE)
# data$mkr3      <- 2*(data$SEX=="M")-1
# data$mkr4      <- 2*(data$CARRIER=="CARRIER")-1
data$mkr3      <- scale(data$SEX=="M")
data$mkr4      <- scale(data$CARRIER=="CARRIER")

data$y         <- data$CHANGE
# Create ID Variable
data$TrueTrtEff <- 0
data$ID <- 1:41

# Table for reference of new names
vars <- data.frame(var.new = c("treatment", paste0("mkr",1:4)),
                   var.old = c("treatment", "severity",
                               "age", "sex", "carrier"))
vars

# Create input list with n_biom, this is needed for our function
input <- list(n_biom=4)
n_biom <- input$n_biom  # n_biom
n_biom_extra<-6
N <- sample.size <- nrow(data)     # Total number of subjects
K <- n_biom*2+2     # Number of coefficients in the model
# Add 6 more biomarkers to the dataset
mkrCols <- data.frame(mvrnorm(N,rep(0,n_biom_extra),diag(1,n_biom_extra)))
names(mkrCols) <- paste0("mkr",(4+1):(4+n_biom_extra))
input$n_biom <- n_biom <- 4 + n_biom_extra
prog <- paste0(" mkr",1:n_biom, collapse = " +")
pred <- paste0(" treatment*mkr",1:n_biom, collapse = " +")
formula <- as.formula(paste0("y ~ treatment + ", paste(prog, pred, sep = " +")))
dataset <- dataset.copy <- data
dataset <- cbind(data,mkrCols)

levels = c("treatment", paste0("treatment:mkr",1:n_biom))
levels.orig <- c("treatment", "treatment:severity",
                 "treatment:age", "treatment:sex",
                 "treatment:carrier", paste0("treatment:mkr",5:n_biom))

## Vectorized gsub function to change names of variables back to original
gsub2 <- function(pattern, replacement, x, ...) {
  for(i in 1:length(pattern))
    x <- gsub(pattern[i], replacement[i], x, ...)
  x
}

###############################################################################-
## Analysis using Scores -------------------------------------------------------

## **Obtain score for ML -------------------------------------------------------
ML.results <- score.lm(dataset, input, verbose = F, alpha = 0.05)
#  Print results for ML
knitr::kable(ML.results$ML.output,
             digits=2, row.names = F,
             caption = 'Model Parameters estimates from Maximum Likelihood')

## **Obtain score for Lasso ----------------------------------------------------
set.seed(47436)
lasso.results <- score.lasso(dataset, input, parameters=NULL,
                             alpha = 0.05, verbose = T,
                             lambda = "lagrange",
                             lam_frac = .5)
knitr::kable(lasso.results$Lasso.output,
             digits=2, row.names = F,
             caption = 'Model Parameters estimates from Lasso')
lambda.min <- score.lasso(dataset = dataset, input = input, alpha = 0.05, verbose = T,
                          lambda = "lambda.min")$bestlam
lambda.1se <- score.lasso(dataset = dataset, input = input, alpha = 0.05, verbose = T,
                          lambda = "lambda.1se")$bestlam
lambda.1se
lasso.results$bestlam
lam <- lasso.results$bestlam*sample.size
#  Print results for Lasso
knitr::kable(lasso.results$Lasso.output,
             digits=2, row.names = F,
             caption = 'Model Parameters estimates from Lasso')
## **Obtain score for Additive Noise Lasso  0.2 ------------------
set.seed(47436)
lasso.an.results1 <- score.lasso.added.noise(dataset, input,
                                             lambda = lam,
                                             alpha = 0.05, verbose = T,
                                             perturb_frac = 0.2,
                                             pite.ci = T)
#  Print results for Additive Noise Lasso
knitr::kable(lasso.an.results1$Lasso.output,
             digits=2, row.names = F,
             caption = 'Model Parameters estimates from Additive Noise Lasso')

## **Obtain score for Additive Noise Lasso 0.5 -----------------
set.seed(47436)
lasso.an.results2 <- score.lasso.added.noise(dataset, input,
                                             lambda = lam,
                                             alpha = 0.05, verbose = T,
                                             perturb_frac = 0.5,
                                             pite.ci = T)
#  Print results for Additive Noise Lasso
knitr::kable(lasso.an.results2$Lasso.output,
             digits=2, row.names = F,
             caption = 'Model Parameters estimates from Additive Noise Lasso')
## **Obtain score for Additive Noise Lasso  0.8 ----------------
set.seed(47436)
lasso.an.results3 <- score.lasso.added.noise(dataset, input,
                                             lambda = lam,
                                             alpha = 0.05, verbose = T,
                                             perturb_frac = 0.8,
                                             pite.ci = T)
#  Print results for Additive Noise Lasso
knitr::kable(lasso.an.results3$Lasso.output,
             digits=2, row.names = F,
             caption = 'Model Parameters estimates from Additive Noise Lasso')

##### Plot coefficients --------------------------------------------------------
myPal <- c("#90ca3d","#82b6cc","#0f5e7f","#024460","#f07318")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
Cols <- c("term", "estimate", "LowConfPt", "UpConfPt")
myPal <- c("black","#EE4035", "#F3A530", "#56B949", "#82b6cc","#0f5e7f","#024460")
methods
names(myPal) <- c("ATE", "full", "Lasso","reduced", "reduced-Scheffe", "rLasso", "rLasso-2")
myShape <- c(NA,16,17,15,3,12,13)
names(myShape) <- c("ATE", "full", "Lasso","reduced", "reduced-Scheffe", "rLasso", "rLasso-2")

rbind(cbind(lasso.results$Lasso.output,method="Lasso"),
      cbind(lasso.results$ML.output.M[-1,c(1,2,4,5,6,7)],method="reduced"),
      cbind(lasso.an.results1$Lasso.output,method="rLasso"),
      cbind(ML.results$ML.output[-1,c(1,2,4,5,6,7)], method="full")) -> results.coef
results.coef %>%
  mutate(method = factor(method, levels=c("full", "Lasso", "reduced", "reduced-Scheffe", "rLasso"))) -> results.coef

results.coef %>%
  filter(grepl("treatment",x = term)) %>%
  mutate(estimate = ifelse(estimate != 0, estimate, NA)) %>%
  mutate(term = gsub2(levels, levels.orig, term)) %>%
  mutate(term = gsub("severity0", "mkr10", term)) %>%
  mutate(term = factor(term, levels= levels.orig)) %>%
  ggplot(aes(ymin=LowConfPt, ymax=UpConfPt, y=estimate, x=term, color=method, group = method)) +
  geom_hline(yintercept = 0, color="gray") +
  geom_errorbar(width = 0.5, size = 0.15, position = position_dodge(width=0.5)) +
  geom_point(aes(shape=method), size = 1, position = position_dodge(width=0.5)) +
  scale_color_manual(name = "Method", values = myPal, drop=FALSE) +
  scale_shape_manual(name = "Method", values = myShape, drop=FALSE) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        strip.text   = element_text(size = 7),
        strip.text.y = element_text(size = 7),
        strip.text.x = element_text(size = 9),
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.title.align = 0.5,
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.text.x    = element_text(size = 7, #color = "black",
                                    angle = 45, hjust = 1,
                                    color = c(rep("black",5),rep("gray40", 6))),
        # axis.text.x  = element_text(size = 7, color = "black"),
        # axis.text.x  = element_blank(),
        # axis.ticks.x  = element_blank(),
        strip.placement = "outside") +
  labs(x="Term",y="Estimate") -> plot.coefficients


##### Plot PITE --------------------------------------------------------
Cols <- c("ID", "estimate", "LowConfPt", "UpConfPt", "method")

d0<-(data.frame(ML.results$scores[c("Dx", "Dx.ll", "Dx.ul")],        method="full"))
d1<-(data.frame(lasso.results$scores[c("Dx", "Dx.ll", "Dx.ul")],     method="Lasso"))
d3<-(data.frame(lasso.an.results1$scores[c("Dx", "Dx.ll", "Dx.ul")], method="rLasso"))
d5<-(data.frame(lasso.results$scoresML[c("Dx", "Dx.ll", "Dx.ul")],   method="reduced"))
d6<-(data.frame(lasso.results$scoresSch[c("Dx", "Dx.ll", "Dx.ul")],  method="reduced-Scheffe"))
nrow(d1)

rbind(d0,d1,#d2,
      d3,#d4,
      d5,d6) -> results.all
results.all %>%
  mutate(method = factor(method, levels=c("full", "Lasso", "reduced", "reduced-Scheffe", "rLasso"))) -> results.all

data.frame(results.all, ID= rep(1:41, nrow(results.all)/41)) -> results.all

results.all %>%
  filter(ID<11) %>%
  ggplot(aes(ymin=Dx.ll, ymax=Dx.ul, x=ID, y=Dx, color=method, group = method)) +
  geom_hline(yintercept = 0, color="gray") +
  geom_point(aes(shape=method), size = 0.5, position = position_dodge(width=0.5)) +
  geom_errorbar(width = 0.5, position = position_dodge(width=0.5)) +
  scale_color_manual(name = "Method", values = myPal) +
  scale_shape_manual(name = "Method", values = myShape) +
  theme_bw() +
  labs(x="Subject Id",y="Estimate") +
  scale_x_continuous(breaks=1:10) +
  theme_classic() +
  theme(legend.position = "right",
        text = element_text(size=12),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())  -> plot.subjects

# plot.subjects

# ## *Print plot in PDF ----
# ggsave(filename = "paper/FigureE1-subjects.pdf",
#        plot = plot.subjects,
#        width=7, height=4)


# Selected patients (grid) -------
head(dataset)
dataset.test <- data.frame(carrier = "NON-CARRIER",
                           severity = mean(data$SEVERITY),
                           age = c(70,85,70,85),
                           sex = c("F","F","M","M"),
                           mkr5 = 0,
                           mkr6 = 0,
                           mkr7 = 0,
                           mkr8 = 0,
                           mkr9 = 0,
                           mkr10 = 0,
                           TrueTrtEff = 0)
table(data$AGE)
table(data$SEX)
mean(data$SEVERITY)
# Recategorize predictors to the 'PARAM' parametrization
# And rename them to use in our function
dataset.test$carrier   <- relevel(dataset.test$carrier,   ref="NON-CARRIER")
dataset.test$mkr1      <- scale(dataset.test$severity, center = attr(data$mkr1, "scaled:center"), scale = attr(data$mkr1, "scaled:scale"))
dataset.test$mkr2      <- scale(dataset.test$age, center = attr(data$mkr2, "scaled:center"), scale = attr(data$mkr2, "scaled:scale"))
dataset.test$mkr3      <- scale(dataset.test$sex=="M", center = attr(data$mkr3, "scaled:center"), scale = attr(data$mkr3, "scaled:scale"))
dataset.test$mkr4      <- scale(dataset.test$carrier=="CARRIER", center = attr(data$mkr4, "scaled:center"), scale = attr(data$mkr4, "scaled:scale"))
# dataset.test$mkr3      <- 2*(dataset.test$sex=="M")-1
# dataset.test$mkr4      <- 2*(dataset.test$carrier=="CARRIER")-1

# dataset.test
head(dataset)
alpha = 0.05
dataset.test.ml <- confidence.intervals.ml.test(dataset = dataset.test, input = input,
                                                parameters = NULL,
                                                ML.results = ML.results,
                                                alpha = alpha)
dataset.test.lasso <- confidence.intervals.lasso.test(dataset = dataset.test,
                                                      input = input,
                                                      parameters = NULL,
                                                      lasso.results = lasso.results,
                                                      alpha = alpha)
dataset.test.mlm <- list(scores = dataset.test.lasso$scoresML)
dataset.test.sch <- list(scores = dataset.test.lasso$scoresSch)
dataset.test.lasso.an1 <- confidence.intervals.lassoan.test(dataset = dataset.test,
                                                            input = input,
                                                            parameters = NULL,
                                                            lasso.results = lasso.an.results1,
                                                            alpha = alpha)

results.test <- rbind(
  cbind(method = "full",   dataset.test[c("age","sex")], ID = 1:4,
        dataset.test.ml$scores[c("Dx", "Dx.ll", "Dx.ul")]),
  cbind(method = "Lasso",dataset.test[c("age","sex")],ID = 1:4,
        dataset.test.lasso$scores[c("Dx", "Dx.ll", "Dx.ul")]),
  cbind(method = "reduced",  dataset.test[c("age","sex")],ID = 1:4,
        dataset.test.mlm$scores[c("Dx", "Dx.ll", "Dx.ul")]),
  cbind(method = "reduced-Scheffe",  dataset.test[c("age","sex")],ID = 1:4,
        dataset.test.sch$scores[c("Dx", "Dx.ll", "Dx.ul")]),
  cbind(method = "rLasso", dataset.test[c("age","sex")],ID = 1:4,
        dataset.test.lasso.an1$scores[c("Dx", "Dx.ll", "Dx.ul")]))
results.test %>%
  mutate(method = factor(method, levels=c("full", "Lasso", "reduced", "reduced-Scheffe", "rLasso"))) -> results.test

##### Plot PITE --------------------------------------------------------
Cols <- c("ID", "estimate", "LowConfPt", "UpConfPt", "method")

results.test %>%
  mutate(subject = paste0(ifelse(sex=="F", "Female, ", "Male, "),paste0("Age: ",age))) %>%
  ggplot(aes(ymin=Dx.ll, ymax=Dx.ul, x=1, y=Dx, color=method, group = method)) +
  geom_hline(yintercept = 0, color="gray") +
  geom_point(aes(shape = method), size = 1, position = position_dodge(width=0.5)) +
  geom_errorbar(width = 0.15, size = 0.20, position = position_dodge(width=0.5)) +
  scale_color_manual(name = "Method", values = myPal) +
  scale_shape_manual(name = "Method", values = myShape) +
  labs(x="Subject Id", y="Estimate") +
  scale_x_continuous(breaks=1:10) +
  facet_wrap(~ subject, nrow = 1, strip.position = "bottom") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        strip.text   = element_text(size = 7),
        strip.text.y = element_text(size = 7),
        strip.text.x = element_text(size = 9),
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.title.align = 0.5,
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.text    = element_text(size = 7, color = "black"),
        # axis.text.x  = element_text(size = 7, color = "black"),
        axis.text.x  = element_blank(),
        axis.ticks.x  = element_blank(),
        strip.placement = "outside") -> plot.grid
# plot.grid

ggarrange(plot.coefficients,
          plot.grid,
          heights = c(5,3.5),
          labels = c("A", "B"),
          common.legend = TRUE, legend = "right",
          ncol=1, nrow = 2) -> p2

rproj <- rprojroot::find_package_root_file()
path.figures <- paste0(rproj, "/paper/figures/")

ggsave(filename = paste0(path.figures,"Fig1.pdf"),
       plot = p2,
       width=17, height=12.97, units = "cm")
ggsave(filename = paste0(path.figures,"Fig1.eps"),
       plot = p2,
       width=17, height=12.97, units = "cm")
