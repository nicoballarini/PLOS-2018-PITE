################################################################################
# Erase everything in environment for reproducible results and set seed
cat("\014"); rm(list = ls())
################################################################################
# Load required packages -------------------------------------------------------
library(PITE)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(glmnet)
library(parallel)
library(sampling)
library(selectiveInference)
library(stats)
library(survminer)
library(survival)

################################################################################
# We use the dataset from Rosenkranz (2016) https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201500147
# to illustrate the methods proposed in this work.
# The data comes from a clinical trial of an prostate cancer treatment
# Data is loaded from the article's supplementary material
# Load Dataset and data transformation -----------------------------------------
library(haven)
data_url = "https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fbimj.201500147&attachmentId=2173117089"
temp <- tempfile()
download.file(data_url,temp)
prca0 <- prca <- read_sas(unz(temp, "adv_prostate_ca.sas7bdat"))

# Find package root and create path to 'figures' folder to save figures for manuscript
rproj <- rprojroot::find_package_root_file()
path.figures <- paste0(rproj, "/paper/figures/")

head(prca0)
# Select the variables that we use for the analysis
prca <- prca.0 <- prca0[, c("SURVTIME","CENS","RX","BM","HX","STAGE","PF", "AGE", "WT" )]
names(prca.0) <- c("survtime","cens","rx","bm",
                   "hx","stage","pf","age","wt")
prca$STAGE <- prca$STAGE - 3
# Name the variables in lower case to simplify in R
names(prca)<- c("survtime","cens","rx","bm",
                "hx","stage","pf","age","wt")
prca.orig <- prca
# cor(prca)

prca$treatment  <- -1*(prca$rx==0)+1*(prca$rx==1)
prca$agegroup <- 1 * (prca0$AGE > 65)
prca$wtgroup <- 1 * (prca0$WT > 100)
# head(prca)
std.vars <- c("bm","hx","stage","pf","age","wt")
prca.std <- scale(prca[std.vars])
x.center <- attr(prca.std, "scaled:center")
x.scale <- attr(prca.std, "scaled:scale")
# Standardize all biomarkers
prca[std.vars] <- prca.std
p.val <- survdiff(Surv(survtime, cens)~rx, data=prca.0)
lrt.p.val <- 1 - pchisq(p.val$chisq, length(p.val$n) - 1)

######################################################################-
## Fit the kaplan meier curves to the data and get plot and estimates -----
##  of P(T>50|G=k) k=0,1
fit<-survfit(Surv(survtime, cens)~rx, data=prca.0)
font.size = 8
# pdf("KM_all.pdf", width = 5, height = 3.5)
# ggsurvplot(
#   fit,                     # survfit object with calculated statistics.
#   pval     = TRUE,             # show p-value of log-rank test.
#   conf.int = TRUE,         # show confidence intervals for
#   palette = c('#0571b0', '#ca0020'),
#   # main  = "Kaplan-Meier curves by treatment",
#   ggtheme = theme_classic(),
#   legend  = c(0.8,0.8),
#   break.time.by = 15, xlab = "Time (days)",
#   font.main = font.size, font.x = font.size,
#   font.y = font.size, font.tickslab =  font.size, font.legend = font.size
# )
# dev.off()
surv<-summary(fit, time=50)
diff<-surv$surv[2]-surv$surv[1]


######################################################################-
## Fit proportional hazard model to create a score-----
## Formula includes all variables and its interaction with treatment
reg.1<-survival::coxph(Surv(survtime, cens) ~
                         treatment + age + bm + stage + pf + hx + wt +
                         treatment*age + treatment*bm + treatment*stage + treatment*pf + treatment*hx +
                         treatment*wt,
                       data=prca,ties="breslow", x=TRUE)
cox.1<-survfit(reg.1)
# diag(cov(reg.1$x))

######################################################################-
## Estimate the score-----
# # head(prca)
# cdf<-ecdf(prca$Dx)
# plot(cdf)
# prca$hDx<-cdf(prca$Dx)
# hcdf<-ecdf(prca$hDx)
# plot(hcdf)

# head(prca)
coef <- reg.1$coefficients
cov.matrix <- reg.1$var
n_biom = 6
all.vars <- names(coef)
var.score <- all.vars[2:(n_biom+1)]
coef.dx <- coef[c(1,(n_biom+2):(2*n_biom+1))]
X.dx <- 2*cbind(1,prca[var.score])
X.dx.full <- 2*cbind(1,matrix(0, nrow(prca), length(var.score)), prca[var.score])
# coef.dx
X.m <- model.matrix(~-1+treatment + age + bm + stage + pf + hx + wt +
                      treatment*age + treatment*bm + treatment*stage + treatment*pf + treatment*hx +
                      treatment*wt, data=prca)

# head(X.dx)
# head(X.m)
alpha. <- 0.05
scores <- data.frame(Dx=numeric(nrow(prca)))
# head(as.matrix(X.m)%*%coef)
scores$Dx <- Dx.m <- as.matrix(X.dx)%*%coef.dx
scores$Dx.v <- vDx.m <- diag(as.matrix(X.dx.full)%*%cov.matrix%*%t(as.matrix(X.dx.full)))

zsigma <- qnorm(1-alpha./2)*sqrt(vDx.m)
scores$Dx.ll <- drop(Dx.m - zsigma)
scores$Dx.ul <- drop(Dx.m + zsigma)
# head(scores)
ml.results <- broom::tidy(reg.1)[,c(1,2,5,6,7)]  #%>%
# ml.results %>%
# knitr::kable(digits = 4)

myPal <- c("#90ca3d","#82b6cc","#0f5e7f","#024460","#f07318")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

reg.0 <- survival::coxph(Surv(survtime, cens) ~
                           treatment ,
                         data=prca,ties="breslow", x=TRUE)
step(object = reg.0,
     scope = Surv(survtime, cens) ~
       treatment + age + bm + stage + pf + hx + wt +
       treatment*age + treatment*bm + treatment*stage + treatment*pf + treatment*hx +
       treatment*wt, trace = 0,
     direction = "forward") -> fw.results
# fw.results
fw.results2 <- broom::tidy(fw.results)
fw.results1 <- data.frame(term = ml.results$term)
fw.results <- left_join(fw.results1, fw.results2)[,-c(3,4)]

step(object = reg.1, scope = Surv(survtime, cens) ~
       treatment + age + bm + stage + pf + hx + wt +
       treatment*age + treatment*bm + treatment*stage + treatment*pf + treatment*hx +
       treatment*wt, trace = 0,
     direction = "backward") -> bw.results
# bw.results
bw.results2 <- broom::tidy(bw.results)
bw.results1 <- data.frame(term = ml.results$term)
bw.results <- left_join(bw.results1, bw.results2)[,-c(3,4)]
# bw.results %>%
#   knitr::kable(digits = 4)

######################################################################-
## Fit proportional hazard model to create a score-----
## Formula includes all variables and its interaction with treatment
y = Surv(prca$survtime, prca$cens)
# X.m
n=nrow(prca)
set.seed(10)
standard_normal <- matrix(rnorm(n*500, 0,1),nrow = n, ncol = 500)
XT.e <- abs(t(X.m)%*%standard_normal)
lam_frac = .5
lam = lam_frac * mean(apply(XT.e,2,max))
lam
bestlam <- lam/n
lambda = bestlam
lambda.min <- glmnet::cv.glmnet(X.m,y,family="cox")$lambda.min
lambda.1se <- glmnet::cv.glmnet(X.m,y,family="cox")$lambda.1se
lambda.min
lambda.1se
lambda
# reg.2 <- glmnet::cv.glmnet(X.m,y,family="cox")
# lambda = reg.2$lambda.1se
# lambda
reg.2 <- glmnet(X.m,y,family="cox")
coef.lasso <- coef(reg.2, s=lambda)[,1]
coef.lasso
# compute fixed lambda p-values and selection intervals
n <- nrow(prca)
coef.lasso

coef.lasso[which(abs(coef.lasso) < 0.1/sqrt(colSums(X.m^2)))] <- 0
out = selectiveInference::fixedLassoInf(x = X.m,
                                        y = prca$survtime,
                                        beta = coef.lasso,
                                        lambda = lambda*n,
                                        status=prca$cens,
                                        family="cox", alpha = 0.05)
fixedLassoInfprint <- function(x, coef) {
  tab = data.frame(names(coef)[x$vars],
                   round(x$coef0,3),
                   round(x$zscore0,3),
                   round(x$pv,3),round(x$ci,3))
  colnames(tab) = c("Var", "Coef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
  tab = cbind(tab,round(x$tailarea,3))
  colnames(tab)[(ncol(tab)-1):ncol(tab)] = c("LowTailArea","UpTailArea")
  rownames(tab) = x$vars
  tab
  tab %>%
    mutate(tmpLow = ifelse(Coef < 0,-UpConfPt, LowConfPt),
           tmpUp  = ifelse(Coef < 0,-LowConfPt, UpConfPt),
           LowConfPt = tmpLow,
           UpConfPt = tmpUp) %>%
    dplyr::select(-tmpLow,-tmpUp)
}
lasso.results1 <- data.frame(term = names(coef.lasso),
                             estimate = coef.lasso)
lasso.results2 <- fixedLassoInfprint(out, coef.lasso)
lasso.results2 <- lasso.results2[c("Var", "P-value", "LowConfPt", "UpConfPt")]
names(lasso.results2) <- c("term", "p.value", "conf.low",   "conf.high")
lasso.results <- left_join(lasso.results1, lasso.results2)
# lasso.results

# lasso.results %>%
#   knitr::kable(digits = 4)

######################################################################-
## Fit proportional hazard model to create a score-----
## Formula includes all variables and its interaction with treatment
vars.lasso <- names(coef(reg.2, s=lambda)[,1])[out$vars]

# score.reduced2(prca, vars.lasso)$reg.out

formula <- paste("Surv(survtime, cens) ~", paste(vars.lasso, collapse = " + "))
reg.m <- survival::coxph(eval(parse(text=formula)),
                         data=prca,ties="breslow")
ML.output.M <- broom::tidy(reg.m)[,c(1,2,5,6,7)]
int <- ML.output.M$term[grepl(":treatment", ML.output.M$term)]
if (length(int)!=0) {
  sub(":treatment", "", int)
  ML.output.M$term[grepl("treatment", ML.output.M$term)] <- paste0("treatment:", sub(":treatment", "", int))
}
ml.results.m <- ml.results
ml.results.m[,-1] <- NA
ml.results.m[,2] <- 0
# ml.results.m$term

ind <- match(ML.output.M$term, ml.results.m$term)
ml.results.m[ind, 2:5] <- ML.output.M[2:5]
coef.m <- ml.results.m$estimate


covb.M <- vcov(reg.m)        # Store covariance matrix for the estimates
# head(X.dx.full)
names(X.dx.full) <- colnames(reg.1$x)

Dx.M.m <- as.matrix(X.dx.full)%*%coef.m #Create score or PITE
covb.M.dx <- covb.M[grepl( "treatment", colnames( covb.M ) ),
                    grepl( "treatment", colnames( covb.M ) )]
# Calculate variance of the score
vDx.M.m   <- diag(as.matrix(X.dx.full[, coef.m!=0]) %*%
                    covb.M %*%
                    t(X.dx.full[, coef.m!=0]))

scores.MLM <- data.frame(Dx=Dx.M.m)
# head(as.matrix(X.m)%*%coef)
# length(vDx.M.m)
scores.MLM$Dx.v <- vDx.M.m
zsigma <- qnorm(1-alpha./2)*sqrt(vDx.M.m)
scores.MLM$Dx.ll <- drop(Dx.M.m - zsigma)
scores.MLM$Dx.ul <- drop(Dx.M.m + zsigma)
# head(scores)

# ml.results.m %>%
# knitr::kable(digits = 4)

lasso.results %>%
  dplyr::select(term, estimate, p.value, conf.low, conf.high) %>%
  mutate(method="Lasso") %>%
  rbind(data.frame(ml.results, method="full"),
        # data.frame(fw.results, method="cox-forward"),
        # data.frame(bw.results, method="cox-backward"),
        data.frame(ml.results.m, method="reduced")) %>%
  mutate(term = factor(term, levels=(all.vars))) %>%
  mutate(estimate = ifelse(estimate==0,NA, estimate)) %>%
  mutate(method = factor(method, levels=c("full",
                                          "cox-forward",
                                          "cox-backward",
                                          "Lasso",
                                          "reduced")))-> datplot

myPal <- c("black","#EE4035", "#F3A530", "#56B949", "#82b6cc","#0f5e7f","#024460")
names(myPal) <- c("ATE", "full", "Lasso","reduced", "reduced-Scheffe", "rLasso", "rLasso-2")
myShape <- c(NA,16,17,15,3,12,13)
names(myShape) <- c("ATE", "full", "Lasso","reduced", "reduced-Scheffe", "rLasso", "rLasso-2")

## Plot Coefficients -----
datplot  %>%
  ggplot(aes(ymin=conf.low, ymax=conf.high,
             x=term, y=estimate, color=method, group = method)) +
  geom_hline(yintercept = 0, color="gray") +
  geom_point(aes(shape = method), position = position_dodge(width=0.6), size=1) +
  geom_errorbar(width = 0.5, size=0.15, position = position_dodge(width=0.6)) +
  scale_color_manual(name = "Method", values = myPal) +
  scale_shape_manual(name = "Method", values = myShape) +
  theme_bw() +
  labs(x="Term",y="Estimate") +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1,
                                   color = c("red", rep("black",6), rep("red", 6))),
        axis.text.y = element_text(size = 7, colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),

        strip.text   = element_text(size = 7),
        strip.text.y = element_text(size = 7),
        strip.text.x = element_text(size = 9),
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.title.align = 0.5,
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        strip.placement = "outside") -> plot_coefficients



out.score <- PITE::fixedCoxLassoInf_eta(x = X.m,
                                                 y = prca$survtime,
                                                 beta = coef.lasso,
                                                 lambda = lam,
                                                 status=prca$cens,
                                                 alpha = 0.05, contrast = X.dx.full)
# print(out)
# head(X.dx.full)
scoresLasso <- data.frame(Dx = as.matrix(X.dx.full)%*%coef.lasso)
# scoresLasso$Dx.2 <- drop(out$coef0)
scoresLasso$Dx.v <- out.score$sd
scoresLasso$Dx.ll <- out.score$ci[,1]
scoresLasso$Dx.ul <- out.score$ci[,2]

scores_all <- rbind(data.frame(scoresLasso, method="Lasso"),
                    data.frame(scores.MLM,  method="reduced"),
                    data.frame(scores,      method="full"))

##### Plot PITE --------------------------------------------------------
Cols <- c("ID", "estimate", "LowConfPt", "UpConfPt", "method")
data.frame(scores_all, ID= rep(1:(nrow(scores_all)/3), 3)) -> results.all

results.all %>%
  filter(ID<11) %>%
  mutate(method = factor(method, levels=c("full",
                                          "forward",
                                          "backward",
                                          "Lasso",
                                          "reduced"))) %>%
  ggplot(aes(ymin=Dx.ll, ymax=Dx.ul, x=ID, y=Dx, color=method, group = method)) +
  geom_hline(yintercept = 0, color="gray") +
  geom_point(aes(shape = method), size = 1, position = position_dodge(width=0.6)) +
  geom_errorbar(width = 0.5, size=0.15, position = position_dodge(width=0.6)) +
  scale_color_manual(name = "Method", values = myPal) +
  scale_shape_manual(name = "Method", values = myShape) +
  theme_bw() +
  labs(x="Subject Id",y="Estimate") +
  scale_x_continuous(breaks=1:10) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 7, colour = "black"),
        axis.text.y = element_text(size = 7, colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        strip.text   = element_text(size = 7),
        strip.text.y = element_text(size = 7),
        strip.text.x = element_text(size = 9),
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.title.align = 0.5,
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        strip.placement = "outside") -> plot_subjects

# plot ----------
ggarrange(plot_coefficients,
          plot_subjects,
          heights = c(5,4),
          labels = c("A", "B"),
          common.legend = TRUE, legend = "right",
          ncol=1, nrow = 2) -> plot4

ggsave(paste0(path.figures,"Fig4.pdf"), plot = plot4, width = 17, height = 15, units = "cm")
ggsave(paste0(path.figures,"Fig4.eps"), plot = plot4, width = 17, height = 15, units = "cm")



# prca$pf
# attr(prca$age, which = "scaled:center")
# attr(prca$age, which = "scaled:scale")
# prca.orig
Age <- prca.orig$age
contrast = expand.grid(age = seq(min(Age), max(Age), 1),
                       bm = c(0,1))
contrast = expand.grid(age = seq(min(Age), 90, 1),
                       bm = c(0,1))
# head(contrast)
# names(X.dx.full)

X.dx.grid <- X.dx.full[1:nrow(contrast),]
X.dx.grid[paste0("treatment:",names(contrast))] <- contrast
X.dx.grid[-which(names(X.dx.grid) %in% paste0("treatment:",names(contrast)))] <- 0
X.dx.grid$`treatment:age` <- scale(X.dx.grid$`treatment:age`,
                             x.center["age"],
                             x.scale["age"])
X.dx.grid$`treatment:bm` <- scale(X.dx.grid$`treatment:bm`,
                            x.center["bm"],
                            x.scale["bm"])
out.score <- PITE::fixedCoxLassoInf_eta(x = X.m,
                                                 y = prca$survtime,
                                                 beta = coef.lasso,
                                                 lambda = lam,
                                                 gridrange = c(-100,10),
                                                 status=prca$cens,
                                                 alpha = 0.050,
                                                 contrast = 2*X.dx.grid)
# print(out)
# head(X.dx.full)
scoresLasso.grid <- data.frame(Dx = 2*as.matrix(X.dx.grid)%*%coef.lasso)
# scoresLasso$Dx.2 <- drop(out$coef0)
scoresLasso.grid$Dx.v <- out.score$sd
scoresLasso.grid$Dx.ll <- out.score$ci[,1]
scoresLasso.grid$Dx.ul <- out.score$ci[,2]
# head(scoresLasso.grid)
scoresLasso.grid %>%
  cbind(contrast) %>%
  mutate(bm = factor(ifelse(bm==1,"Yes", "No"), levels=c("No", "Yes"))) %>%
  ggplot(aes(x=age, y=Dx, color=as.factor(bm), fill=as.factor(bm))) +
  geom_hline(aes(yintercept = 0), color="gray", linetype=2) +
  geom_ribbon(aes(ymin=Dx.ll,ymax=Dx.ul, color=NULL), alpha=0.1) +
  geom_line() +
  theme_bw() +
  scale_color_manual(values = c("No" = "#f1a340", "Yes" = "#998ec3")) +
  scale_fill_manual(values = c("No" = "#f1a340", "Yes" = "#998ec3")) +
  labs(x="Age",y="Estimate", color="Bone\nmetastasis", fill="Bone\nmetastasis") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 7, colour = "black"),
        axis.text.y = element_text(size = 7, colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        strip.text   = element_text(size = 7),
        strip.text.y = element_text(size = 7),
        strip.text.x = element_text(size = 9),
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.title.align = 0.5,
        legend.position = c(0.81, 0.15),
        legend.key.height=unit(0.5,"line"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        strip.placement = "outside") -> grid_pite


scoresLasso.grid %>%
  mutate(subgroup = ifelse(Dx.ul < 0, "hat(B)[u]",
                           ifelse(Dx < 0, "hat(B)",
                                  ifelse(Dx.ll < 0, "hat(B)[l]",
                                         "hat(B)[l]^c")))) %>%
  mutate(subgroup = factor(subgroup,
                           levels = c("hat(B)[l]^c", "hat(B)[l]", "hat(B)", "hat(B)[u]"))) %>%
  cbind(contrast) %>%
  mutate(bm = factor(ifelse(bm==1,"Yes", "No"), levels=c("No", "Yes"))) %>%
  ggplot(aes(x=age, y=as.factor(bm), color=subgroup, fill=subgroup)) +
  geom_tile(height=0.7, colour = "black", size=0.05) +
  scale_fill_brewer(palette = 1,
                    labels = scales::parse_format()(c("hat(B)[l]^c", "hat(B)[l]", "hat(B)", "hat(B)[u]"))) +
  # scale_color_brewer(palette = 1) +
  theme_bw() +
  labs(x="Age",y="Bone\nmetastasis", color="Subgroup", fill="Subgroup") +
  theme_classic() +
  scale_x_continuous(expand=c(0,0), #sec.axis = dup_axis(name = " "),
                     # breaks = seq(min(Age), max(Age),2)) +
                     breaks = c(50,60,70,80,90)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(legend.position = "top",
        axis.line = element_blank(),
        axis.text.x = element_text(size = 7, colour = "black"),
        axis.text.y = element_text(size = 7, colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        strip.text   = element_text(size = 7),
        strip.text.y = element_text(size = 7),
        strip.text.x = element_text(size = 9),
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.title.align = 0.5,
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        strip.placement = "outside") -> grid_subgr


# plot ----------
ggarrange(grid_pite,
          grid_subgr,
          # heights = c(5,4),
          labels = c("A", "B"),
          common.legend = FALSE,
          ncol=2, nrow = 1) -> plot5

ggsave(paste0(path.figures,"Fig5.pdf"), plot = plot5, width = 17, height = 8, units = "cm")
ggsave(paste0(path.figures,"Fig5.eps"), plot = plot5, width = 17, height = 8, units = "cm")
