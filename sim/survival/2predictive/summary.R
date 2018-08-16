# Make sure we are working in a clean workspace for reproducibility
rm(list=ls())
cat("\014")
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
mc.cores =     25 # Number of cores for parallel computing

img.width  = 3.14
img.height = 3.14
path = NULL
# # If running interactive, probide package path
rproj <- rprojroot::find_package_root_file()
path       <- "/sim/survival/2predictive/"
path <- paste0(rproj, path)

# if (!dir.exists(paste0(path, "results/"))) dir.create(paste0(path, "results/"))
# if (!dir.exists(paste0(path, "results/data/"))) dir.create(paste0(path, "results/data/"))
if (!dir.exists(paste0(path, "results/plots/"))) dir.create(paste0(path, "results/plots/"))

path.data <- paste0(path, "results/data/")
path.plot <- paste0(path, "results/plots/")

files.2 <- list.files(path=path.data)
## Coefficients DATA ----
n = c(100,400,700,1000)
files.2.null    <- files.2[grep(files.2, pattern = "result.null")]
files.2.cox     <- files.2[grep(files.2, pattern = "result.cox")]
files.2.lasso   <- files.2[grep(files.2, pattern = "result.lasso")]
files.2.reduced <- files.2[grep(files.2, pattern = "result.reduced")]


i = files.2.cox[1]

results.coeficients <- data.frame()
results.score <- data.frame()
for(i in n){
  load(paste0(path.data, "result.null", i, ".Rda"))
  a <- do.call(rbind,
               lapply(1:1000,
                      function(x)
                        cbind(result.null[[x]]$reg.out, method="null", n = i, nsim = x)
               ))
  b <- do.call(rbind,
               lapply(1:1000,
                      function(x)
                        cbind(result.null[[x]]$scores, method="null", n = i, nsim = x)
               ))
  results.coeficients <- rbind(results.coeficients,
                               a)
  results.score <- rbind(results.score,
                         b)
}
for(i in n){
  load(paste0(path.data, "result.cox", i, ".Rda"))
  a <- do.call(rbind,
               lapply(1:1000,
                      function(x)
                        cbind(result.cox[[x]]$reg.out, method="cox", n = i, nsim = x)
               ))
  b <- do.call(rbind,
               lapply(1:1000,
                      function(x)
                        cbind(result.cox[[x]]$scores, method="cox", n = i, nsim = x)
               ))
  results.coeficients <- rbind(results.coeficients,
                               a)
  results.score <- rbind(results.score,
                         b)
}
for(i in n){
  load(paste0(path.data, "result.lasso", i, ".Rda"))
  a <- do.call(rbind,
               lapply(1:1000,
                      function(x)
                        cbind(result.lasso[[x]]$reg.out, method="cox-lasso", n = i, nsim = x)
               ))
  b <- do.call(rbind,
               lapply(1:1000,
                      function(x)
                        cbind(result.lasso[[x]]$scores, method="cox-lasso", n = i, nsim = x)
               ))
  results.coeficients <- rbind(results.coeficients,
                               a)
  results.score <- rbind(results.score,
                         b)
}
for(i in n){
  load(paste0(path.data, "result.reduced", i, ".Rda"))
  a <- do.call(rbind,
               lapply(1:1000,
                      function(x)
                        cbind(result.reduced[[x]]$reg.out, method="cox-reduced", n = i, nsim = x)
               ))
  b <- do.call(rbind,
               lapply(1:1000,
                      function(x)
                        cbind(result.reduced[[x]]$scores, method="cox-reduced", n = i, nsim = x)
               ))
  results.coeficients <- rbind(results.coeficients,
                               a)
  results.score <- rbind(results.score,
                         b)
}

results.coeficients$method <- factor(results.coeficients$method, levels = c("null", "cox","cox-reduced","cox-lasso"))
results.score$method       <- factor(results.score$method, levels = c("null", "cox","cox-reduced","cox-lasso"))

# Rename levels of factor
levels(results.coeficients$method) = c("null", "full","reduced", "Lasso")
levels(results.score$method) = c("null", "full","reduced", "Lasso")

# Re-order levels of factor
results.coeficients$method <- factor(results.coeficients$method, levels = c("null", "full", "Lasso","reduced"))
results.score$method       <- factor(results.score$method,       levels = c("null", "full", "Lasso","reduced"))



###############################################################################-
## Plot coefficients -----------------------------------------------------------
# do.call(rbind,
#         lapply(1:length(result.lasso), function(x){
#           data.frame(result.lasso[[x]]$reg.out, sim = x, method="lasso")})
# ) -> l.o
# do.call(rbind,lapply(1:length(result.cox), function(x){
#   data.frame(result.cox[[x]]$reg.out, sim = x, method="cox")
# })) -> c.o
# a.o <- rbind(l.o, c.o)
# head(a.o)
# a.o[1:13, "term"]
delta  <- -0.08
delta1 <-  0.18
beta   <- -0.10
true   <- c(beta,
            0.17, 0.18, 0.11, 0.14, 0.21, -0.12,
            delta, delta1,
            rep(0,4))

myPal <- c("black","#EE4035", "#F3A530", "#56B949", "#82b6cc","#0f5e7f","#024460")
names(myPal) <- c("null", "full", "Lasso","reduced", "reduced-Scheffe", "rLasso", "rLasso-2")
myShape <- c(NA,16,17,15,3,12,13)
names(myShape) <- c("null", "full", "Lasso","reduced", "reduced-Scheffe", "rLasso", "rLasso-2")

## ***Coverage of confidence intervals for coefficients -------
pdf(file = paste0(path.plot,"coefficients-coverage.pdf"),
    width = img.width, height = img.height)
results.coeficients %>%
  filter(method != "null") %>%
  mutate(cover = 1*((conf.low < true)&(true < conf.high))) %>%
  # mutate(cover = 1*((conf.low < 0)&(0 < conf.high))) %>%
  group_by(term, method, n) %>%
  summarize(cover = mean(cover, na.rm=T),
            mean.est = mean(estimate, na.rm=T)) %>%
  ggplot() +
  ylim(c(0,1)) +
  geom_hline(yintercept=1-alpha) +
  facet_wrap(~n, nrow = 2)+
  scale_color_manual(name = "Method", values = myPal) +
  scale_fill_manual(name = "Method", values = myPal) +
  guides(shape = guide_legend(title = "Method")) +
  geom_bar(aes(x= term, y=cover, fill=method), stat="identity", position="dodge") +
  ggtitle("Coverage of confidence intervals for coefficients")

dev.off(which = dev.cur())

## ***Estimates for coefficients  -------
pdf(file = paste0(path.plot,"coefficients-estimates.pdf"),
    width = img.width, height = img.height)

results.coeficients %>%
  group_by(term, method, n) %>%
  ggplot() +
  geom_hline(yintercept=0) +
  facet_wrap(~n, nrow = 2)+
  scale_color_manual(name = "Method", values = myPal) +
  scale_fill_manual(name = "Method", values = myPal) +
  # scale_shape_manual(name = "Method", values = myShape) +
  geom_boxplot(aes(x= term, y=estimate, fill=method),outlier.size = 0.1,
               position="dodge") +
  ggtitle("Estimates for coefficients") +
  coord_cartesian(ylim = c(-0.75,0.75))
dev.off()

## ***Percentage of inclusion by terms when using the Lasso  -------
pdf(file = paste0(path.plot,"coefficients-inclusion.pdf"),
    width = img.width, height = img.height)

results.coeficients %>%
  mutate(included = 1*(estimate != 0 )) %>%
  mutate(N = as.factor(n)) %>%
  filter(method=="Lasso") %>%
  filter(grepl(pattern = "z", term)) %>%
  group_by(term, method, N) %>%
  summarize(included = mean(included, na.rm=T)) %>%
  ggplot() +
  geom_hline(yintercept=0) +
  geom_bar(aes(x= term, y = 100*included, fill=N), color="#08519c",
           stat="identity", position="dodge") +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  labs(y="Percentage", x= "Term") +
  # ggtitle("Percentage of inclusion by terms when using the Lasso") +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 7),
        strip.text   = element_text(size = 7),
        strip.text.y = element_text(size = 7),
        strip.text.x = element_text(size = 7),
        legend.text  = element_text(size = 7),
        legend.text.align = 1,
        legend.title = element_text(size = 7),
        legend.title.align = 0.5,
        legend.position = c(0.81, 0.85),
        legend.key.height=unit(0.5,"line"),
        # legend.background = element_rect(color = "black", size = 0.1, linetype = "solid"),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.text    = element_text(size = 6, color = "black"),
        axis.text.x  = element_text(size = 6, color = "black"),
        strip.placement = "outside") -> vars_inclusion
vars_inclusion
dev.off()


## Plot PITE ------------------------------------------------------------------

pdf(file = paste0(path.plot,"pite-coverage-1.pdf"),
    width = img.width, height = img.height)

nsim=1000
p=alpha
q=1-p
e=1.96*sqrt(p*q/nsim)

results.score %>%
  group_by(method, n) %>%
  summarize(cover = mean(cover, na.rm=T)) %>%
  ggplot() +
  ylim(c(0,1)) +
  facet_wrap(~n, nrow = 1)+
  scale_color_manual(name = "Method", values = myPal) +
  scale_fill_manual(name = "Method", values = myPal) +
  guides(shape = guide_legend(title = "Method")) +
  geom_bar(aes(x= 1, y=cover, fill=method), stat="identity", position="dodge") +
  geom_hline(yintercept=1-alpha+e, linetype=2) +
  geom_hline(yintercept=1-alpha) +
  geom_hline(yintercept=1-alpha-e, linetype=2) +
  ggtitle("Coverage for PITE")

dev.off()


## ***Coverage for PITE  -------

pdf(file = paste0(path.plot,"pite-coverage.pdf"),
    width = img.width, height = img.height)

nsim=1000
p=alpha
q=1-p
e=1.96*sqrt(p*q/nsim)

results.score %>%
  mutate(N = (n)) %>%
  group_by(method, N) %>%
  summarize(cover = mean(cover, na.rm=T)) %>%
  ggplot() +
  ylim(c(0.7,1)) +
  geom_line(aes(x = N, y=cover, color=method,  group=method), linetype = 2, size=0.15) +
  geom_point(aes(x = N, y=cover, color=method, shape = method, #fill = method,
                 group=method), size = 1) +
  geom_hline(yintercept=1-alpha+e, linetype=3, size=0.1) +
  geom_hline(yintercept=1-alpha, size=0.1) +
  geom_hline(yintercept=1-alpha-e, linetype=3, size=0.1) +
  scale_x_continuous(breaks = c(0,n), limits = c(0,1000)) +
  labs(y="Coverage", x= "N") +
  scale_color_manual(name = "Method", values = myPal) +
  scale_fill_manual(name = "Method", values = myPal) +
  scale_shape_manual(name = "Method", values = myShape) +
  guides(shape = guide_legend(title = "Method")) +
  # ggtitle("Coverage for PITE")+
  theme_classic()  +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 7),
        strip.text   = element_text(size = 7),
        strip.text.y = element_text(size = 7),
        strip.text.x = element_text(size = 7),
        legend.text  = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.title.align = 0.5,
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        legend.position = c(0.81, 0.15),
        legend.key.height=unit(0.5,"line"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        axis.text    = element_text(size = 6, color = "black"),
        axis.text.x  = element_text(size = 6, color = "black"),
        strip.placement = "outside") -> pite_cover
pite_cover
dev.off()


## PITE estimate ---------------------------------------------------------------
pdf(file = paste0(path.plot,"pite-estimate.pdf"),
    width = 11, height = 8,
    paper = "a4r")
results.score %>%
  group_by(method, n, nsim) %>%
  summarise(Dx = mean(Dx)) %>%
  ggplot() +
  geom_hline(yintercept=0) +
  facet_wrap(~n, nrow = 1)+
  scale_color_manual(name = "Method", values = myPal) +
  scale_fill_manual(name = "Method", values = myPal) +
  geom_boxplot(aes(x= 1, y=Dx, fill=method),outlier.size = 0.1,
               position="dodge") +
  ggtitle("Estimate for PITE")
dev.off()




pdf(file = paste0(path.plot,"pite-bias.pdf"),
    width = img.width, height = img.height)

results.score %>%
  group_by(method, n, nsim) %>%
  summarise(bias = mean(Dx-pite.full)) %>%
  ggplot() +
  geom_hline(yintercept=0) +
  facet_wrap(~n, nrow = 1)+
  scale_color_manual(name = "Method", values = myPal) +
  scale_fill_manual(name = "Method", values = myPal) +
  geom_boxplot(aes(x= 1, y=bias, fill=method),outlier.size = 0.1,
               position="dodge") +
  ggtitle("BIAS for PITE") -> pite_bias
pite_bias
dev.off()

pdf(file = paste0(path.plot,"pite-mse.pdf"),
    width = img.width, height = img.height)

results.score %>%
  group_by(method, n, nsim) %>%
  summarise(mse = mean((Dx-pite.full)^2)) %>%
  ggplot() +
  geom_hline(yintercept=0) +
  facet_wrap(~n, nrow = 1) +
  scale_color_manual(name = "Method", values = myPal) +
  scale_fill_manual(name = "Method", values = myPal) +
  geom_boxplot(aes(x= 1, y=sqrt(mse), fill=method),outlier.size = 0.1,
               position="dodge") +
  ggtitle("sqrt(MSE) for PITE") -> pite_mse
pite_mse
dev.off()




pdf(file = paste0(path.plot,"pite-width.pdf"),
    width = img.width, height = img.height)

results.score %>%
  group_by(method, n, nsim) %>%
  summarise(width = mean(width, na.rm = T)) %>%
  ggplot() +
  geom_hline(yintercept=0) +
  facet_wrap(~n, nrow = 1)+
  scale_color_manual(name = "Method", values = myPal) +
  scale_fill_manual(name = "Method", values = myPal) +
  geom_boxplot(aes(x= method, y=width, fill=method),outlier.size = -1,
               position="dodge") +
  ggtitle("Width of confidence interval for PITE") -> pite_width
pite_width
dev.off()




# head(results.score)

results.score %>%
  group_by(method, n, nsim) %>%
  summarize(sensitivity.ll = sum(Dx.ll<0    & pite.full<0)/sum(pite.full<0),
            sensitivity    = sum(Dx   <0    & pite.full<0)/sum(pite.full<0),
            sensitivity.ul = sum(Dx.ul<0    & pite.full<0)/sum(pite.full<0),
            specificity.ll = sum(!(Dx.ll<0) & !(pite.full<0))/sum(!(pite.full<0)),
            specificity    = sum(!(Dx   <0) & !(pite.full<0))/sum(!(pite.full<0)),
            specificity.ul = sum(!(Dx.ul<0) & !(pite.full<0))/sum(!(pite.full<0))) %>%
  ungroup() %>%
  reshape2::melt(id.vars = c("method", "n", "nsim")) %>%
  group_by(method, n, variable) %>%
  summarize(value = mean(value, na.rm=T)) -> ss.data2
ss.data2$what = rep(c("Sensitivity", "Specificity"), each=3)
ss.data2$who  = factor(rep(c("bold(hat(B)[l])", "bold(hat(B))", "bold(hat(B)[u])"), 2),
                       levels = c("bold(hat(B)[l])", "bold(hat(B))", "bold(hat(B)[u])"))
head(ss.data2)



ss.data2 %>%
  filter(method != "reduced")%>%
  ggplot() +
  # facet_wrap(~n, nrow = 1)+
  # geom_hline(yintercept=0) +
  geom_point(aes(x = n, y=value, color=method, group=method, shape = method), size=1) +
  geom_line(aes(x = n, y=value, color=method, group=method), linetype = 2, size=0.15) +
  facet_grid(what ~ who, switch = "y", labeller = label_parsed) +
  scale_x_continuous(breaks = c(0,n), limits = c(0,1000)) +
  scale_y_continuous(name = "", limits=c(0,1.01)) + #, expand = c(0,0)) +
  labs(color = "Method", shape = "Method", x = "Sample Size") +
  scale_color_manual(name = "Method", values = myPal) +
  scale_fill_manual(name = "Method", values = myPal) +
  scale_shape_manual(name = "Method", values = myShape) +
  theme_bw() +
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
        strip.placement = "outside") -> ss.plot
# ss.plot

ggsave(filename = paste0(path.plot,"pite-sensitivity.eps"),
       plot = ss.plot,
       width = 17, height = 10, unit = "cm")
ggsave(filename = paste0(path.plot,"pite-sensitivity.pdf"),
       plot = ss.plot,
       width = 17, height = 10, unit = "cm")

pdf(file = paste0(path.plot,"all-plots.pdf"),
    width = 5.5, height = 4)
ss.plot
pite_cover + ggtitle("Coverage of Confidence Intervals")
pite_width
pite_bias
pite_mse
vars_inclusion + ggtitle("Percent of inclusion by variable in the Lasso")
dev.off()
