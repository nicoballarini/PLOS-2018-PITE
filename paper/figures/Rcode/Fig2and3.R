## 10 biomarkers -----
library(dplyr)
library(ggplot2)
library(PITE)
library(ggpubr)
rm(list = ls())
cat("\014")
getwd()
nbiom=10
rproj <- rprojroot::find_package_root_file()
path.figures <- paste0(rproj, "/paper/figures/")

path <- sprintf("/sim/normal/nbiom%s/",nbiom)
path <- paste0(rproj, path)
# if (!dir.exists(paste0(path, "results/"))) dir.create(paste0(path, "results/"))
# if (!dir.exists(paste0(path, "results/data/"))) dir.create(paste0(path, "results/data/"))
if (!dir.exists(paste0(path, "results/plots/"))) dir.create(paste0(path, "results/plots/"))
path.data <- paste0(path, "results/data/")
path.plot <- paste0(path, "results/plots/")
cat("Getting data from: ", path.data, "\n")
cat("Saving plots to: ", path.plot, "\n")
files.2 <- list.files(path=path.data)

## PLOTS! ----------------------------------------------------------------------
##### Plot options --------------------------------------------------------
myPal <- c("black","#EE4035", "#F3A530", "#56B949", "#82b6cc","#0f5e7f","#024460")
myShape <- c(NA,16,17,15,3,12,13)
methods <- data.frame(method = c("null","lm","lasso","mlm","sch","an1", 'an2'),
                      Method = factor(c("ATE", "full", "Lasso","reduced", "reduced-Scheffe", "rLasso-1", "rLasso-2"),
                                      levels=c("ATE", "full", "Lasso","reduced", "reduced-Scheffe", "rLasso-1", "rLasso-2")))
methods
names(myPal) <- methods$Method
names(myShape) <- methods$Method
colorize.fill <- scale_fill_manual(name = "Method", values = myPal)
colorize      <- scale_color_manual(name = "Method", values = myPal)

outlier.size = -1

plotsens <- function(results, case.=4){
  sample.size <- c(unique(results$n)*2)
  results %>%
    group_by(case,n, method) %>%
    filter(case==case.) %>%
    summarise(Sensitivity.ll = mean(sensitivity.ll, na.rm = T),
              Sensitivity.dx = mean(sensitivity.dx, na.rm = T),
              Sensitivity.ul = mean(sensitivity.ul, na.rm = T),
              Specificity.ll = mean(specificity.ll, na.rm = T),
              Specificity.dx = mean(specificity.dx, na.rm = T),
              Specificity.ul = mean(specificity.ul, na.rm = T)) %>%
    reshape2::melt(c("case","n","method")) %>%
    mutate(which = factor(substr(variable,start = 13,stop = 14),
                          levels = c("ll","dx","ul")),
           measure = factor(substr(variable,start = 1,stop = 11),
                            levels = c("Sensitivity", "Specificity"))) %>%
    rename(sensitivity = value) -> datplot

  levels(datplot$which) = c("bold(hat(B)[l])", "bold(hat(B))", "bold(hat(B)[u])")

  if (case.==1){
    datplot %>% filter(measure !="Sensitivity") -> datplot
  }
  datplot %>%
    filter(method %in% c("ATE", "full", "Lasso", "reduced-Scheffe", "rLasso-1", "rLasso-2")) %>%
    ggplot() +
    geom_line(aes(y=sensitivity,  color = method, x=2*n), linetype  = 2, size=0.15) +
    geom_point(aes(y=sensitivity,  color = method, shape=method, x=2*n), size=1) +
    colorize.fill + colorize +
    scale_shape_manual(name = "Method", values = myShape) +
    scale_y_continuous(name = "", limits = c(0,1.01)) +
    scale_x_continuous(name = "Sample Size (n)",
                       limits = c(0,350), breaks = c(0,sample.size)) -> p1
  p1 + facet_grid(measure ~ which, switch = "y", labeller = label_parsed) +
    labs(color = "Method", shape = "Method", x = "Sample Size (n)") +
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
          strip.placement = "outside")
}

## ----
## ----
## TEST Data ----

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

results.test$n <- as.numeric(substr(results.test$scenario,
                                    start = pos[,2] + 10, stop = pos[,3] - 1))

results.test %>%
  left_join(methods, by = "method") %>%
  mutate(method = Method)-> results.test

results.test %>% head()

results.test %>%
  group_by(n,case, method) %>%
  summarize(cover = mean(cover)) %>%
  arrange(case,n) -> results.test2
results.test2


## Figures for publication --------
## Figure 2 -----------------------
## ***Coverage for PITE  ----------
nsim=1000
alpha=0.05
p=alpha
q=1-p
e=1.96*sqrt(p*q/nsim)
results.test %>%
  group_by(n,case, method) %>%
  summarize(coverage = mean(cover, na.rm=T)) %>%
  arrange(case,n) -> data.coverage.pite
head(data.coverage.pite)

data.coverage.pite %>%
  mutate(N = 2*n) %>%
  filter(case %in% c(1)) %>%
  ggplot() +
  geom_line(aes(x = N, y=coverage, color=method,  group=method), linetype = 2, size=0.15) +
  geom_point(aes(x = N, y=coverage, color=method, shape = method, #fill = method,
                 group=method), size = 1) +
  geom_hline(yintercept=1-alpha+e, linetype=3, size=0.1) +
  geom_hline(yintercept=1-alpha, size=0.1) +
  geom_hline(yintercept=1-alpha-e, linetype=3, size=0.1) +
  scale_x_continuous(breaks = c(0,40,100,220,350), limits = c(0,350)) +
  scale_y_continuous(limits = c(0.65,1)) +
  scale_color_manual(name = "Method", values = myPal) +
  scale_shape_manual(name = "Method", values = myShape) +
  guides(shape = guide_legend(title = "Method")) +
  labs(y="Coverage", x= "Sample size (n)") +
  theme_bw()  +
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
        legend.position = "none",
        legend.key.height=unit(0.5,"line"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        axis.text    = element_text(size = 6, color = "black"),
        axis.text.x  = element_text(size = 6, color = "black"),
        strip.placement = "outside") -> pite_cover1
# pite_cover1

data.coverage.pite %>%
  mutate(N = 2*n) %>%
  filter(case %in% c(4)) %>%
  ggplot() +
  geom_line(aes(x = N, y=coverage, color=method,  group=method), linetype = 2, size=0.15) +
  geom_point(aes(x = N, y=coverage, color=method, shape = method, #fill = method,
                 group=method), size = 1) +
  geom_hline(yintercept=1-alpha+e, linetype=3, size=0.1) +
  geom_hline(yintercept=1-alpha, size=0.1) +
  geom_hline(yintercept=1-alpha-e, linetype=3, size=0.1) +
  scale_x_continuous(breaks = c(0,40,100,220,350), limits = c(0,350)) +
  scale_y_continuous(limits = c(0.65,1)) +
  scale_color_manual(name = "Method", values = myPal) +
  scale_shape_manual(name = "Method", values = myShape) +
  guides(shape = guide_legend(title = "Method")) +
  labs(y="Coverage", x= "Sample size (n)") +
  theme_bw()  +
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
        legend.position = c(0.81, 0.20),
        legend.key.height=unit(0.5,"line"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        axis.text    = element_text(size = 6, color = "black"),
        axis.text.x  = element_text(size = 6, color = "black"),
        strip.placement = "outside") -> pite_cover4
# pite_cover4

ggarrange(pite_cover1,
          pite_cover4,
          heights = c(1,1),
          labels = c("A", "B"),
          ncol=2, nrow = 1) -> plot2
img.width  = 3.14
img.height = 3.14

ggsave(paste0(path.figures,"Fig2.pdf"),plot = plot2,
       width = 2*img.width, height = img.height, units = "in")
ggsave(paste0(path.figures,"Fig2.eps"),plot = plot2,
       width = 2*img.width, height = img.height, units = "in")


## Figure 3  --------------------------
###***Sensitivity and specificity -----
plotsens(results.test, case. = 1) -> plotA
plotsens(results.test, case. = 4) -> plotB
ggarrange(plotA,
          plotB,
          heights = c(5,8),
          labels = c("A", "B"),
          common.legend = TRUE, legend = "right",
          ncol=1, nrow = 2) -> plot3
ggsave(paste0(path.figures,"Fig3.pdf"), plot = plot3, width = 17, height = 15, units = "cm")
ggsave(paste0(path.figures,"Fig3.eps"), plot = plot3, width = 17, height = 15, units = "cm")

