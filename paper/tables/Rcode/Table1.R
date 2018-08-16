library(dplyr)
rproj <- rprojroot::find_package_root_file()
path.tables <- paste0(rproj, "/paper/tables/")

# Calculate effecs correspondence with Schnell's paper
effects <- expand.grid(b= c(0,1/4), d1=c(0,1/2), d2=c(0,1/3))[c(1,5,3,7,2,8),]
all.ref <- cbind(case=1:6, a=0,effects["b"], g1=0, g2=0, effects[c('d1','d2')])
knitr::kable(all.ref, row.names = F)
all.eff <- data.frame(case=1:6,
                      a = all.ref["a"]+all.ref["b"]/2+all.ref["g1"]/2+all.ref["d1"]/4,
                      b = all.ref["b"]/2+all.ref["d1"]/4,
                      g1 = all.ref["g1"]/2+all.ref["d1"]/4,
                      g2 = sqrt(3)*(all.ref["g2"]+all.ref["d2"]/2),
                      d1 = all.ref["d1"]/4,
                      d2 = sqrt(3)*all.ref["d2"]/2)
knitr::kable(all.eff, row.names = F)
# rownames(all.eff)<-1:6


xtable::xtable(all.eff, row.names = F,
               caption = "Scenarios considered for simulations",
               sanitize.text.function = function(x){x}) %>%
  print(include.rownames=F)

expl <- c("No predictive markers",
          "$X_2$ is predictive",
          "$X_1$ is predictive",
          "$X_1$ and $X_2$ are predictive",
          "Overall treatment effect",
          "Overall treatment effect and $X_1$ and $X_2$ are predictive")

colnames. <- c("Case",
               "$\\alpha$", "$\\beta$",
               "$\\gamma_1$", "$\\gamma_2$",
               "$\\delta_1$", "$\\delta_2$")
cbind(all.eff) -> eff.table
colnames(eff.table) <- colnames.
eff.table$Case <- paste0(eff.table$Case, ": ", expl)


xtable::xtable(eff.table,
               row.names = F,
               caption = "Scenarios considered for simulations",
        digits = c(0,0,rep(2,6))) -> xtable.eff

xtable::align(xtable.eff) <- "ll|cccccc"

print(xtable.eff,
      include.rownames=F,
      type = "latex",
      table.placement=NULL,
      floating.environment = "table*",
      sanitize.text.function = function(x){x},
      hline.after = c(0, nrow(xtable.eff)),
      caption.placement = "top",
      add.to.row = list(pos = list(-1),
                        command = c("\\Hline \n"))) ->Table1

capture.output(cat(Table1), file = paste0(path.tables,"Table1.tex"))
