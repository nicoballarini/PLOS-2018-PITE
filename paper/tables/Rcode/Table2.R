################################################################################
rm(list=ls())
cat("\014")
################################################################################
library(reshape2)
library(ggplot2)
library(dplyr)

# These tables may contain a command addrowspace that is defined as follows
# \newcommand{\addrowspace}{\vspace{-0.25cm}\\}

N = c(40, 100, 220, 350)
sel.cases = c(1,4)
nbiom=10
save = TRUE
latex = TRUE

################################################################################
# BIAS Table
################################################################################

  rproj <- rprojroot::find_package_root_file()
  path <- sprintf("/sim/normal/nbiom%s/",nbiom)
  path <- paste0(rproj, path)
  path
  path.data <- paste0(path, "results/data/")
  path.plot <- paste0(path, "results/plots/")

  path.tables <- paste0(rproj, "/paper/tables/")
  path.tables

  files.2 <- list.files(path=path.data)
  ## test DATA ----
  files.2.test <- files.2[grep(files.2, pattern = "-test")]
  results.test <- data.frame()
  for(i in files.2.test){
    listdata2 = readRDS(paste0(path.data,i))
    results.test <- rbind(results.test, data.frame(listdata2, scenario=i))
  }
  # head(results.test)
  methods <- data.frame(method = c("null","lm","lasso","mlm","sch","an1", 'an2'),
                        Method = c("ATE", "full", "Lasso","reduced", "reduced-Scheffe", "rLasso-1", "rLasso-2"))
  # results.test
  pos = do.call(rbind,gregexpr('-', results.test$scenario))
  results.test$case <- as.numeric(substr(results.test$scenario,
                                          start = pos[,2] - 1, stop = pos[,2] - 1))

  results.test$n <- as.numeric(substr(results.test$scenario,
                                       start = pos[,2] + 10, stop = pos[,3] - 1))

  # Calculate the median width for the confidence intervals
  results.test %>%
    filter(n %in% (N / 2)) %>%
    filter(case %in% sel.cases) %>%
    group_by(case, method, n) %>%
    summarize(w = median(w, na.rm = T)) %>%
    left_join(methods, by = "method") %>%
    arrange(case, n) %>%
    ungroup() -> results.test1_medians

  results.test %>%
    filter(n %in% (N/2)) %>%
    filter(case %in% sel.cases) %>%
    group_by(case, method, n) %>%
    summarize_all(.funs = function(x) mean(x, na.rm=T)) %>%
    left_join(methods, by = "method") %>%
      arrange(case, n) %>%
    mutate(prop.dx.ll = prop.dx.ll*100,
           prop.dx    = prop.dx*100,
           prop.dx.ul = prop.dx.ul*100,
           prop       = prop*100,
           N = n*2) %>%
    ungroup() %>%
    mutate(sqrt.mse = sqrt(mse),
           sqrt.mse.true = sqrt(mse.true)) -> results.test1

  # Replace mean width with median width
  results.test1$w = results.test1_medians$w

  results.test1 %>%
    dplyr::select(case, N, Method,
                  bias.true, sqrt.mse.true, w,
                  prop.dx.ll,
                  prop.dx,
                  prop.dx.ul,
                  prop) -> results.test2
  results.test2

  Case.lookup = data.frame(Case = c(1, 4),
                           case = c("Null", "Predictive"))
  results.test2 -> Table3
  Table3 %>%
    group_by(case) %>%
    mutate(Case = ifelse(row_number()==1, case, NA)) %>%
    ungroup() %>%
    dplyr::select(-case) %>%
    left_join(Case.lookup) %>%
    dplyr::select(case, N, Method,
                  bias.true, sqrt.mse.true, w,
                  prop.dx.ll,
                  prop.dx,
                  prop.dx.ul,
                  prop) -> Table3

  colnames(Table3) <- c("Case", "N","Method",
                        "Bias", "$\\sqrt{MSE}$", "Width",
                        "\\% in $\\hat{B}_{l}$",
                        "\\% in $\\hat{B}$",
                        "\\% in $\\hat{B}_{u}$",
                        "\\% in B")
  Table3
  # cases = max(unique(Table3$Case))
  nmethods = length(levels(unique(Table3$Method)))
  caption  = "Diagnostic measures for 1000 simulations in each null and predictive case.
  Columns 4 and 5 show the bias and the $sqrt(MSE)$ for the point estimate of the PITE.
  The sixth column shows the width of the confidence intervals for the PITE, and the last columns
show the proportion of subjects in the identified subgroup when considering the
using the limits of the confidence intervals and the point estimates.
  Since methods reduced and reduced-Scheffe have the same point estimate, bias and MSE are equal."

  if (latex) {
    Table3.xtable <- xtable::xtable(Table3, digits = c(0,0,0,0,2,2,2,1,1,1,1),
                                    caption = caption)
    xtable::align(Table3.xtable) <- "clrlrrrrrrr"
    xtable::label(Table3.xtable) = "table2"
    print(Table3.xtable,
          floating.environment = "table",
          include.colnames = TRUE,
          include.rownames = FALSE,
          table.placement  = NULL,
          sanitize.text.function = function(x){x},
          size = "\\small",
          hline.after = nrow(Table3.xtable),
          caption.placement = "top",
          add.to.row = list(pos = as.list(c(-1,-1, 0, (1:(2*length(N)-1)*nmethods))),
                            command = c("\\hline \n \\multicolumn{6}{c}{} & \\multicolumn{4}{c}{Proportion of subjects in subgroup} \\\\ \n \\cline{7-10} \n",
                                        " \\addrowspace ", " \\hline\\addrowspace ",
                                        rep("\\addrowspace ", 2*length(N)-1))))  -> Table3.tex
  }
  if (save) {
    # savepath <- sprintf("%s/sim/nbiom%s/TableSS-nbiom%s-n%s.tex", rproj, nbiom, nbiom, N)
    capture.output(cat(Table3.tex), file = sprintf("%sTable2.tex", path.tables))
    cat("\n\n", paste0("Saving table in ", path.tables), "\n\n")
  }

  # knitr::kable(Table3, digits= 2)
  # knitr::kable(Table3 %>% filter(Method=="Lasso"))
  # knitr::kable(Table3 %>% filter(Method=="rLasso-1"))
  # knitr::kable(Table3 %>% filter(Method=="rLasso-2"))
  # knitr::kable(results.test1 %>% filter(Method=="Lasso"))
