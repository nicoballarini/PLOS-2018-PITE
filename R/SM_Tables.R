# Create function to use in the supplementary material
#
# The function gathers the results from the simulations and creates a
# summary table.
do.table <- function(results.test, case.){
  require(dplyr)
  # Calculate the median width for the confidence intervals
  results.test %>%
    filter(n %in% (N / 2)) %>%
    filter(case %in% case.) %>%
    group_by(case, method, n) %>%
    summarize(w = median(w, na.rm = T)) %>%
    left_join(methods, by = "method") %>%
    arrange(case, n) %>%
    mutate(N = n * 2) %>%
    ungroup() -> results.test1_medians


  # Calculate the mean values for the bias, mse, etc
  results.test %>%
    filter(n %in% (N / 2)) %>%
    filter(case %in% case.) %>%
    group_by(case, method, n) %>%
    summarize_all(.funs = function(x) mean(x, na.rm = TRUE)) %>%
    left_join(methods, by = "method") %>%
    arrange(case, n) %>%
    mutate(prop.dx.ll = prop.dx.ll * 100,
           prop.dx    = prop.dx * 100,
           prop.dx.ul = prop.dx.ul * 100,
           prop       = prop * 100,
           N = n * 2) %>%
    ungroup() %>%
    mutate(sqrt.mse = sqrt(mse),
           sqrt.mse.true = sqrt(mse.true)) -> results.test1

  # Replace mean width with median width
  results.test1$w = results.test1_medians$w

  # Select only variables of interest
  results.test1 %>%
    group_by(case) %>%
    mutate(Case = case) %>%
    ungroup() %>%
    dplyr::select(-case) %>%
    dplyr::select(Case, N, Method,
                  bias.true, sqrt.mse.true, w,
                  prop.dx.ll,
                  prop.dx,
                  prop.dx.ul,
                  prop) -> Table3

  # Put columnames
  colnames(Table3) <- c("Case", "N","Method",
                        "Bias", "$\\sqrt{MSE}$", "Width",
                        "\\% in $\\hat{B}_{l}$",
                        "\\% in $\\hat{B}$",
                        "\\% in $\\hat{B}_{u}$",
                        "\\% in B")
  # Calculate number of different methods used
  nmethods = length(levels(unique(Table3$Method)))

  # Create caption for latex table
  caption  = paste0("\\small Diagnostic measures for case ", case.,
  " with ", n_biom, " biomarkers.
  Columns 4 and 5 show the bias and the $\\sqrt{MSE}$ for the point estimate of the PITE. The sixth column shows the median width of the confidence intervals for the PITE, and the last columns show the proportion of subjects in the identified subgroup when considering the using the limits of the confidence intervals and the point estimates. Since methods reduced and reduced-Scheffe have the same point estimate, bias and MSE are equal.")

  # Print Latex table
  knitr::kable(Table3, digits = 2, caption = caption)
}
