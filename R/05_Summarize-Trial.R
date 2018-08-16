#' Summarize the results from the scoring funcitons.
#'
#' Calculates the width of the confidence intervals for pite, their coverage,
#' and the bias and mse of the point estimates.
#' Also calculates the sensitivity and specifity of subgroups.
#'
#' @param scores The dataset with the PITE and their confidence intervals
#' @param dataset The dataset with the clinical trial data
#' @param n.pite number of subjects to consider to summarize
#'
#' @examples
#' input <- MakeInput()
#' parameters <- MakeParameters(input = input)
#' dataset <-OneData(input,parameters)
#' score.lm(dataset, input)
#'
#' @export
summarize_scores <- function(scores, dataset, n.pite = 40) {
  n <- n.pite
  true = dataset[1:n, "TrueTrtEff"]
  scores = scores[1:n.pite, ]

  cover      <- mean(scores$Dx.cover)
  cover.true <- mean((scores$Dx.ll <= true)&(true <= scores$Dx.ul))

  w     <- mean(scores$Dx.w)
  bias  <- mean(scores$Dx - scores$PITE)
  mse   <- mean((scores$Dx - scores$PITE)^2)

  bias.true  <- mean(scores$Dx - true)
  mse.true   <- mean((scores$Dx - true)^2)

  true.positive <- dataset[1:n, "TrueTrtEff"]>0
  true.negative <- !true.positive
  sensitivity.dx <- sum(scores[1:n, "Dx"] > 0      & (true.positive)) / sum(true.positive)
  sensitivity.ll <- sum(scores[1:n, "Dx.ll"] > 0   & (true.positive)) / sum(true.positive)
  sensitivity.ul <- sum(scores[1:n, "Dx.ul"] > 0   & (true.positive)) / sum(true.positive)
  specificity.dx <- sum(scores[1:n, "Dx"] <= 0     & (true.negative)) / sum(true.negative)
  specificity.ll <- sum(scores[1:n, "Dx.ll"] <= 0  & (true.negative)) / sum(true.negative)
  specificity.ul <- sum(scores[1:n, " Dx.ul"] <= 0 & (true.negative)) / sum(true.negative)

  prop       <- sum(true.positive) / n
  prop.dx    <- sum(scores[1:n, "Dx"] > 0) / n
  prop.dx.ll <- sum(scores[1:n, "Dx.ll"] > 0) / n
  prop.dx.ul <- sum(scores[1:n, "Dx.ul"] > 0) / n

  out   <- data.frame(cover = cover,
                      w = w,
                      bias = bias,
                      mse = mse,
                      bias.true = bias.true,
                      mse.true = mse.true,
                      sensitivity.dx, sensitivity.ll, sensitivity.ul,
                      specificity.dx, specificity.ll, specificity.ul,
                      prop = prop,
                      prop.dx    = prop.dx,
                      prop.dx.ll = prop.dx.ll,
                      prop.dx.ul = prop.dx.ul,
                      Dx = scores[1:n, "Dx"],
                      Dx.ll = scores[1:n, "Dx.ll"],
                      Dx.ul = scores[1:n, "Dx.ul"],
                      PITE = scores$PITE,
                      true = true)
  out
}
