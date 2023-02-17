#' Define function to correct the effect of distance
#' on loop enrichment
#' @param x GInteractions Object.
#' @param scores Numeric vector of enrichment scores.
#' @param k Number of observations for rolling window.
#' @param plot Boolean (default=FALSE), of whether to
#'  show diagnostic plot.
#' @importFrom stats smooth.spline
#' @importFrom InteractionSet pairdist
#' @returns Vector of corrected loop enrichment scores.
#' @noRd
.correctEnrich <- function(x, scores, k=25, nknots=10, plot=FALSE) {
  ## Calculate rolling enrichment
  re <- .rollEnrich(x, scores=scores, k=k)
  
  ## Fit smoothed spline
  sp <- smooth.spline(na.omit(re), nknots=nknots)
  m <- median(na.omit(re$rollMedScore))
  
  ## Calculate correction factor
  ratio <- m / predict(sp, pairdist(x))$y
  
  ## Corrected scores
  corrected <- scores*ratio
  
  ## Optional diagnostic plot
  if (plot) {
    plot(re$rollMedSize, re$rollMedScore, type='l',
         main=paste0("k=", k, ", nknots=", nknots), col='lightgrey')
    lines(sp$x, sp$y, col="forestgreen")
    abline(h=m, lty=2)
    cre <- .rollEnrich(x, scores=corrected, k=k)
    lines(cre$rollMedSize, cre$rollMedScore, col="blue")
  }
  
  ## Return corrected scores
  return(corrected)
}