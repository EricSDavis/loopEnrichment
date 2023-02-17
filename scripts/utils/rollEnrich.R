#' Calculate rolling windows of loop enrichment
#' @param x GInteractions Object.
#' @param scores Numeric vector of enrichment scores.
#' @param k Number of observations for rolling window.
#' @param thresh Numeric loop size cutoff (keeps loops)
#'  less than this size.
#' @returns A data table with median loop size and enrichment
#'  after calculating the rolling median within a window.
#' @noRd
.rollEnrich <- function(x, scores, k=200, thresh=Inf) {
  
  ## Data table of loop size and score
  dat <- data.table(
    size = pairdist(x),
    scores = scores
  )
  
  ## Order then filter out interactions >thresh
  dat <- dat[order(size)]
  dat <- dat[size <= thresh]
  
  ## Calculate rollMedians for size and score
  FUN <- \(x) median(x)
  dat[, rollMedSize := frollapply(x=size, n=n, FUN=FUN)]
  dat[, rollMedScore := frollapply(x=scores, n=n, FUN=FUN)]
  
  ## Take the median per group
  ans <- dat[, .(rollMedScore = median(rollMedScore)), by=rollMedSize]
  ans
}
