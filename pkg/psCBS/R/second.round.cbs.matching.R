second.round.cbs.matching <- function(diff.data, alpha) {
  n <- length(diff.data)
  new.cna <- CNA(diff.data, chrom=rep(1, times=n), maploc=1:n)
  new.segment <- segment(new.cna, alpha=alpha)$output

  if (nrow(new.segment) == 1) {
    output <- NULL
  } else {
    output <- cumsum(new.segment[,5])[-nrow(new.segment)]
  }

  output
} # second.round.cbs.matching()

