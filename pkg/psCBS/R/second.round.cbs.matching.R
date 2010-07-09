second.round.cbs.matching <- function(diff.data, alpha) {
  nbrOfLoci <- length(diff.data);
  chrom <- rep(1L, times=nbrOfLoci);
  maploc <- seq(length=nbrOfLoci);
  new.cna <- CNA(diff.data, chrom=chrom, maploc=maploc);
  new.segment <- segment(new.cna, alpha=alpha)$output;

  nbrOfRegions <- nrow(new.segment);
  if (nbrOfRegions == 1) {
    output <- NULL;
  } else {
    output <- cumsum(new.segment[,"num.mark"])[-nbrOfRegions];
  }

  output;
} # second.round.cbs.matching()


############################################################################
# HISTORY:
# 2010-07-08
# o ROBUSTNESS: Segmentation results are now subsetted by names not indices.
# o ROBUSTNESS: Replaced all 1:n with seq(length=n) to deal with n == 0.
############################################################################
