#a.vector - vector of A alleles
#b.vector - vector of B alleles
#lower.bound - lower bound of optimization
#upper.bound=0.25 - upper bound of optimization
#minimum.window - minimum window size
#stepsize.window - how much to increment the overall window size
#stepsize.within - how much to increment the current window size
#Find the optimal window to separate homozygous from heterozygous based on the minimum of the A and B alleles.  Start with window of lower.bound to lower.bound + minimum.window.  Shift window by stepsize.within.  Repeat process with bigger windows starting from lower.bound to lower.bound + stepsize.window
findheterozygous <- function(a.vector, b.vector, lower.bound=0, upper.bound=0.5, minimum.window=0.05, stepsize.window=0.01, stepsize.within=0.01) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'minimum.window', 'upper.bound' & 'lower.bound':
  if (minimum.window > (upper.bound-lower.bound)) {
    stop("Minimum window size is too big for the current lower.bound and upper.bound")
  }


  min.vector <- apply(cbind(a.vector, b.vector), MARGIN=1, FUN=min)
  current.fixed.lower <- lower.bound
  current.lower <- lower.bound
  current.upper <- lower.bound+minimum.window
  current.window.size <- minimum.window
  current.stat <- Inf
  current.bottom <- NA
  current.top <- NA
  while (current.upper <= upper.bound) {
    while (current.upper <= upper.bound) {
      new.stat <- length(which(min.vector >= current.lower & min.vector <= current.upper))/current.window.size
      if (new.stat < current.stat) {
        current.stat <- new.stat
        current.bottom <- current.lower
        current.top <- current.upper
      }

      current.lower <- current.lower+stepsize.within
      current.upper <- current.upper+stepsize.within
    } # while (...)

    current.lower <- current.fixed.lower
    current.window.size <- current.window.size+stepsize.window
    current.upper <- current.lower+current.window.size
  } # while (...)

  is.heterozygous <- rep(NA, times=length(min.vector))
  which.ambiguous <- which(min.vector >= current.bottom & min.vector <= current.top)
  length.ambiguous <- length(which.ambiguous)
  if (length.ambiguous > 0) {
    min.ambigous <- sort(min.vector[which.ambiguous])
    proportion.heterozygous <- length(which(min.vector > current.top)) / (length(which(min.vector > current.top)) + length(which(min.vector < current.bottom)))
    heterozygous.cutoff <- min.ambigous[round(proportion.heterozygous*length.ambiguous)]
  } else {
    heterozygous.cutoff <- (current.top+current.bottom)/2
  }

  is.heterozygous[which(min.vector > heterozygous.cutoff)] <- TRUE
  is.heterozygous[which(min.vector <= heterozygous.cutoff)] <- FALSE
  print(paste("Heterozygous cutoff is", round(heterozygous.cutoff, digits=3)))

  list(is.heterozygous=is.heterozygous, heterozygous.cutoff=heterozygous.cutoff)
} # findheterozygous()

