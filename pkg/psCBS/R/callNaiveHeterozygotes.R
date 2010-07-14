###########################################################################/**
# @RdocDefault callNaiveHeterozygotes
#
# @title "Calls heterozygous SNPs from allele-specific signals"
#
# \description{
#  @get "title" (yA,yB) using the method described in [1].
# }
#
# @synopsis
#
# \arguments{
#   \item{yA,yB}{@numeric @vectors of length J containing 
#     allele A and allele B signals.}
#   \item{range}{A @numeric @vector of length 2 specifying the range 
#     of possible values of the estimated cutoff parameter L (see below).}
#   \item{stepsize.within}{A @numeric specifying the "precision" of L.}
#   \item{minimum.window}{A @numeric specifying the minimum window size
#     used for the objective function.}
#   \item{stepsize.window}{A @numeric specifying the size of each window
#     increment.}
#   \item{...}{Not used.}
#   \item{verbose}{A @logical specifying if verbose output should be
#     given or not.}
# }
#
# \value{
#   Returns a @logical @vector of J calls.
#   Attribute \code{fit} contains the estimated parameters.
# }
# 
# \details{
#   This function calls SNPs that are heterozygous, here heterzygous
#   mean non-homozygous.  A homozygous SNP is a SNP where one of the 
#   two alleles has zero signal/copy.  Note that with this definition,
#   a heterozygous SNP can be any SNP where both alleles have non-zero
#   signal, regardless if they have the same signal or not.  
#   For diploid SNPs, this means that AA and BB SNPs are called 
#   non-heterozygotes (homozygotes) and AB SNPs are called heterozygotes.
#   For non-diploid SNPs, say, AAA, AAB, ABB, and BBB, the 2nd and 3rd
#   are called heterzygotes, the other two not.
#   This function was designed to be able to call SNPs that are 
#   heterozygotes also in tumors.  
#   Note that the function does \emph{not} call genotypes.
#
#   This function finds the optimal window to separate homozygous from
#   heterozygous based on the minimum of the A and B alleles.  
#   It starts with window [range[1], range[1]+minimum.window].
#   It shifts window by stepsize.within.  The process is repeated with bigger
#   and bigger windows where the size is incremented by stepsize.window.
# }
#
# \section{Missing values}{
#   TO DO.
#   Ideally we should return @NA for SNPs that contain a missing value.
# }
#
# @examples "../incl/callNaiveHeterozygotes.Rex"
#
# \references{
#   [1] AB Olshen, RA Olshen, H Bengtsson, P Neuvial, P Spellman, 
#       and VE. Seshan, \emph{Parent-specific copy number in paired 
#       tumor-normal studies using circular binary segmentation},
#       (submitted), 2010.
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("callNaiveHeterozygotes", "default", function(yA, yB, range=c(0,0.5), stepsize.within=0.01, minimum.window=0.05, stepsize.window=0.01, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfLoci <- length(yA);
  if (length(yB) != nbrOfLoci) {
    stop("Arguments 'yA' and 'yB' are of different lengths: ", nbrOfLoci, " != ", length(yB));
  }

  # Argument 'range':
  if (range[2] < range[1]) {
    rangeStr <- paste("[", range[1], ",", range[2], "]", sep="");
    stop("Argument 'range' has an upper bound that is smaller than the lower bound: ", rangeStr);
  }

  # Argument 'minimum.window', 'upper.bound' & 'lower.bound':
  if (minimum.window > range[2]-range[1]) {
    rangeStr <- paste("[", range[1], ",", range[2], "]", sep="");
    stop("Argument 'minimum.window' is too big for the current range ", 
          rangeStr, ": ", minimum.window);
  }

  # Argument 'verbose':
  verbose <- as.logical(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate the minimum of yA and yB for each SNP
  #
  # PSCBS paper:
  # M_i = min(A_i,B_i)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## SLOW: yMin <- apply(cbind(yA, yB), MARGIN=1, FUN=min);
  yMin <- yA;
  idxs <- which(yB < yA);
  yMin[idxs] <- yB[idxs];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find cutoff
  #
  # PSCBS paper: 
  # Find L that minimize Q = (#i: u <= M_i <= v)/(v-u), where
  # [u,v] is a window with v-u >= 0.05
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Values for optimal Q and [u,v]
  uOpt <- vOpt <- as.double(NA);
  Qopt <- Inf;

  # Intital window size (v-u)
  windowSize <- minimum.window;

  # Start with the window "to the very left".
  u <- range[1];
  v <- u + windowSize;
  while (v <= range[2]) {
    if (verbose) {
      printf("Trying with window [%.3f,%.3f]\n", u,v);
    }

    # For the given windows size (v-u), find optimal [u,v] window
    while (v <= range[2]) {
      if (verbose) {
        printf("Current window [%.3f,%.3f]\n", u,v);
      }

      # Calculate (#i: u <= M_i <= v)/(v-u)
      Q <- length(which(u <= yMin & yMin <= v))/windowSize;

      # New optimal [u,v] window found?
      if (Q < Qopt) {
        Qopt <- Q;
        uOpt <- u;
        vOpt <- v;
      }

      # Shift window [u,v] "to the right".
      u <- u + stepsize.within;
      v <- v + stepsize.within;
    } # while (...)

    # Increase window size
    windowSize <- windowSize + stepsize.window;

    # Restart with the window "to the very left".
    u <- range[1];
    v <- u + windowSize;
  } # while (...)

  # Not needed anymore
  rm(u, v, Q, Qopt, windowSize);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate proportion of heterozygotes outside the optimal window [u,v]:
  # (#i: M_i > v)/[(#i:  M_i > v) + (#i: M_i < u)]
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  isAmbiguous <- (uOpt <= yMin & yMin <= vOpt);
  if (any(isAmbiguous)) {
    nbrOfLarge <- sum(yMin > vOpt, na.rm=TRUE);
    nbrOfSmall <- sum(yMin < uOpt, na.rm=TRUE);
    proportionHets <- nbrOfLarge / (nbrOfLarge + nbrOfSmall);

    if (verbose) {
      printf("Proportion of heterozygous SNPs outside window: %.3f\n", proportionHets);
    }

    min.ambigous <- sort(yMin[isAmbiguous]);
    idx <- round(proportionHets * length(min.ambigous));
    L <- min.ambigous[idx];
  } else {
    L <- (vOpt+uOpt) / 2;
  }

  if (verbose) {
    printf("Estimated heterozygous cutoff: %.3f\n", L);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call heterozygous SNPs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # PSCBS paper: Call SNP i heterozygous if M_i > L.
  isHet <- rep(NA, times=nbrOfLoci);
  isHet[yMin  > L] <- TRUE;
  isHet[yMin <= L] <- FALSE;

  # Return results
  fit <- list(heterozygous.cutoff=L);
  attr(isHet, "fit") <- fit;
  
  isHet;
}) # callNaiveHeterozygotes()


############################################################################
# HISTORY:
# 2010-07-14
# o Created from findheterozygous.R.
############################################################################
