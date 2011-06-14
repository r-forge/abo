###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod estimateKappa
#
# @title "Estimate global background in segmented copy numbers"
#
# \description{
#  @get "title".
#  The global background, here called \eqn{\kappa},
#  may have multiple origins where normal contamination is one,
#  but not necessarily the only one.
# }
#
# @synopsis
#
# \arguments{
#   \item{flavor}{A @character string specifying which type of
#    estimator to use.}
#   \item{...}{Additional arguments passed to the estimator.}
# }
#
# \value{
#   Returns the background estimate as a @numeric scalar.
# }
#
# @author
#
# \seealso{
#   Internally, one of the following methods are used:
#   @seemethod "estimateKappaByC1Density".
# }
#
#*/###########################################################################  
setMethodS3("estimateKappa", "PairedPSCBS", function(this, flavor=c("density(C1)"), ...) {
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  if (flavor == "density(C1)") {
    estimateKappaByC1Density(this, ...);
  } else {
    throw("Cannot estimate background. Unsupported flavor: ", flavor);
  }
})



###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod estimateKappaByC1Density
#
# @title "Estimate global background in segmented copy numbers"
#
# \description{
#  @get "title" based on the location of peaks in a weighted
#  density estimator of the minor copy number mean levels. 
#  The global background, here called \eqn{\kappa},
#  may have multiple origins where normal contamination is one,
#  but not necessarily the only one.
# }
#
# @synopsis
#
# \arguments{
#   \item{adjust}{A @numeric scale factor specifying the size of
#    the bandwidth parameter used by the density estimator.}
#   \item{minDensity}{A non-negative @numeric threshold specifying
#    the minimum density a peak should have in order to consider
#    it a peak.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the background estimate as a @numeric scalar.
# }
# 
# \section{Algorithm}{
#  \itemize{
#  \item Retrieve segment-level minor copy numbers and corresponding weights:
#   \enumerate{
#    \item Grabs the segment-level C1 estimates.
#    \item Calculate segment weights proportional to the number of heterozygous SNPs.
#   }
#
#  \item Identify subset of regions with C1=0:
#   \enumerate{
#    \item Estimates the weighted empirical density function
#          (truncated at zero below).  Tuning parameter 'adjust'.
#    \item Find the first two peaks
#          (with a density greater than tuning parameter 'minDensity').
#    \item Assumes that the two peaks corresponds to C1=0 and C1=1.
#    \item Defines threshold Delta0.5 as the center location between
#          these two peaks.
#   }
#
#  \item Estimate the normal contamination:
#   \enumerate{
#    \item For all segments with C1 < Delta0.5, calculate the weighted
#          median of their C1:s.
#    \item Let kappa be the above weighted median. 
#          This is the estimated background.
#   }
#  }
# }
#
# @author
#
# \seealso{
#   Instead of calling this method explicitly, it is recommended
#   to use the @seemethod "estimateKappa" method.
# }
#
# @keyword internal
#*/###########################################################################  
setMethodS3("estimateKappaByC1Density", "PairedPSCBS", function(this, adjust=1, minDensity=0.2, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'adjust':
  adjust <- Arguments$getDouble(adjust, range=c(0,Inf));

  # Argument 'minDensity':
  minDensity <- Arguments$getDouble(minDensity, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Estimate normal contamination");
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the region-level estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  segs <- this$output;
  c1 <- segs$c1Mean;
  stopifnot(!is.null(c1));
  n <- segs$dhNbrOfLoci;

  # Drop missing values
  keep <- (!is.na(c1) & !is.na(n));
  c1 <- c1[keep];
  n <- n[keep];

  verbose && cat(verbose, "Number of segments: ", length(c1));

  # Calculate region weights
  weights <- n / sum(n);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Identify subset of regions with C1=0
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Estimating threshold Delta0.5 from the empirical density of C1:s");
  verbose && cat(verbose, "adjust: ", adjust);
  verbose && cat(verbose, "minDensity: ", minDensity);

  d <- density(c1, weights=weights, adjust=adjust, from=0, na.rm=FALSE);
  fit <- findPeaksAndValleys(d);

  type <- NULL; rm(type); # To please R CMD check
  fit <- subset(fit, type == "peak");
  stopifnot(nrow(fit) >= 2);

  fit <- subset(fit, density >= minDensity);
  stopifnot(nrow(fit) >= 2);
  verbose && cat(verbose, "All peaks:");
  verbose && print(verbose, fit);

  # Keep the first two peaks
  fit <- fit[1:2,,drop=FALSE];
  verbose && cat(verbose, "C1=0 and C1=1 peaks:");
  verbose && print(verbose, fit);

  peaks <- fit$x;
  Delta0.5 <- mean(peaks);
  verbose && cat(verbose, "Estimate of Delta0.5: ", Delta0.5);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Estimate the normal contamination
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  keep <- which(c1 < Delta0.5);
  verbose && cat(verbose, "Number of segments with C1 < Delta0.5: ", length(keep));
  kappa <- weightedMedian(c1[keep], w=weights[keep]);

  verbose && cat(verbose, "Estimate of kappa: ", kappa);

  verbose && exit(verbose);

  kappa;
}, protected=TRUE) # estimateKappaByC1Density()




#############################################################################
# HISTORY:
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-04-08
# o Added Rdoc for estimateKappaByC1Density().
# 2011-02-03
# o Added estimateKappa().
# o Added estimateKappaByC1Density().
# o Created.
#############################################################################