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
  c1 <- segs$c1.mean;
  stopifnot(!is.null(c1));
  n <- segs$dh.num.mark;

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
}) # estimateKappaByC1Density()


setMethodS3("estimateKappa", "PairedPSCBS", function(this, flavor=c("density(C1)"), ...) {
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  estimateKappaByC1Density(this, ...);
})


#############################################################################
# HISTORY:
# 2011-02-03
# o Added estimateKappa().
# o Added estimateKappaByC1Density().
# o Created.
#############################################################################
