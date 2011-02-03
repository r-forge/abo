setMethodS3("estimateStddevByMeanForCACB", "PairedPSCBS", function(fit, ..., tauDH=c(0,1), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'tauDH':
  tauDH <- Arguments$getDoubles(tauDH, range=c(0,1), length=c(2,2));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 

  verbose && enter(verbose, "Estimating linear relationship between std.dev. and mean levels for (CA,CB)");

  # Estimate (CA,CB) means and variances per segment?
  if (is.null(fit$output$ca.avg)) {
    verbose && enter(verbose, "Calculating (CA,CB) statistics per segment");
    fit <- applyByRegion(fit, FUN=.addCACBWithStatitics, ..., verbose=less(verbose, 5));
    verbose && exit(verbose);
  }

  # Stratify by DH?
  if (tauDH[1] > 0 || tauDH[2] < 1) {
    verbose && enter(verbose, "Stratifying by DH");
    verbose && cat(verbose, "tauDH:");
    verbose && print(verbose, tauDH);
    keep <- with(fit$output, dh.mean >= tauDH[1] & dh.mean <= tauDH[2]);
    keep <- is.na(keep) | keep;
    fit <- extractByRegions(fit, which(keep));
    verbose && exit(verbose);
  }

  # Extract signals of interest
  data <- with(fit$output, data.frame(
    mu    = c(ca.avg, cb.avg),
    sigma = c(ca.sd, cb.sd),
    rho   = rep(cacb.cor, times=2),
    dh    = rep(cb.avg/(ca.avg+cb.avg), times=2),
    n     = rep(dh.num.mark, times=2),
    type  = rep(c("ca","cb"), each=length(ca.avg))
  ));


  n <- dh <- NULL; rm(n, dh); # To please R CMD check
  attachLocally(data);

  verbose && enter(verbose, "Fitting a linear relationship (with zero offset) robustly");
  # NOTE: Intercept seems to be unnecessary (very small if estimated)!
  #       /HB 2011-01-27
  w <- sqrt(n);  # Weight by size of region
  w <- w / dh;   # Downweight by DH
  w <- w / sum(w, na.rm=TRUE);
  fit <- MASS::rlm(sigma ~ mu - 1, weights=w);
  verbose && print(verbose, fit);
  verbose && exit(verbose);

  params <- c(a=0, b=coef(fit)[1]);
  attr(params, "modelFit") <- fit;

  verbose && exit(verbose);

  params;
}) # estimateStddevByMeanForCACB()


#############################################################################
# HISTORY:
# 2011-01-27
# o Added estimateStddevByMeanForCACB() for PairedPSCBS.
# o Created.
#############################################################################
