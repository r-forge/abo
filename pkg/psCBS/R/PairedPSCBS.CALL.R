setMethodS3("callAllelicBalanceByDH", "PairedPSCBS", function(fit, tau=0.10, alpha=0.05, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  # Argument 'tau':
  tau <- Arguments$getDouble(tau, range=c(0,Inf));

  # Argument 'alpha':
  alpha <- Arguments$getDouble(alpha, range=c(0,1));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling segments of allelic balance from one-sided DH bootstrap confidence intervals");
  verbose && cat(verbose, "tau (offset adjusting for bias in DH): ", tau);
  verbose && cat(verbose, "alpha (CI quantile; significance level): ", alpha);

  # Calculate DH confidence intervals, if not already done
  statsFcn <- function(x) quantile(x, probs=alpha, na.rm=TRUE);
  fit <- bootstrapTCNandDHByRegion(fit, statsFcn=statsFcn, ..., verbose=less(verbose, 2));

  segs <- as.data.frame(fit);

  # Extract confidence interval
  alphaTag <- sprintf("%g%%", 100*alpha);
  column <- sprintf("dh_%s", alphaTag);
  # Sanity checks
  stopifnot(is.element(column, colnames(segs)));

  # One-sided test
  verbose && enter(verbose, "Calling segments");
  value <- segs[,column, drop=TRUE];
  call <- (value < tau);
  verbose && printf(verbose, "Number of segments called allelic balance (AB): %d (%.2f%%) of %d\n", sum(call), 100*sum(call)/nrow(segs), nrow(segs));
  verbose && exit(verbose);

  segs <- cbind(segs, ab.call=call);
  fit$output <- segs;

  verbose && exit(verbose);

  fit;
}, protected=TRUE) # callAllelicBalanceByDH()


setMethodS3("callExtremeAllelicImbalanceByDH", "PairedPSCBS", function(fit, tau=0.60, alpha=0.05, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  # Argument 'tau':
  tau <- Arguments$getDouble(tau, range=c(0,Inf));

  # Argument 'alpha':
  alpha <- Arguments$getDouble(alpha, range=c(0,1));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling segments of extreme allelic imbalance (AI) from one-sided DH bootstrap confidence intervals");
  verbose && cat(verbose, "tau (offset adjusting for normal contamination and other biases): ", tau);
  verbose && cat(verbose, "alpha (CI quantile; significance level): ", alpha);


  # Calculate DH confidence intervalls, if not already done
  statsFcn <- function(x) quantile(x, probs=alpha, na.rm=TRUE);
  fit <- bootstrapTCNandDHByRegion(fit, statsFcn=statsFcn, ..., verbose=less(verbose, 2));

  segs <- as.data.frame(fit);

  # Extract confidence interval
  alphaTag <- sprintf("%g%%", 100*alpha);
  column <- sprintf("dh_%s", alphaTag);
  # Sanity checks
  stopifnot(is.element(column, colnames(segs)));

  # One-sided test
  value <- segs[,column, drop=TRUE];
  call <- (value >= tau);
  verbose && printf(verbose, "Number of segments called allelic balance (AB): %d (%.2f%%) of %d\n", sum(call), 100*sum(call)/nrow(segs), nrow(segs));
  verbose && exit(verbose);

  segs <- cbind(segs, ai.high.call=call);
  fit$output <- segs;

  fit;
}, protected=TRUE) # callExtremeAllelicImbalanceByDH()



setMethodS3("callABandHighAI", "PairedPSCBS", function(fit, tauAB=0.10, alphaAB=0.05, tauHighAI=0.60, alphaHighAI=0.05, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling segments to be in allelic balance (AB) or extreme allelic imbalance (AI)");

  # Calculate DH confidence intervals, if not already done
  probs <- sort(unique(c(alphaAB, alphaHighAI)));
  statsFcn <- function(x) quantile(x, probs=probs, na.rm=TRUE);
  fit <- bootstrapTCNandDHByRegion(fit, statsFcn=statsFcn, ..., verbose=less(verbose, 1));

  # Call allelic balance
  fit <- callAllelicBalanceByDH(fit, tau=tauAB, alpha=alphaAB, ..., verbose=less(verbose, 1));

  # Call high allelic imbalance
  fit <- callExtremeAllelicImbalanceByDH(fit, tau=tauHighAI, alpha=alphaHighAI, ..., verbose=less(verbose, 1));

  verbose && exit(verbose);

  fit;
}) # callABandHighAI()


##############################################################################
# HISTORY
# 2010-11-26 [HB]
# o BUG FIX: callABandHighAI() for PairedPSCBS used the old DH-only
#   bootstrap method.
# o BUG FIX: The call functions, for instance callABandHighAI(), would throw
#   'Error in quantile.default(x, probs = alpha) : missing values and NaN's
#   not allowed if 'na.rm' is FALSE' unless bootstrapTCNandDHByRegion() was
#   run before.
# 2010-11-22 [HB]
# o Added more verbose output to callABandHighAI().
# o Updated callAllelicBalanceByDH() and callExtremeAllelicImbalanceByDH()
#   to utilize bootstrapTCNandDHByRegion().
# 2010-10-25 [HB]
# o Relaced argument 'ciRange' with 'alpha' for callAllelicBalanceByDH() and
#   callExtremeAllelicImbalanceByDH().
# o Renamed callAllelicBalance() to callAllelicBalanceByDH() and
#   callExtremeAllelicImbalanceByDH() to callExtremeAllelicImbalance().
# o Added arguments 'alphaAB' and 'alphaHighAI' to callABandHighAI().
# o Added sanity checks to the call methods.
# o Now arguments '...' to callABandHighAI() are passed down.
# o Now also arguments '...' to callAllelicBalance() and
#   callExtremeAllelicImbalance() are passed to bootstrapDHByRegion().
# o Added argument 'ciRange' to callAllelicBalance() and
#   callExtremeAllelicImbalance().
# 2010-09-16 [HB]
# o Added callABandHighAI().
# o Added callAllelicBalance() and callExtremeAllelicImbalance().
# o Created.
##############################################################################
