setMethodS3("callAllelicBalanceByDH", "PairedPSCBS", function(fit, tau=0.10, alpha=0.05, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  # Argument 'tau':
  tau <- Arguments$getDouble(tau, range=c(0,Inf));

  # Argument 'alpha':
  alpha <- Arguments$getDouble(alpha, range=c(0,1));


  # Calculate DH confidence intervals, if not already done
  statsFcn <- function(x) quantile(x, probs=alpha);
  fit <- bootstrapDHByRegion(fit, statsFcn=statsFcn, ...);

  segs <- as.data.frame(fit);

  # Extract confidence interval
  alphaTag <- sprintf("%g%%", 100*alpha);
  column <- sprintf("dh.%s", alphaTag);
  # Sanity checks
  stopifnot(is.element(column, colnames(segs)));

  # One-sided test
  value <- segs[,column, drop=TRUE];
  call <- (value < tau);
  segs <- cbind(segs, ab.call=call);
  fit$output <- segs;

  fit;
}) # callAllelicBalanceByDH()


setMethodS3("callExtremeAllelicImbalanceByDH", "PairedPSCBS", function(fit, tau=0.60, alpha=0.05, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  # Argument 'tau':
  tau <- Arguments$getDouble(tau, range=c(0,Inf));

  # Argument 'alpha':
  alpha <- Arguments$getDouble(alpha, range=c(0,1));


  # Calculate DH confidence intervalls, if not already done
  statsFcn <- function(x) quantile(x, probs=alpha);
  fit <- bootstrapDHByRegion(fit, statsFcn=statsFcn, ...);

  segs <- as.data.frame(fit);

  # Extract confidence interval
  alphaTag <- sprintf("%g%%", 100*alpha);
  column <- sprintf("dh.%s", alphaTag);
  # Sanity checks
  stopifnot(is.element(column, colnames(segs)));

  # One-sided test
  value <- segs[,column, drop=TRUE];
  call <- (value >= tau);
  segs <- cbind(segs, ai.high.call=call);
  fit$output <- segs;

  fit;
}) # callExtremeAllelicImbalanceByDH()


setMethodS3("callABandHighAI", "PairedPSCBS", function(fit, tauAB=0.10, alphaAB=0.05, tauHighAI=0.60, alphaHighAI=0.05, ...) {
  # Calculate DH confidence intervals, if not already done
  probs <- sort(unique(c(alphaAB, alphaHighAI)));
  statsFcn <- function(x) quantile(x, probs=probs);
  fit <- bootstrapDHByRegion(fit, statsFcn=statsFcn, ...);

  # Call allelic balance
  fit <- callAllelicBalanceByDH(fit, tau=tauAB, alpha=alphaAB, ...);

  # Call high allelic imbalance
  fit <- callExtremeAllelicImbalanceByDH(fit, tau=tauHighAI, alpha=alphaHighAI, ...);

  fit;
})


##############################################################################
# HISTORY
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
