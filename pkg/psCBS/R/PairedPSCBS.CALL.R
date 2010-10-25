setMethodS3("callAllelicBalance", "PairedPSCBS", function(fit, tau=0.10, ciRange=c(0.05, 0.95), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  # Argument 'tau':
  tau <- Arguments$getDouble(tau, range=c(0,Inf));

  # Argument 'ciRange':
  ciRange <- Arguments$getDoubles(ciRange, range=c(0,1), length=c(2,2));


  # Calculate DH confidence intervalls, if not already done
  statsFcn <- function(x) quantile(x, probs=ciRange);
  fit <- bootstrapDHByRegion(fit, statsFcn=statsFcn, ...);

  segs <- as.data.frame(fit);

  # Extract confidence interval
  ciRangeTags <- sprintf("%g%%", 100*ciRange);
  columns <- sprintf("dh.ci.%s", ciRangeTags);
  stopifnot(all(is.element(columns, colnames(segs))));

  ci <- segs[,columns, drop=FALSE];

  call <- (ci[,1] < tau);
  segs <- cbind(segs, ab.call=call);
  fit$output <- segs;

  fit;
})


setMethodS3("callExtremeAllelicImbalance", "PairedPSCBS", function(fit, tau=0.60, ciRange=c(0.05, 0.95), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  # Argument 'tau':
  tau <- Arguments$getDouble(tau, range=c(0,Inf));

  # Argument 'ciRange':
  ciRange <- Arguments$getDoubles(ciRange, range=c(0,1), length=c(2,2));


  # Calculate DH confidence intervalls, if not already done
  statsFcn <- function(x) quantile(x, probs=ciRange);
  fit <- bootstrapDHByRegion(fit, statsFcn=statsFcn, ...);

  segs <- as.data.frame(fit);

  # Extract confidence interval
  ciRangeTags <- sprintf("%g%%", 100*ciRange);
  columns <- sprintf("dh.ci.%s", ciRangeTags);
  stopifnot(all(is.element(columns, colnames(segs))));

  ci <- segs[,columns, drop=FALSE];

  call <- (ci[,1] >= tau);
  segs <- cbind(segs, ai.high.call=call);
  fit$output <- segs;

  fit;
})


setMethodS3("callABandHighAI", "PairedPSCBS", function(fit, tauAB=0.10, tauHighAI=0.60, ...) {
  fit <- callAllelicBalance(fit, tau=tauAB);
  fit <- callExtremeAllelicImbalance(fit, tau=tauHighAI);
  fit;
})


##############################################################################
# HISTORY
# 2010-10-25 [HB]
# o Now also arguments '...' to callAllelicBalance() and
#   callExtremeAllelicImbalance() are passed to bootstrapDHByRegion().
# o Added argument 'ciRange' to callAllelicBalance() and
#   callExtremeAllelicImbalance().
# 2010-09-16 [HB]
# o Added callABandHighAI().
# o Added callAllelicBalance() and callExtremeAllelicImbalance().
# o Created.
##############################################################################
