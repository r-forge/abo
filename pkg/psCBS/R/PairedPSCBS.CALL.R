setMethodS3("callAllelicBalance", "PairedPSCBS", function(fit, tau=0.10, ...) {
  # Argument 'tau':
  tau <- Arguments$getDouble(tau, range=c(0,Inf));

  # Calculate DH confidence intervalls, if not already done
  fit <- bootstrapDHByRegion(fit);

  segs <- as.data.frame(fit);
  ci <- segs[,c("dh.ci.5%", "dh.ci.95%"), drop=FALSE];
  call <- (ci[,1] < tau);
  segs <- cbind(segs, ab.call=call);
  fit$output <- segs;

  fit;
})


setMethodS3("callExtremeAllelicImbalance", "PairedPSCBS", function(fit, tau=0.60, ...) {
  # Argument 'tau':
  tau <- Arguments$getDouble(tau, range=c(0,Inf));

  # Calculate DH confidence intervalls, if not already done
  fit <- bootstrapDHByRegion(fit);

  segs <- as.data.frame(fit);
  ci <- segs[,c("dh.ci.5%", "dh.ci.95%"), drop=FALSE];
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
# 2010-09-16 [HB]
# o Added callABandHighAI().
# o Added callAllelicBalance() and callExtremeAllelicImbalance().
# o Created.
##############################################################################
