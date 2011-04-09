###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod callAllelicBalance
# @aliasmethod callAB
#
# @title "Calls segments that are in allelic balance"
#
# \description{
#  @get "title", i.e. that have equal minor and major copy numbers.
# }
#
# @synopsis
#
# \arguments{
#   \item{flavor}{A @character string specifying which type of
#    call to use.}
#   \item{...}{Additional arguments passed to the caller.}
# }
#
# \value{
#   Returns a @see "PairedPSCBS" object with allelic-balance calls.
# }
#
# @author
#
# \seealso{
#   Internally, one of the following methods are used:
#   @seemethod "callAllelicBalanceByDH".
# }
#
#*/###########################################################################
setMethodS3("callAllelicBalance", "PairedPSCBS", function(fit, flavor=c("DeltaAB*"), ...) {
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  if (flavor == "DeltaAB*") {
    callAllelicBalanceByDH(fit, ...);
  } else {
    throw("Cannot call allelic balance. Unsupported flavor: ", flavor);
  }
})




###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod callAllelicBalanceByDH
#
# @title "Calls segments that are in allelic balance"
#
# \description{
#  @get "title" by thresholding on DH using a predetermined threshold.
#  The variability of the DH mean levels is taken into account via a
#  bootstrap estimator.
# }
#
# @synopsis
#
# \arguments{
#   \item{flavor}{A @character string specifying which type of
#    call to use.}
#   \item{tau}{(Tuning parameter) A non-negative @numeric threshold.}
#   \item{alpha}{A @numeric in [0,1] specifying the upper and lower
#     quantiles calculated by the bootstrap estimator.}
#   \item{...}{Additional arguments passed to the caller.}
# }
#
# \value{
#   Returns a @see "PairedPSCBS" object with allelic-balance calls.
# }
#
# @author
#
# \section{Algorithm}{
#  \itemize{
#    \item Foo
#    \item Bar
#  }
# }
#
# \seealso{
#   Instead of calling this method explicitly, it is recommended
#   to use the @seemethod "callAllelicBalance" method.
# }
#
# @keyword internal
#*/###########################################################################
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
  probs <- c(alpha, 1-alpha);
  statsFcn <- function(x) quantile(x, probs=probs, na.rm=TRUE);
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
  nbrOfCalls <- sum(call, na.rm=TRUE);
  verbose && printf(verbose, "Number of segments called allelic balance (AB): %d (%.2f%%) of %d\n", nbrOfCalls, 100*nbrOfCalls/nrow(segs), nrow(segs));
  verbose && exit(verbose);

  segs <- cbind(segs, ab.call=call);
  fit$output <- segs;

  # Append 'tau' and 'alpha' to parameters
  params <- fit$params;
  params$tauAB <- tau;
  params$alphaAB <- alpha;
  fit$params <- params;

  verbose && exit(verbose);

  fit;
}, protected=TRUE) # callAllelicBalanceByDH()






##############################################################################
# HISTORY
# 2011-04-08
# o Added Rdoc for callAllelicBalance() and callAllelicBalanceByDH().
# o Extracted from PairedPSCBS.CALL.R.
# 2011-02-03
# o Updated default for 'tauAB' of callABandHighAI() and callABandLowC1()
#   to be estimated from data using estimateTauAB().
# 2010-12-07
# o Added callLowC1ByC1() and callABandLowC1().
# 2010-11-27
# o Corrected verbose output to call results.
# 2010-11-26 [HB]
# o Now all call functions estimate symmetric bootstrap quantiles for
#   convenince of plotting confidence intervals.
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
