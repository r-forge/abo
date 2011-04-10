###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod callLOH
#
# @title "Calls segments that are in LOH"
#
# \description{
#  @get "title", i.e. that have "zero" minor copy number.
# }
#
# @synopsis
#
# \arguments{
#   \item{flavor}{A @character string specifying which type of
#    call to use.}
#   \item{...}{Additional arguments passed to the caller.}
#   \item{force}{If @FALSE, and allelic-balance calls already exits,
#    then nothing is done, otherwise the calls are done.}
# }
#
# \value{
#   Returns a @see "PairedPSCBS" object with LOH calls.
# }
#
# @author
#
# \seealso{
#   Internally, one of the following methods are used:
#   @seemethod "callLowC1ByC1",
#   @seemethod "callExtremeAllelicImbalanceByDH".
# }
#
#*/###########################################################################
setMethodS3("callLOH", "PairedPSCBS", function(fit, flavor=c("SmallC1", "LargeDH"), ..., force=FALSE) {
  # Argument 'flavor':
  flavor <- match.arg(flavor);


  # Already done?
  segs <- as.data.frame(fit);
  calls <- segs$ab.call;
  if (!force && !is.null(calls)) {
    return(invisible(fit));
  }
  

  if (flavor == "SmallC1") {
    fit <- callLowC1ByC1(fit, ...);
  } else if (flavor == "LargeDH") {
    fit <- callExtremeAllelicImbalanceByDH(fit, ...);
  } else {
    throw("Cannot call LOH. Unsupported flavor: ", flavor);
  }

  return(invisible(fit));
})



setMethodS3("callLowC1ByC1", "PairedPSCBS", function(fit, tau=estimateTauLOH(fit, flavor="minC1|nonAB"), alpha=0.05, ..., verbose=FALSE) {
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

  verbose && cat(verbose, "tau (offset adjusting for bias in C1): ", tau);
  verbose && cat(verbose, "alpha (CI quantile; significance level): ", alpha);

  # Calculate C1 confidence intervals, if not already done
  probs <- c(alpha, 1-alpha);
  statsFcn <- function(x) quantile(x, probs=probs, na.rm=TRUE);
  fit <- bootstrapTCNandDHByRegion(fit, statsFcn=statsFcn, ..., verbose=less(verbose, 2));

  segs <- as.data.frame(fit);

  # Extract confidence interval
  alphaTag <- sprintf("%g%%", 100*alpha);
  column <- sprintf("c1_%s", alphaTag);
  # Sanity checks
  stopifnot(is.element(column, colnames(segs)));

  # One-sided test
  verbose && enter(verbose, "Calling segments");
  value <- segs[,column, drop=TRUE];
  call <- (value < tau);
  nbrOfCalls <- sum(call, na.rm=TRUE);
  verbose && printf(verbose, "Number of segments called low C1 (LowC1, \"LOH_C1\"): %d (%.2f%%) of %d\n", nbrOfCalls, 100*nbrOfCalls/nrow(segs), nrow(segs));
  verbose && exit(verbose);

  segs <- cbind(segs, lowc1.call=call);
  fit$output <- segs;

  # Append 'tau' and 'alpha' to parameters
  params <- fit$params;
  params$tauLowC1 <- tau;
  params$alphaLowC1 <- alpha;
  fit$params <- params;

  verbose && exit(verbose);

  fit;
}, private=TRUE) # callLowC1ByC1()




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
  call <- (value >= tau);
  nbrOfCalls <- sum(call, na.rm=TRUE);
  verbose && printf(verbose, "Number of segments called high allelic imbalance (AI/\"LOH_AI\"): %d (%.2f%%) of %d\n", nbrOfCalls, 100*nbrOfCalls/nrow(segs), nrow(segs));
  verbose && exit(verbose);

  segs <- cbind(segs, ai.high.call=call);
  fit$output <- segs;

  # Append 'tau' and 'alpha' to parameters
  params <- fit$params;
  params$tauExtremeDH <- tau;
  params$alphaExtremeDH <- alpha;
  fit$params <- params;

  verbose && exit(verbose);

  fit;
}, private=TRUE) # callExtremeAllelicImbalanceByDH() 



##############################################################################
# HISTORY
# 2011-04-10
# o Added callLOH().
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
