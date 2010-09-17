setMethodS3("bootstrapDHByRegion", "PairedPSCBS", function(fit, B=100, statsFcn=function(x) quantile(x, probs=c(0.05, 0.95)), by=c("betaTN", "betaT"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'B':
  B <- Arguments$getInteger(B, range=c(1,Inf));

  # Argument 'by':
  by <- match.arg(by);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Resample DH signals and reestimate DH mean levels");

  data <- fit$data;
  segs <- fit$output;

  # Already done?
  stats <- statsFcn(1);
  stopifnot(!is.null(names(stats)));
  statsNames <- sprintf("dh.ci.%s", names(stats));
  if (all(is.element(statsNames, names(segs)))) {
    verbose && cat(verbose, "Already done.");
    verbose && exit(verbose);
    return(fit);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Precalculate signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (x,BAF) data
  x <- data$x;
  beta <- data[[by]];

  # Extract subset for heterozygous SNPs
  muN <- data$muN;
  isHet <- (muN == 1/2);
  x <- x[isHet];
  beta <- beta[isHet];

  # Calculate DHs for heterozygous SNPs
  rho <- 2*abs(beta - 1/2);

  # Drop missing values  (Is that ok for bootstrapping? /HB 2010-09-16)
  keep <- is.finite(x) & is.finite(rho);
  x <- x[keep];
  rho <- rho[keep];

  # Not needed anymore
  rm(beta, muN, isHet, keep);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Resample DH within segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nrow(segs);

  rhoMeanList <- vector("list", nbrOfSegments);
  for (jj in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("DH segment #%d of %d", jj, nbrOfSegments));
    segsJJ <- segs[jj,,drop=FALSE];

    # Identify loci in segment
    start <- segsJJ$dh.loc.start[1];
    stop <- segsJJ$dh.loc.end[1];
    units <- whichVector(start <= x & x <= stop);

    if (length(units) >= 1) {
      # Bootstrap B times
      naValue <- as.double(NA);
      rhoMean <- rep(naValue, times=B);
      for (bb in seq(length=B)) {
        # Resample indices
        unitsS <- sample(units, replace=TRUE);
  
        # Resample data
        rhoB <- rho[unitsS];

        # Calculate new mean level (no NAs here)
        rhoMean[bb] <- mean(rhoB, na.rm=FALSE);
      } # for (bb ...)
    } else {
      rhoMean <- double(0);
    }
    rhoMeanList[[jj]] <- rhoMean;

    verbose && exit(verbose);
  } # for (jj ...)

  # Not needed anymore
  rm(x, rho);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Calculate bootstrap statistics
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Calculating bootstrap statistics");
  statsList <- vector("list", length(rhoMeanList));
  for (jj in seq(along=rhoMeanList)) {
    verbose && enter(verbose, sprintf("DH segment #%d of %d", 
                                             jj, length(rhoMeanList)));
    # Get bootstrap sample
    rhoMean <- rhoMeanList[[jj]];

    # Calculate statistics
    stats <- statsFcn(rhoMean);
    # Store
    statsList[[jj]] <- stats;

    verbose && exit(verbose);
  } # for (jj ...)
  verbose && exit(verbose);

  # Flatten statistics, if possible
  stats <- sapply(statsList, FUN=function(x) x);

  # Sanity check
  stopifnot(is.matrix(stats));

  # Column matrix
  stats <- t(stats);
  colnames(stats) <- statsNames;

  # Store results
  segs <- cbind(segs, stats);
  
  fitB <- fit;
  fitB$output <- segs;

  verbose && exit(verbose);

  fitB;
}) # bootstrapDHByRegion()



##############################################################################
# HISTORY
# 2010-09-16 [HB]
# o Added bootstrapDHByRegion(), which is what is used by paired PSCBS.
# o Created.
##############################################################################
