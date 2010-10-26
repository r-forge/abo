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
  chromosome <- data$chromosome;
  x <- data$x;
  muN <- data$muN;
  beta <- data[[by]];

  # SNPs are identifies as those loci that have non-missing 'betaTN' & 'muN'
  isSnp <- (!is.na(beta) & !is.na(muN));
  nbrOfSnps <- sum(isSnp);
  verbose && cat(verbose, "Number of SNPs: ", nbrOfSnps);

  # DH is by definition only defined for heterozygous SNPs.  
  # For simplicity, we set it to be NA for non-heterozygous loci.
  isHet <- isSnp & (muN == 1/2);
  verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                     sum(isHet), 100*sum(isHet)/nbrOfSnps);

  # Extract subset for heterozygous SNPs
  chromosome <- chromosome[isHet];
  x <- x[isHet];
  beta <- beta[isHet];

  # Calculate DHs for heterozygous SNPs
  rho <- 2*abs(beta - 1/2);

  # Drop missing values  (Is that ok for bootstrapping? /HB 2010-09-16)
  keep <- (is.finite(chromosome) & is.finite(x));
  chromosome <- chromosome[keep];
  x <- x[keep];
  rho <- rho[keep];

  # Not needed anymore
  rm(beta, muN, isHet, keep);

  # Sanity checks
  stopifnot(all(is.finite(chromosome)));
  stopifnot(all(is.finite(x)));
  stopifnot(all(is.finite(rho)));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Resample DH within segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nrow(segs);
  naValue <- as.double(NA);

  rhoMeanList <- vector("list", nbrOfSegments);
  for (jj in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("DH segment #%d of %d", jj, nbrOfSegments));
    segsJJ <- segs[jj,,drop=FALSE];

    # Identify loci in segment
    chr <- segsJJ$chromosome[1];
    start <- segsJJ$dh.loc.start[1];
    stop <- segsJJ$dh.loc.end[1];
    units <- whichVector(chr == chromosome & start <= x & x <= stop);
    nbrOfUnits <- length(units);

    if (nbrOfUnits >= 1) {
      # Sanity checks
      dups <- whichVector(duplicated(x[units]));
      nbrOfDups <- length(dups);

      # Sanity check
      mu <- mean(rho[units], na.rm=FALSE);

      # Special case /HB 2010-10-25
      if (nbrOfUnits != segsJJ$dh.num.mark) {
#        print(c(nbrOfUnits=nbrOfUnits, nbrOfDups=nbrOfDups, dh.num.mark=segsJJ$dh.num.mark));
#        stopifnot(nbrOfUnits-nbrOfDups == segsJJ$dh.num.mark);
        if (nbrOfDups == 1) {
          xT <- x[units];
          xDups <- xT[dups];
          if (xT[nbrOfUnits] == xDups) {
            toDrop <- which(xT == xDups)[1];
            mu <- mean(rho[units][-toDrop], na.rm=FALSE);
          }
        }
      }

      dMu <- (mu - segsJJ$dh.mean);
      tol <- 0.0005;
      if (abs(dMu) > tol) {
        str(list(nbrOfUnits=nbrOfUnits, nbrOfDups=nbrOfDups, dh.num.mark=segsJJ$dh.num.mark, mu=mu, dh.mean=segsJJ$dh.mean, dMu=dMu, "abs(dMu)"=abs(dMu), "min(x[units])"=min(x[units]), "x[units][dups]"=x[units][dups], "min(x[units]) == x[units][dups]"=(min(x[units]) == x[units][dups])));
        # Try to find a mean estimate that is correct by dropping 
        # on of the values. /HB 2010-10-25
        for (cc in 1:nbrOfUnits) {
          muT <- mean(rho[units][-cc], na.rm=FALSE);
          dMuT <- (muT - segsJJ$dh.mean);
          if (abs(dMuT) <= tol) {
            print(list(cc=cc, "x[cc]"=x[units][cc], "rho[cc]"=rho[units][cc], muT=muT, dh.mean=segsJJ$dh.mean, "x[units]"=x[units], "rho[units]"=rho[units]));
          }
        }
        # Discrepancies are only observed when there exist exactly
        # one duplicate! /HB 2010-10-25
#        stopifnot(nbrOfDups == 1);
#        stop("INTERNAL ERROR");
      }

      # Bootstrap B times
      rhoMean <- rep(naValue, times=B);
      for (bb in seq(length=B)) {
        # Resample indices
        unitsS <- sample(units, size=nbrOfUnits, replace=TRUE);
  
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
# 2010-10-25 [HB]
# o BUG FIX: bootstrapDHByRegion() for PairedPSCBS would bootstrap from the
#   incorrect set of loci when the DH region contained only one locus.
# o BUG FIX: bootstrapDHByRegion() for PairedPSCBS would bootstrap from the
#   incorrect set of loci if more than one chromosome was available.
# 2010-09-16 [HB]
# o Added bootstrapDHByRegion(), which is what is used by paired PSCBS.
# o Created.
##############################################################################
