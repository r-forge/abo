setMethodS3("bootstrapTCNandDHByRegion", "PairedPSCBS", function(fit, B=1000, statsFcn=function(x) quantile(x, probs=c(0.025, 0.050, 0.95, 0.975), na.rm=TRUE), by=c("betaTN", "betaT"), ..., verbose=FALSE) {
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


  verbose && enter(verbose, "Resample (TCN,DH) signals and re-estimate mean levels");
  data <- fit$data;
  segs <- fit$output;
  params <- fit$params;
  listOfDhLociNotPartOfSegment <- fit$listOfDhLociNotPartOfSegment;

  # Sanity check
  if (!params$joinSegments) {
    throw("Cannot bootstrap TCN and DH by regions unless PSCNs are segmented using joinSegments=TRUE.");
  }
  # Sanity check (same as above, but just in case)
  stopifnot(segs$tcn.loc.start == segs$dh.loc.start);
  stopifnot(segs$tcn.loc.end == segs$dh.loc.end);


  # Find estimates to be done
  stats <- statsFcn(1);
  stopifnot(!is.null(names(stats)));
  tcnStatsNames <- sprintf("tcn.%s", names(stats));
  dhStatsNames <- sprintf("dh.%s", names(stats));
  statsNames <- c(tcnStatsNames, dhStatsNames);
  isDone <- is.element(statsNames, names(segs));

  # Already done?
  if (all(isDone)) {
    verbose && cat(verbose, "Already done.");
    verbose && exit(verbose);
    return(fit);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Reorder signals?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##  index <- data$index;
##  if (!is.null(index)) {
##    verbose && enter(verbose, "Reordering data");
##    o <- order(index);
##    data <- data[o,,drop=FALSE];
##    rm(o);
##    verbose && exit(verbose);
##  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Drop missing values
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##  ### TO COMPLICATED FOR NOW. /HB 2010-11-21
##  # Drop missing values (Is that ok for bootstrapping? /HB 2010-09-16)
##  chromosome <- data$chromosome;
##  x <- data$x;
##  keep <- (is.finite(chromosome) & is.finite(x));
##  data <- data[keep,,drop=FALSE];
##
##  # Update loci to drop
##  listOfDhLociNotPartOfSegment <- fit$listOfDhLociNotPartOfSegment;
##  if (!is.null(listOfDhLociNotPartOfSegment)) {
##    for (ii in seq(along=listOfDhLociNotPartOfSegment)) {
##      listOfDhLociNotPartOfSegment
##    }
##    lociToExclude <- listOfDhLociNotPartOfSegment[[tcnId]][[dhId]];
##  }
##
##  # Not needed anymore
##  rm(keep, chromosome, x);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get (x,TCN,BAF) data
  chromosome <- data$chromosome;
  x <- data$x;
  CT <- data$CT;
  betaT <- data[[by]];
  muN <- data$muN;

  # Sanity checks
  stopifnot(all(is.finite(chromosome)));
  stopifnot(all(is.finite(x)));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Classify each locus as (i) heterozygous SNP, (ii) homozygous SNP,
  # or (iii) non-polymorphic loci
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Identifying heterozygous & homozygous SNPs and non-polymorphic loci");
  nbrOfLoci <- length(muN);
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);

  # SNPs are identifies as those loci that have non-missing 'muN' (& betaTN')
  isSNP <- (!is.na(muN) & !is.na(betaT));
  snps <- which(isSNP);
  nonSNPs <- which(!isSNP);
  nbrOfSNPs <- sum(isSNP);
  nbrOfNonSNPs <- sum(isSNP);
  verbose && cat(verbose, "Number of SNPs: ", nbrOfSNPs);
  verbose && cat(verbose, "Number of non-SNPs: ", nbrOfNonSNPs);
  stopifnot(length(intersect(snps, nonSNPs)) == 0);

  # Heterozygous SNPs
  isHet <- isSNP & (muN == 1/2);
  hets <- which(isHet);
  homs <- which(!isHet);
  nbrOfHets <- length(hets);
  nbrOfHoms <- length(homs);
  verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                      nbrOfHets, 100*nbrOfHets/nbrOfSNPs);
  verbose && printf(verbose, "Number of homozygous SNPs: %d (%.2f%%)\n",
                                      nbrOfHoms, 100*nbrOfHoms/nbrOfSNPs);
  stopifnot(length(intersect(hets, homs)) == 0);

  stopifnot(length(isSNP) == nbrOfLoci);
  stopifnot(length(isHet) == nbrOfLoci);

  # Not needed anymore
  rm(muN);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Precalculate DH signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate DHs for heterozygous SNPs
  rho <- 2*abs(betaT - 1/2);

  # DH is by definition only defined for heterozygous SNPs.  
  # For simplicity, we set it to be NA for non-heterozygous loci.
  rho[!isHet] <- NA;

  # Not needed anymore
  rm(betaT);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Resample (TCN,DH) within each segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nrow(segs);
  naValue <- as.double(NA);

  tcnMeanList <- vector("list", nbrOfSegments);
  dhMeanList <- vector("list", nbrOfSegments);
  for (jj in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("DH segment #%d of %d", jj, nbrOfSegments));
    segJJ <- segs[jj,,drop=FALSE];

    tcnId <- segJJ[,"tcn.id"];
    dhId <- segJJ[,"dh.id"];

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identify loci in segment
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    chr <- segJJ$chromosome[1];
    start <- segJJ$dh.loc.start[1];
    stop <- segJJ$dh.loc.end[1];

    nbrOfTCNs <- segJJ[,"tcn.num.mark"];
    nbrOfDHs <- segJJ[,"dh.num.mark"];
    if (is.na(nbrOfTCNs)) nbrOfTCNs <- 0L;
    if (is.na(nbrOfDHs)) nbrOfDHs <- 0L;

    unitsJJ <- whichVector(chr == chromosome & start <= x & x <= stop);
    lociToExclude <- listOfDhLociNotPartOfSegment[[tcnId]][[dhId]];
    verbose && cat(verbose, "Excluding loci that belongs to a flanking segment: ", length(lociToExclude));
    unitsJJ <- setdiff(unitsJJ, lociToExclude);
    nbrOfUnitsJJ <- length(unitsJJ);
    verbose && cat(verbose, "Number of loci in segment: ", nbrOfUnitsJJ);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identify (i) heterozygous SNPs, (ii) homozygous SNPs, and
    # non-polymorphic loci
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    hetsJJ <- intersect(hets, unitsJJ);
    nbrOfHetsJJ <- length(hetsJJ);
    verbose && cat(verbose, "Number of heterozygous SNPs in segment: ", nbrOfHetsJJ);

    homsJJ <- intersect(homs, unitsJJ);
    nbrOfHomsJJ <- length(homsJJ);
    verbose && cat(verbose, "Number of homozygous SNPs in segment: ", nbrOfHomsJJ);

    # Sanity check
    stopifnot(length(intersect(hetsJJ, homsJJ)) == 0);

    nonSNPsJJ <- intersect(nonSNPs, unitsJJ);
    nbrOfNonSNPsJJ <- length(nonSNPsJJ);
    snpsJJ <- sort(c(hetsJJ, homsJJ));
    stopifnot(length(intersect(snpsJJ, nonSNPsJJ)) == 0);

    # Sanity check
    stopifnot(nbrOfHetsJJ + nbrOfHomsJJ + nbrOfNonSNPsJJ == nbrOfUnitsJJ);
    unitsJJ2 <- sort(c(snpsJJ, nonSNPsJJ));
    stopifnot(all(unitsJJ2 == sort(unitsJJ)));

    # These numbers should be preserved when the resampling
    verbose && printf(verbose, "Number of (#hets, #homs, #nonSNPs): (%d,%d,%d)\n",
                      nbrOfHetsJJ, nbrOfHomsJJ, nbrOfNonSNPsJJ);
                                              

    # Sanity checks
    tol <- 0.0005;
    ys <- rho[hetsJJ];
    mu <- mean(ys, na.rm=TRUE);
    dMu <- (mu - segJJ$dh.mean);
    if (abs(dMu) > tol) {
      str(list(nbrOfUnits=nbrOfUnitsJJ, dh.num.mark=segJJ$dh.num.mark, mu=mu, dh.mean=segJJ$dh.mean, dMu=dMu, "abs(dMu)"=abs(dMu), "min(x[units])"=min(x[unitsJJ])));
      stop("INTERNAL ERROR: Incorrect DH mean!");
    }

    ys <- CT[unitsJJ];
    mu <- mean(ys, na.rm=TRUE);
    dMu <- (mu - segJJ$tcn.mean);
    if (abs(dMu) > tol) {
      str(list(nbrOfUnits=nbrOfUnitsJJ, dh.num.mark=segJJ$tcn.num.mark, mu=mu, tcn.mean=segJJ$tcn.mean, dMu=dMu, "abs(dMu)"=abs(dMu), "min(x[units])"=min(x[unitsJJ])));
      stop("INTERNAL ERROR: Incorrect TCN mean!");
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Bootstrap while preserving (#hets, #homs, #nonSNPs)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Bootstrap B times
    naValue <- as.double(NA);
    tcnMean <- dhMean <- rep(naValue, times=B);
    for (bb in seq(length=B)) {
      # (1) Bootstrap DHs
      if (nbrOfHetsJJ > 0) {
        # (a) Resample heterozygous SNPs (=> resampled DH units)
        hetsBB <- resample(hetsJJ, size=nbrOfHetsJJ, replace=TRUE);

        # Extract signals
        rhoBB <- rho[hetsBB];

        # Calculate bootstrap mean
        dhMean[bb] <- mean(rhoBB, na.rm=TRUE);
      } else {
        hetsBB <- NULL;
      } # if (nbrOfHets > 0)
  
      # (2) Bootstrap TCNs
      if (nbrOfUnitsJJ > 0) {
        # (a) Resample homozygous SNPs
        homsBB <- resample(homsJJ, size=nbrOfHomsJJ, replace=TRUE);
  
        # (b) Resample non-SNPs
        nonSNPsBB <- resample(nonSNPsJJ, size=nbrOfNonSNPsJJ, replace=TRUE);
      
        # (c) Resampled TCN units
        unitsBB <- c(hetsBB, homsBB, nonSNPsBB);

        # Extract signals
        tcnBB <- CT[unitsBB];

        # Calculate bootstrap mean
        tcnMean[bb] <- mean(tcnBB, na.rm=TRUE);
      } # if (nbrOfUnitsJJ > 0)
    } # (for bb ...)

    dhMeanList[[jj]] <- dhMean;
    tcnMeanList[[jj]] <- tcnMean;

    verbose && exit(verbose);
  } # for (jj ...)

  verbose && str(verbose, dhMeanList);
  verbose && str(verbose, tcnMeanList);

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Calculate bootstrap statistics
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Calculating (TCN,DH) bootstrap statistics");
  tcnStatsList <- vector("list", nbrOfSegments);
  dhStatsList <- vector("list", nbrOfSegments);
  for (jj in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("Segment #%d of %d", jj, nbrOfSegments));

    # (1) Calculate DH bootstrap statistics
    dhMean <- dhMeanList[[jj]];
    dhStats <- statsFcn(dhMean);
    dhStatsList[[jj]] <- dhStats;

    # (2) Calculate TCN bootstrap statistics
    tcnMean <- tcnMeanList[[jj]];
    tcnStats <- statsFcn(tcnMean);
    tcnStatsList[[jj]] <- tcnStats;

    verbose && exit(verbose);
  } # for (jj ...)
  verbose && exit(verbose);

  # Flatten statistics, if possible
  tcnStats <- sapply(tcnStatsList, FUN=function(x) x);
  if (!is.matrix(tcnStats)) {
    tcnStats <- as.matrix(tcnStats);
  } else {
    tcnStats <- t(tcnStats);
  }
  # Sanity check
  stopifnot(is.matrix(tcnStats));
  # Column matrix
  colnames(tcnStats) <- tcnStatsNames;

  dhStats <- sapply(dhStatsList, FUN=function(x) x);
  if (!is.matrix(dhStats)) {
    dhStats <- as.matrix(dhStats);
  } else {
    dhStats <- t(dhStats);
  }
  # Sanity check
  stopifnot(is.matrix(dhStats));
  # Column matrix
  colnames(dhStats) <- dhStatsNames;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Store
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Store results
  segs <- cbind(segs, tcnStats, dhStats);

  # Drop previously estimated values
  dups <- duplicated(colnames(segs), fromLast=TRUE);
  if (any(dups)) {
    stats <- stats[,!dups, drop=FALSE];
  }

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Statistical sanity checks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (B >= 100) {
    # Sanity check
    stopifnot(all(tcnStats[,1] <= segs$tcn.mean, na.rm=TRUE));
    stopifnot(all(segs$tcn.mean <= tcnStats[,ncol(tcnStats)], na.rm=TRUE));
  
    # Sanity check
    stopifnot(all(dhStats[,1] <= segs$dh.mean, na.rm=TRUE));
    stopifnot(all(segs$dh.mean <= dhStats[,ncol(dhStats)], na.rm=TRUE));
  }

  
  fitB <- fit;
  fitB$output <- segs;

  verbose && exit(verbose);

  fitB;
}) # bootstrapDHByRegion()



##############################################################################
# HISTORY
# 2010-11-21
# o Added bootstrapTCNandDHByRegion() for PairedPSCBS.
# o Created from PairedPSCBS.BOOT.R.
##############################################################################
