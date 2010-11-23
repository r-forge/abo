setMethodS3("bootstrapTCNandDHByRegion", "PairedPSCBS", function(fit, B=1000, statsFcn=function(x) quantile(x, probs=c(0.025, 0.050, 0.95, 0.975)), by=c("betaTN", "betaT"), ..., seed=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'B':
  B <- Arguments$getInteger(B, range=c(1,Inf));

  # Argument 'by':
  by <- match.arg(by);

  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getInteger(seed);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Resample (TCN,DH) signals and re-estimate mean levels");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Set the random seed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(seed)) {
    verbose && enter(verbose, "Setting (temporary) random seed");
    oldRandomSeed <- NULL;
    if (exists(".Random.seed", mode="integer")) {
      oldRandomSeed <- get(".Random.seed", mode="integer");
    }
    on.exit({
      if (!is.null(oldRandomSeed)) {
        .Random.seed <<- oldRandomSeed;
      }
    }, add=TRUE);
    verbose && cat(verbose, "The random seed will be reset to its original state afterward.");
    verbose && cat(verbose, "Seed: ", seed);
    set.seed(seed);
    verbose && exit(verbose);
  }


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
  nbrOfStats <- length(stats);
  statsNames <- names(stats);

  # Already done?
  tcnStatsNames <- sprintf("tcn_%s", names(stats));
  dhStatsNames <- sprintf("dh_%s", names(stats));
  c1StatsNames <- sprintf("c1_%s", names(stats));
  c2StatsNames <- sprintf("c2_%s", names(stats));
  allStatsNames <- c(tcnStatsNames, dhStatsNames, c1StatsNames, c2StatsNames);
  isDone <- is.element(allStatsNames, names(segs));
  names(isDone) <- allStatsNames;
  verbose && cat(verbose, "Already done?");
  verbose && print(verbose, isDone);

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

  # Allocate JxBx4 matrix M of bootstrap means
  naValue <- as.double(NA);
  dim <- c(nbrOfSegments, B, 4);
  dimnames <- list(NULL, NULL, c("tcn", "dh", "c1", "c2"));
  M <- array(naValue, dim=dim, dimnames=dimnames);
  verbose && str(verbose, M);

  for (jj in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("Segment #%d of %d", jj, nbrOfSegments));
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
    for (bb in seq(length=B)) {
      # (1) Bootstrap DHs
      if (nbrOfHetsJJ > 0) {
        # (a) Resample heterozygous SNPs (=> resampled DH units)
        hetsBB <- resample(hetsJJ, size=nbrOfHetsJJ, replace=TRUE);

        # Extract signals
        rhoBB <- rho[hetsBB];

        # Calculate bootstrap mean
        M[jj,bb,"dh"] <- mean(rhoBB, na.rm=TRUE);
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
        M[jj,bb,"tcn"] <- mean(tcnBB, na.rm=TRUE);
      } # if (nbrOfUnitsJJ > 0)
    } # (for bb ...)

    verbose && exit(verbose);
  } # for (jj ...)

  verbose && cat(verbose, "Bootstrap means");
  verbose && str(verbose, M);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Calculate (C1,C2) bootstrap statistics
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating (C1,C2) from (TCN,DH) bootstraps");
  C1 <- (1-M[,,"dh"]) * M[,,"tcn"] / 2;
  C2 <- M[,,"tcn"] - C1;
  M[,,"c1"] <- C1;
  M[,,"c2"] <- C2;
  verbose && str(verbose, M);
  verbose && exit(verbose);

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Calculate bootstrap statistics
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Calculating (TCN,DH) bootstrap statistics");
  # Allocate JxQx4 matrix S
  naValue <- as.double(NA);
  dim <- dim(M);
  dimnames <- dimnames(M);
  dim[2] <- nbrOfStats;
  dimnames[[2]] <- statsNames;
  S <- array(naValue, dim=dim, dimnames=dimnames);
  verbose && str(verbose, S);

  fields <- dimnames(M)[[3]];
  for (kk in seq(along=fields)) {
    field <- fields[kk];
    verbose && enter(verbose, sprintf("Field #%d ('%s') of %d", kk, field, length(fields)));

    Mkk <- M[,,kk,drop=TRUE];  # An JxB matrix
    # Sanity check
    stopifnot(nrow(Mkk) == nbrOfSegments);
    stopifnot(ncol(Mkk) == B);

    for (jj in seq(length=nbrOfSegments)) {
      verbose && enter(verbose, sprintf("Segment #%d of %d", jj, nbrOfSegments));

      Mkkjj <- Mkk[jj,,drop=TRUE]; # A vector of length B
      S[jj,,kk] <- statsFcn(Mkkjj);

      verbose && exit(verbose);
    } # for (jj ...)

    verbose && exit(verbose);
  } # for (jj ...)
  verbose && exit(verbose);
  verbose && cat(verbose, "Bootstrap statistics");
  verbose && str(verbose, S);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Store
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Reshape JxQx4 array to Jx(4*Q) matrix
  T <- wrap(S, map=list(1,NA), sep="_");
  colnames(T) <- gsub("(.*)_(.*)", "\\2_\\1", colnames(T));

  # Append
  segs <- cbind(segs, T);

  # Drop previously estimated values
  dups <- duplicated(colnames(segs), fromLast=TRUE);
  if (any(dups)) {
    stats <- stats[,!dups, drop=FALSE];
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Statistical sanity checks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (B >= 100) {
    verbose && enter(verbose, "Statistical sanity checks (iff B >= 100)");
    tryCatch({
      fields <- dimnames(M)[[3]];
      for (kk in seq(along=fields)) {
        field <- fields[kk];
        verbose && enter(verbose, sprintf("Field #%d ('%s') of %d", kk, field, length(fields)));
  
        # Bootstrap statistics
        Skk <- S[,,kk, drop=TRUE];
        range <- Skk[,c(1,ncol(Skk))];
        verbose && str(verbose, range);
  
        # Segmentation means
        key <- sprintf("%s.mean", field);
        segMean <- segs[[key]];
        verbose && str(verbose, segMean);
  
        # Sanity check
        stopifnot(all(range[,1] <= segMean, na.rm=TRUE));
        stopifnot(all(segMean <= range[,2], na.rm=TRUE));
  
        verbose && exit(verbose);
      } # for (kk ...)
    }, error = function(ex) {
      # If an error, display the data, then throw the exception
      verbose && print(verbose, segs);
      throw(ex);
    })
    verbose && exit(verbose);
  } # if (B >= 100)
  
  
  fitB <- fit;
  fitB$output <- segs;

  verbose && exit(verbose);

  fitB;
}) # bootstrapTCNandDHByRegion()



##############################################################################
# HISTORY
# 2010-11-22
# o BUG FIX: bootstrapTCNandDHByRegion() for PairedPSCBS would not correctly
#   detect if bootstrap results are already available.
# 2010-11-21
# o Added argument 'seed'.
# o Added bootstrapTCNandDHByRegion() for PairedPSCBS.
# o Created from PairedPSCBS.BOOT.R.
##############################################################################
