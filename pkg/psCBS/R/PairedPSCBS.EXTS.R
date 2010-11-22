setMethodS3("as.data.frame", "PairedPSCBS", function(x, ...) {
  # To please R CMD check
  fit <- x;

  fit$output;
})


setMethodS3("print", "PairedPSCBS", function(x, ...) {
  # To please R CMD check
  fit <- x;

  segs <- as.data.frame(fit, ...);
  print(segs);
})


setMethodS3("subsetByDhSegments", "PairedPSCBS", function(fit, idxs, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  segs <- fit$output;
  stopifnot(!is.null(segs)); 
  nbrOfSegments <- nrow(segs);

  # Argument 'idxs':
  idxs <- Arguments$getIndices(idxs, max=nbrOfSegments);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 



  verbose && enter(verbose, "Extracting a subset of the DH segments");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- fit$data;
  stopifnot(!is.null(data));

  listOfDhLociNotPartOfSegment <- fit$listOfDhLociNotPartOfSegment;

  chromosomes <- getChromosomes(fit);
  nbrOfChromosomes <- length(chromosomes);
  verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes);
  verbose && print(verbose, chromosomes); 

  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  chromosome <- data$chromosome;
  x <- data$x;
  nbrOfLoci <- length(x);

  verbose && enter(verbose, "Subsetting segments");
  # Subset the region-level data
  segs <- segs[idxs,,drop=FALSE];
  verbose && print(verbose, segs);
  verbose && cat(verbose, "Number of TCN markers: ", sum(segs[["tcn.num.mark"]], na.rm=TRUE));
  verbose && exit(verbose);

  verbose && enter(verbose, "Subsetting data");
  # Subset the locus-level data
  keep <- logical(nbrOfLoci);
  for (kk in seq(length=nbrOfSegments)) {
    segKK <- segs[kk,];
    chrKK <- as.numeric(segKK[,"chromosome"]);
    xRange <- as.numeric(segKK[,c("dh.loc.start", "dh.loc.end")]);
    # Skip NA divider?
    if (all(is.na(xRange))) {
      next;
    }
    # Identify all SNPs in the region
    keepKK <- (xRange[1] <= x & x <= xRange[2]);
    if (!is.na(chrKK)) keepKK <- keepKK & (chromosome == chrKK);
    nKK <- sum(keepKK, na.rm=TRUE);

    # Special case?
    if (nKK > segKK[,"dh.num.mark"]) {
      verbose && cat(verbose, "Number of loci in DH segment: ", nKK);
      verbose && cat(verbose, "Number of loci in DH segment: ", segKK[,"dh.num.mark"]);

      # Sanity check
      stopifnot(!is.null(listOfDhLociNotPartOfSegment));

      tcnId <- segKK[,"tcn.id"];
      dhId <- segKK[,"dh.id"];
      dhLociNotPartOfSegment <- listOfDhLociNotPartOfSegment[[tcnId]];
      # Sanity check
      stopifnot(!is.null(dhLociNotPartOfSegment));

      lociToExclude <- dhLociNotPartOfSegment[[dhId]];
      verbose && cat(verbose, "Excluding loci that belongs to a flanking segment: ", length(lociToExclude));
      keepKK[lociToExclude] <- FALSE;
      nKK <- sum(keepKK, na.rm=TRUE);
    }

    verbose && cat(verbose, "Number of units: ", nKK);
    verbose && cat(verbose, "Number of TCN markers: ", segKK[,"tcn.num.mark"]);

    # Sanity check
    stopifnot(nKK == segKK[,"dh.num.mark"]);

    keep <- keep | keepKK;
  } # for (kk ...)
  keep <- whichVector(keep);
  verbose && cat(verbose, "Units to keep:");
  verbose && str(verbose, keep);

  data <- data[keep,,drop=FALSE];
  rm(keep); # Not needed anymore
  verbose && cat(verbose, "Number of units: ", nrow(data));
  verbose && exit(verbose);

  # Sanity check
#  stopifnot(nrow(data) == sum(segs[["tcn.num.mark"]], na.rm=TRUE));

  fitS <- fit;
  fitS$data <- data;
  fitS$output <- segs;

  verbose && exit(verbose);

  fitS;
})




setMethodS3("bootstrapCIs", "PairedPSCBS", function(fit, ...) {
  # ...
})


setMethodS3("callSegments", "PairedPSCBS", function(fit, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # 2. Calling regions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Call region
  # a. Classifying non-LOH regions as balanced or not balanced
  # b. Testing for LOH in the Tumor
}) # callSegments()



setMethodS3("extractTCNAndDHs", "PairedPSCBS", function(fit, ...) {
  segs <- fit$output;
  stopifnot(!is.null(segs));

  data <- segs[,c("tcn.mean", "dh.mean", "tcn.num.mark", "dh.num.mark"), drop=FALSE];
  data;
})


setMethodS3("extractMinorMajorCNs", "PairedPSCBS", function(fit, ...) {
  data <- extractTCNAndDHs(fit, ...);

  gamma <- data[,1];
  rho <- data[,2];
  C1 <- 1/2*(1-rho)*gamma;
  C2 <- gamma - C1;

  data[,1] <- C1;
  data[,2] <- C2;
  colnames(data)[1:2] <- c("C1", "C2");

  data;
})



setMethodS3("extractC1C2", "PairedPSCBS", function(...) {
  extractMinorMajorCNs(...);
})

setMethodS3("extractDeltaC1C2", "PairedPSCBS", function(...) {
  xy <- extractC1C2(...);
  X <- xy[,1:2,drop=FALSE];
  dX <- matrixStats::colDiffs(X);
  dX;
})



setMethodS3("postsegmentTCN", "PairedPSCBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
 
  verbose && enter(verbose, "Post-segmenting TCNs");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- fit$data;
  stopifnot(!is.null(data));

  segs <- fit$output;
  stopifnot(!is.null(segs));

  chromosomes <- getChromosomes(fit);
  nbrOfChromosomes <- length(chromosomes);
  verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes);
  verbose && print(verbose, chromosomes);

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  chromosome <- data$chromosome;
  x <- data$x;
  CT <- data$CT;
  muN <- data$muN;
  isSnp <- !is.na(muN);
  isHet <- isSnp & (muN == 1/2);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update the TCN segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  for (cc in seq(length=nbrOfChromosomes)) {
    chr <- chromosomes[cc];
    chrTag <- sprintf("chr%02d", chr);
    verbose && enter(verbose, sprintf("Chromosome %d ('%s') of %d", chr, chrTag, nbrOfChromosomes));
    rows <- which(is.element(segs[["chromosome"]], chr));
    segsCC <- segs[rows,,drop=FALSE];
    nbrOfSegmentsCC <- nrow(segsCC);
    verbose && cat(verbose, "Number of segments: ", nbrOfSegmentsCC);

    tcnIds <- sort(unique(segsCC[["tcn.id"]]));
    I <- length(tcnIds);
    for (ii in seq(length=I)) {
      tcnId <- tcnIds[ii];
      verbose && enter(verbose, sprintf("TCN segment %d ('%s') of %d", ii, tcnId, I));

      rowsII <- which(segsCC[["tcn.id"]] == tcnId);
      segsII <- segsCC[rowsII,,drop=FALSE];
  
      J <- nrow(segsII);
      for (jj in seq(length=J)) {
        verbose && enter(verbose, sprintf("DH segment %d of %d", jj, J));
        seg <- segsII[jj,,drop=FALSE];
  
        # (a) Start and end of TCN segment is (somewhere) in 
        #     the middle of the end of previous DH segment and
        #     the start of the current DH segment.
        if (jj == 1) {
          start <- seg[["tcn.loc.start"]];
        } else if (jj > 1) {
          starts <- c(segsII[jj-1,"dh.loc.end"], seg[["dh.loc.start"]]);
          start <- mean(starts);
        }
        if (jj == J) {
          end <- seg[["tcn.loc.end"]];
        } else if (jj < J) {
          ends <- c(seg[["dh.loc.end"]], segsII[jj+1,"dh.loc.start"]);
          end <- mean(ends);
        }
  
        # (b) Find the units within this first guess of the TCN segment
        #     HB: (x <= end) or (x < end); is it possible that we
        #         include the same locus in two segments? /HB 2010-09-21
        units <- (start <= x & x <= end);
        if (!is.na(chr)) units <- units & (chromosome == chr);
  
        # (c) Adjust the start and end of the TCN segment
        xJJ <- x[units];
        startSupport <- min(xJJ, na.rm=TRUE);
        endSupport <- max(xJJ, na.rm=TRUE);
        start <- min(start, startSupport, na.rm=TRUE);
        end <- max(end, endSupport, na.rm=TRUE);
  
        # (d) Identify the actual units
        units <- (start <= x & x <= end);
        if (!is.na(chr)) units <- units & (chromosome == chr);
        units <- whichVector(units);
        verbose && cat(verbose, "Units:");
        verbose && str(verbose, units);
  
        # Update the segment, estimates and counts
        seg[["tcn.loc.start"]] <- start;
        seg[["tcn.loc.end"]] <- end;
        seg[["tcn.mean"]] <- mean(CT[units], na.rm=TRUE);
        seg[["tcn.num.mark"]] <- length(units);
        seg[["tcn.num.snps"]] <- sum(isSnp[units]);
        seg[["tcn.num.hets"]] <- sum(isHet[units]);
  
        segsII[jj,] <- seg;
        rm(xJJ, start, end, units, seg);

        verbose && exit(verbose);
      } # for (jj ...)

      segsCC[rowsII,] <- segsII;
      rm(rowsII, segsII);
      verbose && exit(verbose);
    } # for (ii ...)

    segs[rows,] <- segsCC;
    rm(rows, segsCC);
    verbose && exit(verbose);
  } # for (cc ...)

  verbose && enter(verbose, "Update (C1,C2) per segment");
  # Append (C1,C2) estimates
  tcn <- segs$tcn.mean;
  dh <- segs$dh.mean;
  C1 <- 1/2*(1-dh)*tcn;
  C2 <- tcn - C1;
  segs$c1.mean <- C1;
  segs$c2.mean <- C2;
  verbose && exit(verbose);


  # Return results
  fitS <- fit;
  fitS$data <- data;
  fitS$output <- segs;

  verbose && exit(verbose);

  fitS;
}) # postsegmentTCN()



############################################################################
# HISTORY:
# 2010-11-21
# o Adjusted postsegmentTCN() such that the updated TCN segment boundaries 
#   are the maximum of the DH segment and the support by the loci.  This
#   means that postsegmentTCN() will work as expected both when signals
#   where segmented with 'joinSegments' being TRUE or FALSE.
# 2010-10-25
# o Now subsetByDhSegments() for PairedPSCBS handles the rare case when
#   markers with the same positions are split in two different segments.
# o Renamed subsetBySegments() for PairedPSCBS to subsetByDhSegments().
# 2010-09-26
# o Now subsetBySegments() for PairedPSCBS handles multiple chromosomes.
# o Now postsegmentTCN() PairedPSCBS handles multiple chromosomes.
# 2010-09-21
# o Added postsegmentTCN() for PairedPSCBS.
# 2010-09-19
# o BUG FIX: plot() used non-defined nbrOfLoci; now length(x).
# 2010-09-15
# o Added subsetBySegments().
# o Added linesC1C2() and arrowsC1C2().
# o Now the default 'cex' for pointsC1C2() corresponds to 'dh.num.mark'.
# o Now extractTotalAndDH() also returns 'dh.num.mark'.
# 2010-09-08
# o Added argument 'add=FALSE' to plot().
# o Added plotC1C2().
# o Added extractTotalAndDH() and extractMinorMajorCNs().
# 2010-09-04
# o Added drawLevels() for PairedPSCBS.
# o Added as.data.frame() and print() for PairedPSCBS.
# 2010-09-03
# o Added plot() for PairedPSCBS.
# o Created.
############################################################################
