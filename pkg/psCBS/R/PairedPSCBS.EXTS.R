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


setMethodS3("subsetBySegments", "PairedPSCBS", function(fit, idxs, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the segmentation result
  segs <- fit$output;
  stopifnot(!is.null(segs));
  nbrOfSegments <- nrow(segs);

  # Argument 'idxs':
  idxs <- Arguments$getIndices(idxs, max=nbrOfSegments);


  # Extract the data
  data <- fit$data;
  stopifnot(!is.null(data));
  x <- data$x;
  nbrOfLoci <- length(x);

  # Subset the region-level data
  segs <- segs[idxs,,drop=FALSE];

  nbrOfSegments <- nrow(segs);

  # Subset the locus-level data
  keep <- logical(nbrOfLoci);
  for (kk in seq(length=nbrOfSegments)) {
    xRange <- as.numeric(segs[kk,c("dh.loc.start", "dh.loc.end")]);
    # Identify all SNPs in the region
    keep <- keep | (xRange[1] <= x & x <= xRange[2]);
  } # for (kk ...)
  keep <- whichVector(keep);
  data <- data[keep,,drop=FALSE];
  rm(keep); # Not needed anymore

  fitS <- fit;
  fitS$data <- data;
  fitS$output <- segs;

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
  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  chromosome <- segs$chromosome;
  chromosomes <- sort(unique(chromosome), na.last=TRUE);
  nbrOfChromosomes <- length(chromosomes);
  if (nbrOfChromosomes > 1) {
    # Ignore NA chromosomes, because they are dividers
    chromosomes <- chromosomes[is.finite(chromosomes)];
  }
  nbrOfChromosomes <- length(chromosomes);

  verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes);
  verbose && print(verbose, chromosomes);


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
    verbose && enter(verbose, sprintf("Chromosome %d ('%s') of %d", chr, chrTag, 22));
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
        units <- (chromosome == chr & start <= x & x <= end);
  
        # (c) Adjust the start and end of the TCN segment
        xJJ <- x[units];
        start <- min(xJJ, na.rm=TRUE);
        end <- max(xJJ, na.rm=TRUE);
  
        # (d) Identify the actual units
        units <- (chromosome == chr & start <= x & x <= end);
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

  # Return results
  fitS <- fit;
  fitS$data <- data;
  fitS$output <- segs;

  verbose && exit(verbose);

  fitS;
}) # postsegmentTCN()



############################################################################
# HISTORY:
# 2010-09-26
# o Now postsegmentTCN() processes multiple chromosomes.
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
