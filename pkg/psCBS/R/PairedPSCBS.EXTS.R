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



setMethodS3("postsegmentTCN", "PairedPSCBS", function(fit, ..., force=FALSE, verbose=FALSE) {
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

  flavor <- fit$params$flavor;
  if (!force && regexpr("&", flavor, fixed=TRUE) != -1) {
    verbose && cat(verbose, "Nothing to do. Already postsegmentTCN:ed: ", flavor);
    verbose && exit(verbose);
    return(fit);
  }

  joinSegments <- fit$params$joinSegments;
  if (!joinSegments) {
    throw("Postsegmentation of TCNs is only implemented for the case when joinSegments=TRUE: ", joinSegments);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- fit$data;

  segs <- fit$output;
  keep <- is.finite(segs$chromosome);
  segs <- segs[keep,,drop=FALSE];
  tcnSegRows <- fit$tcnSegRows[keep,,drop=FALSE];
  dhSegRows <- fit$dhSegRows[keep,,drop=FALSE];

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
  chromosomes <- getChromosomes(fit);
  nbrOfChromosomes <- length(chromosomes);
  verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes);
  verbose && print(verbose, chromosomes);

  for (cc in seq(length=nbrOfChromosomes)) {
    chr <- chromosomes[cc];
    chrTag <- sprintf("chr%02d", chr);
    verbose && enter(verbose, sprintf("Chromosome %d ('%s') of %d", chr, chrTag, nbrOfChromosomes));
    rows <- which(is.element(segs[["chromosome"]], chr));
    verbose && cat(verbose, "Rows:");
    verbose && print(verbose, rows);

    segsCC <- segs[rows,,drop=FALSE];
    tcnSegRowsCC <- tcnSegRows[rows,,drop=FALSE];
    dhSegRowsCC <- dhSegRows[rows,,drop=FALSE];
    nbrOfSegmentsCC <- nrow(segsCC);
    verbose && cat(verbose, "Number of segments: ", nbrOfSegmentsCC);

    tcnIds <- sort(unique(segsCC[["tcn.id"]]));
    I <- length(tcnIds);
    for (ii in seq(length=I)) {
      tcnId <- tcnIds[ii];
      verbose && enter(verbose, sprintf("TCN segment #%d ('%s') of %d", ii, tcnId, I));

      rowsII <- which(segsCC[["tcn.id"]] == tcnId);
      J <- length(rowsII);
      # Nothing todo?
      if (J == 1) {
        verbose && cat(verbose, "Nothing todo. Only on DH segmentation. Skipping.");    
        verbose && exit(verbose);
        next;    
      }

      verbose && cat(verbose, "Rows:");
      verbose && print(verbose, rowsII);
      segsII <- segsCC[rowsII,,drop=FALSE];

      tcnSegRowsII <- tcnSegRowsCC[rowsII,,drop=FALSE];
      dhSegRowsII <- dhSegRowsCC[rowsII,,drop=FALSE];
##print(tcnSegRowsII);
##print(dhSegRowsII);
      tcnSegRowsIIBefore <- tcnSegRowsII;
      nbrOfTCNsBefore <- segsII[1,"tcn.num.mark"];

      for (jj in seq(length=J)) {
        verbose && enter(verbose, sprintf("DH segment #%d of %d", jj, J));
        seg <- segsII[jj,,drop=FALSE];
        tcnSegRow <- unlist(tcnSegRowsII[jj,,drop=FALSE]);
        dhSegRow <- unlist(dhSegRowsII[jj,,drop=FALSE]);
        # Sanity check
        stopifnot(all(is.na(tcnSegRow)) || (tcnSegRow[1] <= tcnSegRow[2]));
        stopifnot(all(is.na(dhSegRow)) || (dhSegRow[1] <= dhSegRow[2]));

        # Sanity check
        idxsTCN <- tcnSegRow[1]:tcnSegRow[2];
        nbrOfTCNs <- sum(!is.na(CT[idxsTCN]));
        stopifnot(nbrOfTCNs == nbrOfTCNsBefore);
  
        if (joinSegments) {
          # (a) The TCN segment should have identical (start,end) boundaries as the DH region
          xStart <- seg[["dh.loc.start"]];
          xEnd <- seg[["dh.loc.end"]];
          verbose && printf(verbose, "[xStart,xEnd] = [%.0f,%.0f]\n", xStart, xEnd);
          stopifnot(xStart <= xEnd);

          # (b) Identify units
          units <- which(chromosome == chr & xStart <= x & x <= xEnd);
          tcnSegRow <- range(units);
          verbose && printf(verbose, "[idxStart,idxEnd] = [%d,%d]\n", tcnSegRow[1], tcnSegRow[2]);
          verbose && cat(verbose, "Number of TCN loci: ", length(units));

          # (c) Adjust for missing values
          keep <- which(!is.na(CT[units]));
          units <- units[keep];
          tcnSegRow <- range(units);
          verbose && printf(verbose, "[idxStart,idxEnd] = [%d,%d]\n", tcnSegRow[1], tcnSegRow[2]);
          verbose && cat(verbose, "Number of non-missing TCN loci: ", length(units));
        } else {
          throw("Not implemented correctly.")  # /HB 2010-12-02
          # (a) Start and end of TCN segment is (somewhere) in 
          #     the middle of the end of previous DH segment and
          #     the start of the current DH segment.
          if (jj == 1) {
            xStart <- seg[["tcn.loc.start"]];
            idxStart <- tcnSegRow[1];
          } else if (jj > 1) {
            xStarts <- c(segsII[jj-1,"dh.loc.end"], seg[["dh.loc.start"]]);
            xStart <- mean(xStarts);
            idxStart <- max(dhSegRow[1], tcnSegRow[1], na.rm=TRUE);
            idxStart <- min(idxEnd+1L, idxStart);
          }
          if (jj == J) {
            xEnd <- seg[["tcn.loc.end"]];
            idxEnd <- tcnSegRow[2];
          } else if (jj < J) {
            xEnds <- c(seg[["dh.loc.end"]], segsII[jj+1,"dh.loc.start"]);
            xEnd <- mean(xEnds);
            idxEnd <- min(dhSegRow[2], tcnSegRow[2], na.rm=TRUE);
          }
   
          verbose && printf(verbose, "[xStart,xEnd] = [%.0f,%.0f]\n", xStart, xEnd);
          verbose && printf(verbose, "[idxStart,idxEnd] = [%d,%d]\n", idxStart, idxEnd);
          verbose && cat(verbose, "Number of initial loci: ", idxEnd-idxStart+1L);
  
          # Sanity check
          stopifnot(xStart <= xEnd);
          stopifnot(idxStart <= idxEnd);
  
          # (b) Identify units
          units <- which(chromosome == chr & xStart <= x & x <= xEnd);
          verbose && printf(verbose, "[idxStart,idxEnd] = [%d,%d]\n", min(units), max(units));
          verbose && cat(verbose, "Number of initial loci: ", length(units));
  
          # (c) Drop if it belongs to neighboring segment
          keep <- which(idxStart <= units & units <= idxEnd);
          nbrOfDropped <- length(units) - length(keep);
          verbose && cat(verbose, "Number of loci dropped: ", nbrOfDropped);
          units <- units[keep];
          xJJ <- x[units];
          xSupport <- range(xJJ, na.rm=FALSE);
          verbose && printf(verbose, "[xStart,xEnd] = [%.0f,%.0f]\n", xSupport[1], xSupport[2]);
          verbose && printf(verbose, "[idxStart,idxEnd] = [%d,%d]\n", min(units), max(units));
          verbose && cat(verbose, "Number of loci: ", length(units));
  
          # (c) Adjust the start and end of the TCN segment
          keep <- which(!is.na(xJJ));
          xJJ <- xJJ[keep];
          units <- units[keep];
          tcnSegRow <- range(units);
          xSupport <- range(xJJ, na.rm=FALSE);
          verbose && printf(verbose, "[xStart,xEnd] = [%.0f,%.0f]\n", xSupport[1], xSupport[2]);
          verbose && printf(verbose, "[idxStart,idxEnd] = [%d,%d]\n", tcnSegRow[1], tcnSegRow[2]);
          verbose && cat(verbose, "Number of loci: ", length(units));
  
          xRange <- range(c(xStart,xEnd), xSupport, na.rm=TRUE);

          xStart <- xRange[1];
          xEnd <- xRange[2];
        } # if (joinSegments)
  
        gamma <- mean(CT[units]);
        # Sanity check
        stopifnot(length(units) == 0 || !is.na(gamma));

        # Update the segment boundaries, estimates and counts
        seg[["tcn.loc.start"]] <- xStart;
        seg[["tcn.loc.end"]] <- xEnd;
        seg[["tcn.mean"]] <- gamma;
        seg[["tcn.num.mark"]] <- length(units);
        seg[["tcn.num.snps"]] <- sum(isSnp[units]);
        seg[["tcn.num.hets"]] <- sum(isHet[units]);

        # Sanity check
        stopifnot(nrow(seg) == length(jj));

        segsII[jj,] <- seg;
        tcnSegRowsII[jj,] <- tcnSegRow;

        verbose && exit(verbose);
      } # for (jj ...)

      # Sanity check
      stopifnot(nrow(segsII) == length(rowsII));

##print(segsII);
##print(tcnSegRowsIIBefore);
##print(tcnSegRowsII);

      # Sanity check
      nbrOfTCNsAfter <- sum(segsII[,"tcn.num.mark"], na.rm=TRUE);
      verbose && cat(verbose, "Number of TCNs before: ", nbrOfTCNsBefore);
      verbose && cat(verbose, "Number of TCNs after: ", nbrOfTCNsAfter);
      stopifnot(nbrOfTCNsAfter >= nbrOfTCNsBefore);


      segsCC[rowsII,] <- segsII;
      tcnSegRowsCC[rowsII,] <- tcnSegRowsII;

      rm(rowsII, segsII);
      verbose && exit(verbose);
    } # for (ii ...)

    # Sanity check
    stopifnot(nrow(segsCC) == length(rows));

    segs[rows,] <- segsCC;
    tcnSegRows[rows,] <- tcnSegRowsCC;

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
  keep <- which(is.finite(fit$output$chromosome));
  fitS <- fit;
  fitS$data <- data;
  fitS$output[keep,] <- segs;
  fitS$tcnSegRows[keep,] <- tcnSegRows;

  # Update 'flavor'
  fitS$params$flavor <- gsub(",", "&", flavor);

  verbose && exit(verbose);

  fitS;
}) # postsegmentTCN()



############################################################################
# HISTORY:
# 2010-12-02
# o Now postsegmentTCN() assert that total number of TCN loci before
#   and after is the same.
# o Now postsegmentTCN() assert that joinSegment is TRUE.
# 2010-12-01
# o Now postsegmentTCN() checks if it is already postsegmented.
# 2010-11-30
# o TODO: postsegmentTCN() does not make sure of 'dhLociToExclude'. Why?
# o Now postsegmentTCN() recognizes the new 'tcnLociToExclude'.
# 2010-11-28
# o BUG FIX: postsegmentTCN() did not handle loci with the same positions
#   and that are split in two different segments.  It also did not exclude
#   loci with missing values.
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
