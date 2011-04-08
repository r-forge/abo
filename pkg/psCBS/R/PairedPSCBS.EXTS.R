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

  # Sanity check
  stopifnot(nrow(dhSegRows) == nrow(tcnSegRows));
  stopifnot(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE));
#  stopifnot(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE));
  stopifnot(all(dhSegRows[,1] <= dhSegRows[,2], na.rm=TRUE));
  stopifnot(all(dhSegRows[-nrow(dhSegRows),2] < dhSegRows[-1,1], na.rm=TRUE));


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
      if (!force && J == 1) {
        verbose && cat(verbose, "Nothing todo. Only one DH segmentation. Skipping.");    
        verbose && exit(verbose);
        next;    
      }

      verbose && cat(verbose, "Rows:");
      verbose && print(verbose, rowsII);
      segsII <- segsCC[rowsII,,drop=FALSE];

      tcnSegRowsII <- tcnSegRowsCC[rowsII,,drop=FALSE];
      dhSegRowsII <- dhSegRowsCC[rowsII,,drop=FALSE];

      verbose && cat(verbose, "TCN & DH segRows before:");
      verbose && print(verbose, tcnSegRowsII);
      verbose && print(verbose, dhSegRowsII);

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

          # (d) Adjust for DH boundaries
          if (jj > 1L) {
            minIdx <- tcnSegRowsII[jj-1L,2L, drop=TRUE];
            units <- units[units > minIdx];
          }
          if (jj < J) {
            maxIdx <- dhSegRowsII[jj+1L,1L, drop=TRUE];
            units <- units[units < maxIdx];
          }

          tcnSegRow <- range(units);
          verbose && printf(verbose, "[idxStart,idxEnd] = [%d,%d]\n", tcnSegRow[1], tcnSegRow[2]);
          verbose && cat(verbose, "Number of non-missing TCN loci: ", length(units));
        } else {
          throw("Not implemented yet.")  # /HB 2010-12-02
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

      verbose && cat(verbose, "TCN & DH segRows afterward:");
      verbose && print(verbose, tcnSegRowsII);
      verbose && print(verbose, dhSegRowsII);

##print(segsII);

      # Sanity check
      nbrOfTCNsAfter <- sum(segsII[,"tcn.num.mark"], na.rm=TRUE);
      verbose && cat(verbose, "Number of TCNs before: ", nbrOfTCNsBefore);
      verbose && cat(verbose, "Number of TCNs after: ", nbrOfTCNsAfter);
      stopifnot(nbrOfTCNsAfter >= nbrOfTCNsBefore);

      # Sanity check
      stopifnot(nrow(dhSegRowsII) == nrow(tcnSegRowsII));
      stopifnot(all(tcnSegRowsII[,1] <= tcnSegRowsII[,2], na.rm=TRUE));
      stopifnot(all(tcnSegRowsII[-nrow(tcnSegRowsII),2] < tcnSegRowsII[-1,1], na.rm=TRUE));
      stopifnot(all(dhSegRowsII[,1] <= dhSegRowsII[,2], na.rm=TRUE));
      stopifnot(all(dhSegRowsII[-nrow(dhSegRowsII),2] < dhSegRowsII[-1,1], na.rm=TRUE));

      segsCC[rowsII,] <- segsII;
      tcnSegRowsCC[rowsII,] <- tcnSegRowsII;

      rm(rowsII, segsII);
      verbose && exit(verbose);
    } # for (ii ...)

    # Sanity check
    stopifnot(nrow(segsCC) == length(rows));

    # Sanity check
    stopifnot(nrow(dhSegRowsCC) == nrow(tcnSegRowsCC));
    stopifnot(all(tcnSegRowsCC[,1] <= tcnSegRowsCC[,2], na.rm=TRUE));
    stopifnot(all(tcnSegRowsCC[-nrow(tcnSegRowsCC),2] < tcnSegRowsCC[-1,1], na.rm=TRUE));
    stopifnot(all(dhSegRowsCC[,1] <= dhSegRowsCC[,2], na.rm=TRUE));
    stopifnot(all(dhSegRowsCC[-nrow(dhSegRowsCC),2] < dhSegRowsCC[-1,1], na.rm=TRUE));

    segs[rows,] <- segsCC;
    tcnSegRows[rows,] <- tcnSegRowsCC;

    rm(rows, segsCC);
    verbose && exit(verbose);
  } # for (cc ...)

  # Sanity check
  stopifnot(nrow(dhSegRows) == nrow(tcnSegRows));
  stopifnot(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE));
  stopifnot(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE));
  stopifnot(all(dhSegRows[,1] <= dhSegRows[,2], na.rm=TRUE));
  stopifnot(all(dhSegRows[-nrow(dhSegRows),2] < dhSegRows[-1,1], na.rm=TRUE));

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

  # Sanity check
  tcnSegRows <- fitS$tcnSegRows;
  dhSegRows <- fitS$dhSegRows;
  stopifnot(nrow(dhSegRows) == nrow(tcnSegRows));
  stopifnot(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE));
  stopifnot(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE));
  stopifnot(all(dhSegRows[,1] <= dhSegRows[,2], na.rm=TRUE));
  stopifnot(all(dhSegRows[-nrow(dhSegRows),2] < dhSegRows[-1,1], na.rm=TRUE));

  # Update 'flavor'
  fitS$params$flavor <- gsub(",", "&", flavor, fixed=TRUE);

  verbose && exit(verbose);

  fitS;
}) # postsegmentTCN()




setMethodS3("extractByRegion", "PairedPSCBS", function(this, region, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'region':
  region <- Arguments$getIndex(region, max=nbrOfSegments(this));

  extractByRegions(this, regions=region, ...);
}) # extractByRegion()



setMethodS3("extractByRegions", "PairedPSCBS", function(this, regions, ..., verbose=FALSE) {
  fit <- this;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  updateSegRows <- function(segRows, regions) {
    segRows <- segRows[regions,,drop=FALSE];
    ns <- segRows[,2] - segRows[,1] + 1L;
    from <- c(1L, cumsum(ns)[-length(ns)]);
    to <- from + (ns - 1L);
    segRows[,1] <- from;
    segRows[,2] <- to;
    verbose && str(verbose, segRows);
    segRows;
  } # updateSegRows()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'regions':
  regions <- Arguments$getIndices(regions, max=nbrOfSegments(fit));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

 
  verbose && enter(verbose, "Extracting subset by regions");

  verbose && cat(verbose, "Number of segments: ", length(regions));
  verbose && str(verbose, regions);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- fit$data;
  tcnSegRows <- fit$tcnSegRows;
  dhSegRows <- fit$dhSegRows;
  segs <- fit$output;
  params <- fit$params;

  # Sanity checks
  stopifnot(all(!is.na(data$chromosome) & !is.na(data$x)));
  stopifnot(length(tcnSegRows) == length(dhSegRows));

  # Sanity checks
  if (!params$joinSegments) {
    throw("Cannot bootstrap TCN and DH by segments unless PSCNs are segmented using joinSegments=TRUE.");
  } 

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Subset segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Update 'output'");
  segsT <- segs[regions,,drop=FALSE];
  verbose && str(verbose, segsT);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Subset data accordingly
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Update 'data'");

  segRows <- tcnSegRows;
  segRows <- segRows[regions,,drop=FALSE];
  from <- segRows[[1]];
  to <- segRows[[2]];
  ok <- (!is.na(from) & !is.na(to));
  from <- from[ok];
  to <- to[ok];
  keep <- logical(nrow(data));
  for (rr in seq(along=from)) {
    keep[from[rr]:to[rr]] <- TRUE;
  }
  keep <- which(keep);
  verbose && printf(verbose, "Identified %d (%.2f%%) of %d data rows:\n", length(keep), 100*length(keep)/nrow(data), nrow(data));
  verbose && str(verbose, keep);

  dataT <- data[keep,,drop=FALSE];
  verbose && str(verbose, dataT);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update 'segRows'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Update 'segRows'");

  segRows <- updateSegRows(tcnSegRows, regions=regions);
  d <- tcnSegRows[regions,] - segRows;
  # Sanity check
  stopifnot(identical(d[,1], d[,2]));
  d <- d[,1];
  verbose && cat(verbose, "Row deltas:");
  verbose && str(verbose, d);

  tcnSegRows <- tcnSegRows[regions,,drop=FALSE] - d;
  verbose && str(verbose, tcnSegRows);
  # Sanity checks
  stopifnot(max(tcnSegRows, na.rm=TRUE) <= nrow(dataT));

  dhSegRows <- dhSegRows[regions,,drop=FALSE] - d;
  verbose && str(verbose, dhSegRows);
  # Sanity checks
  stopifnot(max(dhSegRows, na.rm=TRUE) <= nrow(dataT));

  verbose && exit(verbose);


  # Create new object
  res <- fit;
  res$data <- dataT;
  res$output <- segsT;
  res$tcnSegRows <- tcnSegRows;
  res$dhSegRows <- dhSegRows;

  verbose && exit(verbose);

  res;
}) # extractByRegions()



setMethodS3("estimateStdDevForHeterozygousBAF", "PairedPSCBS", function(this, tauDH=0.20, tauTCN=5, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'tauDH':
  tauDH <- Arguments$getDouble(tauDH, range=c(0,1));

  # Argument 'tauTCN':
  tauTCN <- Arguments$getDouble(tauTCN, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Estimating standard deviation of tumor BAFs for heterozygous SNPs");
  verbose && cat(verbose, "DH threshold: ", tauDH);
  verbose && cat(verbose, "TCN threshold: ", tauTCN);

  segs <- as.data.frame(this);

  verbose && cat(verbose, "Number of segments: ", nrow(segs));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find segments to be used for the estimation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find segments that have low DHs
  idxsDH <- which(segs$dh.mean <= tauDH);
  verbose && cat(verbose, "Identified segments with small DH levels: ", length(idxsDH));
  verbose && str(verbose, idxsDH);

  # Sanity check
  if (length(idxsDH) == 0) {
    throw("Cannot estimate standard deviation.  There exist no segments with DH less or equal to the given threshold: ", tauDH); 
  }

  # Find segments that have low TCNs
  idxsTCN <- which(segs$tcn.mean <= tauTCN);
  verbose && cat(verbose, "Identified segments with small TCN levels: ", length(idxsTCN));
  verbose && str(verbose, idxsTCN);

  # Sanity check
  if (length(idxsTCN) == 0) {
    throw("Cannot estimate standard deviation.  There exist no segments with TCN less or equal to the given threshold: ", tauTCN); 
  }

  # Segments with small DH and small TCN
  idxs <- intersect(idxsDH, idxsTCN);
  verbose && cat(verbose, "Identified segments with small DH and small TCN levels: ", length(idxs));
  verbose && str(verbose, idxs);

  # Sanity check
  if (length(idxs) == 0) {
    throw("Cannot estimate standard deviation.  There exist no segments with small DH and small TCN.");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and estimate parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract those segments
  verbose && enter(verbose, "Extracting identified segments");
  fitT <- extractByRegions(this, idxs);
  verbose && exit(verbose);

  # Get the tumor BAFs for the heterozygous SNPs
  verbose && enter(verbose, "Extracting BAFs for the heterozygous SNPs");
  beta <- with(fitT$data, betaTN[muN == 1/2]);
  verbose && str(verbose, beta);
  verbose && exit(verbose);

  # Estimate the standard deviation for those
  sd <- mad(beta, na.rm=TRUE);
  verbose && cat(verbose, "Estimated standard deviation: ", sd);


  verbose && exit(verbose);

  sd;
}) # estimateStdDevForHeterozygousBAF()




setMethodS3("estimateMeanForDH", "PairedPSCBS", function(this, tauDH=0.20, tauTCN=5, robust=TRUE, trim=0.05, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'tauDH':
  tauDH <- Arguments$getDouble(tauDH, range=c(0,1));

  # Argument 'tauTCN':
  tauTCN <- Arguments$getDouble(tauTCN, range=c(0,Inf));

  # Argument 'robust':
  robust <- Arguments$getLogical(robust);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Estimating mean of tumor DHs for heterozygous SNPs");
  verbose && cat(verbose, "DH threshold: ", tauDH);
  verbose && cat(verbose, "TCN threshold: ", tauTCN);
  verbose && cat(verbose, "Robust estimator: ", robust);
  verbose && cat(verbose, "Trim: ", trim);

  segs <- as.data.frame(this);

  verbose && cat(verbose, "Number of segments: ", nrow(segs));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find segments to be used for the estimation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find segments that have low DHs
  idxsDH <- which(segs$dh.mean <= tauDH);
  verbose && cat(verbose, "Identified segments with small DH levels: ", length(idxsDH));
  verbose && str(verbose, idxsDH);

  # Sanity check
  if (length(idxsDH) == 0) {
    throw("Cannot estimate standard deviation.  There exist no segments with DH less or equal to the given threshold: ", tauDH); 
  }

  # Find segments that have low TCNs
  idxsTCN <- which(segs$tcn.mean <= tauTCN);
  verbose && cat(verbose, "Identified segments with small TCN levels: ", length(idxsTCN));
  verbose && str(verbose, idxsTCN);

  # Sanity check
  if (length(idxsTCN) == 0) {
    throw("Cannot estimate standard deviation.  There exist no segments with TCN less or equal to the given threshold: ", tauTCN); 
  }

  # Segments with small DH and small TCN
  idxs <- intersect(idxsDH, idxsTCN);
  verbose && cat(verbose, "Identified segments with small DH and small TCN levels: ", length(idxs));
  verbose && str(verbose, idxs);

  # Sanity check
  if (length(idxs) == 0) {
    throw("Cannot estimate standard deviation.  There exist no segments with small DH and small TCN.");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and estimate parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract those segments
  verbose && enter(verbose, "Extracting identified segments");
  fitT <- extractByRegions(this, idxs);
  verbose && exit(verbose);

  # Get the tumor DHs for the heterozygous SNPs
  verbose && enter(verbose, "Extracting DHs for the heterozygous SNPs");
  rho <- with(fitT$data, rho[muN == 1/2]);
  verbose && str(verbose, rho);
  verbose && exit(verbose);

  # Estimate the average for those
  rho <- rho[is.finite(rho)];
  if (robust) {
    mu <- median(rho, na.rm=FALSE);
    qlow <- quantile(rho, probs=0.05, na.rm=FALSE);
    delta <- mu-qlow;
    print(list(qlow=qlow, mu=mu, delta=delta, "mu+delta"=mu+delta));
  } else {
    mu <- mean(rho, trim=trim, na.rm=FALSE);
  }
  verbose && cat(verbose, "Estimated mean: ", mu);


  verbose && exit(verbose);

  mu;
}) # estimateMeanForDH()



setMethodS3("estimateHighDHQuantileAtAB", "PairedPSCBS", function(this, quantile=0.99, scale=1, tauDH=0.20, tauTCN=5, robust=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'quantile':
  quantile <- Arguments$getDouble(quantile, range=c(0.5,1));

  # Argument 'scale':
  scale <- Arguments$getDouble(scale, range=c(0,Inf));

  # Argument 'tauDH':
  tauDH <- Arguments$getDouble(tauDH, range=c(0,1));

  # Argument 'tauTCN':
  tauTCN <- Arguments$getDouble(tauTCN, range=c(0,Inf));

  # Argument 'robust':
  robust <- Arguments$getLogical(robust);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Estimating DH quantile of tumor DHs for heterozygous SNPs");
  verbose && cat(verbose, "DH threshold: ", tauDH);
  verbose && cat(verbose, "TCN threshold: ", tauTCN);
  verbose && cat(verbose, "Robust estimator: ", robust);
  verbose && cat(verbose, "Scale factor: ", scale);

  segs <- as.data.frame(this);

  verbose && cat(verbose, "Number of segments: ", nrow(segs));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find segments to be used for the estimation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Finding some segments that are likely to in allelic balance (AB)");

  # Find some segments that have low DHs
  idxsDH <- which(segs$dh.mean <= tauDH);
  verbose && cat(verbose, "Identified segments with small DH levels: ", length(idxsDH));
  verbose && str(verbose, idxsDH);

  # Sanity check
  if (length(idxsDH) == 0) {
    throw("Cannot estimate standard deviation.  There exist no segments with DH less or equal to the given threshold: ", tauDH); 
  }

  # Find segments that have low TCNs
  idxsTCN <- which(segs$tcn.mean <= tauTCN);
  verbose && cat(verbose, "Identified segments with small TCN levels: ", length(idxsTCN));
  verbose && str(verbose, idxsTCN);

  # Sanity check
  if (length(idxsTCN) == 0) {
    throw("Cannot estimate standard deviation.  There exist no segments with TCN less or equal to the given threshold: ", tauTCN); 
  }

  # Segments with small DH and small TCN
  idxs <- intersect(idxsDH, idxsTCN);
  verbose && cat(verbose, "Identified segments with small DH and small TCN levels: ", length(idxs));
  verbose && str(verbose, idxs);

  # Sanity check
  if (length(idxs) == 0) {
    throw("Cannot estimate standard deviation.  There exist no segments with small DH and small TCN.");
  }

  # Extract those segments
  verbose && enter(verbose, "Extracting identified segments");
  fitT <- extractByRegions(this, idxs);
  verbose && exit(verbose);

  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and estimate parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the tumor DHs for the heterozygous SNPs
  verbose && enter(verbose, "Extracting DHs for the heterozygous SNPs");
  rho <- with(fitT$data, rho[muN == 1/2]);
  verbose && str(verbose, rho);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimating the DH quantile
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Estimating the quantile of interest");
  verbose && cat(verbose, "Quantile: ", quantile);

  # Drop missing values
  rho <- rho[is.finite(rho)];

  if (robust) {
    lq <- quantile(rho, probs=1-quantile, na.rm=FALSE);
    verbose && printf(verbose, "Estimated lower quantile (%.3f): %f\n", 1-quantile, lq);
    mu <- median(rho, na.rm=FALSE);
    verbose && cat(verbose, "Estimated median: ", mu);
    delta <- mu-lq;
    verbose && printf(verbose, "Estimated \"spread\": %f\n", delta);
    uq <- mu + scale*delta;
    verbose && printf(verbose, "Scale parameter: %f\n", scale);
    qs <- c(lq, mu, mu+delta, uq);
    names(qs) <- sprintf("%.1f%%", 100*c(1-quantile, 0.5, quantile, 0.5+scale*(quantile-0.5)));
    names(qs)[3:4] <- sprintf("%s*", names(qs)[3:4]);
    attr(uq, "quantiles") <- qs;
  } else {
    uq <- quantile(rho, probs=quantile, na.rm=FALSE);
  }

  names(uq) <- uq;
  verbose && printf(verbose, "Estimated upper quantile (%.3f): %f\n", quantile, uq);

  verbose && exit(verbose);

  verbose && exit(verbose);

  uq;
}) # estimateHighDHQuantileAtAB()



setMethodS3("estimateTauAB", "PairedPSCBS", function(this, scale=NULL, flavor=c("q(DH)", "mad(hBAF)", "median(DH)"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Estimating DH threshold for calling allelic imbalances");
  verbose && cat(verbose, "flavor: ", flavor);

  if (flavor == "mad(hBAF)") {
    if (is.null(scale)) scale <- 3;
    verbose && cat(verbose, "scale: ", scale);
    # sigma = mad(hBAF) = 1.4826*median(|hBAF-m|),
    # where m = median(hBAF) ~= 1/2
    sd <- estimateStdDevForHeterozygousBAF(this, ..., verbose=verbose);
    verbose && printf(verbose, "sd: %.3g\n", sd);
    tau <- scale * sd;
  } else if (flavor == "median(DH)") {
    if (is.null(scale)) scale <- 3;
    verbose && cat(verbose, "scale: ", scale);
    # sigma = 1/2*1.4826*median(|hBAF-1/2|), 
    # because DH = 2*|hBAF-1/2|
    mu <- estimateMeanForDH(this, ..., verbose=verbose);
    verbose && printf(verbose, "mu: %.3g\n", mu);
    sd <- 1/2 * 1.4826 * mu;
    verbose && printf(verbose, "sd: %.3g\n", sd);
    tau <- scale * sd;
  } else if (flavor == "q(DH)") {
    if (is.null(scale)) scale <- 1;
    verbose && cat(verbose, "scale: ", scale);
    tau <- estimateHighDHQuantileAtAB(this, scale=scale, ..., verbose=verbose);
  } else {
    throw("Unkown flavor: ", flavor);
  }

##   } else if (flavor == "DHskew") {
##     fit <- this;
##     if (is.null(fit$output$dh.skew)) {
##       verbose && enter(verbose, "Estimating DH skewness for each segment");
##       fit <- applyByRegion(fit, FUN=.addTcnDhStatitics, verbose=less(verbose, 5));
##       verbose && exit(verbose);
##     }
##     mu <- fit$output$dh.mean;
##     skew <- fit$output$dh.skew;
## 
##     tauSkew <- -0.55;
##     keep <- which(skew < tauSkew);
##     verbose && printf(verbose, "Number of segments heavily skewed (< %.3f): %d\n", tauSkew, length(keep));
##     # Sanity check
##     if (length(keep) == 0) {
##       throw("Cannot estimate DH threshold for AB. No segments with strong skewness exists.");
##     }
##     tauDH <- median(mu[keep], na.rm=TRUE);
##     verbose && printf(verbose, "tauDH: %.3g\n", tauDH);
##     tauDH <- 1.10*tauDH;
##     verbose && printf(verbose, "Adjusted +10%% tauDH: %.3g\n", tauDH);
## 
##     # sigma = 1/2*1.4826*median(|hBAF-1/2|), 
##     # because DH = 2*|hBAF-1/2|
##     mu <- estimateMeanForDH(this, tau=tauDH, ...);
##     verbose && printf(verbose, "mu: %.3g\n", mu);
##     sd <- 1/2 * 1.4826 * mu;
##     verbose && printf(verbose, "sd: %.3g\n", sd);
##  }

  verbose && printf(verbose, "tau: %.3g\n", tau);

  verbose && exit(verbose);

  tau;
}) # estimateTauAB()


############################################################################
# HISTORY:
# 2011-04-08
# o BUG FIX: postsegmentTCN() for PairedPSCBS could generate an invalid
#   'tcnSegRows' matrix, where the indices for two consecutive segments
#   would overlap, which is invalid.
# 2011-04-05
# o BUG FIX: estimateHighDHQuantileAtAB() for PairedPSCBS would throw
#   an error on an undefined 'trim' if verbose output was used.
# 2011-02-17
# o Added arguments 'robust' and 'trim' to estimateMeanForDH().
# 2011-02-03
# o Added argument 'tauTCN' to estimateMeanForDH().
# 2011-01-27
# o Added flavor="DHskew" to estimateTauAB().
# o Added flavor="DH" to estimateTauAB() to estimate from DH instead 
#   of hBAF.  As argued by the equations in the comments, these two
#   approaches gives virtually the same results.  The advantage with the
#   DH approach is that it requires one less degree of freedom.
# o Added estimateMeanForDH().
# 2011-01-18
# o BUG FIX: 'tcnSegRows' and 'dhSegRows' where not updated by
#   extractByRegions() for PairedPSCBS.
# 2011-01-14
# o Added estimateTauAB() for estimating the DeltaAB parameter.
# o Added estimateStdDevForHeterozygousBAF() for PairedPSCBS.
# o BUG FIX: extractByRegions() did not handle the case where multiple loci
#   at the same position are split up in two different segments.
# 2011-01-12
# o Added extractByRegions() and extractByRegion() for PairedPSCBS.
# o Now postsegmentTCN(..., force=TRUE) for PairedPSCBS also updates
#   the TCN estimates even for segments where the DH segmentation did
#   not find any additional change points.
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
