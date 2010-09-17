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

  for (ff in seq(along=data)) {
    values <- data[[ff]];
    if (length(values) == nbrOfLoci) {
      values <- values[keep];
      data[[ff]] <- values;
    }
  }
  rm(keep, values); # Not needed anymore

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



setMethodS3("plot", "PairedPSCBS", function(x, pch=".", Clim=c(0,6), Blim=c(0,1), ..., add=FALSE) {
  # To please R CMD check
  fit <- x;

  # Extract the input data
  data <- fit$data;
  if (is.null(data)) {
    throw("Cannot plot segmentation results. No input data available.");
  }

  x <- data$x;
  CT <- data$CT;
  betaT <- data$betaT;
  betaN <- data$betaN;
  betaTN <- data$betaTN;
  muN <- data$muN;

  # Extract the segmentation
  segs <- fit$output;

  if (!add) {
    subplots(5, ncol=1);
    par(mar=c(1,4,1,2)+1);
  }

  plot(x, CT, pch=pch, col="gray", ylim=Clim);
  drawLevels(fit, what="tcn");

  col <- c("gray", "black")[(muN == 1/2) + 1];
  plot(x, betaN, pch=pch, col=col, ylim=Blim);
  plot(x, betaT, pch=pch, col=col, ylim=Blim);
  plot(x, betaTN, pch=pch, col=col, ylim=Blim);

  isSnp <- (!is.na(betaTN) & !is.na(muN));
  isHet <- isSnp & (muN == 1/2);
	  naValue <- as.double(NA);
  rho <- rep(naValue, length=nbrOfLoci);
  rho[isHet] <- 2*abs(betaTN[isHet]-1/2);
  plot(x, rho, pch=pch, col=col, ylim=Blim);
  drawLevels(fit, what="dh");
})


setMethodS3("drawLevels", "PairedPSCBS", function(fit, what=c("tcn", "dh"), ...) {
  # Argument 'what':
  what <- match.arg(what);

  # Get segmentation results
  segs <- as.data.frame(fit);

  # Extract subset of segments
  fields <- c("loc.start", "loc.end", "mean");
  fields <- sprintf("%s.%s", what, fields);
  segs <- segs[,fields, drop=FALSE];
  segs <- unique(segs);

  colnames(segs) <- c("loc.start", "loc.end", "seg.mean");

  dummy <- list(output=segs);
  class(dummy) <- "DNAcopy";
  drawLevels(dummy, ...);
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


setMethodS3("plotC1C2", "PairedPSCBS", function(fit, ..., xlab=expression(C[1]), ylab=expression(C[2]), Clim=c(0,4)) {
  plot(NA, xlim=Clim, ylim=Clim, xlab=xlab, ylab=ylab);
  abline(a=0, b=1, lty=3);
  pointsC1C2(fit, ...);
})


setMethodS3("pointsC1C2", "PairedPSCBS", function(fit, cex=NULL, ...) {
  data <- extractC1C2(fit);
  X <- data[,1:2,drop=FALSE];
  n <- data[,4,drop=TRUE];
  n <- sqrt(n);
  w <- n / sum(n, na.rm=TRUE);

  if (is.null(cex)) {
    cex <- w;
    cex <- cex / mean(cex, na.rm=TRUE);
    cex <- cex + 1/2;
  }

  points(X, cex=cex, ...);
})


setMethodS3("linesC1C2", "PairedPSCBS", function(fit, ...) {
  xy <- extractMinorMajorCNs(fit);
  xy <- xy[,1:2,drop=FALSE];
  lines(xy, ...);
}) # linessC1C2()


setMethodS3("arrowsC1C2", "PairedPSCBS", function(fit, length=0.05, ...) {
  xy <- extractMinorMajorCNs(fit);
  xy <- xy[,1:2,drop=FALSE];
  x <- xy[,1,drop=TRUE];
  y <- xy[,2,drop=TRUE];
  s <- seq(length=length(x)-1);
  arrows(x0=x[s],y=y[s], x1=x[s+1],y1=y[s+1], code=2, length=length, ...);
}) # arrowsC1C2()




setMethodS3("plotDeltaC1C2", "PairedPSCBS", function(fit, ..., xlab=expression(Delta*C[1]), ylab=expression(Delta*C[2]), Clim=c(-2,2)) {
  plot(NA, xlim=Clim, ylim=Clim, xlab=xlab, ylab=ylab);
  abline(h=0, lty=3);
  abline(v=0, lty=3);
  pointsDeltaC1C2(fit, ...);
})


setMethodS3("pointsDeltaC1C2", "PairedPSCBS", function(fit, ...) {
  data <- extractDeltaC1C2(fit);
  X <- data[,1:2,drop=FALSE];
  points(X, ...);
})


setMethodS3("linesDeltaC1C2", "PairedPSCBS", function(fit, ...) {
  xy <- extractDeltaC1C2(fit);
  xy <- xy[,1:2,drop=FALSE];
  lines(xy, ...);
})


setMethodS3("arrowsDeltaC1C2", "PairedPSCBS", function(fit, length=0.05, ...) {
  xy <- extractDeltaC1C2(fit);
  xy <- xy[,1:2,drop=FALSE];
  x <- xy[,1,drop=TRUE];
  y <- xy[,2,drop=TRUE];
  s <- seq(length=length(x)-1);
  arrows(x0=x[s],y=y[s], x1=x[s+1],y1=y[s+1], code=2, length=length, ...);
})






############################################################################
# HISTORY:
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
