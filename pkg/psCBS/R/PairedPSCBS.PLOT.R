###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod plotTracks
# @alias plotTracks
#
# @title "Plots parental specific copy numbers along the genome"
#
# \description{
#  @get "title" for one or more chromosomes.
#  It is possible to specify what type of tracks to plot.
#  Each type of track is plotted in its own panel.
# }
#
# @synopsis
#
# \arguments{
#   \item{x}{A @see "PairedPSCBS" result object.}
#   \item{tracks}{A @character @vector specifying what types of tracks to plot.}
#   \item{pch}{The type of points to use.}
#   \item{Clim}{The range of copy numbers.}
#   \item{Blim}{The range of allele B fractions (BAFs) and 
#     decrease of heterozygosity (DHs).}
#   \item{xScale}{The scale factor used for genomic positions.}
#   \item{...}{Not used.}
#   \item{add}{If @TRUE, the panels plotted are added to the existing plot,
#     otherwise a new plot is created.}
# }
#
# \value{
#   Returns nothing.
# }
# 
# @author
#
# @keyword IO
# @keyword internal
#*/########################################################################### 
setMethodS3("plotTracks", "PairedPSCBS", function(x, tracks=c("tcn", "rho", "tcn,c1,c2", "betaN", "betaT", "betaTN")[1:3], pch=".", Clim=c(0,6), Blim=c(0,1), xScale=1e-6, ..., add=FALSE) {
  # To please R CMD check
  fit <- x;
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'fit':
  if (nbrOfChromosomes(fit) > 1) {
    return(plotTracksManyChromosomes(fit, tracks=tracks, pch=pch, Clim=Clim, Blim=Blim, xScale=xScale, ..., add=add));
  }

  # Argument 'tracks':
  tracks <- match.arg(tracks, several.ok=TRUE);
  tracks <- unique(tracks);

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf));

  # Extract the input data
  data <- fit$data;
  if (is.null(data)) {
    throw("Cannot plot segmentation results. No input data available.");
  }

  chromosomes <- getChromosomes(fit);
  chromosome <- chromosomes[1];
  x <- data$x;
  CT <- data$CT;
  betaT <- data$betaT;
  betaN <- data$betaN;
  betaTN <- data$betaTN;
  muN <- data$muN;
  nbrOfLoci <- length(x);

  # Extract the segmentation
  segs <- fit$output;


  if (chromosome != 0) {
    chrTag <- sprintf("Chr%02d", chromosome);
  } else {
    chrTag <- "";
  }

  if (xScale != 1) {
    x <- xScale * x;
  }

  if (!add) {
    subplots(length(tracks), ncol=1);
    par(mar=c(1,4,1,2)+1);
  }

  for (track in tracks) {
    if (track == "tcn") {
      plot(x, CT, pch=pch, col="gray", ylim=Clim, ylab="TCN");
      stext(side=3, pos=1, chrTag);
      drawLevels(fit, what="tcn", xScale=xScale);
    }
  
    if (track == "tcn,c1,c2") {
      plot(x, CT, pch=pch, col="gray", ylim=Clim, ylab="C1, C2, TCN");
      stext(side=3, pos=1, chrTag);
      drawLevels(fit, what="tcn", xScale=xScale);
      drawLevels(fit, what="c2", col="purple", xScale=xScale);
      drawLevels(fit, what="c1", col="blue", xScale=xScale);
    }
  
    col <- c("gray", "black")[(muN == 1/2) + 1];
    if (track == "betaN") {
      plot(x, betaN, pch=pch, col=col, ylim=Blim, ylab="BAF_N");
      stext(side=3, pos=1, chrTag);
    }
  
    if (track == "betaT") {
      plot(x, betaT, pch=pch, col=col, ylim=Blim, ylab="BAF_T");
      stext(side=3, pos=1, chrTag);
    }
  
    if (track == "betaTN") {
      plot(x, betaTN, pch=pch, col=col, ylim=Blim, ylab="BAF_TN");
      stext(side=3, pos=1, chrTag);
    }
  
    if (track == "rho") {
      isSnp <- (!is.na(betaTN) & !is.na(muN));
      isHet <- isSnp & (muN == 1/2);
      naValue <- as.double(NA);
      rho <- rep(naValue, length=nbrOfLoci);
      rho[isHet] <- 2*abs(betaTN[isHet]-1/2);
      plot(x, rho, pch=pch, col="gray", ylim=Blim, ylab="DH");
      stext(side=3, pos=1, chrTag);
      drawLevels(fit, what="dh", xScale=xScale);
    }
  } # for (track ...)
})


setMethodS3("plot", "PairedPSCBS", function(x, ...) {
  plotTracks(x, ...);
})


setMethodS3("drawLevels", "PairedPSCBS", function(fit, what=c("tcn", "dh", "c1", "c2"), xScale=1e-6, ...) {
  # Argument 'what':
  what <- match.arg(what);

  # Get segmentation results
  segs <- as.data.frame(fit);


  if (is.element(what, c("c1", "c2"))) {
    c1c2 <- extractC1C2(fit);
    c1c2 <- c1c2[,1:2,drop=FALSE];
    # Sanity check
    stopifnot(nrow(c1c2) == nrow(segs));

    # AD HOC
    if (what == "c1") {
      segs[,"dh.mean"] <- c1c2[,1];
    } else if (what == "c2") {
      segs[,"dh.mean"] <- c1c2[,2];
    }
    what <- "dh";
  }

  # Extract subset of segments
  fields <- c("loc.start", "loc.end", "mean");
  fields <- sprintf("%s.%s", what, fields);
  segs <- segs[,fields, drop=FALSE];
  segs <- unique(segs);

  # Reuse drawLevels() for the DNAcopy class
  colnames(segs) <- c("loc.start", "loc.end", "seg.mean");
  dummy <- list(output=segs);
  class(dummy) <- "DNAcopy";
  drawLevels(dummy, xScale=xScale, ...);
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



setMethodS3("arrowsC1C2", "PairedPSCBS", function(fit, length=0.05, ...) {
  xy <- extractMinorMajorCNs(fit);
  xy <- xy[,1:2,drop=FALSE];
  x <- xy[,1,drop=TRUE];
  y <- xy[,2,drop=TRUE];
  s <- seq(length=length(x)-1);
  arrows(x0=x[s],y=y[s], x1=x[s+1],y1=y[s+1], code=2, length=length, ...);
}) # arrowsC1C2()


setMethodS3("arrowsDeltaC1C2", "PairedPSCBS", function(fit, length=0.05, ...) {
  xy <- extractDeltaC1C2(fit);
  xy <- xy[,1:2,drop=FALSE];
  x <- xy[,1,drop=TRUE];
  y <- xy[,2,drop=TRUE];
  s <- seq(length=length(x)-1);
  arrows(x0=x[s],y=y[s], x1=x[s+1],y1=y[s+1], code=2, length=length, ...);
})




setMethodS3("tileChromosomes", "PairedPSCBS", function(fit, chrStarts=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'chrStarts':
  if (!is.null(chrStarts)) {
    chrStarts <- Arguments$getDoubles(chrStarts);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Tile chromosomes");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- fit$data;
  segs <- fit$output;

  # Identify all chromosome
  chromosomes <- getChromosomes(fit);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Additional chromosome annotations
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(chrStarts)) {
    xRange <- matrix(0, nrow=length(chromosomes), ncol=2);
    for (kk in seq(along=chromosomes)) {
      chromosome <- chromosomes[kk];
      idxs <- which(data$chromosome == chromosome);
      x <- data$x[idxs];
      r <- range(x, na.rm=TRUE);
      r <- r / 1e6;
      r[1] <- floor(r[1]);
      r[2] <- ceiling(r[2]);
      r <- 1e6 * r;
      xRange[kk,] <- r;
    } # for (kk ...)

    chrLength <- xRange[,2];
    chrStarts <- c(0, cumsum(chrLength)[-length(chrLength)]);
    chrEnds <- chrStarts + chrLength;

    # Not needed anymore
    rm(x, idxs);
  } # if (is.null(chrStarts))

  verbose && cat(verbose, "Chromosome starts:");
  chromosomeStats <- cbind(start=chrStarts, end=chrEnds, length=chrEnds-chrStarts);
  verbose && print(chromosomeStats);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Offset...
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  segFields <- grep("(start|end)$", colnames(segs), value=TRUE);
  for (kk in seq(along=chromosomes)) {
    chromosome <- chromosomes[kk];
    chrTag <- sprintf("Chr%02d", chromosome);
    verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", 
                                         kk, chrTag, length(chromosomes)));

    # Get offset for this chromosome
    offset <- chrStarts[kk];
    verbose && cat(verbose, "Offset: ", offset);

    # Offset data
    idxs <- which(data$chromosome == chromosome);
    data$x[idxs] <- offset + data$x[idxs];

    # Offset segmentation
    idxs <- which(segs$chromosome == chromosome);
    segs[idxs,segFields] <- offset + segs[idxs,segFields];

    verbose && exit(verbose);
  } # for (kk ...)

  # Update results
  fit$data <- data;
  fit$output <- segs;
  fit$chromosomeStats <- chromosomeStats;

  verbose && exit(verbose);

  fit;
}) # tileChromosomes()



setMethodS3("plotTracksManyChromosomes", "PairedPSCBS", function(x, tracks=c("tcn", "rho", "tcn,c1,c2", "betaN", "betaT", "betaTN")[1:3], pch=".", Clim=c(0,6), Blim=c(0,1), xScale=1e-6, ..., subset=0.1, add=FALSE, verbose=FALSE) {
  # To please R CMD check
  fit <- x;
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'fit':

  # Argument 'tracks':
  tracks <- match.arg(tracks, several.ok=TRUE);
  tracks <- unique(tracks);

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf));

  # Argument 'subset':
  if (!is.null(subset)) {
    subset <- Arguments$getDouble(subset);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Tile chromosomes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit <- tileChromosomes(fit, verbose=verbose);
  verbose && str(verbose, fit);

  # Extract the input data
  data <- fit$data;
  if (is.null(data)) {
    throw("Cannot plot segmentation results. No input data available.");
  }

  # Subset of the loci?
  if (!is.null(subset)) {
    n <- nrow(data);
    keep <- sample(n, size=subset*n);
    data <- data[keep,];
  }

  attachLocally(data);
  x <- xScale * x;
  vs <- xScale * fit$chromosomeStats[,1:2];
  mids <- (vs[,1]+vs[,2])/2;

  chromosomes <- getChromosomes(fit);
  chrTags <- sprintf("Chr%02d", chromosomes);

  if (!add) {
    subplots(length(tracks), ncol=1);
    par(mar=c(1,4,1,2)+1);
  }

  for (track in tracks) {
    if (track == "tcn") {
      plot(x, CT, pch=pch, col="gray", ylim=Clim, ylab="TCN", axes=FALSE);
      mtext(text=chrTags, side=rep(c(1,3), length.out=length(chrTags)), at=mids, line=0.1, cex=0.7);
      abline(v=vs, lty=3);
      axis(side=2); box();
      drawLevels(fit, what="tcn", xScale=xScale);
    }
  
    if (track == "tcn,c1,c2") {
      plot(x, CT, pch=pch, col="gray", ylim=Clim, ylab="C1, C2, TCN", axes=FALSE);
      mtext(text=chrTags, side=rep(c(1,3), length.out=length(chrTags)), at=mids, line=0.1, cex=0.7);
      abline(v=vs, lty=3);
      axis(side=2); box();
      drawLevels(fit, what="tcn", xScale=xScale);
      drawLevels(fit, what="c2", col="purple", xScale=xScale);
      drawLevels(fit, what="c1", col="blue", xScale=xScale);
    }
  
    col <- c("gray", "black")[(muN == 1/2) + 1];
    if (track == "betaN") {
      plot(x, betaN, pch=pch, col=col, ylim=Blim, ylab="BAF_N", axes=FALSE);
      mtext(text=chrTags, side=rep(c(1,3), length.out=length(chrTags)), at=mids, line=0.1, cex=0.7);
      abline(v=vs, lty=3);
      axis(side=2); box();
    }
  
    if (track == "betaT") {
      plot(x, betaT, pch=pch, col=col, ylim=Blim, ylab="BAF_T", axes=FALSE);
      mtext(text=chrTags, side=rep(c(1,3), length.out=length(chrTags)), at=mids, line=0.1, cex=0.7);
      abline(v=vs, lty=3);
      axis(side=2); box();
    }
  
    if (track == "betaTN") {
      plot(x, betaTN, pch=pch, col=col, ylim=Blim, ylab="BAF_TN", axes=FALSE);
      mtext(text=chrTags, side=rep(c(1,3), length.out=length(chrTags)), at=mids, line=0.1, cex=0.7);
      abline(v=vs, lty=3);
      axis(side=2); box();
    }
  
    if (track == "rho") {
      isSnp <- (!is.na(betaTN) & !is.na(muN));
      isHet <- isSnp & (muN == 1/2);
      naValue <- as.double(NA);
      rho <- rep(naValue, length=nbrOfLoci);
      rho[isHet] <- 2*abs(betaTN[isHet]-1/2);
      plot(x, rho, pch=pch, col="gray", ylim=Blim, ylab="DH", axes=FALSE);
      mtext(text=chrTags, side=rep(c(1,3), length.out=length(chrTags)), at=mids, line=0.1, cex=0.7);
      abline(v=vs, lty=3);
      axis(side=2); box();
      drawLevels(fit, what="dh", xScale=xScale);
    }
  } # for (track ...)
}) # plotTracksManyChromosomes()




############################################################################
# HISTORY:
# 2010-10-18
# o Now plotTracks() can plot whole-genome data.
# o Now plotTracks() utilized plotTracksManyChromosomes() if there is
#   more than one chromosome.
# o Added internal plotTracksManyChromosomes().
# o Added internal tileChromosomes().
# 2010-10-03
# o Now the default is that plotTracks() for PairedPSCBS generated three
#   panels: (1) TCN, (2) DH, and (3) C1+C2+TCN.
# o Added plotTracks() to be more explicit than just plot().
# o Added argument 'xScale' to plot() for PairedPSCBS.
# o Now plot() for PairedPSCBS adds a chromosome tag.
# 2010-09-21
# o Added argument 'what' to plot() for PairedPSCBS.
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
