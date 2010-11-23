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
#   \item{x}{A result object returned by @see "segmentByPairedPSCBS".}
#   \item{tracks}{A @character @vector specifying what types of tracks to plot.}
#   \item{calls}{A @character @vector of regular expression identifying
#     call labels to be highlighted in the panels.}
#   \item{pch}{The type of points to use.}
#   \item{xlim}{(Optional) The genomic range to plot.}
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
setMethodS3("plotTracks", "PairedPSCBS", function(x, tracks=c("tcn", "dh", "tcn,c1,c2", "tcn,c1", "tcn,c2", "c1,c2", "betaN", "betaT", "betaTN")[1:3], calls=".*", scatter=TRUE, pch=".", quantiles=c(0.025,0.975), cex=1, grid=FALSE, xlim=NULL, Clim=c(0,6), Blim=c(0,1), xScale=1e-6, ..., add=FALSE) {
  # To please R CMD check
  fit <- x;
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'fit':
  if (nbrOfChromosomes(fit) > 1) {
    return(plotTracksManyChromosomes(fit, tracks=tracks, scatter=scatter, pch=pch, Clim=Clim, Blim=Blim, xScale=xScale, ..., add=add));
  }

  # Argument 'tracks':
  knownTracks <- c("tcn", "dh", "tcn,c1,c2", "tcn,c1", "tcn,c2", "c1,c2", "betaN", "betaT", "betaTN");
  tracks <- match.arg(tracks, choices=knownTracks, several.ok=TRUE);
  tracks <- unique(tracks);

  # Argument 'calls':
  if (!is.null(calls)) {
    calls <- sapply(calls, FUN=Arguments$getRegularExpression);
  }

  # Argument 'xlim':
  if (!is.null(xlim)) {
    xlim <- Arguments$getNumerics(xlim, length=c(2,2));
  }

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
  segs <- as.data.frame(fit);

  # Identify available calls
  if (!is.null(calls)) {
    verbose && enter(verbose, "Identifying calls");

    pattern <- "[.]call$";
    callColumns <- grep(pattern, colnames(segs), value=TRUE);
    if (length(callColumns) > 0) {
      keep <- sapply(calls, FUN=function(pattern) {
        (regexpr(pattern, callColumns) != -1);
      });
      if (is.matrix(keep)) {
        keep <- apply(keep, MARGIN=1, FUN=any);
      }
      callColumns <- callColumns[keep];
      callLabels <- gsub(pattern, "", callColumns);
      callLabels <- toupper(callLabels);
    }
    verbose && cat(verbose, "Call columns:");
    verbose && print(verbose, callColumns);

    verbose && exit(verbose);
  } else {
    callColumns <- NULL;
  }

  if (chromosome != 0) {
    chrTag <- sprintf("Chr%02d", chromosome);
  } else {
    chrTag <- "";
  }

  if (xScale != 1) {
    x <- xScale * x;
    if (!is.null(xlim)) {
      xlim <- xScale * xlim;
    }
  }

  if (!add) {
    subplots(length(tracks), ncol=1);
    par(mar=c(1,4,1,2)+1);
  }

  # Color loci by genotype
  col <- c("gray", "black")[(muN == 1/2) + 1];

  for (track in tracks) {
    pchT <- if (scatter) { pch } else { NA };

    if (track == "tcn") {
      plot(x, CT, pch=pchT, cex=cex, col="black", xlim=xlim, ylim=Clim, ylab="TCN");
      stext(side=3, pos=1, chrTag);
      if (grid) {
        abline(h=seq(from=0, to=Clim[2], by=2), lty=3, col="gray");
        abline(h=0, lty=1, col="black");
      }
      drawConfidenceBands(fit, what="tcn", quantiles=quantiles, col="purple", xScale=xScale);
      drawLevels(fit, what="tcn", col="purple", xScale=xScale);
    }
  
    if (is.element(track, c("tcn,c1,c2", "tcn,c1", "tcn,c2", "c1,c2"))) {
      subtracks <- strsplit(track, split=",", fixed=TRUE)[[1]];
      ylab <- paste(toupper(subtracks), collapse=", ");
      plot(x, CT, pch=pchT, cex=cex, col="black", xlim=xlim, ylim=Clim, ylab=ylab);
      stext(side=3, pos=1, chrTag);
      if (grid) {
        abline(h=seq(from=0, to=Clim[2], by=2), lty=3, col="gray");
        abline(h=0, lty=1, col="black");
      }
      if (is.element("tcn", subtracks)) {
        drawConfidenceBands(fit, what="tcn", quantiles=quantiles, col="purple", xScale=xScale);
      }
      if (is.element("c2", subtracks)) {
        drawConfidenceBands(fit, what="c2", quantiles=quantiles, col="red", xScale=xScale);
      }
      if (is.element("c1", subtracks)) {
        drawConfidenceBands(fit, what="c1", quantiles=quantiles, col="blue", xScale=xScale);
      }
      if (is.element("tcn", subtracks)) {
        drawLevels(fit, what="tcn", col="purple", xScale=xScale);
      }
      if (is.element("c2", subtracks)) {
        drawLevels(fit, what="c2", col="red", xScale=xScale);
      }
      if (is.element("c1", subtracks)) {
        drawLevels(fit, what="c1", col="blue", xScale=xScale);
      }
    }
  
    if (track == "betaN") {
      plot(x, betaN, pch=pchT, cex=cex, col=col, xlim=xlim, ylim=Blim, ylab=expression(BAF[N]));
      stext(side=3, pos=1, chrTag);
    }
  
    if (track == "betaT") {
      plot(x, betaT, pch=pchT, cex=cex, col=col, xlim=xlim, ylim=Blim, ylab=expression(BAF[T]));
      stext(side=3, pos=1, chrTag);
    }
  
    if (track == "betaTN") {
      plot(x, betaTN, pch=pchT, cex=cex, col=col, xlim=xlim, ylim=Blim, ylab=expression(BAF[T]^"*"));
      stext(side=3, pos=1, chrTag);
    }
  
    if (track == "dh") {
      isSnp <- (!is.na(betaTN) & !is.na(muN));
      isHet <- isSnp & (muN == 1/2);
      naValue <- as.double(NA);
      rho <- rep(naValue, length=nbrOfLoci);
      rho[isHet] <- 2*abs(betaTN[isHet]-1/2);
      plot(x, rho, pch=pchT, cex=cex, col=col[isHet], xlim=xlim, ylim=Blim, ylab="DH");
      stext(side=3, pos=1, chrTag);
      drawConfidenceBands(fit, what="dh", quantiles=quantiles, col="orange", xScale=xScale);
      drawLevels(fit, what="dh", col="orange", xScale=xScale);
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # For each panel of tracks, annotate will calls?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (length(callColumns) > 0) {
      for (cc in seq(along=callColumns)) {
        callColumn <- callColumns[cc];
        callLabel <- callLabels[cc];

        segsT <- segs[,c("dh.loc.start", "dh.loc.end", callColumn)];
        segsT <- segsT[segsT[,callColumn],1:2,drop=FALSE];
        segsT <- xScale * segsT;

        side <- 2*((cc+1) %% 2) + 1;
        # For each segment called...
        for (ss in seq(length=nrow(segsT))) {
          x0 <- segsT[ss,1,drop=TRUE];
          x1 <- segsT[ss,2,drop=TRUE];
          abline(v=c(x0,x1), lty=3, col="gray");
          xMid <- (x0+x1)/2;
          mtext(side=side, at=xMid, line=-1, cex=0.7, col="#666666", callLabel);
        } # for (ss in ...)
      } # for (cc in ...)
    } # if (length(callColumns) > 0)

  } # for (track ...)

  invisible();  
}) # plotTracks()


setMethodS3("plot", "PairedPSCBS", function(x, ...) {
  plotTracks(x, ...);
})


setMethodS3("drawLevels", "PairedPSCBS", function(fit, what=c("tcn", "dh", "c1", "c2"), xScale=1e-6, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'what':
  what <- match.arg(what);

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf));


  # Get segmentation results
  segs <- as.data.frame(fit);

  # Extract subset of segments
  fields <- c("loc.start", "loc.end");
  fields <- sprintf("%s.%s", ifelse(what == "tcn", what, "dh"), fields);
  fields <- c(fields, sprintf("%s.mean", what));
  segsT <- segs[,fields, drop=FALSE];
  segsT <- unique(segsT);

  # Reuse drawLevels() for the DNAcopy class
  colnames(segsT) <- c("loc.start", "loc.end", "seg.mean");
  dummy <- list(output=segsT);
  class(dummy) <- "DNAcopy";

  drawLevels(dummy, xScale=xScale, ...);
})


setMethodS3("drawConfidenceBands", "PairedPSCBS", function(fit, what=c("tcn", "dh", "c1", "c2"), quantiles=c(0.025,0.975), col=col, alpha=0.4, xScale=1e-6, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'what':
  what <- match.arg(what);

  # Argument 'quantiles':
  if (!is.null(quantiles)) {
    quantiles <- Arguments$getNumerics(quantiles, range=c(0,1), length=c(2,2));
  }

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf));


  # Nothing todo?
  if (is.null(quantiles)) {
    return();
  }


  # Get segmentation results
  segs <- as.data.frame(fit);

  # Extract subset of segments
  fields <- c("loc.start", "loc.end");
  fields <- sprintf("%s.%s", ifelse(what == "tcn", what, "dh"), fields);

  tags <- sprintf("%g%%", 100*quantiles);
  qFields <- sprintf("%s_%s", what, tags);

  # Nothing todo?
  if (!all(is.element(qFields, colnames(segs)))) {
    return();
  }

  fields <- c(fields, qFields);

  segsT <- segs[,fields, drop=FALSE];
  segsT <- unique(segsT);

  # Rescale x-axis
  segsT[,1:2] <- xScale * segsT[,1:2];

  colQ <- col2rgb(col, alpha=TRUE);
  colQ["alpha",] <- alpha*colQ["alpha",];
  colQ <- rgb(red=colQ["red",], green=colQ["green",], blue=colQ["blue",], alpha=colQ["alpha",], maxColorValue=255);

  for (kk in seq(length=nrow(segsT))) {
    rect(xleft=segsT[kk,1], xright=segsT[kk,2], ybottom=segsT[kk,3], ytop=segsT[kk,4], col=colQ, border=FALSE);
  }
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



setMethodS3("plotTracksManyChromosomes", "PairedPSCBS", function(x, tracks=c("tcn", "dh", "tcn,c1,c2", "betaN", "betaT", "betaTN")[1:3], calls=".*", quantiles=c(0.025,0.975), pch=".", Clim=c(0,6), Blim=c(0,1), xScale=1e-6, ..., subset=0.1, add=FALSE, onBegin=NULL, onEnd=NULL, verbose=FALSE) {
  # To please R CMD check
  fit <- x;
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'fit':

  # Argument 'tracks':
  tracks <- match.arg(tracks, several.ok=TRUE);
  tracks <- unique(tracks);

  # Argument 'calls':
  if (!is.null(calls)) {
    calls <- sapply(calls, FUN=Arguments$getRegularExpression);
  }

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

  # Extract the segmentation
  segs <- as.data.frame(fit);

  # Identify available calls
  if (!is.null(calls)) {
    verbose && enter(verbose, "Identifying calls");

    pattern <- "[.]call$";
    callColumns <- grep(pattern, colnames(segs), value=TRUE);
    if (length(callColumns) > 0) {
      keep <- sapply(calls, FUN=function(pattern) {
        (regexpr(pattern, callColumns) != -1);
      });
      if (is.matrix(keep)) {
        keep <- apply(keep, MARGIN=1, FUN=any);
      }
      callColumns <- callColumns[keep];
      callLabels <- gsub(pattern, "", callColumns);
      callLabels <- toupper(callLabels);
    }
    verbose && cat(verbose, "Call columns:");
    verbose && print(verbose, callColumns);

    verbose && exit(verbose);
  } else {
    callColumns <- NULL;
  }

  # Subset of the loci?
  if (!is.null(subset)) {
    n <- nrow(data);
    keep <- sample(n, size=subset*n);
    data <- data[keep,];
  }

  # To please R CMD check
  CT <- muN <- betaT <- betaN <- betaTN <- NULL;
  rm(CT, muN, betaT, betaN, betaTN);
  attachLocally(data);
  x <- xScale * x;
  vs <- xScale * fit$chromosomeStats[,1:2,drop=FALSE];
  mids <- (vs[,1]+vs[,2])/2;

  nbrOfLoci <- length(x);
  chromosomes <- getChromosomes(fit);
  chrTags <- sprintf("Chr%02d", chromosomes);

  if (!add) {
    subplots(length(tracks), ncol=1);
    par(mar=c(1,4,1,2)+1);
  }

  gh <- fit;
  gh$xScale <- xScale;

  xlim <- range(x, na.rm=TRUE);
  xlab <- "Genomic position";

  for (track in tracks) {
    if (track == "tcn") {
      plot(NA, xlim=xlim, ylim=Clim, xlab=xlab, ylab="TCN", axes=FALSE);
      if (!is.null(onBegin)) onBegin(gh=gh);
      points(x, CT, pch=pch, col="black");
      mtext(text=chrTags, side=rep(c(1,3), length.out=length(chrTags)), at=mids, line=0.1, cex=0.7);
      abline(v=vs, lty=1, lwd=2);
      axis(side=2); box();
      drawConfidenceBands(fit, what="tcn", quantiles=quantiles, col="purple", xScale=xScale);
      drawLevels(fit, what="tcn", col="purple", xScale=xScale);
      if (!is.null(onEnd)) onEnd(gh=gh);
    }
  
    if (track == "tcn,c1,c2") {
      plot(NA, xlim=xlim, ylim=Clim, xlab=xlab, ylab="C1, C2, TCN", axes=FALSE);
      if (!is.null(onBegin)) onBegin(gh=gh);
      points(x, CT, pch=pch, col="black");
      mtext(text=chrTags, side=rep(c(1,3), length.out=length(chrTags)), at=mids, line=0.1, cex=0.7);
      abline(v=vs, lty=1, lwd=2);
      axis(side=2); box();
      drawConfidenceBands(fit, what="tcn", quantiles=quantiles, col="purple", xScale=xScale);
      drawConfidenceBands(fit, what="c2", quantiles=quantiles, col="red", xScale=xScale);
      drawConfidenceBands(fit, what="c1", quantiles=quantiles, col="blue", xScale=xScale);
      drawLevels(fit, what="tcn", col="purple", xScale=xScale);
      drawLevels(fit, what="c2", col="red", xScale=xScale);
      drawLevels(fit, what="c1", col="blue", xScale=xScale);
      if (!is.null(onEnd)) onEnd(gh=gh);
    }
  
    col <- c("gray", "black")[(muN == 1/2) + 1];
    if (track == "betaN") {
      plot(NA, xlim=xlim, ylim=Blim, xlab=xlab, ylab="BAF_N", axes=FALSE);
      if (!is.null(onBegin)) onBegin(gh=gh);
      points(x, betaN, pch=pch, col=col);
      mtext(text=chrTags, side=rep(c(1,3), length.out=length(chrTags)), at=mids, line=0.1, cex=0.7);
      abline(v=vs, lty=1, lwd=2);
      axis(side=2); box();
      if (!is.null(onEnd)) onEnd(gh=gh);
    }
  
    if (track == "betaT") {
      plot(NA, xlim=xlim, ylim=Blim, xlab=xlab, ylab="BAF_T", axes=FALSE);
      if (!is.null(onBegin)) onBegin(gh=gh);
      points(x, betaT, pch=pch, col=col);
      mtext(text=chrTags, side=rep(c(1,3), length.out=length(chrTags)), at=mids, line=0.1, cex=0.7);
      abline(v=vs, lty=1, lwd=2);
      axis(side=2); box();
      if (!is.null(onEnd)) onEnd(gh=gh);
    }
  
    if (track == "betaTN") {
      plot(NA, xlim=xlim, ylim=Blim, xlab=xlab, ylab="BAF_TN", axes=FALSE);
      if (!is.null(onBegin)) onBegin(gh=gh);
      points(x, betaTN, pch=pch, col=col);
      mtext(text=chrTags, side=rep(c(1,3), length.out=length(chrTags)), at=mids, line=0.1, cex=0.7);
      abline(v=vs, lty=1, lwd=2);
      axis(side=2); box();
      if (!is.null(onEnd)) onEnd(gh=gh);
    }
  
    if (track == "dh") {
      isSnp <- (!is.na(betaTN) & !is.na(muN));
      isHet <- isSnp & (muN == 1/2);
      naValue <- as.double(NA);
      rho <- rep(naValue, length=nbrOfLoci);
      rho[isHet] <- 2*abs(betaTN[isHet]-1/2);
      plot(NA, xlim=xlim, ylim=Blim, xlab=xlab, ylab="DH", axes=FALSE);
      if (!is.null(onBegin)) onBegin(gh=gh);
      points(x, rho, pch=pch, col="black");
      mtext(text=chrTags, side=rep(c(1,3), length.out=length(chrTags)), at=mids, line=0.1, cex=0.7);
      abline(v=vs, lty=1, lwd=2);
      axis(side=2); box();
      drawConfidenceBands(fit, what="dh", quantiles=quantiles, col="orange", xScale=xScale);
      drawLevels(fit, what="dh", col="orange", xScale=xScale);
      if (!is.null(onEnd)) onEnd(gh=gh);
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # For each panel of tracks, annotate will calls?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (length(callColumns) > 0) {
      for (cc in seq(along=callColumns)) {
        callColumn <- callColumns[cc];
        callLabel <- callLabels[cc];

        segsT <- segs[,c("dh.loc.start", "dh.loc.end", callColumn)];
        segsT <- segsT[segsT[,callColumn],1:2,drop=FALSE];
        segsT <- xScale * segsT;

        side <- 2*((cc+1) %% 2) + 1;
        # For each segment called...
        for (ss in seq(length=nrow(segsT))) {
          x0 <- segsT[ss,1,drop=TRUE];
          x1 <- segsT[ss,2,drop=TRUE];
          abline(v=c(x0,x1), lty=3, col="gray");
          xMid <- (x0+x1)/2;
          mtext(side=side, at=xMid, line=-1, cex=0.7, col="#666666", callLabel);
        } # for (ss in ...)
      } # for (cc in ...)
    } # if (length(callColumns) > 0)
  } # for (track ...)

  invisible(gh);
}) # plotTracksManyChromosomes()




############################################################################
# HISTORY:
# 2010-11-22
# o ROBUSTNESS: Now drawConfidenceBands() of PairedPSCBS silently does
#   nothing if the requested bootstrap quantiles are available.
# o Added argument 'calls' to plotTracks() and plotTracksManyChromosomes()
#   for highlighing called regions.
# 2010-11-21
# o Now plotTracks() supports tracks "tcn,c1", "tcn,c2" and "c1,c2" too.
# o Added argument 'xlim' to plotTracks() making it possible to zoom in.
# o Now plotTracks() and plotTracksManyChromosomes() draws confidence
#   bands, iff argument quantiles is given.
# o Added drawConfidenceBands() for PairedPSCBS.
# 2010-11-09
# o Added argument 'cex=1' to plotTracks().
# o BUG FIX: It was not possible to plot BAF tracks with plotTracks().
# 2010-10-20
# o Added arguments 'onBegin' and 'onEnd' to plotTracksManyChromosomes().
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
