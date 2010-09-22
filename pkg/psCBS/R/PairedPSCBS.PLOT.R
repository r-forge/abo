setMethodS3("plot", "PairedPSCBS", function(x, what=c("tcn", "betaN", "betaT", "betaTN", "rho"), pch=".", Clim=c(0,6), Blim=c(0,1), ..., add=FALSE) {
  # To please R CMD check
  fit <- x;
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'what':
  what <- match.arg(what, several.ok=TRUE);
  what <- unique(what);

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
  nbrOfLoci <- length(x);

  # Extract the segmentation
  segs <- fit$output;

  if (!add) {
    subplots(length(what), ncol=1);
    par(mar=c(1,4,1,2)+1);
  }

  if (is.element("tcn", what)) {
    plot(x, CT, pch=pch, col="gray", ylim=Clim);
    drawLevels(fit, what="tcn");
  }

  col <- c("gray", "black")[(muN == 1/2) + 1];
  if (is.element("betaN", what)) {
    plot(x, betaN, pch=pch, col=col, ylim=Blim);
  }

  if (is.element("betaT", what)) {
    plot(x, betaT, pch=pch, col=col, ylim=Blim);
  }

  if (is.element("betaTN", what)) {
    plot(x, betaTN, pch=pch, col=col, ylim=Blim);
  }

  if (is.element("rho", what)) {
    isSnp <- (!is.na(betaTN) & !is.na(muN));
    isHet <- isSnp & (muN == 1/2);
    naValue <- as.double(NA);
    rho <- rep(naValue, length=nbrOfLoci);
    rho[isHet] <- 2*abs(betaTN[isHet]-1/2);
    plot(x, rho, pch=pch, col=col, ylim=Blim);
    drawLevels(fit, what="dh");
  }
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
