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

  data <- segs[,c("tcn.mean", "dh.mean", "tcn.num.mark"), drop=FALSE];
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




setMethodS3("plotC1C2", "PairedPSCBS", function(fit, ..., xlab=expression(C[1]), ylab=expression(C[2]), Clim=c(0,4)) {
  plot(NA, xlim=Clim, ylim=Clim, xlab=xlab, ylab=ylab);
  abline(a=0, b=1, lty=3);
  pointsC1C2(fit, ...);
})


setMethodS3("pointsC1C2", "PairedPSCBS", function(fit, cex=NULL, ...) {
  data <- extractMinorMajorCNs(fit);
  X <- data[,1:2,drop=FALSE];
  n <- data[,3,drop=TRUE];
  w <- n / sum(n, na.rm=TRUE);

  if (is.null(cex)) {
    cex <- w;
    cex <- cex / mean(cex, na.rm=TRUE);
    cex <- cex + 1/2;
  }

  points(X, cex=cex, ...);
})




############################################################################
# HISTORY:
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
