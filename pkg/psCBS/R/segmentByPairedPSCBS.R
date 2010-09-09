###########################################################################/**
# @RdocDefault segmentByPairedPSCBS
# @alias segmentByPairedPSCBS
#
# @title "Segment copy numbers using the paired PSCBS method"
#
# \description{
#  @get "title".
#  This is a low-level segmentation method.
#  It is intended to be applied to one sample and one chromosome at the time.
# }
#
# @synopsis
#
# \arguments{
#   \item{CT}{A @numeric @vector of J tumor total copy numbers to be segmented.}
#   \item{betaT}{A @numeric @vector of J tumor allele B fractions (BAFs).}
#   \item{betaN}{A @numeric @vector of J matched normal BAFs.}
#   \item{muN}{An optional @numeric @vector of J genotype calls.
#       If not given, they are estimated from the normal BAFs.}
#   \item{x}{Optional @numeric @vector of J genomic locations.
#            If @NULL, index locations \code{1:J} are used.}
#   \item{...}{Not used.}
#   \item{flavor}{A @character specifying what type of segmentation and 
#     calling algorithm to be used.}
#   \item{seed}{An (optional) @integer specifying the random seed to be 
#     set before calling the segmentation method.  The random seed is
#     set to its original state when exiting.  If @NULL, it is not set.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the fit object.
# }
# 
# \details{
#   Internally @see "segmentByCBS" is used for segmentation.
#   This segmentation method does \emph{not} support weights.
#
#   The "DNAcopy::segment" implementation of CBS uses approximation
#   through random sampling for some estimates.  Because of this,
#   repeated calls using the same signals may result in slightly 
#   different results, unless the random seed is fixed.  
# }
#
# \section{Missing and non-finite values}{
#   The total copy number signals as well as any optional positions
#   must not contain missing values, i.e. @NAs or @NaNs.
#   If there are any, an informative error is thrown.
#   Allele B fractions may contain missing values, because such are
#   interpreted as representing non-polymorphic loci.
#
#   None of the input signals can have infinite values, i.e. -@Inf or @Inf.
#   If so, an informative error is thrown.
# }
#
# @examples "../incl/segmentByPairedPSCBS.Rex"
#
# @author
#
# \seealso{
#   Internally @see "segmentByCBS" is used for segmentation.
# }
#
# @keyword IO
#*/########################################################################### 
setMethodS3("segmentByPairedPSCBS", "default", function(CT, betaT, betaN, muN=NULL, x=NULL, ..., flavor=c("dh|tcn", "tcn|dh", "tcn&dh"), seed=NULL, verbose=FALSE) {
  require("R.utils") || throw("Package not loaded: R.utils");
  require("aroma.light") || throw("Package not loaded: aroma.light");
  ver <- packageDescription("aroma.light")$Version;
  stopifnot(compareVersion(ver, "1.17.2") >= 0);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'CT':
  disallow <- c("NA", "NaN", "Inf");
  CT <- Arguments$getDoubles(CT, disallow=disallow);
  nbrOfLoci <- length(CT);
  length2 <- rep(nbrOfLoci, 2);

  # Argument 'betaT':
  betaT <- Arguments$getDoubles(betaT, length=length2, disallow="Inf");

  # Argument 'betaN':
  betaN <- Arguments$getDoubles(betaN, length=length2, disallow="Inf");

  # Argument 'muN':
  if (!is.null(muN)) {
    muN <- Arguments$getDoubles(muN, length=length2, range=c(0,1), disallow="Inf");
  }

  # Argument 'x':
  if (!is.null(x)) {
    x <- Arguments$getDoubles(x, length=length2, disallow=disallow);
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);
  if (flavor != "dh|tcn") {
    throw("Segmentation flavor not supported. Currently only \"dh|tcn\" is implemented: ", flavor);
  }

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


  verbose && enter(verbose, "Segmenting paired tumor-normal signals using PSCBS");
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);


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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call genotypes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # If muN is missing, call genotypes from betaN
  if (is.null(muN)) {
    verbose && enter(verbose, "Calling genotypes from normal allele B fractions");
    verbose && str(verbose, betaN);
    muN <- aroma.light::callNaiveGenotypes(betaN);
    verbose && cat(verbose, "Called genotypes:");
    verbose && str(verbose, muN);
    verbose && print(verbose, table(muN));
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize betaT using betaN (TumorBoost normalization)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Normalizing betaT using betaN (TumorBoost)");
  betaTN <- aroma.light::normalizeTumorBoost(betaT=betaT, betaN=betaN, muN=muN, preserveScale=TRUE);
  verbose && cat(verbose, "Normalized BAFs:");
  verbose && str(verbose, betaTN);
  verbose && exit(verbose);

  # Assert that no missing values where introduced
  keep <- (is.finite(betaT) & is.finite(betaN) & is.finite(muN));
  if (any(is.na(betaTN[keep]))) {
    throw("Internal error: normalizeTumorBoost() introduced missing values.");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Record input data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Physical positions of loci
  if (is.null(x)) {
    x <- seq(length=nbrOfLoci);
  }

  # SNPs are identifies as those loci that have non-missing 'betaTN' & 'muN'
  isSnp <- (!is.na(betaTN) & !is.na(muN));
  nbrOfSnps <- sum(isSnp);
  verbose && cat(verbose, "Number of SNPs: ", nbrOfSnps);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate decrease-of-heterozygosity signals (DHs)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating DHs");
  # DH is by definition only defined for heterozygous SNPs.  
  # For simplicity, we set it to be NA for non-heterozygous loci.
  isHet <- isSnp & (muN == 1/2);
  verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                     sum(isHet), 100*sum(isHet)/nbrOfSnps);
  naValue <- as.double(NA);
  rho <- rep(naValue, length=nbrOfLoci);
  rho[isHet] <- 2*abs(betaTN[isHet]-1/2);
  verbose && cat(verbose, "Normalized DHs:");
  verbose && str(verbose, rho);
  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1a. Identification of change points in total copy numbers
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # TO DO: Deal with negative values
  h <- function(x) { 
    y <- x;
    y[x >= 0] <- sqrt(x[x >= 0]);
    y;
  }
  hinv <- function(y) {
    x <- y;
    x[y >= 0] <- y[y >= 0]^2;
    x;
  }

  verbose && enter(verbose, "Identification of change points by total copy numbers");
  if (!is.null(h)) {
    verbose && enter(verbose, "Transforming signals");
    verbose && cat(verbose, "Input signals:");
    verbose && str(verbose, CT);

    # Segment total CN ratios (CT) using CBS.
    # We actually segment sqrt(CT).
    # All loci are utilized, i.e. all SNPs and all non-polymorphic markers
    yT <- h(CT);

    verbose && cat(verbose, "Signals to be segmented:");
    verbose && str(verbose, yT);
    verbose && exit(verbose);
  } else {
    yT <- CT;
  } # if (!is.null(h))

  fit <- segmentByCBS(yT, x=x, verbose=verbose);
  verbose && str(verbose, fit);
  tcnSegments <- fit$output;
  rm(yT, fit);

  if (!is.null(hinv)) {
    verbose && enter(verbose, "Backtransforming segmented mean levels");
    tcnSegments[,"seg.mean"] <- hinv(tcnSegments[,"seg.mean"]);
    verbose && exit(verbose);
  }

  # Drop dummy columns
  keep <- setdiff(colnames(tcnSegments), c("ID", "chrom"));
  tcnSegments <- tcnSegments[,keep,drop=FALSE];

  # Tag fields by TCN
  names <- names(tcnSegments);
  names <- gsub("seg.mean", "mean", names, fixed=TRUE);
  names <- sprintf("tcn.%s", names);
  names(tcnSegments) <- names;
  rm(names);
  verbose && print(verbose, tcnSegments);

  nbrOfSegs <- nrow(tcnSegments);
  verbose && cat(verbose, "Number of TCN segments: ", nbrOfSegs);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1b. Identification of additional change points using DH
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For each segment independently, segment decrease of heterozygousity (DH)
  # using CBS. By definition, only heterozygous SNPs are used.

  # For each TCN segment...
  segs <- vector("list", length=nbrOfSegs);
  for (kk in seq(length=nbrOfSegs)) {
    # Extract the region
    xStart <- tcnSegments[kk,"tcn.loc.start"];
    xEnd <- tcnSegments[kk,"tcn.loc.end"];
    gammaT <- tcnSegments[kk,"tcn.mean"];

    regionTag <- sprintf("[%d,%d]", xStart, xEnd);
    verbose && enter(verbose, sprintf("Total CN segment #%d (%s) of %d", kk, regionTag, nbrOfSegs));

    # Identify subset of loci
    keep <- (xStart <= x & x <= xEnd);
    xKK <- x[keep];
    nbrOfLociKK <- length(xKK);
    rhoKK <- rho[keep];
    isSnpKK <- isSnp[keep];
    nbrOfSnpsKK <- sum(isSnpKK);
    isHetKK <- isHet[keep];
    verbose && cat(verbose, "Number of loci in region: ", nbrOfLociKK);
    verbose && cat(verbose, "Number of SNPs in region: ", nbrOfSnpsKK);
    rm(keep);

    # Identify heterozygous SNPs
    keep <- isHetKK;
    xKKHet <- xKK[keep];
    rhoKKHet <- rhoKK[keep];
    nbrOfHetsKK <- length(xKKHet);
    verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                  nbrOfHetsKK, 100*nbrOfHetsKK/nbrOfSnpsKK);
    rm(keep);

    verbose && enter(verbose, "Segmenting DH signals");
    fit <- segmentByCBS(rhoKKHet, x=xKKHet, verbose=less(verbose, 10));
    verbose && str(verbose, fit);
    dhSegments <- fit$output;
    verbose && exit(verbose);
    rm(rhoKKHet, xKKHet, fit);

    # Drop dummy columns
    keep <- setdiff(colnames(dhSegments), c("ID", "chrom"));
    dhSegments <- dhSegments[,keep,drop=FALSE];

    # Tag fields by DH
    names <- names(dhSegments);
    names <- gsub("seg.mean", "mean", names, fixed=TRUE);
    names <- sprintf("dh.%s", names);
    names(dhSegments) <- names;
    rm(names);

    # Special case: If there where not enough data to segment DH...
    if (nrow(dhSegments) == 0) {
      dhSegments <- dhSegments[as.integer(NA),];
    }

    # Insert that last column
    at <- which(colnames(dhSegments) == "dh.num.mark");

    verbose && cat(verbose, "DH segmentation table:");
    verbose && print(verbose, dhSegments);


    # Expand the TCN segmentation result data frame
    if (nrow(dhSegments) > 0) {
      tcnSegmentsKK <- rep(list(tcnSegments[kk,]), times=nrow(dhSegments));
      tcnSegmentsKK <- Reduce(rbind, tcnSegmentsKK);
    } else {
      tcnSegmentsKK <- tcnSegments[integer(0),,drop=FALSE];
    }
    # Append information on number of SNPs and hets in CN region
    tcnSegmentsKK <- cbind(
      tcnSegmentsKK, 
      tcn.num.snps=nbrOfSnpsKK,
      tcn.num.hets=nbrOfHetsKK
    );
    verbose && cat(verbose, "Total CN segmentation table (expanded):");
    verbose && print(verbose, tcnSegmentsKK);

    # Combine TCN and DH segmentation results
    tcndhSegments <- cbind(
      tcn.id=rep(kk, times=nrow(dhSegments)),
      dh.id=seq(length=nrow(dhSegments)),
      tcnSegmentsKK,
      dhSegments
    );

    segs[[kk]] <- tcndhSegments;

    verbose && exit(verbose);    
  } # for (kk ...)

  segs <- Reduce(rbind, segs);
  rownames(segs) <- NULL;

  # Reorder columns
#  fields <- c("loc.start", "loc.end", "num.mark", "num.het.mark", "dh.mean", "tcn.mean");
#  segs <- segs[,fields,drop=FALSE];
  verbose && print(verbose, segs);

  nbrOfSegs <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegs);

  verbose && exit(verbose);

  # Create result object
  data <- list(CT=CT, betaT=betaT, betaTN=betaTN, betaN=betaN, muN=muN, x=x);

  fit <- list(
    data = data,
    output = segs
  );
  class(fit) <- c("PairedPSCBS", "PSCBS");

  # Return segments found
  fit;
}) # segmentByPairedPSCBS()




############################################################################
# HISTORY:
# 2010-09-08
# o Now segmentByPairedPSCBS() also returns the TumorBoost normalized data.
#   This also means that plot() for PairedPSCBS no longer has to 
#   recalculate them.
# 2010-09-04
# o Added drawLevels() for PairedPSCBS.
# o Added as.data.frame() and print() for PairedPSCBS.
# 2010-09-03
# o Added plot().
# 2010-07-09
# o The segmentByPairedPSCBS() method was written completely from scratch.
# o Created.
############################################################################
