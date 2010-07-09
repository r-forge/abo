###########################################################################/**
# @RdocDefault segmentByPairedPSCBS
#
# @title "Segment copy numbers using the paired PSCBS method"
#
# \description{
#  @get "title".
#  This is a low-level method.
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
#   different results.  
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
setMethodS3("segmentByPairedPSCBS", "default", function(CT, betaT, betaN, muN=NULL, x=NULL, ..., seed=NULL, verbose=FALSE) {
  require("R.utils") || throw("Package not loaded: R.utils");
  require("aroma.light") || throw("Package not loaded: aroma.light");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'CT':
  CT <- Arguments$getDoubles(CT);
  nbrOfLoci <- length(CT);
  length2 <- rep(nbrOfLoci, 2);

  # Argument 'betaT':
  betaT <- Arguments$getDoubles(betaT, length=length2);

  # Argument 'betaN':
  betaN <- Arguments$getDoubles(betaN, length=length2);

  # Argument 'muN':
  if (!is.null(muN)) {
    muN <- Arguments$getDoubles(muN, length=length2, range=c(0,1));
  }

  # Argument 'x':
  if (!is.null(x)) {
    x <- Arguments$getDoubles(x, length=length2);
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
  betaTN <- aroma.light::normalizeTumorBoost(betaT=betaT, betaN=betaN, muN=muN);
  verbose && cat(verbose, "Normalized betaT:");
  verbose && str(verbose, betaTN);
  verbose && exit(verbose);


  if (is.null(x)) {
    x <- seq(length=nbrOfLoci);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1a. Identification of change points in total copy numbers
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identification of change points by total copy numbers");
  h <- function(x) { sqrt(x) };
  hinv <- function(y) { y^2 };
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

  fit <- segmentByCBS(yT, x=x, verbose=verbose);
  tcnSegments <- fit$output;
  rm(yT);

  verbose && enter(verbose, "Backtransforming segmented mean levels");
  tcnSegments[,"seg.mean"] <- hinv(tcnSegments[,"seg.mean"]);
  verbose && exit(verbose);

  verbose && print(verbose, tcnSegments);
  nbrOfSegs <- nrow(tcnSegments);
  verbose && cat(verbose, "Number of TCN segments: ", nbrOfSegs);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1b. Identification of additional change points using DH
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For each segment independently, segment decrease of heterozygousity (DH)
  # using CBS. DH is by definition only defined for heterozygous SNPs.
  # Only heterozygous SNPs are used.

  # Calculate DH, which is NA for non-heterozygous loci.
  isHet <- (muN == 1/2);
  naValue <- as.double(NA);
  rho <- rep(naValue, length=nbrOfLoci);
  rho[isHet] <- 2*abs(betaTN[isHet]-1/2);

  # For each TCN segment...
  segs <- vector("list", length=nbrOfSegs);
  for (kk in seq(length=nbrOfSegs)) {
    verbose && enter(verbose, sprintf("Total CN segment #%d of %d", kk, nbrOfSegs));

    xStart <- tcnSegments[kk,"loc.start"];
    xEnd <- tcnSegments[kk,"loc.end"];
    gammaT <- tcnSegments[kk,"seg.mean"];

    keep <- (xStart <= x & x <= xEnd);
    xKK <- x[keep];
    verbose && cat(verbose, "Number of loci: ", length(xKK));

    keep <- keep & isHet;
    xKKHet <- x[keep];
    rhoKKHet <- rho[keep];
    verbose && cat(verbose, "Number of heterozygous SNPs: ", length(xKKHet));

    # Drop non-finite data points
    keep <- (is.finite(xKKHet) & is.finite(rhoKKHet));
    xKKHet <- xKKHet[keep];
    rhoKKHet <- rhoKKHet[keep];

    # Nothing to do?
    if (length(xKKHet) == 0) {
      verbose && cat(verbose, "No heterozygous SNPs with (finite) DHs. Skipping DH segmentation.");
      # Need to return something for the total CNs...
      verbose && exit(verbose);    
      next;
    }

    fit <- segmentByCBS(rhoKKHet, x=xKKHet, verbose=less(verbose, 5));

    # Extract regions
    dhSegments <- fit$output;
    names <- names(dhSegments);
    names <- gsub("seg.mean", "dh.mean", names, fixed=TRUE);
    names <- gsub("num.mark", "num.het.mark", names, fixed=TRUE);
    names(dhSegments) <- names;

    # For each DH segment, count number of loci and update TCN mean
    dhSegments$num.mark <- -1L;
    dhSegments$tcn.mean <- gammaT;
    for (ll in seq(length=nrow(dhSegments))) {
      xStart <- dhSegments[ll,"loc.start"];
      xEnd <- dhSegments[ll,"loc.end"];
      keep <- (xStart <= x & x <= xEnd);
      dhSegments[ll,"num.mark"] <- sum(keep);
      dhSegments[ll,"tcn.mean"] <- mean(CT[keep], na.rm=TRUE);
    }

    verbose && print(verbose, dhSegments);

    segs[[kk]] <- dhSegments;

    verbose && exit(verbose);    
  } # for (kk ...)

  segs <- Reduce(rbind, segs);

  # Reorder columns
  fields <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "num.het.mark", "dh.mean", "tcn.mean");
  segs <- segs[,fields];
  verbose && print(verbose, segs);

  nbrOfSegs <- nrow(segs);
  verbose && cat(verbose, "Number of TCN and DH segments: ", nbrOfSegs);

  verbose && exit(verbose);

  # Create result object
  data <- list(CT=CT, betaT=betaT, betaN=betaN, muN=muN, x=x);

  fit <- list(
    data = data,
    output = segs
  );
  class(fit) <- c("PSCBS", "DNAcopy");

  # Return segments found
  fit;
}) # segmentByPairedPSCBS()




setMethodS3("callSegmentsByPSCBS", "default", function(data, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # 2. Calling regions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Call region
  # a. Testing for LOH in the Tumor
  # b. Classifying non-LOH regions as balanced or not balanced
}) # callSegmentsByPSCBS()


############################################################################
# HISTORY:
# 2010-07-09
# o The segmentByPairedPSCBS() method was written completely from scratch.
# o Created.
############################################################################
