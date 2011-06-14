###########################################################################/**
# @RdocDefault callSegmentationOutliers
#
# @title "Calls single-locus outliers along the genome"
#
# \description{
#  @get "title" that have a signal that differ significantly from the
#  neighboring loci.
# }
#
# @synopsis
#
# \arguments{
#   \item{y}{A @numeric @vector of J genomic signals to be segmented.}
#   \item{chromosome}{(Optional) An @integer scalar 
#       (or a @vector of length J contain a unique value).
#       Only used for annotation purposes.}
#   \item{x}{Optional @numeric @vector of J genomic locations.
#            If @NULL, index locations \code{1:J} are used.}
#   \item{method}{A @character string specifying the Method
#        used for calling outliers.}
#   \item{...}{Additional arguments passed to internal outlier
#        detection method.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns @logical @vector of length J.
# }
# 
# \section{Missing and non-finite values}{
#   Signals as well as genomic positions may contain missing
#   values, i.e. @NAs or @NaNs.  By definition, these cannot
#   be outliers.
# }
#
# @author
#
# \seealso{
#   Internally @see "DNAcopy::smooth.CNA" is utilized to identify
#   the outliers.
# }
#
# @keyword IO
#*/########################################################################### 
setMethodS3("callSegmentationOutliers", "default", function(y, chromosome=0, x=NULL, method="DNAcopy::smooth.CNA", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  disallow <- c("Inf");
  y <- Arguments$getDoubles(y, disallow=disallow);
  nbrOfLoci <- length(y);

  length2 <- rep(nbrOfLoci, times=2);

  # Argument 'chromosome':
  disallow <- c("Inf");
  chromosome <- Arguments$getIntegers(chromosome, range=c(0,Inf), disallow=disallow);
  if (length(chromosome) == 1) {
    chromosome <- rep(chromosome, times=nbrOfLoci);
  } else {
    chromosome <- Arguments$getVector(chromosome, length=length2);
  }

  # Argument 'x':
  if (!is.null(x)) {
    disallow <- c("Inf");
    x <- Arguments$getDoubles(x, length=length2, disallow=disallow);
  }

  # Argument 'method':
  method <- match.arg(method);

  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }




  verbose && enter(verbose, "Identifying outliers");
  uChromosomes <- sort(unique(chromosome));
  nbrOfChromosomes <- length(uChromosomes);
  verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes);
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);
  verbose && cat(verbose, "Detection method: ", method);

  # Allocate result vector
  isOutlier <- logical(nbrOfLoci);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Filter missing data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying loci with non-missing data");
  keep <- (!is.na(x) & !is.na(y));
  if (!is.null(chromosome)) {
    keep <- (keep & !is.na(chromosome));
  }
  keep <- which(keep);
  chromosome <- chromosome[keep];
  x <- x[keep];
  y <- y[keep];
  nbrOfLoci <- length(x);
  verbose && cat(verbose, "Number of loci with non-missing data: ", nbrOfLoci);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For each chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  isOutlierT <- logical(nbrOfLoci);
  for (kk in seq(along=uChromosomes)) {
    chr <- uChromosomes[kk];
    verbose && enter(verbose, sprintf("Chromosome #%d ('Chr%02d') of %d", kk, chr, length(uChromosomes)));
    keepKK <- which(chromosome == chr);
    nbrOfLociKK <- length(keepKK);
    verbose && cat(verbose, "Number of loci on chromosome: ", nbrOfLociKK);

    # Extract data
    yKK <- y[keepKK];
    xKK <- x[keepKK];
    chromosomeKK <- chromosome[keepKK];

    # Order loci along chromosome
    o <- order(xKK);
    xKK <- xKK[o];
    yKK <- yKK[o];
    chromosomeKK <- chromosomeKK[o];
    keepKK <- keepKK[o];
    
    dataKK <- DNAcopy::CNA(genomdat=yKK, chrom=chromosomeKK, maploc=xKK, sampleid="y", presorted=TRUE);
    yKKs <- DNAcopy::smooth.CNA(dataKK, ...)$y;

    # Sanity check
    stopifnot(length(yKKs) == nbrOfLociKK);

    outliersKK <- which(yKKs != yKK);
    nbrOfOutliers <- length(outliersKK);
    verbose && cat(verbose, "Number of outliers: ", nbrOfOutliers);

    outliers <- keepKK[outliersKK];
    isOutlierT[outliers] <- TRUE;

    verbose && exit(verbose);
  } # for (kk ...)

  isOutlier[keep] <- isOutlierT;

  nbrOfOutliers <- sum(isOutlier, na.rm=TRUE);
  verbose && cat(verbose, "Total number of outliers: ", nbrOfOutliers);

  verbose && exit(verbose);

  isOutlier;
}) # callSegmentationOutliers()


setMethodS3("dropSegmentationOutliers", "default", function(y, ...) {
  isOutlier <- callSegmentationOutliers(y, ...);
  y[isOutlier] <- as.double(NA);
  y;
})


############################################################################
# HISTORY:
# 2011-05-31
# o Now explicitly using DNAcopy::nnn() to call DNAcopy functions.
# 2010-11-27
# o Added dropSegmentationOutliers() which sets outliers to missing values.
# o Added callSegmentationOutliers(), which utilizes the detection method
#   of DNAcopy::smooth.CNA() as suggested by ABO.
# o Created.
############################################################################