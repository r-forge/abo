###########################################################################/**
# @RdocDefault segmentByCBS
#
# @title "Segment genomic signals using the CBS method"
#
# \description{
#  @get "title" of the \pkg{DNAcopy} package.
#  This is a convenient low-level wrapper for the \code{DNAcopy::segment()}
#  method.  It is intended to be applied to one sample and one chromosome
#  at the time.
#  For more details on the Circular Binary Segmentation (CBS) method 
#  see [1,2].
# }
#
# @synopsis
#
# \arguments{
#   \item{y}{A @numeric @vector of J genomic signals to be segmented.}
#   \item{x}{Optional @numeric @vector of J genomic locations.
#            If @NULL, index locations \code{1:J} are used.}
#   \item{w}{Optional @numeric @vector in [0,1] of J weights.}
#   \item{...}{Additional arguments passed to the segmentation function.}
#   \item{knownCPs}{Optional @numeric @vector of known 
#     change point locations.}
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
#   Internally @see "DNAcopy::segment" is used to segment the signals.
#   This segmentation method support weighted segmentation.
#
#   The "DNAcopy::segment" implementation of CBS uses approximation
#   through random sampling for some estimates.  Because of this,
#   repeated calls using the same signals may result in slightly 
#   different results.  
# }
#
# \section{Missing and non-finite values}{
#   The signals to be segmented as well as any optional positions
#   must not contain missing values, i.e. @NAs or @NaNs.
#   If there are any, an informative error is thrown.
#   Furthermore, none of the input signals may have infinite values,
#   i.e. -@Inf or @Inf. If so, an informative error is thrown.
# }
#
# @examples "../incl/segmentByCBS.Rex"
#
# @author
#
# \references{
#  [1] A.B. Olshen, E.S. Venkatraman (aka Venkatraman E. Seshan), 
#      R. Lucito and M. Wigler, \emph{Circular binary segmentation for 
#      the analysis of array-based DNA copy number data},
#      Biostatistics, 2004.\cr
#  [2] E.S. Venkatraman and A.B. Olshen, \emph{A faster circular binary
#      segmentation algorithm for the analysis of array CGH data}. 
#      Bioinformatics, 2007.\cr
# }
#
# @keyword IO
#*/########################################################################### 
setMethodS3("segmentByCBS", "default", function(y, x=NULL, w=NULL, ..., knownCPs=NULL, seed=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  disallow <- c("NA", "NaN", "Inf");
  y <- Arguments$getDoubles(y, disallow=disallow);
  nbrOfLoci <- length(y);

  length2 <- rep(nbrOfLoci, 2);

  # Argument 'x':
  if (!is.null(x)) {
    x <- Arguments$getDoubles(x, length=length2, disallow=disallow);
  }

  # Argument 'w':
  hasWeights <- !is.null(w);
  if (hasWeights) {
    w <- Arguments$getDoubles(w, range=c(0,1), length=length2, disallow=disallow);
  }

  # Argument 'knownChangePoints':
  if (!is.null(knownCPs)) {
    if (is.null(x)) {
      knownCPs <- Arguments$getIndices(knownCPs, max=nbrOfLoci);
    } else {
      knownCPs <- Arguments$getDoubles(knownCPs);
    }
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


  verbose && enter(verbose, "Segmenting by CBS");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving the fit function");
  pkgName <- "DNAcopy";
  # Assert that package is installed
  isPackageInstalled(pkgName) || throw("Package is not installed: ", pkgName);
  pkg <- packageDescription(pkgName);
  pkgVer <- pkg$Version;
  pkgDetails <- sprintf("%s v%s", pkgName, pkgVer);

  methodName <- "segment";
  verbose && cat(verbose, "Method: ", methodName);
  verbose && cat(verbose, "Package: ", pkgDetails);

  # We need to load package
  require(pkgName, character.only=TRUE) || throw("Package not loaded: ", pkgName);

  # Get the fit function for the segmentation method
#  fitFcn <- getExportedValue(pkgName, methodName);
  fitFcn <- getFromNamespace(methodName, pkgName);
  verbose && str(verbose, "Function: ", fitFcn);
  formals <- formals(fitFcn);
  verbose && cat(verbose, "Formals:");
  verbose && str(verbose, formals);
  verbose && exit(verbose);
 

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up arguments to pass to segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Setting up method arguments");

  verbose && enter(verbose, "Setting up ", pkgName, " data structure"); 

  sampleName <- "Unnamed sample";
  chromosome <- 0;

  if (is.null(x)) {
    x <- seq(length=nbrOfLoci);
  }

  cnData <- DNAcopy::CNA(
    genomdat  = y,
    chrom     = rep(chromosome, times=nbrOfLoci),
    data.type = "logratio",
    maploc    = x,
    sampleid  = sampleName
  );
  verbose && str(verbose, cnData);
  names(cnData)[3] <- sampleName;
  verbose && str(verbose, cnData);
  verbose && exit(verbose);


  params <- list();
  if (hasWeights) {
    params$weights <- w;
    verbose && cat(verbose, "Segmentation parameters:");
    verbose && str(verbose, params);
  }

  userArgs <- list(...);
  if (length(userArgs) > 0) {
    verbose && cat(verbose, "User and segmentation arguments:");
    verbose && str(verbose, userArgs);
    # Assign/overwrite by user arguments
    for (ff in names(userArgs)) {
      params[[ff]] <- userArgs[[ff]];
    }
  }

  # Cleaning out unknown parameters
  keep <- (names(params) %in% names(formals));
  params <- params[keep];

  args <- c(list(cnData), params, verbose=as.logical(verbose));
  verbose && cat(verbose, "Final arguments:");
  verbose && str(verbose, args);

  verbose && exit(verbose);
 

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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calling segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, sprintf("Calling %s() of %s", methodName, pkgName));

  # WORKAROUND for the case when there are no data points.
  if (nbrOfLoci == 0) {
    args[[1]] <- CNA(genomdat=0, chrom=0, maploc=0);
  }

  # In case the method writes to stdout, we capture it
  # Note: DNAcopy::segment() *does* this.
  stdout <- capture.output({
    # Does not work, because some internal function of the fit function
    # may only be accessible from within the namespace
    # How to do this for DNAcopy::segment()? /HB
##    fit <- do.call(fitFcn, args);
    # This works, but requires that one loads the package and that the
    # function is not masked in the search() path.
    t <- system.time({
      fit <- do.call(methodName, args);
    }); 
    # Drop the 'call' (because it will be huge due to the do.call() call)
    fit$call <- NULL;
  });
  attr(fit, "processingTime") <- t;
  attr(fit, "pkgDetails") <- pkgDetails;
  attr(fit, "randomSeed") <- seed;

  # WORKAROUND for the case when there are no data points.
  if (nbrOfLoci == 0) {
    # Drop dummy data point...
    fit$data <- fit$data[-1,,drop=FALSE];
    # ...dummy region found
    fit$output <- fit$output[-1,,drop=FALSE];
  }

  verbose && cat(verbose, "Captured output that was sent to stdout:");
  stdout <- paste(stdout, collapse="\n");
  verbose && cat(verbose, stdout);

  verbose && cat(verbose, "Fitting time (in seconds):");
  verbose && print(verbose, t);

  verbose && cat(verbose, "Fitting time per 1000 loci (in seconds):");
  verbose && print(verbose, 1000*t/nbrOfLoci);

  # Coerce
  fit$output$num.mark <- as.integer(fit$output$num.mark);

  verbose && cat(verbose, "Results object:");
  verbose && str(verbose, fit);

  verbose && exit(verbose); 


  verbose && exit(verbose);
  
  fit; 
}) # segmentByCBS()



############################################################################
# HISTORY:
# 2010-09-02
# o ROBUSTNESS: Now segmentByCBS() also works if there are no data points.
# 2010-07-09
# o Created from segmentByCBS() for RawGenomicSignals in aroma.core.
#   The latter will eventually call this method.
############################################################################
