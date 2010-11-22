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
#   \item{chromosome}{(Optional) An @integer scalar 
#       (or a @vector of length J contain a unique value).
#       Only used for annotation purposes.}
#   \item{x}{Optional @numeric @vector of J genomic locations.
#            If @NULL, index locations \code{1:J} are used.}
#   \item{w}{Optional @numeric @vector in [0,1] of J weights.}
#   \item{undo}{A non-negative @numeric.  If less than +@Inf, then
#       arguments \code{undo.splits="sdundo"} and \code{undo.SD=undo}
#       are passed to \code{DNAcopy::segment()}.}
#   \item{...}{Additional arguments passed to the segmentation function.}
#   \item{preserveOrder}{If @TRUE, the returned \code{data} field is
#     ordered as \code{x}, otherwise along the genome.}
#   \item{joinSegments}{If @TRUE, there are no gaps between neighboring
#     segments.
#     If @FALSE, the boundaries of a segment are defined by the support
#     that the loci in the segments provides, i.e. there exist a locus
#     at each end point of each segment.  This also means that there
#     is a gap between any neighboring segments, unless the change point
#     is in the middle of multiple loci with the same position.
#     The latter is what \code{DNAcopy::segment()} returns.
#   }
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
#   Signals as well as genomic positions may contain missing
#   values, i.e. @NAs or @NaNs.  If so, they are silently excluded
#   before performing the segmentation.
#
#   None of the input signals may have infinite values,
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
setMethodS3("segmentByCBS", "default", function(y, chromosome=0, x=NULL, w=NULL, undo=Inf, ..., preserveOrder=FALSE, joinSegments=TRUE, knownCPs=NULL, seed=NULL, verbose=FALSE) {
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
  if (length(chromosome) > 1) {
    chromosome <- Arguments$getIntegers(chromosome, length=length2, disallow=disallow);
    # If 'chromosome' is a vector of length J, then it must contain
    # a unique chromosome.
    chromosomes <- sort(unique(chromosome));
    if (length(chromosomes) > 1) {
      throw("Argument 'chromosome' specifies more than one unique chromosome: ", paste(seqToHumanReadable(chromosomes), collapse=", "));
    }
    chromosome <- chromosomes;
  }

  # For future usage
  chrom <- rep(chromosome, times=nbrOfLoci);

  # Argument 'x':
  if (!is.null(x)) {
    disallow <- c("Inf");
    x <- Arguments$getDoubles(x, length=length2, disallow=disallow);
  }

  # Argument 'w':
  hasWeights <- !is.null(w);
  if (hasWeights) {
    disallow <- c("NA", "NaN", "Inf");
    w <- Arguments$getDoubles(w, range=c(0,1), length=length2, disallow=disallow);
  }

  # Argument 'undo':
  undo <- Arguments$getDouble(undo, range=c(0,Inf));

  # Argument 'preserveOrder':
  preserveOrder <- Arguments$getLogical(preserveOrder);

  # Argument 'cpFlavor':
  joinSegments <- Arguments$getLogical(joinSegments);

  # Argument 'knownCPs':
  if (!is.null(knownCPs)) {
    if (is.null(x)) {
      knownCPs <- Arguments$getIndices(knownCPs, max=nbrOfLoci);
    } else {
      knownCPs <- Arguments$getDoubles(knownCPs);
    }
    if (length(knownCPs) != 2) {
      throw("Currently argument 'knownCPs' can be used to specify the boundaries of the region to be segmented: ", length(knownCPs));
      throw("Support for specifying known change points (argument 'knownCPs') is not yet implemented as of 2010-10-02.");
    }
    if (!joinSegments) {
      throw("Argument 'knownCPs' should only be specified if argument 'joinSegments' is TRUE.");
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
  # Reorder data points along the genome, because that is what
  # DNAcopy::segment() will return.  At the end, we will undo
  # the sort such that the returned 'data' object is always in
  # the same order and number of loci as the input data.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  wasSorted <- FALSE;
  if (is.null(x)) {
    x <- seq(length=nbrOfLoci);
    index <- NULL;
  } else {
    x0 <- x;
    # Sort the data along the chromosome
    o <- order(x, decreasing=FALSE, na.last=TRUE);

    # Any change?
    wasSorted <- any(o != seq(length=nbrOfLoci));

    # Then reorder.
    if (wasSorted) {
      verbose && enter(verbose, "Ordering loci along genome");
      chrom <- chrom[o];
      x <- x[o];
      y <- y[o];
      if (hasWeights) {
        w <- w[o];
      }
      # Record the ordering
      index <- o;
      verbose && exit(verbose);
    }
    rm(o); # Not needed anymore
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop data points without known genomic positions, because that
  # is what DNAcopy::CNA() will do otherwise.  At the end, we will 
  # undo this such that the returned 'data' object is complete.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nok <- is.na(x);
  hadMissing <- any(nok);
  if (hadMissing) {
    verbose && enter(verbose, "Dropping loci with unknown locations");
    keep <- which(!nok);
    verbose && cat(verbose, "Number of loci dropped: ", sum(nok));
    chrom <- chrom[keep];
    x <- x[keep];
    y <- y[keep];
    if (hasWeights) {
      w <- w[keep];
    }
    rm(keep);
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up arguments to pass to segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Setting up method arguments");

  verbose && enter(verbose, "Setting up ", pkgName, " data structure"); 

  sampleName <- "y";  # This is going to be the name of the data field

  cnData <- DNAcopy::CNA(
    genomdat  = y,
    chrom     = chrom,
    data.type = "logratio",
    maploc    = x,
    sampleid  = sampleName,
    presorted = TRUE
  );
  verbose && str(verbose, cnData);
  names(cnData)[3] <- sampleName;
  verbose && str(verbose, cnData);
  verbose && exit(verbose);
  rm(x, y, chromosome, sampleName);  # Not needed anymore


  params <- list();
  if (hasWeights) {
    params$weights <- w;
  }

  if (undo < Inf) {
    params$undo.splits <- "sdundo";
    params$undo.SD <- undo;
  }

  verbose && cat(verbose, "Segmentation parameters:");
  verbose && str(verbose, params);

  userArgs <- list(...);
  if (length(userArgs) > 0) {
    verbose && cat(verbose, "User arguments:");
    verbose && str(verbose, userArgs);
    # Assign/overwrite by user arguments
    for (ff in names(userArgs)) {
      params[[ff]] <- userArgs[[ff]];
    }
  }

  verbose && cat(verbose, "Segmentation and user parameters:");
  verbose && str(verbose, params);

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

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Undo the dropping of loci with unknown locations
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (hadMissing) {
    verbose && enter(verbose, "Undo dropping of loci with unknown locations");
    idxs <- rep(as.double(NA), times=nbrOfLoci);
    keep <- which(!nok);
    idxs[keep] <- seq(along=keep);
    fit$data <- fit$data[idxs,];
    rm(keep);
    verbose && exit(verbose);
  }

  # Sanity check (does DNAcopy::segment() drop data points?!?)
  stopifnot(nrow(fit$data) == nbrOfLoci);

  verbose && cat(verbose, "Captured output that was sent to stdout:");
  stdout <- paste(stdout, collapse="\n");
  verbose && cat(verbose, stdout);

  verbose && cat(verbose, "Fitting time (in seconds):");
  verbose && print(verbose, t);

  verbose && cat(verbose, "Fitting time per 1000 loci (in seconds):");
  verbose && print(verbose, 1000*t/nbrOfLoci);

  # Coerce
  fit$output$num.mark <- as.integer(fit$output$num.mark);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Undo the ordering along the genome, or add the index map?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (wasSorted) {
    if (preserveOrder) {
      verbose && enter(verbose, "Undoing reordering along genome");

      o <- order(index);
      fit$data <- fit$data[o,];
      rm(o);
      verbose && cat(verbose, "Reordering of data undone:");
      verbose && str(verbose, fit$data);

      # Sanity check
      stopifnot(all(fit$data$maploc == x0));

      verbose && exit(verbose);
    } else {
      fit$data$index <- index;

      # Sanity check
      stopifnot(all(fit$data$maploc == x0[fit$data$index], na.rm=TRUE));
    }
  }
  rm(x0); # Not needed anymore

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Special case: Neighboring segments with "overlapping" end points
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Details: In case there are multiple loci at the same genomic
  # positions ('maploc') and CBS happens to detect a change points
  # in the middle of these loci, we here will infer, by "reverse
  # engineering" which of these loci belongs to which segment.
  segs <- fit$output;
  nbrOfSegs <- nrow(segs);
  if (nbrOfSegs > 1) {
    # Identify change points of such cases
    cpIdxs <- which(segs[-nbrOfSegs,"loc.end"] == segs[-1,"loc.start"]);
    if (length(cpIdxs) > 0) {
      verbose && enter(verbose, "Detected segments that share the same end points");
      verbose && cat(verbose, "Number of shared end points: ", length(cpIdxs));

      # The genomic positions for which there exist loci with
      # the sample 'maploc' and that have been split between 
      # two segments.
      maplocs <- segs[cpIdxs,"loc.end"];
      verbose && cat(verbose, "Locations of shared end points:");
      verbose && print(verbose, maplocs);

      segIdxs <- unique(sort(c(cpIdxs, cpIdxs+1L)));
      verbose && cat(verbose, "Number of such segments: ", length(segIdxs));
      segsT <- segs[segIdxs,];
      verbose && print(verbose, segsT);

      # WARNING: Above we ordered the units along the chromosome.  This means
      # that the 'data' object used here is ordered that way, which in turn
      # means that the unit indices identified below match 'data' and not the
      # input data, e.g. 'x'.  This is also why we record the inverse index map.
      # /HB 2010-11-16
      x <- fit$data$maploc;

      lociNotPartOfSegment <- vector("list", length=nbrOfSegs);
      names(lociNotPartOfSegment) <- rownames(segs);

      prevSeg <- segsT[1L,];
      for (ss in 2:nrow(segsT)) {
        currSeg <- segsT[ss,];
        currStart <- currSeg[,"loc.start"];
        prevEnd <- prevSeg[,"loc.end"];
        if (currStart == prevEnd) {
          currCount <- currSeg[,"num.mark"];
          currEnd <- currSeg[,"loc.end"];
          units <- which(currStart <= x & x <= currEnd);
          nbrOfUnits <- length(units);
          nbrToMove <- nbrOfUnits - currCount;

          # It should always be the case that nbrToMove > 0
          verbose && cat(verbose, "Number of loci to move: ", nbrToMove);

          # Identify units to move
          units <- units[x[units] == currStart];
          verbose && cat(verbose, "Number of loci with the sample genomic position: ", length(units));

          # Sanity check
          stopifnot(length(units) > nbrToMove);

          # Alt 1. Move as they occur in the data vector
          toMove <- 1:nbrToMove;

          # Alt 2. ...or by some other rule?!?

          unitsToMove <- units[toMove];
          unitsToKeep <- units[-toMove];

          segIdx <- segIdxs[ss];
          lociNotPartOfSegment[[segIdx]] <- c(lociNotPartOfSegment[[segIdx]], unitsToMove);
          lociNotPartOfSegment[[segIdx-1L]] <- c(lociNotPartOfSegment[[segIdx-1L]], unitsToKeep);
        }

        prevSeg <- currSeg;
      } # for (ss ...)
      verbose && cat(verbose, "Loci not part of segment:");
      verbose && str(verbose, lociNotPartOfSegment);

      verbose && exit(verbose);

      fit$lociNotPartOfSegment <- lociNotPartOfSegment;
    } # if (length(cpIdxs) > 0)
  } # if (nbrOfSegs > 1)

  params <- list(
    preserveOrder = preserveOrder,
    joinSegments = joinSegments,
    knownCPs = knownCPs
  );

  fit$params <- params;

  class(fit) <- c("CBS", class(fit));

  # Join segments?
  if (joinSegments) {
    fit <- joinSegments(fit, range=knownCPs, verbose=verbose);
  }

  verbose && cat(verbose, "Results object:");
  verbose && str(verbose, fit);

  verbose && exit(verbose); 


  verbose && exit(verbose);
  
  fit; 
}) # segmentByCBS()



############################################################################
# HISTORY:
# 2010-11-21
# o Now segmentByCBS(..., joinSegments=TRUE) utilizes joinSegments().
# 2010-11-20
# o Now it is possible to specify the boundaries of the regions to be
#   segmented as known change points via argument 'knownCPs'.
# 2010-11-19
# o Added argument 'joinSegments' to segmentByCBS() in order to specify
#   if neighboring segments should be joined or not.
# o Now segmentByCBS() returns an object of class CBS.
# o Now segmentByCBS() allows for unknown genomic positions.
# o Now segmentByCBS() allows for missing signals.
# o Added argument 'preserveOrder' to segmentByCBS().  If TRUE, then
#   the loci in the returned 'data' object are ordered as the input
#   data, otherwise it is ordered along the genome.
# 2010-11-16
# o Now the 'data' object returned by segmentByCBS() contains field
#   'index' if and only if the loci had to be reorder along the genome.
# 2010-11-02
# o Added argument 'undo' to segmentByCBS(), which corresponds to
#   undo.splits="sdundo" and undo.SD=undo, if undo < Inf.
# 2010-10-25
# o Now segmentByCBS() also returns element 'lociNotPartOfSegment',
#   if there are segments that share end points, which can happen if
#   a change point is called in middle of a set of loci that have the
#   same genomic positions.  In such cases, 'lociNotPartOfSegment'
#   specifies which loci are *not* part of which segment.  Then by
#   identifying the loci that are within a segment by their positions
#   and excluding any of the above, one knows exactly which loci
#   CBS included in each segment.
# 2010-10-02
# o Added argument optional 'chromosome'.
# 2010-09-02
# o ROBUSTNESS: Now segmentByCBS() also works if there are no data points.
# 2010-07-09
# o Created from segmentByCBS() for RawGenomicSignals in aroma.core.
#   The latter will eventually call this method.
############################################################################
