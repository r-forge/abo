###########################################################################/**
# @RdocDefault segmentByPairedPSCBS
# @alias segmentByPairedPSCBS
#
# @title "Segment copy numbers using the paired PSCBS method"
#
# \description{
#  @get "title" [1].
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
#   \item{muN}{An optional @numeric @vector of J genotype calls in 
#        \{0,1/2,1\} for AA, AB, and BB, respectively. If not given,
#        they are estimated from the normal BAFs using
#        @see "aroma.light::callNaiveGenotypes" as described in [2].}
#   \item{chromosome}{(Optional) An @integer scalar 
#       (or a @vector of length J contain a unique value).
#       Only used for annotation purposes.}
#   \item{x}{Optional @numeric @vector of J genomic locations.
#            If @NULL, index locations \code{1:J} are used.}
#   \item{alphaTCN, alphaDH}{The significance levels for segmenting total
#        copy numbers (TCNs) and decrease-in-heterozygosity signals (DHs),
#        respectively.}
#   \item{undoTCN, undoDH}{Non-negative @numerics.  If less than +@Inf, 
#        then a cleanup of segmentions post segmentation is done.
#        See argument \code{undo} of @see "segmentByCBS" for more
#        details.}
#   \item{...}{Additional arguments passed to @see "segmentByCBS".}
#   \item{flavor}{A @character specifying what type of segmentation and 
#     calling algorithm to be used.}
#   \item{tbn}{If @TRUE, \code{betaT} is normalized before segmentation
#     using the TumorBoost method [2], otherwise not.}
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
#   None of the input signals may have infinite values, i.e. -@Inf or @Inf.
#   If so, an informative error is thrown.
# }
#
# @examples "../incl/segmentByPairedPSCBS.Rex"
#
# @author
#
# \seealso{
#   Internally, @see "aroma.light::callNaiveGenotypes" is used to 
#   call naive genotypes, @see "aroma.light::normalizeTumorBoost" is 
#   used for TumorBoost normalization, and @see "segmentByCBS" is used 
#   to segment TCN and DH separately.
# }
#
# \references{
#   [1] A. Olshen, R. Olshen, H. Bengtsson, P. Neuvial, P. Spellman and P. Seshan, \emph{Parent-specific copy number in paired tumor-normal studies using circular binary segmentation}, 2010 (submitted).\cr
#   [2] H. Bengtsson, P. Neuvial and T.P. Speed, \emph{TumorBoost: Normalization of allele-specific tumor copy numbers from a single pair of tumor-normal genotyping microarrays}, BMC Bioinformatics, 2010.
# }
#
# @keyword IO
#*/########################################################################### 
setMethodS3("segmentByPairedPSCBS", "default", function(CT, betaT, betaN, muN=NULL, chromosome=0, x=NULL, alphaTCN=0.009, alphaDH=0.001, undoTCN=Inf, undoDH=Inf, ..., flavor=c("tcn&dh", "tcn,dh", "sqrt(tcn),dh", "sqrt(tcn)&dh"), tbn=TRUE, joinSegments=TRUE, knownCPs=NULL, seed=NULL, verbose=FALSE) {
  require("R.utils") || throw("Package not loaded: R.utils");
  require("aroma.light") || throw("Package not loaded: aroma.light");
  ver <- packageDescription("aroma.light")$Version;
  stopifnot(compareVersion(ver, "1.17.2") >= 0);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'CT':
  disallow <- c("Inf");
  CT <- Arguments$getDoubles(CT, disallow=disallow);
  nbrOfLoci <- length(CT);
  length2 <- rep(nbrOfLoci, times=2);

  # Argument 'betaT':
  betaT <- Arguments$getDoubles(betaT, length=length2, disallow="Inf");

  # Argument 'betaN':
 betaN <- Arguments$getDoubles(betaN, length=length2, disallow="Inf");

  # Argument 'muN':
  if (!is.null(muN)) {
    muN <- Arguments$getDoubles(muN, length=length2, range=c(0,1), disallow="Inf");
  }

  # Argument 'chromosome':
  disallow <- c("Inf");
  chromosome <- Arguments$getIntegers(chromosome, range=c(0,Inf), disallow=disallow);
  if (length(chromosome) > 1) {
    chromosome <- Arguments$getIntegers(chromosome, length=length2, disallow=disallow);
  }

  # Argument 'x':
  if (!is.null(x)) {
    disallow <- c("Inf");
    x <- Arguments$getDoubles(x, length=length2, disallow=disallow);
  }

  # Argument 'alphaTCN':
  alphaTCN <- Arguments$getDouble(alphaTCN, range=c(0,1));

  # Argument 'alphaDH':
  alphaDH <- Arguments$getDouble(alphaDH, range=c(0,1));

  # Argument 'undoTCN':
  undoTCN <- Arguments$getDouble(undoTCN, range=c(0,Inf));

  # Argument 'undoDH':
  undoDH <- Arguments$getDouble(undoDH, range=c(0,Inf));


  # Argument 'flavor':
  flavor <- match.arg(flavor);
  knownFlavors <- c("tcn,dh", "tcn&dh", "sqrt(tcn),dh", "sqrt(tcn)&dh");
  if (!is.element(flavor, knownFlavors)) {
    throw("Segmentation flavor is not among the supported ones (", paste(sprintf("\"%s\"", knownFlavors), collapse=", "), "): ", flavor);
  }

  # Argument 'tbn':
  tbn <- Arguments$getLogical(tbn);

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


  verbose && enter(verbose, "Segmenting paired tumor-normal signals using PSCBS");

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
    muN <- aroma.light::callNaiveGenotypes(betaN, censorAt=c(0,1));
    verbose && cat(verbose, "Called genotypes:");
    verbose && str(verbose, muN);
    verbose && print(verbose, table(muN));
    verbose && exit(verbose);
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize betaT using betaN (TumorBoost normalization)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (tbn) {
    verbose && enter(verbose, "Normalizing betaT using betaN (TumorBoost)");
    betaTN <- aroma.light::normalizeTumorBoost(betaT=betaT, betaN=betaN, muN=muN, preserveScale=TRUE);
    verbose && cat(verbose, "Normalized BAFs:");
    verbose && str(verbose, betaTN);

    # Assert that no missing values where introduced
    keep <- (is.finite(betaT) & is.finite(betaN) & is.finite(muN));
    if (any(is.na(betaTN[keep]))) {
      throw("Internal error: normalizeTumorBoost() introduced missing values.");
    }
    rm(keep);
    verbose && exit(verbose);
  } else {
    betaTN <- betaT;
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Multiple chromosomes?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify all chromosomes, excluding missing values
  chromosomes <- sort(unique(chromosome), na.last=NA);
  nbrOfChromosomes <- length(chromosomes);
  if (nbrOfChromosomes > 1) {
    verbose && enter(verbose, "Segmenting multiple chromosomes");
    verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes);

    fitList <- list();
    for (kk in seq(length=nbrOfChromosomes)) {
      chromosomeKK <- chromosomes[kk];
      chrTag <- sprintf("Chr%02d", chromosomeKK);
      verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", kk, chrTag, nbrOfChromosomes));

      units <- whichVector(chromosome == chromosomeKK);
      verbose && cat(verbose, "Units:");
      verbose && str(verbose, units); 

      if (is.null(x)) {
        xKK <- NULL;
      } else {
        xKK <- x[units];
      }
      fit <- segmentByPairedPSCBS(CT=CT[units], betaT=betaTN[units], 
                betaN=betaN[units], muN=muN[units], chromosome=chromosomeKK,
                x=xKK, tbn=FALSE, joinSegments=joinSegments,
                alphaTCN=alphaTCN, alphaDH=alphaDH,
                undoTCN=undoTCN, undoDH=undoDH,
                ..., verbose=verbose);
      verbose && print(verbose, head(as.data.frame(fit)));
      verbose && print(verbose, tail(as.data.frame(fit)));
      
      fitList[[chrTag]] <- fit;

      # Not needed anymore
      rm(units, fit);
      verbose && exit(verbose);
    } # for (kk ...)

    verbose && enter(verbose, "Merging");
    fit <- Reduce(append, fitList);
    # Not needed anymore
    rm(fitList);
    verbose && exit(verbose);

    verbose && print(verbose, head(as.data.frame(fit)));
    verbose && print(verbose, tail(as.data.frame(fit)));
   
    verbose && exit(verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Return results
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return(fit);
  } # if (nbrOfChromosomes > 1)




  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup input data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "alphaTCN: ", alphaTCN);
  verbose && cat(verbose, "alphaDH: ", alphaDH);
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);

  # SNPs are identifies as those loci that have non-missing 'betaTN' & 'muN'
  isSnp <- (!is.na(betaTN) & !is.na(muN));
  nbrOfSnps <- sum(isSnp);
  verbose && cat(verbose, "Number of SNPs: ", nbrOfSnps);

  # Expand 'chromosome' to be of full length
  if (length(chromosome) == 1) {
    chromosome <- rep(chromosome, times=nbrOfLoci);
  }


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
  # 1a. Transform total copy-number signals?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The default is not to transform signals
  h <- hinv <- NULL;

  if (is.element(flavor, c("sqrt(tcn),dh", "sqrt(tcn)&dh"))) {
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
  }

  if (!is.null(h)) {
    verbose && enter(verbose, "Transforming signals");
    verbose && cat(verbose, "Transform:");
    verbose && print(verbose, h);

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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1b. Identification of change points in total copy numbers
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identification of change points by total copy numbers");

  # Physical positions of loci
  if (is.null(x)) {
    x <- seq(length=nbrOfLoci);
  }

  fit <- segmentByCBS(yT, chromosome=chromosome, x=x, 
                      joinSegments=joinSegments, knownCPs=knownCPs,
                      alpha=alphaTCN, undo=undoTCN, ..., verbose=verbose);
  verbose && str(verbose, fit);

  tcnSegments <- fit$output;

  # Extract loci that should be excluded for each segment
  tcnLociNotPartOfSegment <- fit$lociNotPartOfSegment;
  if (!is.null(tcnLociNotPartOfSegment)) {
    # The specified loci are indexed according to the fit$data object,
    # which may or may not have been reordered and therefore does not
    # necessarily match the ordering of 'x' here.  If reordered,
    # we need to map the indices in 'tcnLociNotPartOfSegment' back
    # accordingly.
    
    # Was input data reordered? If so, then reorder.
    index <- fit$data$index;
    if (!is.null(index)) {
      tcnLociNotPartOfSegment <- lapply(tcnLociNotPartOfSegment, FUN=function(idxs) {
        index[idxs];
      });
    }

    rm(index); # Not needed anymore
  }

  rm(fit);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1c. Backtransform estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(hinv)) {
    verbose && enter(verbose, "Backtransforming segmented mean levels");
    verbose && cat(verbose, "Backtransform:");
    verbose && print(verbose, hinv);

    tcnSegments[,"seg.mean"] <- hinv(tcnSegments[,"seg.mean"]);

    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1d. Restructure data 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop dummy columns
  keep <- setdiff(colnames(tcnSegments), c("ID"));
  tcnSegments <- tcnSegments[,keep,drop=FALSE];

  # Tag fields by TCN
  names <- names(tcnSegments);
  names <- gsub("seg.mean", "mean", names, fixed=TRUE);
  names <- sprintf("tcn.%s", names);
  names <- gsub("tcn.chrom", "chromosome", names, fixed=TRUE);
  names(tcnSegments) <- names;
  rm(names);
  verbose && print(verbose, tcnSegments);

  nbrOfSegs <- nrow(tcnSegments);
  verbose && cat(verbose, "Number of TCN segments: ", nbrOfSegs);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2a. Identification of additional change points using DH
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For each segment independently, segment decrease of heterozygousity (DH)
  # using CBS. By definition, only heterozygous SNPs are used.

  listOfDhLociNotPartOfSegment <- vector("list", length=nbrOfSegs);
  names(listOfDhLociNotPartOfSegment) <- seq(length=nbrOfSegs);

  # Identify all loci with non-missing signals and locations
  ok <- (!is.na(yT) & !is.na(x));

  # For each TCN segment...
  segs <- vector("list", length=nbrOfSegs);
  for (kk in seq(length=nbrOfSegs)) {
    # Extract the region
    xStart <- tcnSegments[kk,"tcn.loc.start"];
    xEnd <- tcnSegments[kk,"tcn.loc.end"];
    gammaT <- tcnSegments[kk,"tcn.mean"];

    regionTag <- sprintf("[%g,%g]", xStart, xEnd);
    verbose && enter(verbose, sprintf("Total CN segment #%d (%s) of %d", kk, regionTag, nbrOfSegs));

    verbose && cat(verbose, "Number of TCN loci in segment: ", tcnSegments[kk,"tcn.num.mark"]);

    # Identify subset of finite loci
    keep <- which(ok & xStart <= x & x <= xEnd);
    
    # Special case?
    if (!is.null(tcnLociNotPartOfSegment)) {
      lociToExclude <- tcnLociNotPartOfSegment[[kk]];
      if (length(lociToExclude) > 0) {
        verbose && cat(verbose, "Identified number of DHs before: ", length(keep));
        verbose && printf(verbose, "Excluding %d loci that belongs to a flanking segment: units (%s) at positions (%s)\n",
               length(lociToExclude), paste(lociToExclude, collapse=", "), paste(x[lociToExclude], collapse=", "));

        # Sanity check
        stopifnot(all(is.element(lociToExclude, keep)));

        keep <- setdiff(keep, lociToExclude);
        verbose && cat(verbose, "Identified number of DHs afterward: ", length(keep));
      }
    }

    # Sanity check
    stopifnot(length(keep) == tcnSegments[kk,"tcn.num.mark"]);

    chromosomeKK <- chromosome[keep];
    xKK <- x[keep];
    nbrOfLociKK <- length(xKK);
    rhoKK <- rho[keep];
    isSnpKK <- isSnp[keep];
    nbrOfSnpsKK <- sum(isSnpKK);
    isHetKK <- isHet[keep];
    verbose && cat(verbose, "Number of loci in segment: ", nbrOfLociKK);
    verbose && cat(verbose, "Number of SNPs in segment: ", nbrOfSnpsKK);
    rm(keep);

    # Identify heterozygous SNPs
    keep <- isHetKK;
    chromosomeKK <- chromosomeKK[keep];
    xKKHet <- xKK[keep];
    rhoKKHet <- rhoKK[keep];
    nbrOfHetsKK <- length(xKKHet);
    verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                  nbrOfHetsKK, 100*nbrOfHetsKK/nbrOfSnpsKK);
    rm(keep);

    verbose && enter(verbose, "Segmenting DH signals");
    if (joinSegments) {
      knownCPsKK <- c(xStart, xEnd);
    } else {
      knownCPsKK <- NULL;
    }
    fit <- segmentByCBS(rhoKKHet, chromosome=chromosomeKK, x=xKKHet, 
                        joinSegments=joinSegments, knownCPs=knownCPsKK,
                        alpha=alphaDH, undo=undoDH, ..., verbose=verbose);
    verbose && str(verbose, fit);
    dhSegments <- fit$output;
    dhLociNotPartOfSegment <- fit$lociNotPartOfSegment;
    if (!is.null(dhLociNotPartOfSegment)) {
      listOfDhLociNotPartOfSegment[[kk]] <- dhLociNotPartOfSegment;
    }
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

    # Sanity check
    stopifnot(nrow(tcnSegmentsKK) == nrow(dhSegments));

    # Combine TCN and DH segmentation results
    tcndhSegments <- cbind(
      tcn.id=rep(kk, times=nrow(dhSegments)),
      dh.id=seq(length=nrow(dhSegments)),
      tcnSegmentsKK,
      dhSegments
    );

    segs[[kk]] <- tcndhSegments;

    verbose && cat(verbose, "(TCN,DH) segmentation for one total CN segment:");
    verbose && print(verbose, segs[[kk]]);

    verbose && exit(verbose);    
  } # for (kk ...)

  segs <- Reduce(rbind, segs);
  rownames(segs) <- NULL;

  # Move 'chrom' column to the first column
  idx <- match("chromosome", names(segs));
  idxs <- c(idx, seq(length=ncol(segs))[-idx]);
  segs <- segs[,idxs,drop=FALSE];
  verbose && print(verbose, segs);

  verbose && enter(verbose, "Calculating (C1,C2) per segment");
  # Append (C1,C2) estimates
  tcn <- segs$tcn.mean;
  dh <- segs$dh.mean;
  C1 <- 1/2*(1-dh)*tcn;
  C2 <- tcn - C1;
  segs <- cbind(segs, c1.mean=C1, c2.mean=C2);
  verbose && exit(verbose);

  nbrOfSegs <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegs);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create result object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  params <- list(
    alphaTCN = alphaTCN,
    alphaDH = alphaDH,
    flavor = flavor,
    tbn = tbn,
    joinSegments = joinSegments,
    knownCPs = knownCPs,
    seed = seed
  );

  data <- data.frame(
    CT=CT,
    betaT=betaT,
    betaTN=betaTN,
    betaN=betaN,
    muN=muN,
    chromosome=chromosome,
    x=x
  );
  # Should we drop attributes? /HB 2010-09-24
  class(data) <- c("PairedPSCNData", class(data));

  class(segs) <- c("PairedPSCNSegments", class(segs));

  fit <- list(
    data = data,
    output = segs,
    params = params
  );

  if (!is.null(tcnLociNotPartOfSegment)) {
    fit$tcnLociNotPartOfSegment <- tcnLociNotPartOfSegment;
  }

  if (any(sapply(listOfDhLociNotPartOfSegment, FUN=length) > 0)) {
    fit$listOfDhLociNotPartOfSegment <- listOfDhLociNotPartOfSegment;
  }

  class(fit) <- c("PairedPSCBS", "PSCBS");

  # Update 
  if (is.element(flavor, c("tcn&dh", "sqrt(tcn)&dh"))) {
    fit <- postsegmentTCN(fit, verbose=verbose);
  }

  verbose && print(verbose, head(as.data.frame(fit)));
  verbose && print(verbose, tail(as.data.frame(fit)));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit;
}) # segmentByPairedPSCBS()


############################################################################
# HISTORY:
# 2010-11-28
# o BUG FIX: Iff argument 'chromosome' to segmentByPairedPSCBS() was of
#   length greater than one and specified exactly one unique chromosome,
#   then exception "Number of elements in argument 'chromosome' should
#   be exactly 8712 not 86209 value(s)" would be thrown.
# 2010-11-27
# o BUG FIX: segmentByPairedPSCBS() would not accept missing values in
#   argument 'chromosome'.
# o Now arguments '...' of segmentByPairedPSCBS() are passed to
#   the two segmentByCBS() calls.
# 2010-11-22
# o BUG FIX: segmentByPairedPSCBS() would not subset the correct set of
#   DH signals if there were some missing values in TCN.
# 2010-11-21
# o Changed the default to flavor="tch&dh".
# o Added support for flavors "tcn&dh", which, contrary to "tcn,dh",
#   enforces TCN and DH to have the same change points.
# o Now segmentByPairedPSCBS() also returns minor and major copy numbers
#   for each segment.
# o Forgot to return arguments 'joinSegments' & 'knownCPs' in 'params'.
# 2010-11-20
# o Now it is possible to specify the boundaries of the regions to be
#   segmented as known change points via argument 'knownCPs'.
# o Added argument 'joinSegments' to segmentByPairedPSCBS() in order to 
#   specify if neighboring segments should be joined or not.
# o Now segmentByCBS() allows for unknown genomic positions.
# o Now segmentByCBS() allows also for missing total CN signals.
# 2010-11-16
# o BUG FIX: In the rare cases where two loci at the same positions are
#   split up into two neighboring segments, then segmentByPairedPSCBS()
#   would fail to infer which they were if and only if the loci were not
#   ordered along the genome.  This could happen with for instance
#   Affymetrix GenomeWideSNP_6 data.
# o DOCUMENTATION: Clarified the form of argument 'muN', and added
#   references to papers and cross links to more internal methods.
# 2010-11-04
# o BUG FIX: There was a stray/debug stop() statement left in  
#   segmentByPairedPSCBS() causing an "error" in the rare case 
#   when loci that have the same physical locations are split
#   into two different segments.
# 2010-11-02
# o Added arguments 'undoTCN' and 'undoDH' to segmentByPairedPSCBS().
# o BUG FIX: Arguments 'alphaTCN' and 'alphaDH' of segmentByPairedPSCBS() 
#   were not used when more than one chromosome were segmented.
# 2010-10-25
# o BUG FIX: Now the correct set of loci are extracted from each TCN
#   segment, in the rare case that two neighboring TCN segments have
#   the same end points.
# 2010-10-18
# o Added arguments 'alphaTCN' and 'alphaDH' to segmentByPairedPSCBS().
# o Now segmentByPairedPSCBS() can segment multiple chromosomes.
# 2010-10-17
# o Added argument 'tbn' to segmentByPairedPSCBS() specifying whether
#   TumorBoostNormalization should be applied or not.
# 2010-10-10
# o The default is now to segment TCN on the original scale, not the sqrt().
# o Added flavor "sqrt(tcn),dh", which is segments sqrt(TCN) and then DH,
#   as original proposed by ABO.
# 2010-10-03
# o CLEAN UP: Now segmentByPairedPSCBS() is making use of argument
#   'chromosome' of segmentByCBS().
# 2010-10-02
# o Argument 'chromosome' default to 0 and have to be a finite integer.
# 2010-09-24
# o Now the 'data' field returned is a data.frame (no longer a list).
# o Now the 'chromosome' field of the data field is expanded to have the
#   same number of elements as the other locus fields.
# 2010-09-18
# o Added argument 'chromosome' to segmentByPairedPSCBS(), which, if given,
#   adds a chromosome column to the data and segmentation results.
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
