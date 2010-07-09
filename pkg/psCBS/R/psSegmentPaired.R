###########################################################################/**
# @set "class=psCNA"
# @RdocMethod "psSegmentPaired"
# 
# @title "Parent specific DNA copy number segmentation algorithm"
# 
# \description{
#   This function segments allele specific copy number data into parent
#   specific DNA copy number segments using CBS (segment function from
#   DNAcopy) for initial segmentation of Total DNA copy numbers and
#   further segmenting when parental copy numbers change.
# }
# 
# @synopsis
# 
# \arguments{
#   \item{x}{an object of class @see "psCNA".}
#   \item{alpha}{overall significance level for the test to accept
#     change-points.}
#   \item{alpha1}{significance level for test on total copy number.  Must
#     be less than alpha.}
#   \item{het.lower}{lowest possible threshold on minor copy number to
#     distinguish homozygotes from heterozygotes.}
#   \item{het.upper}{highest possible threshold on minor copy number to
#     distinguish homozygotes from heterozygotes.}
#   \item{het.minimum.window}{smallest window to find threshold on
#     genotypes.} 
#   \item{het.stepsize.window}{how much to increment the overall window
#     size to find threshold on genotypes.}
#   \item{het.stepsize.within}{how much to increment the current window
#     size to find threshold on genotypes.}
#   \item{trim.mean}{how much to trim when estimating parent-specific
#     means.}
#   \item{min.hetero}{minimum number of heterozygotes to estimate
#     parent-specific means.}
#   \item{alpha.homozygous}{significance level to test whether region is
#     LOH.}
#   \item{alpha.equality}{significance level to test whether
#     parent-specific means are equal}
#   \item{maf}{the assumed minor allele frequency used when estimating
#     whether region is LOH; decreasing it will lead to more LOH regions.}
#   \item{min.homo.region}{the smallest number of SNPs where no attempt
#     made to combine regions when testing heterozygosity.}
#   \item{modeOffset}{a tuning parameter for estimating whether the parental
#     copy numbers are unequal in the second round of segmentation;
#     increasing it will lead to fewer regions estimated to be unequal.}
#   \item{smooth.segmentation}{whether to smooth before segmenting on
#     total copy number.}
#   \item{smooth.region}{number of points to consider on the left and the
#     right of a point to detect it as an outlier.  Note this is a
#     parameter in smooth.CNA.}
#   \item{outlier.SD.scale}{the number of SDs away from the nearest point
#     in the smoothing region to call a point an outlier.  Note this is a
#     parameter in smooth.CNA.}
#   \item{smooth.SD.scale}{the number of SDs from the median in the
#     smoothing region where a smoothed point is positioned.}
#   \item{trim}{proportion of data to be trimmed for variance calculation
#     for smoothing outliers and undoing splits based on SD.}
#   \item{impute.LOH}{whether to combine small regions with testing for
#     LOH.}
#   \item{zero.homo}{whether to make the minor copy number 0 for
#     homozygotes. Default is @FALSE.}
#   \item{verbose}{verbosity of CBS}
#   \item{...}{other arguments to be passed on to @see "DNAcopy::segment"}
# }
#
# \value{
#   A @data.frame with 13 columns that give the following characteristics
#   of each segment: "ID", "chrom", "loc.start", "loc.end", "num.mark",
#   "num.hetero", "mean.diff.mbaf", "mindiff.mean", "maxdiff.mean",
#   "min.mean", "max.mean", "total.mean"
# }
#
# \details{
#   This function implements the parent specific segmentation. It segments
#   the total copy number.  Heterozygotes are found.  Then within each
#   segment the allele specific copy number data for the heterozygotes are
#   used to determine further segmentation is needed to allow for matched
#   parent specific changes occur i.e. a gain in a parental chromosome is
#   offset by a matching loss in the other keeping the total copy number
#   constant.  After the segmentation is complete segments are called for
#   LOH, equal parental copies or unequal copies.  
# }
#
# @examples "../incl/psSegmentPaired.Rex"
#
# @keyword nonparametric
#*/###########################################################################
setMethodS3("psSegmentPaired", "psCNA", function(x, alpha=0.01, alpha1=0.009, het.lower=0.1, het.upper=0.75, het.minimum.window=0.05, het.stepsize.window=0.01, het.stepsize.within=0.01, trim.mean=0.1, min.hetero=5, alpha.homozygous=0.05, alpha.equality=0.05, maf=0.075, min.homo.region=100, modeOffset=0.025, smooth.segmentation=TRUE, smooth.region=2, outlier.SD.scale=4, smooth.SD.scale=2, trim=0.025, impute.LOH, zero.homo=FALSE, verbose=1, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  findHeterozygousSNPs <- function(yTA, yTB) {
    res <- findheterozygous(yTA, yTB, het.lower, het.upper, het.minimum.window, het.stepsize.window, het.stepsize.within);
    res$is.heterozygous;
  } # findHeterozygousSNPs()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSamples <- nbrOfSamples(x);
  nbrOfRefSamples <- nbrOfReferenceSamples(x);
  if (nbrOfRefSamples != nbrOfSamples) {
    stop("The number of test and reference samples must be the same: ", nbrOfSamples, " != ", nbrOfRefSamples)
  }

  # Argument 'alpha' and 'alpha1':
  if (alpha1 >= alpha) {
    stop("alpha1 must be less than alpha: ", alpha1, " >= ", alpha)
  }



  call <- match.call()
  alpha2 <- alpha-alpha1
  sampleid <- getSampleNames(x);
  nbrOfLoci <- nbrOfLoci(x);

  if (is.null(x$genomdat.hetmatrix)) {
    predict.genomdat.heterozygous <- TRUE;
  } else {
    predict.genomdat.heterozygous <- FALSE;
  }

  if (is.null(x$normaldat.hetmatrix)) {
    predict.normal.heterozygous <- TRUE;
  } else {
    predict.normal.heterozygous <- FALSE;
  }

  # Extract (chromosome, position)
  chr <- x$chrom;
  pos <- x$maploc;

  # Find out which loci has a known location
  hasPosition <- (!is.na(chr) & !is.na(pos));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Segment each sample independently
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- NULL;

  for (ii in seq(length=nbrOfSamples)) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Extract tumor signals
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Extract (yA,yB)
    yTA <- x$genomdat.a[,ii];
    yTB <- x$genomdat.b[,ii];

    # Calculate allele B fractions (beta)
    yT <- yTA + yTB;
    betaT <- yTB / yT;

    betaT[abs(betaT) == Inf] <- 1;  # Does this happen? /HB 2010-07-08
     
    # Calculate mirrored allele B fractions (rho)
    ## SLOW: rhoT <- apply(cbind(betaT, 1-betaT), MARGIN=1, FUN=max)
    rhoT <- abs(betaT - 1/2) + 1/2;
    rhoT[abs(rhoT) == Inf] <- 1;  # Does this happen? /HB 2010-07-08

    # Check which loci has valid data
    hasBetaT <- (!is.na(betaT));
    okT <- (hasPosition & hasBetaT);

    isHetsN <- rep(NA, length=nbrOfLoci);

    if (predict.normal.heterozygous) {
      # Extract (yA,yB)
      yNA <- x$normaldat.a[,ii];
      yNB <- x$normaldat.b[,ii];

      # Calculate allele B fractions (beta)
      yN <- yNA + yNB;
      betaN <- yNB/yN;
      hasBetaN <- (!is.na(betaN));

      # Calculate mirrored allele B fractions (rho)
      ## SLOW: rhoN <- apply(cbind(betaN, 1-betaN), MARGIN=1, FUN=max)
      rhoN <- abs(betaN - 1/2) + 1/2;
      rhoN[abs(rhoN) == Inf] <- 1;  # Does this happen? /HB 2010-07-08

      ok <- which(okT & hasBetaN);
      isHetsN[ok] <- findHeterozygousSNPs(yNA[ok], yNB[ok]);
      normal.het.vector.notmissing <- isHetsN
    } else {
      normal.het.vector <- x$normaldat.hetmatrix[,ii]
      ok <- which(okT & !is.na(normal.het.vector))
      # ok <- which(okT & !is.na(normal.het.vector) & !is.na(genomdat.het.vector))                
      normal.het.vector.notmissing <- normal.het.vector[ok]
      # stop("Code not yet written for normal genotyping data; must use normal copy number data")
    } # if (predict.normal.heterozygous)


    if (predict.normal.heterozygous) {
      # TumorBoost normalization
      idxs <- which(betaT > betaN);
      betaT[ok] <- 0.5*betaT[ok]/betaN[ok];
      betaT[ok][idxs] <- (1-0.5*((1-betaT[ok])/(1-betaN[ok])))[idxs];

      diff.mbaf.vector.notmissing <- rhoT[ok]-rhoN[ok];
    } else {
      diff.mbaf.vector.notmissing <- rep(NA, times=length(betaT[ok]));
    }


    yTA.notmissing <- yTA[ok];
    yTB.notmissing <- yTB[ok];

    if (predict.genomdat.heterozygous) {
      isHets <- findHeterozygousSNPs(yTA[ok], yTB[ok]);
      het.vector.notmissing <- isHets
    } else {
      het.vector <- x$genomdat.hetmatrix[,ii]
      het.vector.notmissing <- het.vector[ok]
#     stop("Code not yet written for test genotyping data; must use test copy number data")
    } # if (predict.genomdat.heterozygous)

    # Make the smaller value of all homozygotes be zero if desired
    if (zero.homo) {
      min.indicator <- rep(FALSE, times=length(yTA[ok]));
      min.indicator[yTB[ok] < yTA[ok]] <- TRUE;
      yTA.notmissing[!het.vector.notmissing & !min.indicator] <- 0;
      yTB.notmissing[!het.vector.notmissing &  min.indicator] <- 0;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identification of change points in total copy numbers
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calculate yT = yTA + yTB
    yT.notmissing <- yTA.notmissing + yTB.notmissing;
    # Must make negatives zero else thrown out, but negatives come back later
    yT.notmissing[yT.notmissing < 0] <- 0
    sqrt.yT.notmissing <- sqrt(yT.notmissing)
    CNA.sqrt.yT.notmissing <- CNA(sqrt.yT.notmissing, chr[ok], pos[ok], sampleid=sampleid[ii])

    # By default smoothing is done
    if (smooth.segmentation) {
      smoothed.CNA.sqrt.yT.notmissing <- smooth.CNA(CNA.sqrt.yT.notmissing, smooth.region=smooth.region, outlier.SD.scale=outlier.SD.scale, smooth.SD.scale=smooth.SD.scale, trim=trim)
      # smoothed.CNA.sqrt.yT.notmissing <- smooth.CNA(CNA.sqrt.yT.notmissing, ...)            
      segmented.smoothed.CNA.sqrt.yT.notmissing <- segment(smoothed.CNA.sqrt.yT.notmissing, alpha=alpha1, verbose=verbose, ...)
    } else {
      segmented.smoothed.CNA.sqrt.yT.notmissing <- segment(CNA.sqrt.yT.notmissing, alpha=alpha1, verbose=verbose, ...)
    }

    segmented.output <- segmented.smoothed.CNA.sqrt.yT.notmissing$output
    segment.lengths <- as.numeric(segmented.output[,5])
    segment.ends <- cumsum(segment.lengths)
    n.segment <- length(segment.ends)
    segment.starts <- c(1, segment.ends[-n.segment]+1)
    chr.starts <- chr[ok][segment.starts]
    # first.LOH <- first.call(segment.starts, segment.ends, het.vector.notmissing, chr.starts, maf, homozygous.pvalue.cutoff)
    # print(first.LOH)
    # segmented.output.LOH <- cbind(segmented.output, first.LOH)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identification of additional change points using DH
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Second round of segmentation
    foldedTumorBoostBAFNotmissing <- 2*abs(betaT[ok]-0.5);
    second.output <- second.segment.matching(TRUE, foldedTumorBoostBAFNotmissing, normal.het.vector.notmissing, alpha2, segmented.output, pos[ok]);



    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Build table of PSCN segments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    keep <- c("ID", "chrom", "loc.start", "loc.end", "num.mark");
    headII <- second.output[,keep,drop=FALSE];
    headII <- as.data.frame(headII);
    headII$ID <- as.character(headII$ID);
    headII$chrom <- as.character(headII$chrom);
    headII$loc.start <- as.integer(as.character(headII$loc.start));
    headII$loc.end <- as.integer(as.character(headII$loc.end));
    headII$num.mark <- as.integer(as.character(headII$num.mark));

    segment.lengths <- as.numeric(second.output[,"num.mark"]);
    # print(paste("segment.lengths=", segment.lengths))

    segment.ends <- cumsum(segment.lengths);
    # print(paste("segment.ends=", segment.ends))

    n.segment <- length(segment.ends);
    # print(paste("n.segment=", n.segment))

    segment.starts <- c(1, segment.ends[-n.segment]+1);
    # print(paste("segement.starts=", segment.starts))

    chr.starts <- as.numeric(chr[ok])[segment.starts];
    # print(paste("chr.starts=", chr.starts))

    total.vector.starts <- as.numeric(second.output[,"seg.mean"]);
    # print(paste("total.vector.starts=", total.vector.starts))
    # print("Before called segments")

    called.segments <- callsegments.matching(matching=TRUE, segment.starts=segment.starts, segment.ends=segment.ends, a.vector=yTA.notmissing, b.vector=yTB.notmissing, diff.mbaf.vector=diff.mbaf.vector.notmissing, het.vector=het.vector.notmissing, normal.het.vector=normal.het.vector.notmissing, foldedTumorBoostBAF=foldedTumorBoostBAFNotmissing, chrom.vector.starts=chr.starts, total.vector.starts=total.vector.starts, min.hetero=min.hetero, trim.mean=trim.mean, alpha.equality=alpha.equality, modeOffset=modeOffset, maf=1-maf, alpha.homozygous=alpha.homozygous, min.homo.region=min.homo.region, impute.LOH=impute.LOH)

    keep <- c("n.hetero", "mean.diff.mbaf", "mindiff.output", "maxdiff.output", "confidence95NewModes", "min.output", "max.output", "total.output");
    resII <- called.segments[keep];
    names(resII) <- c("num.hetero", "mean.diff.mbaf", "mindiff.mean", "maxdiff.mean", "confidence95", "min.mean", "max.mean", "total.mean");

    resII <- cbind(headII, resII);

    if (is.null(res)) {
      res <- resII;
    } else {
      res <- rbind(res, resII);
    }
  } # for (ii in ...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Merge PSCN segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (nrow(res) >= 2) {
    kk <- 2;
    while (kk <= nrow(res)) {
      # Is region #kk same as region #(kk-1)...?
      rows <- c(kk-1,kk);
      data <- res[rows,c("ID", "chrom", "min.mean", "total.mean")];
      if (isTRUE(all.equal(data[1,], data[2,]))) {
        # ...then merge #kk into #(kk-1)
        res[kk-1,"num.mark"] <- res[kk-1,"num.mark"] + res[kk,"num.mark"];
        res[kk-1,"num.hetero"] <- res[kk-1,"num.hetero"] + res[kk,"num.hetero"];
        res[kk-1,"loc.start"] <- min(res[rows,"loc.end"]);
        res[kk-1,"loc.end"] <- max(res[rows,"loc.end"]);
        # ...and drop #kk
        res <- res[-kk,,drop=FALSE];
      } else {
        kk <- kk + 1;
      }
    } # for (kk ...)
  } # if (nrow(res) >= 2)

  res
}) # psSegmentPaired()



############################################################################
# HISTORY:
# 2010-07-08
# o Now psSegmentPaired() returns a data frame (no longer a matrix).
# o CLEANUP: Major cleanup, i.e. renaming variables, reordering etc.
# o CLEANUP: Created psSegmentPaired() from psSegment().
# o ROBUSTNESS: Now samples are index by 'ii' and regions by 'kk'.
# o ROBUSTNESS: Replaced all 1:n with seq(length=n) to deal with n == 0.
# o ROBUSTNESS: Now all interator variables i & j are written as ii & jj.
# o ROBUSTNESS: Now all list elements are referenced by name.
# o Made psSegment() an S3 method for the psCNA class.
############################################################################
