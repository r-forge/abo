%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  segmentByPairedPSCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{segmentByPairedPSCBS}
\alias{segmentByPairedPSCBS.default}
\alias{segmentByPairedPSCBS}
 \alias{segmentByPairedPSCBS}

 \title{Segment total copy numbers and allele B fractions using the paired PSCBS method}

 \description{
  Segment total copy numbers and allele B fractions using the paired PSCBS method [1].
  This method requires matched normals.
  This is a low-level segmentation method.
  It is intended to be applied to one tumor-normal sample at the time.
 }

 \usage{\method{segmentByPairedPSCBS}{default}(CT, betaT, betaN, muN=NULL, chromosome=0, x=NULL, alphaTCN=0.009, alphaDH=0.001, undoTCN=Inf, undoDH=Inf, ..., flavor=c("tcn&dh", "tcn,dh", "sqrt(tcn),dh", "sqrt(tcn)&dh"), tbn=TRUE, joinSegments=TRUE, knownCPs=NULL, seed=NULL, verbose=FALSE)}

 \arguments{
   \item{CT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor total copy numbers to be segmented.}
   \item{betaT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor allele B fractions (BAFs).}
   \item{betaN}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J matched normal BAFs.}
   \item{muN}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genotype calls in
        \{0,1/2,1\} for AA, AB, and BB, respectively. If not given,
        they are estimated from the normal BAFs using
        \code{\link[aroma.light]{callNaiveGenotypes}} as described in [2].}
   \item{chromosome}{(Optional) An \code{\link[base]{integer}} scalar
       (or a \code{\link[base]{vector}} of length J contain a unique value).
       Only used for annotation purposes.}
   \item{x}{Optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genomic locations.
            If \code{\link[base]{NULL}}, index locations \code{1:J} are used.}
   \item{alphaTCN, alphaDH}{The significance levels for segmenting total
        copy numbers (TCNs) and decrease-in-heterozygosity signals (DHs),
        respectively.}
   \item{undoTCN, undoDH}{Non-negative \code{\link[base]{numeric}}s.  If less than +\code{\link[base:is.finite]{Inf}},
        then a cleanup of segmentions post segmentation is done.
        See argument \code{undo} of \code{\link{segmentByCBS}}() for more
        details.}
   \item{...}{Additional arguments passed to \code{\link{segmentByCBS}}().}
   \item{flavor}{A \code{\link[base]{character}} specifying what type of segmentation and
     calling algorithm to be used.}
   \item{tbn}{If \code{\link[base:logical]{TRUE}}, \code{betaT} is normalized before segmentation
     using the TumorBoost method [2], otherwise not.}
   \item{joinSegments}{If \code{\link[base:logical]{TRUE}}, there are no gaps between neighboring
     segments.
     If \code{\link[base:logical]{FALSE}}, the boundaries of a segment are defined by the support
     that the loci in the segments provides, i.e. there exist a locus
     at each end point of each segment.  This also means that there
     is a gap between any neighboring segments, unless the change point
     is in the middle of multiple loci with the same position.
     The latter is what \code{DNAcopy::segment()} returns.
   }
   \item{knownCPs}{Optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of known
     change point locations.}
   \item{seed}{An (optional) \code{\link[base]{integer}} specifying the random seed to be
     set before calling the segmentation method.  The random seed is
     set to its original state when exiting.  If \code{\link[base]{NULL}}, it is not set.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns the fit object.
 }

 \details{
   Internally \code{\link{segmentByCBS}}() is used for segmentation.
   This segmentation method does \emph{not} support weights.

   The "DNAcopy::segment" implementation of CBS uses approximation
   through random sampling for some estimates.  Because of this,
   repeated calls using the same signals may result in slightly
   different results, unless the random seed is fixed.
 }

 \section{Missing and non-finite values}{
   The total copy number signals as well as any optional positions
   must not contain missing values, i.e. \code{\link[base]{NA}}s or \code{\link[base:is.finite]{NaN}}s.
   If there are any, an informative error is thrown.
   Allele B fractions may contain missing values, because such are
   interpreted as representing non-polymorphic loci.

   None of the input signals may have infinite values, i.e. -\code{\link[base:is.finite]{Inf}} or \code{\link[base:is.finite]{Inf}}.
   If so, an informative error is thrown.
 }

 \examples{
library("R.utils");
library("aroma.light");


pathname <- system.file("data-ex/pcCBS,exData,chr01.Rbin", package="psCBS");
data <- R.utils::loadObject(pathname);

# Subset to speed up example
#data0 <- data;
#keep <- seq(from=1, to=nrow(data), by=5);
#data <- data[keep,,drop=FALSE];
# Order
o <- order(data$position);
data <- data[o,];
str(data);
R.oo::attachLocally(data);
x <- position;

# Drop single-locus outliers
CTs <- dropSegmentationOutliers(CT, chromosome=1, x=x, verbose=-10);

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CTs, betaT=betaT, betaN=betaN,
                            chromosome=1, x=x, verbose=-10);
print(fit);

# Plot results
plotTracks(fit);

# Estimate the below tauAB tuning parameter for call AB regions
tauAB <- estimateTauAB(fit, scale=2);
print(tauAB);
## [1] 0.1354

# Estimate "normal contamination"
kappa <- estimateKappa(fit);
print(kappa);
## [1] 0.2179

fit2 <- bootstrapTCNandDHByRegion(fit, verbose=-10);
print(fit2);
plotTracks(fit2);

fit3 <- callABandHighAI(fit2, tauAB=tauAB, verbose=-10);
print(fit3);
plotTracks(fit3);


## # Call allelic balance and extreme allelic imbalance ("LOH")
## cfit <- callABandHighAI(fit, tauAB=0.10, alphaAB=0.05,
##                   tauHighAI=0.60, alphaHighAI=0.05, verbose=-10);
##
## # Display segmentation and calls
## print(fit);
}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   Internally, \code{\link[aroma.light]{callNaiveGenotypes}} is used to
   call naive genotypes, \code{\link[aroma.light]{normalizeTumorBoost}} is
   used for TumorBoost normalization, and \code{\link{segmentByCBS}}() is used
   to segment TCN and DH separately.
 }

 \references{
   [1] A. Olshen, R. Olshen, H. Bengtsson, P. Neuvial, P. Spellman and P. Seshan, \emph{Parent-specific copy number in paired tumor-normal studies using circular binary segmentation}, 2010 (submitted).\cr
   [2] H. Bengtsson, P. Neuvial and T.P. Speed, \emph{TumorBoost: Normalization of allele-specific tumor copy numbers from a single pair of tumor-normal genotyping microarrays}, BMC Bioinformatics, 2010.
 }


\keyword{IO}
