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

 \title{Segment copy numbers using the paired PSCBS method}

 \description{
  Segment copy numbers using the paired PSCBS method.
  This is a low-level segmentation method.
  It is intended to be applied to one sample and one chromosome at the time.
 }

 \usage{\method{segmentByPairedPSCBS}{default}(CT, betaT, betaN, muN=NULL, x=NULL, ..., flavor=c("dh|tcn", "tcn|dh", "tcn&dh"), seed=NULL, verbose=FALSE)}

 \arguments{
   \item{CT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor total copy numbers to be segmented.}
   \item{betaT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor allele B fractions (BAFs).}
   \item{betaN}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J matched normal BAFs.}
   \item{muN}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genotype calls.
       If not given, they are estimated from the normal BAFs.}
   \item{x}{Optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genomic locations.
            If \code{\link[base]{NULL}}, index locations \code{1:J} are used.}
   \item{...}{Not used.}
   \item{flavor}{A \code{\link[base]{character}} specifying what type of segmentation and
     calling algorithm to be used.}
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

   None of the input signals can have infinite values, i.e. -\code{\link[base:is.finite]{Inf}} or \code{\link[base:is.finite]{Inf}}.
   If so, an informative error is thrown.
 }

 \examples{
library("R.utils");
library("aroma.light");

pathname <- system.file("data-ex/pcCBS,exData,chr01.Rbin", package="psCBS");
data <- R.utils::loadObject(pathname);

# Subset to speed up example
data0 <- data;
#keep <- seq(from=1, to=nrow(data), by=5);
#data <- data[keep,,drop=FALSE];
# Order
o <- order(data$position);
data <- data[o,];
str(data);
R.oo::attachLocally(data);
x <- position;

# PSCBS segmentation
fit <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN, x=x, verbose=-10);

# Extract segments
segs <- as.data.frame(fit$output);
print(segs);

# Plot results
plot(fit);

}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   Internally \code{\link{segmentByCBS}}() is used for segmentation.
 }


\keyword{IO}
