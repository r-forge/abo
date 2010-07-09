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

 \title{Segment copy numbers using the paired PSCBS method}

 \description{
  Segment copy numbers using the paired PSCBS method.
  This is a low-level method.
  It is intended to be applied to one sample and one chromosome at the time.
 }

 \usage{\method{segmentByPairedPSCBS}{default}(CT, betaT, betaN, muN=NULL, x=NULL, ..., seed=NULL, verbose=FALSE)}

 \arguments{
   \item{CT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor total copy numbers to be segmented.}
   \item{betaT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor allele B fractions (BAFs).}
   \item{betaN}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J matched normal BAFs.}
   \item{muN}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genotype calls.
       If not given, they are estimated from the normal BAFs.}
   \item{x}{Optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genomic locations.
            If \code{\link[base]{NULL}}, index locations \code{1:J} are used.}
   \item{...}{Not used.}
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
   different results.
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
segs <- fit$output;
segs <- segs[,-(1:2)];
print(segs);

# Plotting data
xMb <- x/1e6;
xlab <- "Position (Mb)";
pch <- ".";

muN <- callNaiveGenotypes(betaN);
isHet <- (muN == 1/2);
betaTN <- normalizeTumorBoost(betaT, betaN, muN=muN);
rhoTN <- 2*abs(betaTN-1/2);
cps <- sort(c(segs$loc.start, segs$loc.end));

layout(matrix(1:5, ncol=1));
par(mar=c(2,4,1,1)+0.1);
plot(xMb,CT, pch=pch, col="#000000", ylim=c(0,5), xlab=xlab);
abline(v=cps/1e6, lwd=2);
drawLevels(fit, field="tcn.mean", xScale=1e-6);
col <- c("#aaaaaa", "#000000", "#aaaaaa")[2*muN+1];
plot(xMb,betaN, pch=pch, col=col, ylim=c(0,1), xlab=xlab);
abline(v=cps/1e6, lwd=2);
plot(xMb,betaT, pch=pch, col=col, ylim=c(0,1), xlab=xlab);
abline(v=cps/1e6, lwd=2);
plot(xMb,betaTN, pch=pch, col=col, ylim=c(0,1), xlab=xlab);
abline(v=cps/1e6, lwd=2);
rho <- rhoTN[isHet];
plot(xMb[isHet],rho, pch=pch, col="#000000", ylim=c(0,1), xlab=xlab);
abline(v=cps/1e6, lwd=2);
drawLevels(fit, field="dh.mean", xScale=1e-6);
}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   Internally \code{\link{segmentByCBS}}() is used for segmentation.
 }


\keyword{IO}
