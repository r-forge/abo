%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  segmentByCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{segmentByCBS}
\alias{segmentByCBS.default}
\alias{segmentByCBS}

 \title{Segment genomic signals using the CBS method}

 \description{
  Segment genomic signals using the CBS method of the \pkg{DNAcopy} package.
  This is a convenient low-level wrapper for the \code{DNAcopy::segment()}
  method.  It is intended to be applied to one sample and one chromosome
  at the time.
  For more details on the Circular Binary Segmentation (CBS) method
  see [1,2].
 }

 \usage{\method{segmentByCBS}{default}(y, x=NULL, w=NULL, ..., knownCPs=NULL, seed=NULL, verbose=FALSE)}

 \arguments{
   \item{y}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genomic signals to be segmented.}
   \item{x}{Optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genomic locations.
            If \code{\link[base]{NULL}}, index locations \code{1:J} are used.}
   \item{w}{Optional \code{\link[base]{numeric}} \code{\link[base]{vector}} in [0,1] of J weights.}
   \item{...}{Additional arguments passed to the segmentation function.}
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
   Internally \code{\link[DNAcopy]{segment}} is used to segment the signals.
   This segmentation method support weighted segmentation.

   The "DNAcopy::segment" implementation of CBS uses approximation
   through random sampling for some estimates.  Because of this,
   repeated calls using the same signals may result in slightly
   different results.
 }

 \section{Missing and non-finite values}{
   The signals to be segmented as well as any optional positions
   must not contain missing values, i.e. \code{\link[base]{NA}}s or \code{\link[base:is.finite]{NaN}}s.
   If there are any, an informative error is thrown.
   Furthermore, none of the input signals may have infinite values,
   i.e. -\code{\link[base:is.finite]{Inf}} or \code{\link[base:is.finite]{Inf}}. If so, an informative error is thrown.
 }

 \examples{
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Number of loci
J <- 1000;

mu <- double(J);
mu[200:300] <- mu[200:300] + 1;
mu[650:800] <- mu[650:800] - 1;
eps <- rnorm(J, sd=1/2);
y <- mu + eps;
x <- sort(runif(length(y), max=length(y))) * 1e5;
w <- runif(J);
w[650:800] <- 0.001;


xlab <- "Position (Mb)";
ylim <- c(-3,3);
xMb <- x/1e6;
plot(xMb,y, pch=20, col="#aaaaaa", xlab=xlab, ylim=ylim);

fit <- segmentByCBS(y, x=x);
print(fit);

drawLevels(fit, col="red", lwd=2, xScale=1e-6);
}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
  [1] A.B. Olshen, E.S. Venkatraman (aka Venkatraman E. Seshan),
      R. Lucito and M. Wigler, \emph{Circular binary segmentation for
      the analysis of array-based DNA copy number data},
      Biostatistics, 2004.\cr
  [2] E.S. Venkatraman and A.B. Olshen, \emph{A faster circular binary
      segmentation algorithm for the analysis of array CGH data}.
      Bioinformatics, 2007.\cr
 }


\keyword{IO}
