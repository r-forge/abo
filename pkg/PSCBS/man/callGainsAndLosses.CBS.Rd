%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.CALL.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{callGainsAndLosses.CBS}
\alias{callGainsAndLosses.CBS}
\alias{CBS.callGainsAndLosses}
\alias{callGainsAndLosses.CBS}
\alias{callGainsAndLosses,CBS-method}


 \title{Calls gains and losses}

 \description{
  Calls gains and losses.
 }

 \usage{\method{callGainsAndLosses}{CBS}(fit, adjust=1, method=c("ucsf-mad"), ...)}

 \arguments{
  \item{adjust}{A positive scale factor adjusting the sensitivity of the
    caller, where a value less (greater) than 1.0 makes the caller
    less (more) sensitive.}
  \item{method}{A \code{\link[base]{character}} string specifying the calling algorithm to use.}
  \item{...}{Additional/optional arguments used to override the default
    parameters used by the caller.}
 }

 \value{
  Returns a \code{\link[PSCBS]{CBS}} object where \code{\link[base]{logical}} columns
  'lossCall' and 'gainCall' have been appended to the segmentation table.
 }

 \section{The UCSF caller}{
   "Previously, we defined a segment to be gained/lost if it were at
    least two times the sample MAD away from the median segmented value."
    (D. Albertson 2011-07-18)
 }

 \examples{
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xBEEF)

# Number of loci
J <- 1000

mu <- double(J)
mu[200:300] <- mu[200:300] + 1
mu[650:800] <- mu[650:800] - 1
eps <- rnorm(J, sd=1/2)
y <- mu + eps
x <- sort(runif(length(y), max=length(y))) * 1e5
w <- runif(J)
w[650:800] <- 0.001


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- segmentByCBS(y, x=x)
print(fit)

  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CALLS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Call gains and losses by segments
fitC <- callGainsAndLosses(fit)

# Call amplifications by segments
fitC <- callAmplifications(fitC)

# Call outliers by loci
fitC <- callOutliers(fitC)

 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[PSCBS:callAmplifications.CBS]{*callAmplifications}()}.
   \code{\link[PSCBS:callOutliers.CBS]{*callOutliers}()}.
   For more information see \code{\link{CBS}}.
 }


\keyword{internal}
\keyword{methods}