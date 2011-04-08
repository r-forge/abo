%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.estimateTauAB.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{estimateTauAB.PairedPSCBS}
\alias{estimateTauAB.PairedPSCBS}
\alias{PairedPSCBS.estimateTauAB}
\alias{estimateTauAB.PairedPSCBS}
\alias{estimateTauAB,PairedPSCBS-method}


\title{Estimate a threshold for calling allelic balance from DH}

\description{
 Estimate a threshold for calling allelic balance from DH.
}

\usage{\method{estimateTauAB}{PairedPSCBS}(this, scale=NULL, flavor=c("q(DH)", "mad(hBAF)", "median(DH)"), ..., verbose=FALSE)}

\arguments{
  \item{scale}{An optional \code{\link[base]{numeric}} scale factor.}
  \item{flavor}{A \code{\link[base]{character}} string specifying which type of
   estimator to use.}
  \item{...}{Additional arguments passed to the estimator.}
}

\value{
  Returns the threshold estimate as a \code{\link[base]{numeric}} scalar.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally, one of the following methods are used:
  \code{\link[psCBS:estimateStdDevForHeterozygousBAF.PairedPSCBS]{*estimateStdDevForHeterozygousBAF}()},
  \code{\link[psCBS:estimateMeanForDH.PairedPSCBS]{*estimateMeanForDH}()}, and
  \code{\link[psCBS:estimateHighDHQuantileAtAB.PairedPSCBS]{*estimateHighDHQuantileAtAB}()}.
}
\keyword{internal}
\keyword{methods}
