%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.callLOH.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{callLOH.PairedPSCBS}
\alias{callLOH.PairedPSCBS}
\alias{PairedPSCBS.callLOH}
\alias{callLOH.PairedPSCBS}
\alias{callLOH,PairedPSCBS-method}


\title{Calls segments that are in LOH}

\description{
 Calls segments that are in LOH, i.e. that have "zero" minor copy number.
}

\usage{\method{callLOH}{PairedPSCBS}(fit, flavor=c("SmallC1", "LargeDH"), ..., xorCalls=TRUE, force=FALSE)}

\arguments{
  \item{flavor}{A \code{\link[base]{character}} string specifying which type of
   call to use.}
  \item{...}{Additional arguments passed to the caller.}
  \item{xorCalls}{If \code{\link[base:logical]{TRUE}}, a region already called AB, will
   not be called LOH.}
  \item{force}{If \code{\link[base:logical]{FALSE}}, and allelic-balance calls already exits,
   then nothing is done, otherwise the calls are done.}
}

\value{
  Returns a \code{\link{PairedPSCBS}} object with LOH calls.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally, one of the following methods are used:
  \code{\link[psCBS:callLowC1ByC1.PairedPSCBS]{*callLowC1ByC1}()},
  \code{\link[psCBS:callExtremeAllelicImbalanceByDH.PairedPSCBS]{*callExtremeAllelicImbalanceByDH}()}.
}
\keyword{internal}
\keyword{methods}
