%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.estimateTauLOH.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{estimateTauLOHByMinC1ForNonAB.PairedPSCBS}
\alias{estimateTauLOHByMinC1ForNonAB.PairedPSCBS}
\alias{PairedPSCBS.estimateTauLOHByMinC1ForNonAB}
\alias{estimateTauLOHByMinC1ForNonAB.PairedPSCBS}
\alias{estimateTauLOHByMinC1ForNonAB,PairedPSCBS-method}


\title{Estimate a threshold for calling LOH from DH}

\description{
 Estimate a threshold for calling LOH from DH based on the location of guessed C1=0 and C1=1 peaks.
}

\usage{\method{estimateTauLOHByMinC1ForNonAB}{PairedPSCBS}(this, tauMax=1/3, ..., verbose=FALSE)}

\arguments{
  \item{tauMax}{A \code{\link[base]{numeric}} specifying the maximum estimate returned.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns the estimated LOH treshold as a \code{\link[base]{numeric}} scalar.
}

\details{
  This method requires that calls for allelic balances already have
  been me made, cf. \code{\link[psCBS:callAllelicBalance.PairedPSCBS]{*callAllelicBalance}()}.
}

\section{Algorithm}{
 \itemize{
  \item Grabs the segment-level C1 estimates.
  \item Calculate segment weights proportional to the number of heterozygous SNPs.
  \item Estimate the C1=1 location as the weighted median C1 for segments that have been called to be in allelic balance.
  \item Estimate the C1=0 location as the smallest C1 among segments that are not in allelic balance.
  \item Let the LOH threshold be the midpoint of the estimates C1=0 and C1=1 locations.
  \item If the estimated LOH threshold is greater than \code{tauMax}, set it to \code{tauMax}.
 }
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Instead of calling this method explicitly, it is recommended
  to use the \code{\link[psCBS:estimateTauLOH.PairedPSCBS]{*estimateTauLOH}()} method.
}


\keyword{internal}
\keyword{methods}
