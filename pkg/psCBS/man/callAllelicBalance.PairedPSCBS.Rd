%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.callAllelicBalance.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{callAllelicBalance.PairedPSCBS}
\alias{callAllelicBalance.PairedPSCBS}
\alias{PairedPSCBS.callAllelicBalance}
\alias{callAllelicBalance.PairedPSCBS}
\alias{callAllelicBalance,PairedPSCBS-method}

\alias{PairedPSCBS.callAB}
\alias{callAB.PairedPSCBS}
\alias{callAB,PairedPSCBS-method}


\title{Calls segments that are in allelic balance}

\description{
 Calls segments that are in allelic balance, i.e. that have equal minor and major copy numbers.
}

\usage{\method{callAllelicBalance}{PairedPSCBS}(fit, flavor=c("DeltaAB*"), ...)}

\arguments{
  \item{flavor}{A \code{\link[base]{character}} string specifying which type of
   call to use.}
  \item{...}{Additional arguments passed to the caller.}
}

\value{
  Returns a \code{\link{PairedPSCBS}} object with allelic-balance calls.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally, one of the following methods are used:
  \code{\link[psCBS:callAllelicBalanceByDH.PairedPSCBS]{*callAllelicBalanceByDH}()}.
}
\keyword{internal}
\keyword{methods}
