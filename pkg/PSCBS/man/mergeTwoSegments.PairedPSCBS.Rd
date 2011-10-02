%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.EXTS3.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{mergeTwoSegments.PairedPSCBS}
\alias{mergeTwoSegments.PairedPSCBS}
\alias{PairedPSCBS.mergeTwoSegments}
\alias{mergeTwoSegments.PairedPSCBS}
\alias{mergeTwoSegments,PairedPSCBS-method}


 \title{Merge two neighboring segments}

 \description{
   Merge two neighboring segments by recalculating segment statistics.
 }

 \usage{\method{mergeTwoSegments}{PairedPSCBS}(this, left, verbose=FALSE, ...)}

 \arguments{
  \item{left}{An \code{\link[base]{integer}} specifying the segments (left, left+1)
    to be merged.}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
  \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link{PairedPSCBS}} with one less segment.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   To drop regions (a connected set of segments)
   see \code{\link[PSCBS:dropByRegions.PairedPSCBS]{*dropByRegions}()}.
   For more information see \code{\link{PairedPSCBS}}.
 }
\keyword{internal}
\keyword{methods}