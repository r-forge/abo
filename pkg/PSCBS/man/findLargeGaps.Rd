%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findLargeGaps.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{findLargeGaps}
\alias{findLargeGaps.default}
\alias{findLargeGaps}
 \alias{findLargeGaps.data.frame}

 \title{Identifies gaps of a genome where there exist no observations}

 \description{
  Identifies gaps of a genome where there exist no observations.
 }

 \usage{\method{findLargeGaps}{default}(chromosome=NULL, x, minLength, resolution=1, ...)}

 \arguments{
   \item{chromosome}{(Optional) An \code{\link[base]{integer}} \code{\link[base]{vector}} of length J of
     chromosome indices.}
   \item{x}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J of genomic locations.}
   \item{minLength}{A positive \code{\link[base]{numeric}} scalar specifying the minimum
     length of a gap.}
   \item{resolution}{A non-negative \code{\link[base]{numeric}} specifying the minimum
     length unit, which by default equals one nucleotide/base pair.}
   \item{...}{Not used.}
 }

 \value{
   Returns \code{\link[base]{data.frame}} with columns \code{chromosome} (if given),
   \code{start}, \code{stop} and \code{length}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link{gapsToSegments.data.frame}}().
 }



\keyword{IO}
\keyword{internal}