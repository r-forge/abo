%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{plotTracks.PairedPSCBS}
\alias{plotTracks.PairedPSCBS}
\alias{PairedPSCBS.plotTracks}
\alias{plotTracks.PairedPSCBS}
\alias{plotTracks,PairedPSCBS-method}

 \alias{plotTracks}

 \title{Plots parental specific copy numbers along the genome}

 \description{
  Plots parental specific copy numbers along the genome for one or more chromosomes.
  It is possible to specify what type of tracks to plot.
  Each type of track is plotted in its own panel.
 }

 \usage{\method{plotTracks}{PairedPSCBS}(x, tracks=c("tcn", "dh", "tcn,c1,c2", "tcn,c1", "tcn,c2", "c1,c2", "betaN", "betaT", "betaTN")[1:3], scatter="*", calls=".*", pch=".", quantiles=c(0.05, 0.95), cex=1, grid=FALSE, xlim=NULL, Clim=c(0, 6), Blim=c(0, 1), xScale=1e-06, ..., add=FALSE, verbose=FALSE)}

 \arguments{
   \item{x}{A result object returned by \code{\link{segmentByPairedPSCBS}}().}
   \item{tracks}{A \code{\link[base]{character}} \code{\link[base]{vector}} specifying what types of tracks to plot.}
   \item{scatter}{A \code{\link[base]{character}} \code{\link[base]{vector}} specifying which of the tracks should
     have scatter plot.}
   \item{calls}{A \code{\link[base]{character}} \code{\link[base]{vector}} of regular expression identifying
     call labels to be highlighted in the panels.}
   \item{pch}{The type of points to use.}
   \item{xlim}{(Optional) The genomic range to plot.}
   \item{Clim}{The range of copy numbers.}
   \item{Blim}{The range of allele B fractions (BAFs) and
     decrease of heterozygosity (DHs).}
   \item{xScale}{The scale factor used for genomic positions.}
   \item{...}{Not used.}
   \item{add}{If \code{\link[base:logical]{TRUE}}, the panels plotted are added to the existing plot,
     otherwise a new plot is created.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns nothing.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}



\keyword{internal}
\keyword{methods}
\keyword{IO}
