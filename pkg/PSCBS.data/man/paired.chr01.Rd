%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  paired.chr01.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{paired.chr01}
\alias{paired.chr01}
\docType{data}


\title{Matched tumor-normal SNP microarray data}

\description{
 Matched tumor-normal SNP microarray data on chromosome 1 (73,346 loci) for a single anonymous individual.
}

\usage{paired.chr01}

\format{
 A \code{\link[base]{data.frame}} with 73,346 rows (one per SNP) and columns
 'chromosome', 'x' 'CT', 'betaT', 'CN' and 'betaN'.
}

\examples{
  library("PSCBS.data")
  data("paired.chr01")
  str(paired.chr01)
}


\references{
 [1] A.B. Olshen, H. Bengtsson, P. Neuvial, P.T. Spellman, R.A. Olshen, V.E. Seshan, \emph{Parent-specific copy number in paired tumor-normal studies using circular binary segmentation}, Bioinformatics, 2011

}

\author{Henrik Bengtsson}
\keyword{datasets}
