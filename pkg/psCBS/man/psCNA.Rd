\name{psCNA}
\alias{psCNA}
\alias{print.psCNA}

\title{Create `Parent-Specific Copy Number Array' data object}

\description{
  Creates a `parent-sprcific copy number array' data object used for
  Parent-Specific  DNA copy number analyses by psSegment.
}

\usage{
  psCNA(genomdat.a, genomdat.b, chrom, maploc, normaldat.a, normaldat.b,
      genomdat.hetmatrix=NULL, normaldat.hetmatrix=NULL, sampleid=NULL,
      arraytype=NULL) 
  \method{print}{psCNA}(x, \dots)
}

\arguments{
  \item{genomdat.a}{a vector or matrix of parental copy number data form
    a test SNP array. If it is a matrix the rows correspond to the
    markers and the columns to the samples.}
  \item{genomdat.b}{a second vector or matrix of test parental copy
    number data.}
  \item{chrom}{the chromosomes (or other group identifier) from which
    the markers came.  Vector of length same as the number of rows of
    genomdat.  If one wants the chromosomes to be ordered in the
    natural order, this variable should be numeric or ordered category.} 
  \item{maploc}{the locations of marker on the genome.  Vector of length
    same as the number of rows of genomdat. This has to be numeric.}
  \item{normaldat.a}{a vector or matrix of parental copy number data
    form a reference SNP array. If it is a matrix the rows correspond to
    the markers and the columns to the samples.}
  \item{normaldat.b}{a second vector or matrix of reference parental
    copy number data.}
  \item{genomdat.hetmatrix}{A vector or matrix of TRUE/FALSE as to
    whether a test SNP is heterozygous.  Must be the same size as
    genomdat.a.  Not yet implemented.  Currently, heterozygosity is
    determined in psSegment.}
  \item{normaldat.hetmatrix}{A vector or matrix of TRUE/FALSE as to
    whether a test SNP is heterozygous.  Must be the same size as
    normaldat.a.  Not yet implemented.  Currently, heterozygosity is
    determined in psSegment.} 
  \item{sampleid}{sample identifier.  If missing the samples are named
    by prefixing "Sample" to consecutive integers.}
  \item{arraytype}{The type of SNP array being analyzed.  Could be used
    to set some parameters.  Not yet implemented}
  \item{x}{object returned by CNA}
  \item{...}{arguments to be passed onto print command called within.}
}

\value{
  An object of class \code{psCNA}.  There is a \code{print} method that
  gives the number of samples and probes and the type of data.
}

\keyword{nonparametric}
