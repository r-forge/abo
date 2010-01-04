%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  as.segment.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{as.segment}
\alias{as.segment}


 \title{Allele-specific segmentation and calling using CBS}

 \description{
  Allele-specific segmentation and calling using CBS.
 }

 \usage{as.segment(x, matching.reference=FALSE, alpha=0.01, alpha1=0.009, het.lower=0.1, het.upper=0.75, het.minimum.window=0.05, het.stepsize.window=0.01, het.stepsize.within=0.01, trim.mean=0.1, min.hetero=10, cutoff.equality=0.1, homozygous.pvalue.cutoff=0.05, alpha.equality=0.05, maf=0.025, p0=0.5, min.homo.region=100, het.bias=0.04, smooth.segmentation=TRUE, impute.LOH, zero.homo=FALSE, verbose=1, ...)}

 \arguments{
   \item{x}{An \code{\link{as.CNA}}() object.}
   \item{matching.reference}{...is whether reference data are paired; can use single or matrix of unpaired reference.}
   \item{alpha}{}
   \item{alpha1}{}
   \item{het.lower, het.upper}{The range to look for the cutoff of the minimum value to call a SNP heterozygous; for crma v2 they need to be increased.}
   \item{het.minimum.window}{}
   \item{het.stepsize.window}{}
   \item{het.stepsize.within}{}
   \item{trim.mean}{}
   \item{min.hetero}{The minimum number of heterozyzgous markers to estimate parental copy number.}
   \item{cutoff.equality}{}
   \item{homozygous.pvalue.cutoff}{}
   \item{alpha.equality}{}
   \item{maf}{The minor allele frequency used when calculating whether region is LOH; increasing it will lead to fewer LOH region.}
   \item{p0}{}
   \item{min.homo.region}{}
   \item{het.bias}{Tuning parameter for deciding whether the parental copy numbers are unequal in the second round of segmentation; increasing it will lead to fewer regions called unequal.}
   \item{smooth.segmentation}{}
   \item{impute.LOH}{... has to do with how regions are called LOH.}
   \item{zero.homo}{}
   \item{verbose}{}
   \item{...}{Arguments passed to (CBS) \code{\link[DNAcopy]{segment}}()
              of \pkg{DNAcopy}.}
 }

 \value{
   Returns a \code{\link[base]{data.frame}} with the following variable:
   \item{ID, chrom, loc. start, loc.end, num.mark}{...are the same as CBS.}
   \item{num.hetero}{The number of heterozygous markers in the segment.}
   \item{mean.diff.mbar}{The average difference in the mirrored BAF between the tumor segment and the reference; higher implies unequal parental copy number.}
   \item{mindiff.mean}{The lower parental copy number assuming they parental copy numbers are unequal.}
   \item{maxdiff.mean}{The higher parental copy number assuming they parental copy numbers are unequal.}
   \item{min.mean}{The estimated lower parental copy number.}
   \item{max.mean}{The estimated higher parental copy number.}
   \item{total.mean}{The estimated total copy number.}
 }

 \author{Adam Olshen <olshena@biostat.ucsf.edu>}

 \seealso{
   \code{\link{as.CNA}}().
 }

