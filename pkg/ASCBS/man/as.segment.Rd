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
   \item{x}{}
   \item{matching.reference}{...is whether reference data are paired; can use single or matrix of unpaired reference.}
   \item{alpha}{}
   \item{alpha1}{}
   \item{het.lower, het.upper}{...is the range to look for the cutoff of the minimum value to call a SNP heterozygous; for crma v2 they need to be increased.}
   \item{het.minimum.window}{}
   \item{het.stepsize.window}{}
   \item{het.stepsize.within}{}
   \item{trim.mean}{}
   \item{min.hetero}{...is the minimum number of heterozyzgous markers to estimate parental copy number.}
   \item{cutoff.equality}{}
   \item{homozygous.pvalue.cutoff}{}
   \item{alpha.equality}{}
   \item{maf}{...is the minor allele frequency used when calculating whether region is LOH; increasing it will lead to fewer LOH region.}
   \item{p0}{}
   \item{min.homo.region}{}
   \item{het.bias}{...is a tuning parameter for deciding whether the parental copy numbers are unequal in the second round of segmentation; increasing it will lead to fewer regions called unequal.}
   \item{smooth.segmentation}{}
   \item{impute.LOH}{... has to do with how regions are called LOH.}
   \item{zero.homo}{}
   \item{verbose}{}
   \item{...}{Arguments passed to \code{\link[DNAcopy]{segment}} (CBS).}
 }

 \value{
   Returns a \code{\link[base]{data.frame}} with the following variable:
ID,chrom,loc.start,loc.end,num.mark are the same as CBS
num.hetero is the number of heterozygous markers in the segment
mean.diff.mbar is the average difference in the mirrored BAF between the tumor segment and the reference; higher implies unequal parental copy number
mindiff.mean is the lower parental copy number assuming they parental copy numbers are unequal
maxdiff.mean is the higher parental copy number assuming they parental copy numbers are unequal
min.mean is the estimated lower parental copy number
max.mean is the estimated higher parental copy number
total.mean is the estimated total copy number
 }

 \author{Adam Olshen <olshena@biostat.ucsf.edu>}

 \seealso{
   \code{\link{as.CNA}}().
 }

