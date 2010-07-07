%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  psSegment.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{psSegment}
\alias{psSegment}


 \title{Parent specific DNA copy number segmentation algorithm}

 \description{
   This function segments allele specific copy number data into parent
   specific DNA copy number segments using CBS (segment function from
   DNAcopy) for initial segmentation of Total DNA copy numbers and
   further segmenting when parental copy numbers change.
 }

 \usage{psSegment(x, matching.reference=FALSE, alpha=0.01, alpha1=0.009, het.lower=0.1, het.upper=0.75, het.minimum.window=0.05, het.stepsize.window=0.01, het.stepsize.within=0.01, trim.mean=0.1, min.hetero=5, alpha.homozygous=0.05, alpha.equality=0.05, maf=0.075, min.homo.region=100, modeOffset=0.025, smooth.segmentation=TRUE, smooth.region=2, outlier.SD.scale=4, smooth.SD.scale=2, trim=0.025, impute.LOH, zero.homo=FALSE, verbose=1, ...)}

 \arguments{
   \item{x}{an object of class \code{\link{psCNA}}().}
   \item{matching.reference}{whether matching reference sample is
     available for every test sample.}
   \item{alpha}{overall significance level for the test to accept
     change-points.}
   \item{alpha1}{significance level for test on total copy number.  Must
     be less than alpha.}
   \item{het.lower}{lowest possible threshold on minor copy number to
     distinguish homozygotes from heterozygotes.}
   \item{het.upper}{highest possible threshold on minor copy number to
     distinguish homozygotes from heterozygotes.}
   \item{het.minimum.window}{smallest window to find threshold on
     genotypes.}
   \item{het.stepsize.window}{how much to increment the overall window
     size to find threshold on genotypes.}
   \item{het.stepsize.within}{how much to increment the current window
     size to find threshold on genotypes.}
   \item{trim.mean}{how much to trim when estimating parent-specific
     means.}
   \item{min.hetero}{minimum number of heterozygotes to estimate
     parent-specific means.}
   \item{alpha.homozygous}{significance level to test whether region is
     LOH.}
   \item{alpha.equality}{significance level to test whether
     parent-specific means are equal}
   \item{maf}{the assumed minor allele frequency used when estimating
     whether region is LOH; decreasing it will lead to more LOH regions.}
   \item{min.homo.region}{the smallest number of SNPs where no attempt
     made to combine regions when testing heterozygosity.}
   \item{modeOffset}{a tuning parameter for estimating whether the parental
     copy numbers are unequal in the second round of segmentation;
     increasing it will lead to fewer regions estimated to be unequal.}
   \item{smooth.segmentation}{whether to smooth before segmenting on
     total copy number.}
   \item{smooth.region}{number of points to consider on the left and the
     right of a point to detect it as an outlier.  Note this is a
     parameter in smooth.CNA.}
   \item{outlier.SD.scale}{the number of SDs away from the nearest point
     in the smoothing region to call a point an outlier.  Note this is a
     parameter in smooth.CNA.}
   \item{smooth.SD.scale}{the number of SDs from the median in the
     smoothing region where a smoothed point is positioned.}
   \item{trim}{proportion of data to be trimmed for variance calculation
     for smoothing outliers and undoing splits based on SD.}
   \item{impute.LOH}{whether to combine small regions with testing for
     LOH.}
   \item{zero.homo}{whether to make the minor copy number 0 for
     homozygotes. Default is \code{\link[base:logical]{FALSE}}.}
   \item{verbose}{verbosity of CBS}
   \item{...}{other arguments to be passed on to \code{\link[DNAcopy]{segment}}}
 }

 \value{
   A \code{\link[base]{data.frame}} with 13 columns that give the following characteristics
   of each segment: "ID", "chrom", "loc.start", "loc.end", "num.mark",
   "num.hetero", "mean.diff.mbaf", "mindiff.mean", "maxdiff.mean",
   "min.mean", "max.mean", "total.mean"
 }

 \details{
   This function implements the parent specific segmentation. It segments
   the total copy number.  Heterozygotes are found.  Then within each
   segment the allele specific copy number data for the heterozygotes are
   used to determine further segmentation is needed to allow for matched
   parent specific changes occur i.e. a gain in a parental chromosome is
   offset by a matching loss in the other keeping the total copy number
   constant.  After the segmentation is complete segments are called for
   LOH, equal parental copies or unequal copies.
 }

 \examples{
if (require("R.utils")) {

pathname <- system.file("data-ex/pcCBS,exData,chr01.Rbin", package="psCBS")
data <- R.utils::loadObject(pathname)

# Transform data
data$CTB <- data$betaT * data$CT
data$CTA <- data$CT - data$CTB
data$CNB <- data$betaN * data$CN
data$CNA <- data$CN - data$CNB

# Subset to speed up example
data0 <- data
keep <- seq(from=1, to=nrow(data), by=5)
data <- data[keep,,drop=FALSE]
str(data)
R.oo::attachLocally(data)

# Setup psCNA object
pd <- psCNA(genomdat.a=CTA, genomdat.b=CTB,
            chrom=chromosome, maploc=position,
            normaldat.a=CNA, normaldat.b=CNB)
print(pd)
## Number of Test Samples 1 Number of Reference Samples 1
## Number of Probes  14670

fit <- psSegment(pd, matching.reference=TRUE, impute.LOH=TRUE, verbose=1)
## [1] "Heterozygous cutoff is 0.702"
## [1] "Heterozygous cutoff is 0.697"
## Analyzing: Sample.1
## Analyzing: Sample.1
## Analyzing: Sample.1
## Analyzing: Sample.1
## Analyzing: Sample.1

str(fit)
## 'data.frame':   5 obs. of  13 variables:
##  $ ID            : Factor w/ 1 level "Sample.1": 1 1 1 1 1
##  $ chrom         : Factor w/ 1 level "1": 1 1 1 1 1
##  $ loc.start     : Factor w/ 5 levels "   554484"," 90446342",..: 1 2 3 4 5
##  $ loc.end       : Factor w/ 5 levels " 90437574"," 90564186",..: 1 2 3 4 5
##  $ num.mark      : Factor w/ 5 levels "1942","2675",..: 4 5 1 2 3
##  $ num.hetero    : Factor w/ 5 levels "0","1281","156",..: 3 1 4 5 2
##  $ mean.diff.mbaf: Factor w/ 4 levels "-0.00580908781811386",..: 3 NA 2 1 4
##  $ mindiff.mean  : Factor w/ 4 levels "0.532319235794478",..: 1 NA 2 3 4
##  $ maxdiff.mean  : Factor w/ 4 levels "0.837964234142425",..: 1 NA 2 3 4
##  $ confidence95  : Factor w/ 4 levels "0.00584787059363753",..: 4 NA 3 1 2
##  $ min.mean      : Factor w/ 3 levels "0","0.996010651148487",..: 1 1 1 3 2
##  $ max.mean      : Factor w/ 5 levels "0.746089002916825",..: 3 1 4 2 5
##  $ total.mean    : Factor w/ 5 levels "0.746089002916825",..: 2 1 3 4 5


} # if (require("R.utils"))
}


\keyword{nonparametric}
