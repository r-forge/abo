%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  callNaiveHeterozygotes.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{callNaiveHeterozygotes}
\alias{callNaiveHeterozygotes.default}
\alias{callNaiveHeterozygotes}

 \title{Calls heterozygous SNPs from allele-specific signals}

 \description{
  Calls heterozygous SNPs from allele-specific signals (yA,yB) using the method described in [1].
 }

 \usage{\method{callNaiveHeterozygotes}{default}(yA, yB, range=c(0, 0.5), stepsize.within=0.01, minimum.window=0.05, stepsize.window=0.01, ..., verbose=FALSE)}

 \arguments{
   \item{yA,yB}{\code{\link[base]{numeric}} \code{\link[base]{vector}}s of length J containing
     allele A and allele B signals.}
   \item{range}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length 2 specifying the range
     of possible values of the estimated cutoff parameter L (see below).}
   \item{stepsize.within}{A \code{\link[base]{numeric}} specifying the "precision" of L.}
   \item{minimum.window}{A \code{\link[base]{numeric}} specifying the minimum window size
     used for the objective function.}
   \item{stepsize.window}{A \code{\link[base]{numeric}} specifying the size of each window
     increment.}
   \item{...}{Not used.}
   \item{verbose}{A \code{\link[base]{logical}} specifying if verbose output should be
     given or not.}
 }

 \value{
   Returns a \code{\link[base]{logical}} \code{\link[base]{vector}} of J calls.
   Attribute \code{fit} contains the estimated parameters.
 }

 \details{
   This function calls SNPs that are heterozygous, here heterzygous
   mean non-homozygous.  A homozygous SNP is a SNP where one of the
   two alleles has zero signal/copy.  Note that with this definition,
   a heterozygous SNP can be any SNP where both alleles have non-zero
   signal, regardless if they have the same signal or not.
   For diploid SNPs, this means that AA and BB SNPs are called
   non-heterozygotes (homozygotes) and AB SNPs are called heterozygotes.
   For non-diploid SNPs, say, AAA, AAB, ABB, and BBB, the 2nd and 3rd
   are called heterzygotes, the other two not.
   This function was designed to be able to call SNPs that are
   heterozygotes also in tumors.
   Note that the function does \emph{not} call genotypes.

   This function finds the optimal window to separate homozygous from
   heterozygous based on the minimum of the A and B alleles.
   It starts with window [range[1], range[1]+minimum.window].
   It shifts window by stepsize.within.  The process is repeated with bigger
   and bigger windows where the size is incremented by stepsize.window.
 }

 \section{Missing values}{
   TO DO.
   Ideally we should return \code{\link[base]{NA}} for SNPs that contain a missing value.
 }

 \examples{
library("R.utils");
library("aroma.light");

pathname <- system.file("data-ex/pcCBS,exData,chr01.Rbin", package="psCBS");
data <- R.utils::loadObject(pathname);

# Subset to speed up example
data0 <- data;
#keep <- seq(from=1, to=nrow(data), by=5);
#data <- data[keep,,drop=FALSE];
# Order
o <- order(data$position);
data <- data[o,];
str(data);
R.oo::attachLocally(data);
x <- position;


# Plotting data
xMb <- x/1e6;
xlab <- "Position (Mb)";
pch <- ".";

yB <- betaN * CN;
yA <- CN - yB;

layout(matrix(1:4, ncol=1));
par(mar=c(2,4,1,1)+0.1);

# Default parameters for psSement() that are not the same as
# the defaults of findheterozygous() and callNaiveHeterozygotes().
het.lower=0.1;
het.upper=0.75;

fit <- findheterozygous(yA, yB, lower.bound=het.lower, upper.bound=het.upper);
isHet1 <- isHet <- fit$is.heterozygous;
print(table(isHet));
col <- c("#aaaaaa", "#000000")[isHet+1];
plot(xMb,betaN, pch=pch, col=col, ylim=c(0,1), xlab=xlab);

range <- c(het.lower, het.upper);
isHet2 <- isHet <- callNaiveHeterozygotes(yA, yB, range=range);
print(table(isHet));
col <- c("#aaaaaa", "#000000")[isHet+1];
plot(xMb,betaN, pch=pch, col=col, ylim=c(0,1), xlab=xlab);

# Sanity check
stopifnot(all(isHet2 == isHet1));

muN <- callNaiveGenotypes(betaN);
isHet3 <- isHet <- (muN == 1/2);
print(table(isHet));
col <- c("#aaaaaa", "#000000")[isHet+1];
plot(xMb,betaN, pch=pch, col=col, ylim=c(0,1), xlab=xlab);

isHet4 <- isHet <- callNaiveHeterozygotes(1-betaN, betaN, range=c(0,1/2));
print(table(isHet));
col <- c("#aaaaaa", "#000000")[isHet+1];
plot(xMb,betaN, pch=pch, col=col, ylim=c(0,1), xlab=xlab);
}

 \references{
   [1] AB Olshen, RA Olshen, H Bengtsson, P Neuvial, P Spellman,
       and VE. Seshan, \emph{Parent-specific copy number in paired
       tumor-normal studies using circular binary segmentation},
       (submitted), 2010.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{internal}
