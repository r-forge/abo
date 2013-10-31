\name{babel}
\alias{babel}
\title{Ribosome Profiling Analysis Program}
\description{
From paired mRNA and ribosome protected fragment count data run Babel analysis to detect changes in translation.
}
\usage{
babel(rna,rp,group,nreps,method.adjust="BH",...)
}
\arguments{
  \item{rna}{a matrix or data frame of mRNA count data whose rows
    correspond to genes and whose columns correspond to samples.}
  \item{rp}{a matrix or data frame of ribosome protected fragment count
    data whose rows correspond to genes and whose columns correspond to
    samples.  The row names and column names (if specified) must be the same
    as for rna.}
  \item{group}{class labels corresponding to the samples in rna and rp.
    There must be exactly two unique values.}
  \item{nreps}{number of permutations for within comparison.  Must be
    divisible by 10,000.}
  \item{method.adjust}{method of adjustment for multiple comparisons used
    by p.adjust.}
  \item{...}{Additional arguments.}
}

\value{
  A list with the following named components:
  \item{within}{A list with one component per sample that is named after the sample 
    that contains a data frame consisting of a gene id, direction (1 for
    translation greater than expected, -1 for less than expected), a
    one-sided p-value for the same, a two-sided p-value, and the 
    corresponding FDR.}
  \item{combined}{A list with two components (one per group label) that is named 
  after the group label that contains a data frame consisting of a
  gene id, direction ( 1 for translation greater than expected, -1 for  
  less than expected), a (two-sided) p-value, and the 
  corresponding FDR.}
\item{between}{A data frame consisting of a gene id, p-value, corresponding FDR,
  and direction (1 for translation higher in first group label, -1 for 
  lower in first group label).)  If there is at least two samples per group
  there are three addition variables: rna log fold-change, rna FDR,
  and change type (translation only or both rna and translation).}
}

\details{
  No missing values are allowed.
  
  Making nreps larger gives more precise estimates.  At least 100,000
  is required.  A million would be better, and ten million would be ideal, but will
  execute very slowly in the alpha version of the serial implementation.}  

\examples{
data(ribo.prof)
#Get rna count data
test.rna <- ribo.prof$test.rna
#Read in rp count data, if rownames or colnames differ from rna count
# data, adjust them so that they are the same
test.rp <- ribo.prof$test.rp
#Assign group labels to samples.  All groups will be tested pairwise.
test.group=c("A","B","A","B")
#Set the seed so the result is reproducibile
set.seed(12345)
#Run babel with 100000 repetitions.
test.babel <- babel(test.rna,test.rp,group=test.group,nreps=100000)
#Extract within sample p-values (is translation expected given mRNA level per sample)
within.babel <- test.babel$within
#Extract within group combined p-values (is translation expected given
# mRNA level per group label)
Combined.babel <- test.babel$combined
#Extract between group p-values (is translation given mRNA level the same between group)
between.babel <- test.babel$between
}
