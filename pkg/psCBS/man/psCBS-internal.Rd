\name{psCBS-internal}
\alias{callsegments.matching}
\alias{findheterozygous}
\alias{second.round.cbs.matching}
\alias{second.segment.matching}
\title{Internal psCBS functions}
\description{
  Internal functions of package psCBS.
}
\usage{
callsegments.matching(matching, segment.starts, segment.ends, a.vector,
    b.vector, diff.mbaf.vector, het.vector, normal.het.vector,
    chrom.vector.starts, total.vector.starts, min.hetero, trim.mean, 
    alpha.equality, het.bias, maf, alpha.homozygous, min.homo.region,
    impute.LOH)
findheterozygous(a.vector, b.vector, lower.bound=0, upper.bound=0.5,
    minimum.window=0.05, stepsize.window=0.01, stepsize.within=0.01)
second.round.cbs.matching(diff.data,alpha)
second.segment.matching(matching, diff.mbaf.vector.notmissing,
    normal.data, alpha.level, input.matrix,map.data)
}
\details{These are not to be called directly by the user}
\keyword{internal}
