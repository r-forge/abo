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
keep <- seq(from=1, to=nrow(data), by=10)
data <- data[keep,,drop=FALSE]
str(data)
R.oo::attachLocally(data)

# Setup psCNA object
pd <- psCNA(genomdat.a=CTA, genomdat.b=CTB, 
            chrom=chromosome, maploc=position, 
            normaldat.a=CNA, normaldat.b=CNB)
print(pd)
## Number of Test Samples 1 Number of Reference Samples 1
## Number of Probes  7335

fit <- psSegment(pd, matching.reference=TRUE, impute.LOH=TRUE, verbose=1)
## [1] "Heterozygous cutoff is 0.695"
## [1] "Heterozygous cutoff is 0.686"
## Analyzing: Sample.1
## Analyzing: Sample.1
## Analyzing: Sample.1
## Analyzing: Sample.1

# AD HOC: Turn fit into a data frame
fit <- as.data.frame(fit, stringsAsFactors=FALSE);
for (cc in 3:5) { fit[,cc] <- as.integer(fit[,cc]); }
for (cc in 6:ncol(fit)) { fit[,cc] <- as.double(fit[,cc]); }

str(fit)
## 'data.frame':   3 obs. of  13 variables:
##  $ ID            : chr  "Sample.1" "Sample.1" "Sample.1"
##  $ chrom         : chr  "1" "1" "1"
##  $ loc.start     : int  554484 143726347 185505853
##  $ loc.end       : int  143656390 185438403 246679946
##  $ num.mark      : int  3766 1365 2204
##  $ num.hetero    : num  113 360 667
##  $ mean.diff.mbaf: num  -0.05794 -0.00943 0.04095
##  $ mindiff.mean  : num  0.521 0.898 0.992
##  $ maxdiff.mean  : num  0.85 1.15 1.64
##  $ confidence95  : num  0.471 0.036 0.191
##  $ min.mean      : num  0 0.898 0.992
##  $ max.mean      : num  1.37 1.15 1.64
##  $ total.mean    : num  1.37 2.05 2.63


} # if (require("R.utils"))
