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
