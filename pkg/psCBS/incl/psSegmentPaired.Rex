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

fit <- psSegmentPaired(pd, impute.LOH=TRUE, verbose=1)
## [1] "Heterozygous cutoff is 0.695"
## [1] "Heterozygous cutoff is 0.686"
## Analyzing: Sample.1
## Analyzing: Sample.1
## Analyzing: Sample.1
## Analyzing: Sample.1

str(fit)


} # if (require("R.utils"))
