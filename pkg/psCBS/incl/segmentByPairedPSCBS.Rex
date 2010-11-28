library("R.utils");
library("aroma.light");


pathname <- system.file("data-ex/pcCBS,exData,chr01.Rbin", package="psCBS");
data <- R.utils::loadObject(pathname);

# Subset to speed up example
#data0 <- data;
#keep <- seq(from=1, to=nrow(data), by=5);
#data <- data[keep,,drop=FALSE];
# Order
o <- order(data$position);
data <- data[o,];
str(data);
R.oo::attachLocally(data);
x <- position;

# Drop single-locus outliers
CTs <- dropSegmentationOutliers(CT, chromosome=1, x=x, verbose=-10);

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CTs, betaT=betaT, betaN=betaN,
                            chromosome=1, x=x, verbose=-10);
print(fit);

# Plot results
plotTracks(fit);

fit2 <- bootstrapTCNandDHByRegion(fit, verbose=-10);
print(fit2);
plotTracks(fit2);

fit3 <- callABandHighAI(fit2, verbose=-10);
print(fit3);
plotTracks(fit3);


## # Call allelic balance and extreme allelic imbalance ("LOH")
## cfit <- callABandHighAI(fit, tauAB=0.10, alphaAB=0.05,
##                   tauHighAI=0.60, alphaHighAI=0.05, verbose=-10);
## 
## # Display segmentation and calls
## print(fit);
