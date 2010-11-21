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

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN,
                                 chromosome=1, x=x, verbose=-10);

# Call allelic balance and extreme allelic imbalance ("LOH")
cfit <- callABandHighAI(fit, tauAB=0.10, alphaAB=0.05,
                  tauHighAI=0.60, alphaHighAI=0.05, verbose=-10);

# Display segmentation and calls
print(fit);

# Plot results
plotTracks(fit, xScale=1e-6);

x11();
fit2 <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN,
          chromosome=1, x=x, joinSegments=TRUE, verbose=-10);
plotTracks(fit2, xScale=1e-6);
