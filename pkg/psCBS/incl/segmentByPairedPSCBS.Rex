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
fit <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN, x=x, verbose=-10);

# Call allelic balance and "LOH" (= extreme allelic imbalance)
cfit <- callABandHighAI(fit, verbose=-10);

# Display segmentation and calls
print(fit);

# Plot results
plot(fit, xScale=1e-6);
