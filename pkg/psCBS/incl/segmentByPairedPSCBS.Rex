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

# PSCBS segmentation
fit <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN, x=x, verbose=-10);
segs <- fit$output;
segs <- segs[,-(1:2)];
print(segs);

# Plotting data
xMb <- x/1e6;
xlab <- "Position (Mb)";
pch <- ".";

muN <- callNaiveGenotypes(betaN);
isHet <- (muN == 1/2);
betaTN <- normalizeTumorBoost(betaT, betaN, muN=muN);
rhoTN <- 2*abs(betaTN-1/2);
cps <- sort(c(segs$loc.start, segs$loc.end));

layout(matrix(1:5, ncol=1));
par(mar=c(2,4,1,1)+0.1);
plot(xMb,CT, pch=pch, col="#000000", ylim=c(0,5), xlab=xlab);
abline(v=cps/1e6, lwd=2);
drawLevels(fit, field="tcn.mean", xScale=1e-6);
col <- c("#aaaaaa", "#000000", "#aaaaaa")[2*muN+1];
plot(xMb,betaN, pch=pch, col=col, ylim=c(0,1), xlab=xlab);
abline(v=cps/1e6, lwd=2);
plot(xMb,betaT, pch=pch, col=col, ylim=c(0,1), xlab=xlab);
abline(v=cps/1e6, lwd=2);
plot(xMb,betaTN, pch=pch, col=col, ylim=c(0,1), xlab=xlab);
abline(v=cps/1e6, lwd=2);
rho <- rhoTN[isHet];
plot(xMb[isHet],rho, pch=pch, col="#000000", ylim=c(0,1), xlab=xlab);
abline(v=cps/1e6, lwd=2);
drawLevels(fit, field="dh.mean", xScale=1e-6);
