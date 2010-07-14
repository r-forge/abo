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
