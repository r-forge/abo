callsegments.matching <- function(matching, segment.starts, segment.ends, a.vector, b.vector, diff.mbaf.vector, het.vector, normal.het.vector, foldedTumorBoostBAF, chrom.vector.starts, total.vector.starts, min.hetero, trim.mean, alpha.equality, modeOffset, maf, alpha.homozygous, min.homo.region, impute.LOH, nBootstrapSamples=1000) {
  n.segments <- length(segment.starts)
  n.homo <- rep(NA, times=n.segments)
  total.output <- rep(NA, times=n.segments)
  min.output <- rep(NA, times=n.segments)
  max.output <- rep(NA, times=n.segments)
  mindiff.output <- rep(NA, times=n.segments)
  maxdiff.output <- rep(NA, times=n.segments)
  n.hetero <- rep(NA, times=n.segments)
  mean.diff.mbaf <- rep(NA, times=n.segments)
  mean.hetero <- rep(NA, times=n.segments)
  var.hetero <- rep(NA, times=n.segments)
  meandiff.hetero <- rep(NA, times=n.segments)
  sumsquarediff.hetero <- rep(NA, times=n.segments)
  pvalues.equality <- rep(NA, times=n.segments)
  length.segments <- segment.ends-segment.starts+1
  confidence95NewModes <- rep(NA, times=n.segments)

  # New loop for total output
  mean.index <- rep(NA, times=n.segments)
  mean.index[1] <- 1
  if (n.segments > 1) {
    total.same <- rep(TRUE, times=n.segments-1)
    total.same[(total.vector.starts[2:n.segments] != total.vector.starts[1:(n.segments-1)]) | (chrom.vector.starts[2:n.segments] != chrom.vector.starts[1:(n.segments-1)])] <- FALSE

    for (kk in 2:n.segments) {
      if (total.same[kk-1]) {
        mean.index[kk] <- mean.index[kk-1]
      } else {
        mean.index[kk] <- mean.index[kk-1] + 1
      }
    }
  } # if (n.segments > 1)

  for (kk in seq(length=max(mean.index))) {
    new.i <- which(mean.index == kk)
    which.i <- segment.starts[min(new.i)]:segment.ends[max(new.i)]
    new.a <- a.vector[which.i]
    new.b <- b.vector[which.i]
    total.output[min(new.i):max(new.i)] <- mean(new.a+new.b, trim=trim.mean)
  }

  # Now check for LOH
  for (kk in seq(length=n.segments)) {
    n.homo[kk] <- length.segments[kk]-sum(het.vector[segment.starts[kk]:segment.ends[kk]])
    pvalue.LOH <- binom.test(n.homo[kk], length.segments[kk], p=(maf^2 + (1-maf)^2), alt="greater")$p.value
    if (pvalue.LOH < alpha.homozygous) {
      min.output[kk] <- 0
      max.output[kk] <- total.output[kk]
    }
  }

  n.hetero <- length.segments-n.homo

  # Think about adding code for small regions
  if (impute.LOH) {
    large.segments <- rep(TRUE, times=n.segments)
    large.segments[length.segments < min.homo.region] <- FALSE
    chrom.unique <- unique(chrom.vector.starts)
    expected.homo <- maf^2 + (1-maf)^2
    p.homo <- n.homo / length.segments
    convert.homo <- rep(FALSE, times=n.segments)

    for (chr in chrom.unique) {
      which.candidates <- which(chrom.vector.starts == chr & (is.na(min.output) | min.output != 0) & !large.segments)

      if (length(which.candidates) > 0) {
        merge.candidates <- NULL
        which.i <- which(chrom.vector.starts == chr)
        min.i <- min(which.i)
        max.i <- max(which.i)

        if (min.i < max.i) {
          merge.candidates <- NULL

          for (jj in which.candidates) {
            if (jj == min.i) {
              if (!is.na(min.output[jj+1]) && min.output[jj+1] == 0 && 
                   p.homo[jj] > expected.homo) {
                # If first region and next door is homo and more homo
                # than expected, convert...
                convert.homo[jj] <- TRUE
              } else if ((jj+1) <= n.segments && chrom.vector.starts[jj+1] == chr && !large.segments[jj+1]) {
                # ...else if next door is not a large segment merge
                # candidates
                merge.candidates <- jj
              }
            } else if (jj == max.i) {
              if (!is.na(min.output[jj-1]) && min.output[jj-1] == 0 && p.homo[jj] > expected.homo) {
                # If last region and next door is homo and more homo
                # than expected, convert...
                convert.homo[jj] <- TRUE
              } else if (!is.na(match(jj-1, merge.candidates))) {
                # ...merge everything together and test
                merge.candidates <- c(merge.candidates, jj)
                new.length <- sum(length.segments[merge.candidates])
                new.homo <- sum(n.homo[merge.candidates])
                new.pvalue <- binom.test(new.homo, new.length, p=(maf^2 + (1-maf)^2), alt="greater")$p.value
                if (new.pvalue < alpha.homozygous) {
                  convert.homo[merge.candidates] <- TRUE
                }
              }
            } else {
              # Both sides are homozygous convert                            
              if (!is.na(min.output[jj-1]) && min.output[jj-1] == 0 &&
                  !is.na(min.output[jj+1]) && min.output[jj+1] == 0 &&
                   p.homo[jj] > expected.homo) {
                convert.homo[jj] <- TRUE
              } else if (!is.na(match(jj-1, merge.candidates)) && !is.na(match(jj+1, which.candidates))) {
                # Merge if left and right are candidates
                merge.candidates <- c(merge.candidates, jj)
              } else if (!is.na(match(jj-1, merge.candidates)) && is.na(match(jj+1, which.candidates))) {
                # If candidates just on left, merge and start merging over
                merge.candidates <- c(merge.candidates, jj)
                new.length <- sum(length.segments[merge.candidates])
                new.homo <- sum(n.homo[merge.candidates])
                new.pvalue <- binom.test(new.homo, new.length, p=(maf^2 + (1-maf)^2), alt="greater")$p.value

                if (new.pvalue < alpha.homozygous) {
                  convert.homo[merge.candidates] <- TRUE
                }

                merge.candidates <- NULL
              }
            }
          } # for (jj in which.candidates)
        } # if (min.i < max.i)
      } # if (length(which.candidates) > 0)
    } # for (chr in chrom.unique)

    if (sum(convert.homo) > 0) {
      for (kk in seq(length=n.segments)) {
        if (convert.homo[kk]) {
          min.output[kk] <- 0
          max.output[kk] <- total.output[kk]
        }
      }
    }
  } # if (impute.LOH)

  # print(min.output)
  # print(max.output)

  # Now make estimates for non-LOH segments
  # Big problem if only one segment!!!    

  if (matching) {
    for (kk in seq(length=n.segments)) {
      which.i <- segment.starts[kk]:segment.ends[kk]
      new.a <- a.vector[which.i]
      new.b <- b.vector[which.i]
      newFoldedTumorBoostBAF <- foldedTumorBoostBAF[which.i]
      new.diff.mbaf <- diff.mbaf.vector[which.i] 
      # total.output[kk] <- mean(new.a+new.b, trim=trim.mean)
      new.normal.hetero <- normal.het.vector[which.i]
      new.tumor.hetero <- het.vector[which.i]
      # n.hetero[kk] <- sum(new.normal.hetero)
      sumNewNormalHetero <- sum(new.normal.hetero)
      # if (sumNewNormalHetero >= min.hetero)
      # if (n.hetero[kk] >= min.hetero)          
      if (sum(new.normal.hetero) >= min.hetero &&
          sum(new.tumor.hetero) >= min.hetero) {
        new.a.hetero <- new.a[new.tumor.hetero]
        new.b.hetero <- new.b[new.tumor.hetero]
        meandiff.hetero[kk] <- mean(abs(new.a.hetero-new.b.hetero), trim=trim.mean)
        newFoldedTumorBoostBAFHetero <- newFoldedTumorBoostBAF[new.normal.hetero]
        newModes <- rep(NA, times=nBootstrapSamples)
        for (jj in seq(length=nBootstrapSamples)) {
          newModes[jj] <- hrmode(sample(newFoldedTumorBoostBAFHetero, replace=TRUE))
        }

        confidence95NewModes[kk] <- quantile(newModes, 0.05)-modeOffset
        new.diff.mbaf.hetero <- new.diff.mbaf[new.tumor.hetero]
        mean.diff.mbaf[kk] <- mean(new.diff.mbaf.hetero)

        # Add small amount of noise in case mbaf are constant
        # pvalues.equality[kk] <- t.test(new.diff.mbaf.hetero-modeOffset+rnorm(length(new.diff.mbaf.hetero), 0, 0.000001), alt="greater")$p.value
        # pvalues.equality[kk] <- t.test(new.diff.mbaf.hetero-het.bias, alt="greater")$p.value            
        # print(paste("Segment", kk))
        # print(paste("n=", length(new.diff.mbaf.hetero)))
        # print(paste("mean(new.diff.mbaf.hetero)=", mean(new.diff.mbaf.hetero)))
        # print(paste("mean(new.diff.mbaf.hetero)-het.bias=", mean(new.diff.mbaf.hetero)-het.bias))
        # print(paste("pvalues.equality[i]=", pvalues.equality[kk], sep=""))

        maxdiff.output[kk] <- (total.output[kk]+meandiff.hetero[kk])/2
        mindiff.output[kk] <- (total.output[kk]-meandiff.hetero[kk])/2
        if (is.na(min.output[kk]) || min.output[kk] != 0) {
          # if (pvalues.equality[kk] > alpha.equality) 
          if (confidence95NewModes[kk] < 0) {
            min.output[kk] <- total.output[kk]/2
            max.output[kk] <- total.output[kk]/2
          } else {
            max.output[kk] <- maxdiff.output[kk]
            min.output[kk] <- mindiff.output[kk]
            # max.output[kk] <- (total.output[kk]+meandiff.hetero[kk])/2
            # min.output[kk] <- (total.output[kk]-meandiff.hetero[kk])/2
          }
        }
      }
    } # for (kk in ...)
  } else {  # if (matching) 
    for (kk in seq(length=n.segments)) {
      which.i <- segment.starts[kk]:segment.ends[kk]
      new.a <- a.vector[which.i]
      new.b <- b.vector[which.i]
      new.diff.mbaf <- diff.mbaf.vector[which.i] 
      new.tumor.hetero <- het.vector[which.i]
      sumNewTumorHetero <- sum(new.tumor.hetero)

      if (sumNewTumorHetero >= min.hetero) {
        new.a.hetero <- new.a[new.tumor.hetero]
        new.b.hetero <- new.b[new.tumor.hetero]
        meandiff.hetero[kk] <- mean(abs(new.a.hetero-new.b.hetero), trim=trim.mean)
        new.diff.mbaf.hetero <- new.diff.mbaf[new.tumor.hetero]
        mean.diff.mbaf[kk] <- mean(new.diff.mbaf.hetero)
        pvalues.equality[kk] <- t.test(new.diff.mbaf.hetero-modeOffset + rnorm(length(new.diff.mbaf.hetero), mean=0, sd=0.000001), alt="greater")$p.value
        maxdiff.output[kk] <- (total.output[kk]+meandiff.hetero[kk])/2
        mindiff.output[kk] <- (total.output[kk]-meandiff.hetero[kk])/2
        if (is.na(min.output[kk]) || min.output[kk] != 0) {
          if (pvalues.equality[kk] > alpha.equality) {
            min.output[kk] <- total.output[kk]/2
            max.output[kk] <- total.output[kk]/2
          } else {
            max.output[kk] <- maxdiff.output[kk]
            min.output[kk] <- mindiff.output[kk]
          }
        }
      }
    } # for (kk in ...)
  } # if (matching)
    
  list(total.output=total.output, min.output=min.output, max.output=max.output, n.hetero=n.hetero, mindiff.output=mindiff.output, maxdiff.output=maxdiff.output, mean.diff.mbaf=mean.diff.mbaf, confidence95NewModes=(confidence95NewModes+modeOffset))
} # callsegments.matching()


############################################################################
# HISTORY:
# 2010-07-08
# o ROBUSTNESS: Now regions are indexed by 'kk' and chromosomes by 'chr'.
# o ROBUSTNESS: Replaced all 1:n with seq(length=n) to deal with n == 0.
# o ROBUSTNESS: Now all interator variables i & j are written as ii & jj.
############################################################################
