# \arguments{
#  \item{diff.mbaf.vector.notmissing}{absolute differences of A and B alleles}
#  \item{alpha.level}{alpha-level of segmentation}
#  \item{window.size}{size of windows for#comparing left to right distributions}
#  \item{quantile.values}{quantiles of distribution being compared}
#  \item{frequency.sample.markers}{how often to compute p-values before smoothing}
#  \item{f}{amount of smoothing in lowess}
#  \item{input.matrix}{the results of the first round of segmentation}
#  \item{map.data}{map data}
# }
second.segment.matching <- function(matching, diff.mbaf.vector.notmissing, normal.data, alpha.level, input.matrix, map.data) {
  output.matrix <- as.matrix(input.matrix)
  colnames(output.matrix) <- colnames(input.matrix)
  nbrOfRegions <- nrow(input.matrix)
  new.ends <- cumsum(input.matrix[,"num.mark"])
  new.starts <- c(1, new.ends[-nbrOfRegions]+1)
  counter <- 1

  for (kk in seq(length=nbrOfRegions)) {
    new.row <- input.matrix[kk,]
    # print(new.row)

    new.normal.data <- normal.data[new.starts[kk]:new.ends[kk]]
    # Changed 02/09/10, no reason different between matched and unmatched
    which.het.data <- which(new.normal.data)
    # if (matching) which.het.data <- which(new.normal.data)        
    # else which.het.data <- 1:length(new.normal.data)
    # Dropped criterion of not second round of segmentation if 
    # first round results in LOH
    # if (sum(new.normal.data) >= 10 && !first.LOH[kk])          
    # if (new.row[5] >= minimum.markers.second.round)          
    if (sum(new.normal.data) >= 10) {
      new.diff.mbaf.vector.notmissing <- (diff.mbaf.vector.notmissing[new.starts[kk]:new.ends[kk]])[new.normal.data]
      second.output <- second.round.cbs.matching(new.diff.mbaf.vector.notmissing, alpha.level)
      # print(second.output)
      # second.output <- second.round.cbs(new.diff.mbaf.vector.notmissing, new.normal.data, alpha.level)
      # second.output <- second.round(new.abs.diff.data, alpha.level, window.size, quantile.values, frequency.sample.markers)            
      # second.output <- second.round(new.abs.diff.data, alpha.level, window.size, quantile.values, frequency.sample.markers)

      if (length(second.output) > 0) {
        new.map.data <- map.data[new.starts[kk]:new.ends[kk]]
        new.output <- matrix(rep(as.character(new.row), times=length(second.output)+1), ncol=6, byrow=TRUE)
        # print(new.output)
        ends.output <- c(which.het.data[second.output], length(new.map.data))
        # print(paste("ends.output=", ends.output))

        starts.output <- c(1, ends.output[-length(ends.output)]+1)
        # print(paste("starts.output=", starts.output))

        new.output[,3] <- new.map.data[starts.output]
        new.output[,4] <- new.map.data[ends.output]
        new.output[,5] <- ends.output-starts.output+1
        # print(new.output)

        if (kk > 1 && kk < nbrOfRegions) {
          output.matrix <- rbind(output.matrix[1:(counter-1),], new.output, output.matrix[(counter+1):nrow(output.matrix),])                    
        } else if (kk == 1) {
          if (nbrOfRegions > 1) {
            # print(dim(new.output))
            # print(dim(output.matrix))
            # print("1")
            output.matrix <- rbind(new.output, output.matrix[2:nbrOfRegions,])
            # print("2")
          } else {
            # print("3")
            output.matrix <- new.output
            # print("4")
          }
        } else if (kk == nbrOfRegions) {
          output.matrix <- rbind(output.matrix[1:(counter-1),], new.output)
        } # if (kk > 1 && kk < nbrOfRegions)

        counter <- counter+length(second.output)
      } # if (length(second.output) > 0)
    } # if (sum(new.normal.data) >= 10)

    counter <- counter+1
  } # for (kk in ...)

  # print("end of call")

  output.matrix
} # second.segment.matching()


############################################################################
# HISTORY:
# 2010-07-08
# o ROBUSTNESS: Segmentation results are now subsetted by names not indices.
# o ROBUSTNESS: Now regions are index by 'kk'.
# o ROBUSTNESS: Replaced all 1:n with seq(length=n) to deal with n == 0.
# o ROBUSTNESS: Now all interator variables i & j are written as ii & jj.
############################################################################
