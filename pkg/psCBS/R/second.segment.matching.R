#diff.mbaf.vector.notmissing - absolute differences of A and B alleles alpha.level -
#alpha-level of segmentation window.size - size of windows for
#comparing left to right distributions quantile.values - quantiles of
#distribution being compared frequency.sample.markers - how often to
#compute p-values before smoothing f - amount of smoothing in lowess
#input.matrix - the results of the first round of segmentation
#map.data - map data
second.segment.matching <- function(matching,diff.mbaf.vector.notmissing,normal.data,alpha.level,input.matrix,map.data)
  {
    output.matrix <- as.matrix(input.matrix)
    colnames(output.matrix) <- colnames(input.matrix)
    n <- nrow(input.matrix)
    new.ends <- cumsum(input.matrix[,5])
    new.starts <- c(1,new.ends[-n]+1)
    counter <- 1
    for(i in 1:n)
      {
        new.row <- input.matrix[i,]
#        print(new.row)
        new.normal.data <- normal.data[new.starts[i]:new.ends[i]]
#Changed 02/09/10, no reason different between matched and unmatched
        which.het.data <- which(new.normal.data)
#        if(matching) which.het.data <- which(new.normal.data)        
#        else which.het.data <- 1:length(new.normal.data)
#Dropped criterion of not second round of segmentation if first round results in LOH
        if(sum(new.normal.data)>=10)
#        if(sum(new.normal.data)>=10 & !first.LOH[i])          
#        if(new.row[5]>=minimum.markers.second.round)          
          {
            new.diff.mbaf.vector.notmissing <- (diff.mbaf.vector.notmissing[new.starts[i]:new.ends[i]])[new.normal.data]
            second.output <- second.round.cbs.matching(new.diff.mbaf.vector.notmissing,alpha.level)
#            print(second.output)
#                        second.output <- second.round.cbs(new.diff.mbaf.vector.notmissing,new.normal.data,alpha.level)
#            second.output <- second.round(new.abs.diff.data,alpha.level,window.size,quantile.values,frequency.sample.markers)            
#            second.output <- second.round(new.abs.diff.data,alpha.level,window.size,quantile.values,frequency.sample.markers)            
            if(length(second.output)>0)
              {
                new.map.data <- map.data[new.starts[i]:new.ends[i]]
                new.output <- matrix(rep(as.character(new.row),length(second.output)+1),ncol=6,byrow=TRUE)
#                print(new.output)
                ends.output <- c(which.het.data[second.output],length(new.map.data))
#                print(paste("ends.output=",ends.output))
                starts.output <- c(1,ends.output[-length(ends.output)]+1)
#                print(paste("starts.output=",starts.output))
                new.output[,3] <- new.map.data[starts.output]
                new.output[,4] <- new.map.data[ends.output]
                new.output[,5] <- ends.output-starts.output+1
#                print(new.output)
                if(i>1 & i<n)
                  {
                    output.matrix <- rbind(output.matrix[1:(counter-1),],new.output,output.matrix[(counter+1):nrow(output.matrix),])                    
                  }
                else if(i==1)
                  {
                    if(n>1)
                      {
#                        print(dim(new.output))
#                        print(dim(output.matrix))
#                        print("1")
                        output.matrix <- rbind(new.output,output.matrix[2:n,])
#                        print("2")
                      }
                    else
                      {
#                        print("3")
                        output.matrix <- new.output
#                        print("4")
                      }
                  }
                else if (i==n)
                  {
                    output.matrix <- rbind(output.matrix[1:(counter-1),],new.output)
                  }
                counter <- counter+length(second.output)
              }
          }
        counter <- counter+1
      }
#    print("end of call")
    output.matrix
  }

