callsegments.matching <- function(matching,segment.starts,segment.ends,a.vector,b.vector,diff.mbaf.vector,het.vector,normal.het.vector,foldedTumorBoostBAF,chrom.vector.starts,total.vector.starts,min.hetero,trim.mean,alpha.equality,modeOffset,maf,alpha.homozygous,min.homo.region,impute.LOH,nBootstrapSamples=1000)  
  {
    n.segments <- length(segment.starts)
    n.homo <- rep(NA,n.segments)
    total.output <- rep(NA,n.segments)
    min.output <- rep(NA,n.segments)
    max.output <- rep(NA,n.segments)
    mindiff.output <- rep(NA,n.segments)
    maxdiff.output <- rep(NA,n.segments)
    n.hetero <- rep(NA,n.segments)
    mean.diff.mbaf <- rep(NA,n.segments)
    mean.hetero <- rep(NA,n.segments)
    var.hetero <- rep(NA,n.segments)
    meandiff.hetero <- rep(NA,n.segments)
    sumsquarediff.hetero <- rep(NA,n.segments)
    pvalues.equality <- rep(NA,n.segments)
    length.segments <- segment.ends-segment.starts+1
    confidence95NewModes <- rep(NA,n.segments)
                                        #New loop for total output
    mean.index <- rep(NA,n.segments)
    mean.index[1] <- 1
    if(n.segments>1)
      {
        total.same <- rep(TRUE,n.segments-1)
        total.same[(total.vector.starts[2:n.segments]!=total.vector.starts[1:(n.segments-1)])|(chrom.vector.starts[2:n.segments]!=chrom.vector.starts[1:(n.segments-1)])] <- FALSE
        for(i in 2:n.segments)
          {
            if(total.same[i-1]) mean.index[i] <- mean.index[i-1]
            else mean.index[i] <- mean.index[i-1]+1
          }
      }
    for(i in 1:max(mean.index))
      {
        new.i <- which(mean.index==i)
        which.i <- segment.starts[min(new.i)]:segment.ends[max(new.i)]
        new.a <- a.vector[which.i]
        new.b <- b.vector[which.i]
        total.output[min(new.i):max(new.i)] <- mean(new.a+new.b,trim=trim.mean)
      }
#Now check for LOH

    for(i in 1:n.segments)
      {
        n.homo[i] <- length.segments[i]-sum(het.vector[segment.starts[i]:segment.ends[i]])
        pvalue.LOH <- binom.test(n.homo[i],length.segments[i],p=(maf^2 + (1-maf)^2),alt="greater")$p.value
        if(pvalue.LOH<alpha.homozygous)
          {
            min.output[i] <- 0
            max.output[i] <- total.output[i]
          }
      }
    n.hetero <- length.segments-n.homo
#Think about adding code for small regions
    if(impute.LOH)
      {
        large.segments <- rep(TRUE,n.segments)
        large.segments[length.segments<min.homo.region] <- FALSE
        chrom.unique <- unique(chrom.vector.starts)
        expected.homo <- maf^2 + (1-maf)^2
        p.homo <- n.homo/length.segments
        convert.homo <- rep(FALSE,n.segments)
        for(i in chrom.unique)
          {
            which.candidates <- which(chrom.vector.starts==i & (is.na(min.output)|min.output!=0) & !large.segments)
            if(length(which.candidates)>0)
              {
                merge.candidates <- NULL
                which.i <- which(chrom.vector.starts==i)
                min.i <- min(which.i)
                max.i <- max(which.i)
                if(min.i<max.i)
                  {
                    merge.candidates <- NULL
                    for(j in which.candidates)
                      {

                        if(j==min.i)
                          {
#If first region and next door is homo and more homo than expected, convert
                            if(!is.na(min.output[j+1]) & min.output[j+1]==0 & p.homo[j]>expected.homo)
                              {
                                convert.homo[j] <- TRUE
                              }
#Else if next door is not a large segment merge candidates                            
                            else if((j+1)<=n.segments & chrom.vector.starts[j+1]==i & !large.segments[j+1])
                              {
                                merge.candidates <- j
                              }
                          }
                        else if(j==max.i)
                          {
#If last region and next door is homo and more homo than expected, convert                            
                            if(!is.na(min.output[j-1]) & min.output[j-1]==0 & p.homo[j]>expected.homo)
                              {
                                convert.homo[j] <- TRUE
                              }
#Merge everything together and test                            
                            else if(!is.na(match(j-1,merge.candidates)))
                              {
                                merge.candidates <- c(merge.candidates,j)
                                new.length <- sum(length.segments[merge.candidates])
                                new.homo <- sum(n.homo[merge.candidates])
                                new.pvalue <- binom.test(new.homo,new.length,p=(maf^2 + (1-maf)^2),alt="greater")$p.value
                                if(new.pvalue<alpha.homozygous)
                                  convert.homo[merge.candidates] <- TRUE
                              }
                          }
                        else
                          {
#Both sides are homozygous convert                            
                            if(!is.na(min.output[j-1]) & min.output[j-1]==0 &
                               !is.na(min.output[j+1]) & min.output[j+1]==0 & p.homo[j]>expected.homo)
                              {
                                convert.homo[j] <- TRUE
                              }
#Merge if left and right are candidates                            
                            else if (!is.na(match(j-1,merge.candidates)) & !is.na(match(j+1,which.candidates)))
                              {
                              
                                merge.candidates <- c(merge.candidates,j)
                              }
#If candidates just on left, merge and start merging over                            
                            else if (!is.na(match(j-1,merge.candidates)) & is.na(match(j+1,which.candidates)))
                              {
                                merge.candidates <- c(merge.candidates,j)
                                new.length <- sum(length.segments[merge.candidates])
                                new.homo <- sum(n.homo[merge.candidates])
                                new.pvalue <- binom.test(new.homo,new.length,p=(maf^2 + (1-maf)^2),alt="greater")$p.value
                                if(new.pvalue<alpha.homozygous)
                                  convert.homo[merge.candidates] <- TRUE
                                merge.candidates <- NULL
                              }
                          }
                      }
                  }
              }
          }
        if(sum(convert.homo)>0)
          {
            for(i in 1:n.segments)
              {
                if(convert.homo[i])
                  {
                    min.output[i] <- 0
                    max.output[i] <- total.output[i]
                  }
              }
          }
      }
#    print(min.output)
#    print(max.output)
#Now make estimates for non-LOH segments
#Big problem if only one segment!!!    
#    for(i in which(is.na(min.output)))
    if(matching)
      {
        for(i in 1:n.segments)      
          {
            which.i <- segment.starts[i]:segment.ends[i]
            new.a <- a.vector[which.i]
            new.b <- b.vector[which.i]
            newFoldedTumorBoostBAF <- foldedTumorBoostBAF[which.i]
            new.diff.mbaf <- diff.mbaf.vector[which.i] 
                                        #        total.output[i] <- mean(new.a+new.b,trim=trim.mean)
            new.normal.hetero <- normal.het.vector[which.i]
            new.tumor.hetero <- het.vector[which.i]
                                        #        n.hetero[i] <- sum(new.normal.hetero)
            sumNewNormalHetero <- sum(new.normal.hetero)
#            if(sumNewNormalHetero>=min.hetero)
            if(sum(new.normal.hetero)>=min.hetero & sum(new.tumor.hetero)>=min.hetero)              
                                        #        if(n.hetero[i]>=min.hetero)          
              {
                new.a.hetero <- new.a[new.tumor.hetero]
                new.b.hetero <- new.b[new.tumor.hetero]
                meandiff.hetero[i] <- mean(abs(new.a.hetero-new.b.hetero),trim=trim.mean)
                newFoldedTumorBoostBAFHetero <- newFoldedTumorBoostBAF[new.normal.hetero]
                newModes <- rep(NA,nBootstrapSamples)
                for(j in 1:nBootstrapSamples)
                  {
                    newModes[j] <- hrmode(sample(newFoldedTumorBoostBAFHetero,replace=TRUE))
                  }

                confidence95NewModes[i] <- quantile(newModes,0.05)-modeOffset
                new.diff.mbaf.hetero <- new.diff.mbaf[new.tumor.hetero]
                mean.diff.mbaf[i] <- mean(new.diff.mbaf.hetero)

                                        # Add small amount of noise in case mbaf are constant
                
#                pvalues.equality[i] <- t.test(new.diff.mbaf.hetero-modeOffset+rnorm(length(new.diff.mbaf.hetero),0,0.000001),alt="greater")$p.value
                                        #            pvalues.equality[i] <- t.test(new.diff.mbaf.hetero-het.bias,alt="greater")$p.value            
                                        #            print(paste("Segment",i))
                                        #            print(paste("n=",length(new.diff.mbaf.hetero)))
                                        #            print(paste("mean(new.diff.mbaf.hetero)=",mean(new.diff.mbaf.hetero)))
                                        #            print(paste("mean(new.diff.mbaf.hetero)-het.bias=",mean(new.diff.mbaf.hetero)-het.bias))
                                        #            print(paste("pvalues.equality[i]=",pvalues.equality[i],sep=""))
                maxdiff.output[i] <- (total.output[i]+meandiff.hetero[i])/2
                mindiff.output[i] <- (total.output[i]-meandiff.hetero[i])/2
                if(is.na(min.output[i])|min.output[i]!=0)
                  {
                    if(confidence95NewModes[i]<0)
#                    if(pvalues.equality[i]>alpha.equality)                      
                      {
                        min.output[i] <- total.output[i]/2
                        max.output[i] <- total.output[i]/2
                      }
                    else
                      {
                        max.output[i] <- maxdiff.output[i]
                        min.output[i] <- mindiff.output[i]
                                        #                max.output[i] <- (total.output[i]+meandiff.hetero[i])/2
                                        #                min.output[i] <- (total.output[i]-meandiff.hetero[i])/2
                      }
                  }
              }
          }#End of n.segments loop
      }
    else
      {
        for(i in 1:n.segments)      
          {
            which.i <- segment.starts[i]:segment.ends[i]
            new.a <- a.vector[which.i]
            new.b <- b.vector[which.i]
            new.diff.mbaf <- diff.mbaf.vector[which.i] 
            new.tumor.hetero <- het.vector[which.i]
            sumNewTumorHetero <- sum(new.tumor.hetero)
            if(sumNewTumorHetero>=min.hetero)
              {
                new.a.hetero <- new.a[new.tumor.hetero]
                new.b.hetero <- new.b[new.tumor.hetero]
                meandiff.hetero[i] <- mean(abs(new.a.hetero-new.b.hetero),trim=trim.mean)
                new.diff.mbaf.hetero <- new.diff.mbaf[new.tumor.hetero]
                mean.diff.mbaf[i] <- mean(new.diff.mbaf.hetero)
                pvalues.equality[i] <- t.test(new.diff.mbaf.hetero-modeOffset+rnorm(length(new.diff.mbaf.hetero),0,0.000001),alt="greater")$p.value
                maxdiff.output[i] <- (total.output[i]+meandiff.hetero[i])/2
                mindiff.output[i] <- (total.output[i]-meandiff.hetero[i])/2
                if(is.na(min.output[i])|min.output[i]!=0)
                  {
                    if(pvalues.equality[i]>alpha.equality)
                      {
                        min.output[i] <- total.output[i]/2
                        max.output[i] <- total.output[i]/2
                      }
                    else
                      {
                        max.output[i] <- maxdiff.output[i]
                        min.output[i] <- mindiff.output[i]
                      }
                  }
              }
          }#End of n.segments loop
      }
    
    list(total.output=total.output,min.output=min.output,max.output=max.output,n.hetero=n.hetero,mindiff.output=mindiff.output,maxdiff.output=maxdiff.output,mean.diff.mbaf=mean.diff.mbaf,confidence95NewModes=(confidence95NewModes+modeOffset))
  }
