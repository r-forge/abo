#library(splines)
#
#cubic.splines <- function(y,markers.per.knot)
#  {
#    n <- length(y)
#    p <- round(n/markers.per.knot)
#    knot.distances <- n/(p+1)
#    knot.positions <- (1:p)*knot.distances
#    indices <- 1:n
#    bs.fit <- bs(indices,knots=knot.positions)
#    bs.lm <- lm(y~bs.fit)$fitted.values
#    bs.lm
#  }

compute.chisquare <- function(abs.diff.data,quantile.values,i,n,q,window.size,window.size.minusone,window.size.onefifth,window.size.negativeonefifth,theta.matrix,gamma.matrix)
  {
    f <- rep(NA,q)
    g <- rep(NA,q)
    scaled.theta.matrix <- theta.matrix
    scaled.gamma.matrix <- gamma.matrix
    left.endpoint <- i-window.size.minusone
    if(left.endpoint>0) left.data <- abs.diff.data[left.endpoint:i]
    else left.data <- abs.diff.data[c((n+left.endpoint):n,1:i)]
    left.data.quantiles <- quantile(left.data,c(quantile.values,0.25,0.75))
    right.endpoint <- i+window.size
    if(right.endpoint<=n) right.data <- abs.diff.data[(i+1):right.endpoint]
    else if (i<n) right.data <- abs.diff.data[c((i+1):n,1:(right.endpoint-n))]
    else right.data <- abs.diff.data[1:window.size]
    right.data.quantiles <- quantile(right.data,c(quantile.values,0.25,0.75))
    QF <- 2*(left.data.quantiles[q+2]-left.data.quantiles[q+1])
    for(l in 1:q)
      {
        K <- 1-abs(left.data.quantiles[l]-left.data)/(window.size.negativeonefifth*QF)
        K[K<0] <- 0
        f[l] <- sum(window.size.onefifth*(QF^(-1))*K)/window.size
      }
    for(l in 1:q)
      {
        for(m in 1:q)
          {
            scaled.theta.matrix[l,m] <- theta.matrix[l,m]/(window.size*f[l]*f[m])
          }
      }
    QG <- 2*(right.data.quantiles[q+2]-right.data.quantiles[q+1])
    for(l in 1:q)
      {
        K <- 1-abs(right.data.quantiles[l]-right.data)/(window.size.negativeonefifth*QG)
        K[K<0] <- 0
        g[l] <- sum(window.size.onefifth*(QG^(-1))*K)/window.size
      }
    for(l in 1:q)
      {
        for(m in 1:q)
          {
            scaled.gamma.matrix[l,m] <- gamma.matrix[l,m]/(window.size*g[l]*g[m])
          }
      }
    psi.matrix <- (2*window.size)*(scaled.theta.matrix+scaled.gamma.matrix)
#        write(abs.diff.data,"abs-diff-data.txt",ncol=1)
#        write(i,"i.txt",ncol=1)
#        write(left.endpoint,"left-endpoint.txt",ncol=1)
#        write(right.endpoint,"right-endpoint.txt",ncol=1)
#        write(left.data.quantiles,"left-data-quantiles.txt",ncol=1)
#        write(right.data.quantiles,"right-data-quantiles.txt",ncol=1)
#        write(QF,"QF.txt",ncol=1)
#        write(QG,"QG.txt",ncol=1)
#        write(f,"f.txt",ncol=1)
#        write(g,"g.txt",ncol=1)
#        write.table(psi.matrix,"psi-matrix.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
#        write.table(scaled.theta.matrix,"scaled-theta-matrix.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
#        write.table(scaled.gamma.matrix,"scaled-gamma-matrix.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)    
    d <- sqrt(2*window.size)*(left.data.quantiles-right.data.quantiles)[1:q]    
    output <- d%*%solve(psi.matrix)%*%d
    output
  }
                              
#New code to implement Kosorok's quantile test
#Now only works for on window size, see older code for multiple window sizes
#abs.diff.data - absolute differences of A and B alleles
#quantile.values - quantiles of distribution being compared
#window.size - size of windows for comparing left to right distributions
#frequency.sample.markers - how often to compute p-values before smoothing
#pvalue.calculate - set all statistics above this p-alue to 0

quantile.test.windows <- function(abs.diff.data,quantile.values=seq(.02,0.1,0.02),window.size=200,frequency.sample.markers=10,pvalue.calculate=0.001)
  {
    n <- length(abs.diff.data)
    quantile.values <- sort(quantile.values)
    q <- length(quantile.values)
    window.size.minusone <- window.size-1
    window.size.onefifth <- window.size^(1/5)
    window.size.negativeonefifth <- window.size^(-1/5)
    chisquare.vector <- rep(0,n)
    theta.matrix <- matrix(NA,q,q)
    scaled.theta.matrix <- matrix(NA,q,q)
    gamma.matrix <- matrix(NA,q,q)
    scaled.gamma.matrix <- matrix(NA,q,q)
    f <- rep(NA,q)
    g <- rep(NA,q)
    for(i in 1:q)
      {        
        for(j in 1:q)
          {
            
            theta.matrix[i,j] <- gamma.matrix[i,j] <- min(quantile.values[i],quantile.values[j])-quantile.values[i]*quantile.values[j]
          }
      }
    for(i in seq(2,n-2,frequency.sample.markers))            
      {
        chisquare.vector[i] <- compute.chisquare(abs.diff.data,quantile.values,i,n,q,window.size,window.size.minusone,window.size.onefifth,window.size.negativeonefifth,theta.matrix,gamma.matrix)
      }#End of i in 1:n
#Now fill in chisquare values where pvalues are low
    chisquare.cutoff <- qchisq(1-pvalue.calculate,q)
    chisquare.vector[chisquare.vector<chisquare.cutoff] <- 0
    which.notzero <- which(chisquare.vector>0)
    if(length(which.notzero)>0)
      {
        new.i <- NULL
        frequency.sample.markers.minusone <- frequency.sample.markers-1
        for(i in 1:length(which.notzero))
          {
            new.i <- c(new.i,(which.notzero[i]-frequency.sample.markers.minusone):(which.notzero[i]+frequency.sample.markers.minusone))
          }
#Get rid of duplicates, those outside range, and those already computed        
        new.i <- unique(new.i)
        new.i <- new.i[new.i>=2 & new.i<=(n-2)]
        new.i <- new.i[is.na(match(new.i,which.notzero))]
        for(i in new.i)
          {
            chisquare.vector[i] <- compute.chisquare(abs.diff.data,quantile.values,i,n,q,window.size,window.size.minusone,window.size.onefifth,window.size.negativeonefifth,theta.matrix,gamma.matrix)
          }#End of i in 1:n
      }
    chisquare.vector
  }

simes.significant <- function(chisquare.statistics,degrees.of.freedom,alpha.level,n)
  {
    order.statistics <- order(chisquare.statistics)
    chisquare.statistics.ordered <- chisquare.statistics[order.statistics]
    pvalues.final <- pvalues.chisquare.statistics.ordered <- 1-pchisq(chisquare.statistics.ordered,degrees.of.freedom)
    pvalues.adjusted.ordered <- n*pvalues.chisquare.statistics.ordered/(1:n)
    pvalues.adjusted.ordered.final <- pvalues.adjusted.ordered
    for(i in n:2)
      {
        if(pvalues.adjusted.ordered.final[i-1]< pvalues.adjusted.ordered.final[i])
          {             
            pvalues.adjusted.ordered.final[i-1] <- pvalues.adjusted.ordered.final[i]
          }
      }
    pvalues.final[order.statistics] <- pvalues.adjusted.ordered.final
    output <- rep(FALSE,n)
    output[pvalues.final<alpha.level] <- TRUE
    output
  }

find.left.min <- function(data.vector,simes.indicator)
  {
    n <- length(data.vector)
    diff.vector <- data.vector[1:(n-1)]-data.vector[2:n]
    simes.indicator <- simes.indicator[1:(n-1)]
    max(1,which(diff.vector>=0 & simes.indicator==FALSE))
  }
    
find.right.min <- function(data.vector,simes.indicator)
  {
    n <- length(data.vector)
    diff.vector <- data.vector[2:n]-data.vector[1:(n-1)]
    simes.indicator <- simes.indicator[1:(n-1)]
    min(which(diff.vector>=0 & simes.indicator==FALSE),n)
  }

#abs.diff.data - absolute differences of A and B alleles
#alpha.level - alpha-level of segmentation
#window.size - size of windows for comparing left to right distributions
#quantile.values - quantiles of distribution being compared
#frequency.sample.markers - how often to compute p-values before smoothing
#f - amount of smoothing in lowess

second.round <- function(abs.diff.data,normal.data,alpha.level,window.size,frequency.sample.markers=10)
{
#  chisquare.statistics <- quantile.test.windows(abs.diff.data,quantile.values,window.size,frequency.sample.markers=10,pvalue.calculate=alpha.level)
  n <- length(abs.diff.data)
  chisquare.statistics[1] <- 0
    chisquare.statistics[(n-1):n] <- 0
  degrees.of.freedom <- length(quantile.values)
  spline.statistics <- smooth.spline(chisquare.statistics)$y
#  spline.statistics <- cubic.splines(chisquare.statistics)
#  lowess.statistics <- lowess(chisquare.statistics,markers.per.knot=markers.per.knot,delta=1/n)$y  
  simes.output <- simes.significant(spline.statistics,degrees.of.freedom,alpha.level,n)
#  simes.output <- simes.significant(lowess.statistics,degrees.of.freedom,alpha.level,n)
  output <- NULL
  number.significant <- sum(simes.output)
  while(number.significant>0)
    {
      which.max <- which(spline.statistics==max(spline.statistics) & simes.output)[1]
#      which.max <- which(lowess.statistics==max(lowess.statistics) & simes.output)[1]      
      if(!is.na(which.max))
        {
          output <- c(output,which.max)
          left.min <- find.left.min(spline.statistics[1:which.max],simes.output[1:which.max])
          right.min <- find.right.min(spline.statistics[(which.max+1):n],simes.output[(which.max+1):n])
#          left.min <- find.left.min(lowess.statistics[1:which.max])
#          right.min <- find.right.min(lowess.statistics[(which.max+1):n])
          left.delete <- left.min:which.max
          right.delete <- (which.max+1):(which.max+right.min)
          simes.output[left.delete] <- FALSE
          simes.output[right.delete] <- FALSE
          spline.statistics[left.delete] <- 0
          spline.statistics[right.delete] <- 0
#          lowess.statistics[left.delete] <- 0
#          lowess.statistics[right.delete] <- 0
          number.significant <- sum(simes.output)
        }
      else
        {
          number.significant <- 0
        }
    }
  if(!is.null(output))
    {
      output <- sort(output)
      output <- output[output>1 & output<(n-1)]
    }
  output
}

second.round.cbs <- function(diff.data,normal.data,alpha)
  {
    t.data <- diff.data-normal.data
    n <- length(t.data)
    new.cna <- CNA(t.data,rep(1,n),1:n)
    new.segment <- segment(new.cna,alpha=alpha)$output
    if(nrow(new.segment)==1) output <- NULL
    else output <- cumsum(new.segment[,5])[-nrow(new.segment)]
    output
  }

second.round.cbs.matching <- function(diff.data,alpha)
  {
    n <- length(diff.data)
    new.cna <- CNA(diff.data,rep(1,n),1:n)
    new.segment <- segment(new.cna,alpha=alpha)$output
    if(nrow(new.segment)==1) output <- NULL
    else output <- cumsum(new.segment[,5])[-nrow(new.segment)]
    output
  }

#abs.diff.data - absolute differences of A and B alleles
#alpha.level - alpha-level of segmentation
#window.size - size of windows for comparing left to right distributions
#quantile.values - quantiles of distribution being compared
#frequency.sample.markers - how often to compute p-values before smoothing
#f - amount of smoothing in lowess
#input.matrix - the results of the first round of segmentation
#map.data - map data
secondsegment <- function(abs.diff.data,normal.data,alpha.level,window.size,quantile.values=seq(0.02,.1,0.02),frequency.sample.markers=10,input.matrix,map.data,first.LOH)
  {
    output.matrix <- as.matrix(input.matrix)
    colnames(output.matrix) <- colnames(input.matrix)
    n <- nrow(input.matrix)
    new.ends <- cumsum(input.matrix[,5])
    new.starts <- c(1,new.ends[-n]+1)
    minimum.markers.second.round <- 2*window.size+1
    counter <- 1
    for(i in 1:n)
      {
        new.row <- input.matrix[i,]
        if(new.row[5]>=20 & !first.LOH[i])
#        if(new.row[5]>=minimum.markers.second.round)          
          {
            new.abs.diff.data <- abs.diff.data[new.starts[i]:new.ends[i]]
            new.normal.data <- normal.data[new.starts[i]:new.ends[i]]
            new.map.data <- map.data[new.starts[i]:new.ends[i]]
#            new.f <- f/new.row[5]
            second.output <- second.round.cbs(new.abs.diff.data,new.normal.data,alpha.level)
#            second.output <- second.round(new.abs.diff.data,alpha.level,window.size,quantile.values,frequency.sample.markers)            
#            second.output <- second.round(new.abs.diff.data,alpha.level,window.size,quantile.values,frequency.sample.markers)
            print(second.output)
            if(length(second.output)>0)
              {
                new.output <- matrix(rep(as.character(new.row),length(second.output)+1),ncol=6,byrow=TRUE)
                ends.output <- c(second.output,length(new.abs.diff.data))
                starts.output <- c(1,ends.output[-length(ends.output)]+1)
                new.output[,3] <- new.map.data[starts.output]
                new.output[,4] <- new.map.data[ends.output]
                new.output[,5] <- ends.output-starts.output+1
                if(i>1 & i<n)
                  {
                    output.matrix <- rbind(output.matrix[1:(counter-1),],new.output,output.matrix[(counter+1):nrow(output.matrix),])                    
                  }
                else if(i==1)
                  {
                    if(n>1)
                      {
                        output.matrix <- rbind(new.output,output.matrix[2:n,])
                      }
                    else
                      output.matrix <- new.output
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
    output.matrix
  }

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
        if(matching) which.het.data <- which(new.normal.data)
        else which.het.data <- 1:length(new.normal.data)
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

