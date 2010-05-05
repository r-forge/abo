psSegment <- function(x, matching.reference=FALSE, alpha=0.01, alpha1=0.009, het.lower=0.1, het.upper=0.75, het.minimum.window=0.05, het.stepsize.window=0.01, het.stepsize.within=0.01, trim.mean=0.1, min.hetero=5, alpha.homozygous=0.05, alpha.equality=0.05, maf=0.075, min.homo.region=100, modeOffset=0.025, smooth.segmentation=TRUE, smooth.region=2, outlier.SD.scale=4, smooth.SD.scale=2, trim=0.025, impute.LOH, zero.homo=FALSE,verbose=1, ...)
  {
    if (!inherits(x, 'psCNA')) stop("First arg must be an psCNA object")
    if(alpha1>=alpha) stop ("alpha1 must be less than alpha")
    call <- match.call()
    alpha2 <- alpha-alpha1
    sampleid <- colnames(x[[1]])
    nsample <- ncol(x[[1]])
    nprobe <-  nrow(x[[1]])
    nsample.normal <- ncol(x[[5]])
    if(matching.reference & nsample.normal!=nsample) stop ("If matching normal, number of test and reference samples must be the same")
    if(!matching.reference & is.null(x[[8]]))
      {
        count.heterozygous <- rep(0,nprobe)
        typical.mbaf.heterozygous <- rep(0,nprobe)
        for(i in 1:nsample.normal)
          {
            normal.a.vector <- x[[5]][,i]
            normal.b.vector <- x[[6]][,i]
            normal.baf.vector <- normal.b.vector/(normal.a.vector+normal.b.vector)
            normal.baf.vector[abs(normal.baf.vector)==Inf] <- 1             
            normal.mbaf.vector <- apply(cbind(normal.baf.vector,1-normal.baf.vector),1,max)
            normal.mbaf.vector[abs(normal.mbaf.vector)==Inf] <- 1                
            normal.fit.heterozygous <- findheterozygous(normal.a.vector,normal.b.vector,het.lower,het.upper,het.minimum.window,het.stepsize.window,het.stepsize.within)
            normal.het.vector <- which(normal.fit.heterozygous$is.heterozygous)
            count.heterozygous[normal.het.vector] <- count.heterozygous[normal.het.vector] + 1
            typical.mbaf.heterozygous[normal.het.vector] <- typical.mbaf.heterozygous[normal.het.vector] + normal.mbaf.vector[normal.het.vector]
          }
        reference.mbaf.heterozygous <- typical.mbaf.heterozygous/count.heterozygous
        reference.mbaf.heterozygous[is.na(reference.mbaf.heterozygous)] <- mean(reference.mbaf.heterozygous,na.rm=TRUE)
      }      
    if(is.null(x[[7]])) predict.genomdat.heterozygous <- TRUE
    else predict.genomdat.heterozygous <- FALSE
    if(is.null(x[[8]])) predict.normal.heterozygous <- TRUE
    else predict.normal.heterozygous <- FALSE
    chrom.vector <- x[[3]]
    position.vector <- x[[4]]
    for(i in 1:nsample)
      {
        a.vector <- x[[1]][,i]
        b.vector <- x[[2]][,i]
        baf.vector <- b.vector/(a.vector+b.vector)
        baf.vector[abs(baf.vector)==Inf] <- 1             
        mbaf.vector <- apply(cbind(baf.vector,1-baf.vector),1,max)
        mbaf.vector[abs(mbaf.vector)==Inf] <- 1
        if(matching.reference)
          {
            if(predict.normal.heterozygous)
              {
                normal.a.vector <- x[[5]][,i]
                normal.b.vector <- x[[6]][,i]
                normal.baf.vector <- normal.b.vector/(normal.a.vector+normal.b.vector)
                normal.mbaf.vector <- apply(cbind(normal.baf.vector,1-normal.baf.vector),1,max)
                normal.mbaf.vector[abs(normal.mbaf.vector)==Inf] <- 1                
                which.notmissing <- which(!is.na(a.vector) & !is.na(b.vector) & !is.na(chrom.vector) & !is.na(position.vector) & !is.na(normal.a.vector) & !is.na(normal.b.vector))
                normal.a.vector.notmissing <- normal.a.vector[which.notmissing]
                normal.b.vector.notmissing <- normal.b.vector[which.notmissing]
                normal.baf.vector.notmissing <- normal.baf.vector[which.notmissing]
                normal.mbaf.vector.notmissing <- normal.mbaf.vector[which.notmissing]                
                normal.fit.heterozygous <- findheterozygous(normal.a.vector.notmissing,normal.b.vector.notmissing,het.lower,het.upper,het.minimum.window,het.stepsize.window,het.stepsize.within)
                normal.het.vector.notmissing <- normal.fit.heterozygous$is.heterozygous
              }
            else
              {
                normal.het.vector <- x[[8]][,i]
                which.notmissing <- which(!is.na(a.vector) & !is.na(b.vector) & !is.na(chrom.vector) & !is.na(position.vector) & !is.na(normal.het.vector))
#                which.notmissing <- which(!is.na(a.vector) & !is.na(b.vector) & !is.na(chrom.vector) & !is.na(position.vector) & !is.na(normal.het.vector) & !is.na(genomdat.het.vector))                
                normal.het.vector.notmissing <- normal.het.vector[which.notmissing]
#                stop("Code not yet written for normal genotyping data; must use normal copy number data")
              }
          }
        else
          {
            which.notmissing <- which(!is.na(a.vector) & !is.na(b.vector) & !is.na(chrom.vector) & !is.na(position.vector))
          }
        a.vector.notmissing <- a.vector[which.notmissing]
        b.vector.notmissing <- b.vector[which.notmissing]
        baf.vector.notmissing <- baf.vector[which.notmissing]
        mbaf.vector.notmissing <- mbaf.vector[which.notmissing]        
        chrom.vector.notmissing <- chrom.vector[which.notmissing]
        position.vector.notmissing <- position.vector[which.notmissing]
        if(predict.genomdat.heterozygous)
          {
            fit.heterozygous <- findheterozygous(a.vector.notmissing,b.vector.notmissing,het.lower,het.upper,het.minimum.window,het.stepsize.window,het.stepsize.within)
            het.vector.notmissing <- fit.heterozygous$is.heterozygous
          }
        else
          {
            het.vector <- x[[7]][,i]
            het.vector.notmissing <- het.vector[which.notmissing]
#            stop("Code not yet written for test genotyping data; must use test copy number data")
          }        
#Make the smaller value of all homozygotes be zero if desired
        if(zero.homo)
          {
            min.indicator <- rep(FALSE,length(a.vector.notmissing))
            min.indicator[b.vector.notmissing<a.vector.notmissing] <- TRUE
            a.vector.notmissing[!het.vector.notmissing & !min.indicator] <- 0
            b.vector.notmissing[!het.vector.notmissing & min.indicator] <- 0
          }
        sum.notmissing <- a.vector.notmissing+b.vector.notmissing
#Must make negatives zero else thrown out, but negatives come back later
        sum.notmissing[sum.notmissing<0] <- 0
        sqrt.sum.notmissing <- sqrt(sum.notmissing)
        CNA.sqrt.sum.notmissing<- CNA(sqrt.sum.notmissing,chrom.vector.notmissing,position.vector.notmissing,sampleid=sampleid[i])
#By default smoothing is done
        if(smooth.segmentation)
          {
            smoothed.CNA.sqrt.sum.notmissing <- smooth.CNA(CNA.sqrt.sum.notmissing,smooth.region=smooth.region,outlier.SD.scale=outlier.SD.scale,smooth.SD.scale=smooth.SD.scale,trim=trim)
#            smoothed.CNA.sqrt.sum.notmissing <- smooth.CNA(CNA.sqrt.sum.notmissing,...)            
            segmented.smoothed.CNA.sqrt.sum.notmissing <- segment(smoothed.CNA.sqrt.sum.notmissing,alpha=alpha1,verbose=verbose,...)
          }
        else segmented.smoothed.CNA.sqrt.sum.notmissing <- segment(CNA.sqrt.sum.notmissing,alpha=alpha1,verbose=verbose,...)
        segmented.output <- segmented.smoothed.CNA.sqrt.sum.notmissing$output
        segment.lengths <- as.numeric(segmented.output[,5])
        segment.ends <- cumsum(segment.lengths)
        n.segment <- length(segment.ends)
        segment.starts <- c(1,segment.ends[-n.segment]+1)
        chrom.vector.starts <- chrom.vector.notmissing[segment.starts]
#        first.LOH <- first.call(segment.starts,segment.ends,het.vector.notmissing,chrom.vector.starts,maf,homozygous.pvalue.cutoff)
#        print(first.LOH)
#        segmented.output.LOH <- cbind(segmented.output,first.LOH)
#Second round of segmentation
        if(matching.reference)
          {
            if(predict.normal.heterozygous)
              {
                tumorBoostBAFNotmissing <- 0.5*baf.vector.notmissing/normal.baf.vector.notmissing
                whichTumorGreater <- which(baf.vector.notmissing>normal.baf.vector.notmissing)
                tumorBoostBAFNotmissing[whichTumorGreater] <- (1-0.5*((1-baf.vector.notmissing)/(1-normal.baf.vector.notmissing)))[whichTumorGreater]
                diff.mbaf.vector.notmissing <- mbaf.vector.notmissing-normal.mbaf.vector.notmissing
              }
            else
              {
                tumorBoostBAFNotmissing <- baf.vector.notmissing
                diff.mbaf.vector.notmissing <- rep(NA,length(baf.vector.notmissing))
              }
            foldedTumorBoostBAFNotmissing <- 2*abs(tumorBoostBAFNotmissing-0.5)
            second.output <- second.segment.matching(matching.reference,foldedTumorBoostBAFNotmissing,normal.het.vector.notmissing,alpha2,segmented.output,position.vector.notmissing)
#            second.output <- second.segment.matching(matching.reference,diff.mbaf.vector.notmissing,normal.het.vector.notmissing,alpha2,segmented.output,position.vector.notmissing)            
          }
        else
          {
            reference.mbaf.heterozygous.notmissing <- reference.mbaf.heterozygous[which.notmissing]
            diff.mbaf.vector.notmissing <- mbaf.vector.notmissing-reference.mbaf.heterozygous.notmissing
            second.output <- second.segment.matching(matching.reference,diff.mbaf.vector.notmissing,het.vector.notmissing,alpha2,segmented.output,position.vector.notmissing)
          }

#            print("second.output")
#            print(second.output)
#            write.table(second.output,"second-output.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

#            second.output <- secondsegment(het.vector.notmissing,normal.het.vector.notmissing,alpha2,window.size=window.size,quantile.values=quantile.values,frequency.sample.markers=frequency.sample.markers,segmented.output,position.vector.notmissing,first.LOH)
#        second.output <- secondsegment(abs(a.vector.notmissing-b.vector.notmissing),alpha2,window.size=window.size,quantile.values=quantile.values,frequency.sample.markers=frequency.sample.markers,segmented.output,position.vector.notmissing)        
        segment.lengths <- as.numeric(second.output[,5])
#        print(paste("segment.lengths=",segment.lengths))
        segment.ends <- cumsum(segment.lengths)
#        print(paste("segment.ends=",segment.ends))
        n.segment <- length(segment.ends)
#        print(paste("n.segment=",n.segment))
        segment.starts <- c(1,segment.ends[-n.segment]+1)
#        print(paste("segement.starts=",segment.starts))
        chrom.vector.starts <- as.numeric(chrom.vector.notmissing)[segment.starts]
#        print(paste("chrom.vector.starts=",chrom.vector.starts))
        total.vector.starts <- as.numeric(second.output[,6])
#        print(paste("total.vector.starts=",total.vector.starts))
#        print("Before called segments")
        if(matching.reference)
          {
            called.segments <- callsegments.matching(matching=matching.reference,segment.starts=segment.starts,segment.ends=segment.ends,a.vector=a.vector.notmissing,b.vector=b.vector.notmissing,diff.mbaf.vector=diff.mbaf.vector.notmissing,het.vector=het.vector.notmissing,normal.het.vector=normal.het.vector.notmissing,foldedTumorBoostBAF=foldedTumorBoostBAFNotmissing,chrom.vector.starts=chrom.vector.starts,total.vector.starts=total.vector.starts,min.hetero=min.hetero,trim.mean=trim.mean,alpha.equality=alpha.equality,modeOffset=modeOffset,maf=1-maf,alpha.homozygous=alpha.homozygous,min.homo.region=min.homo.region,impute.LOH=impute.LOH)
          }
        else
          {
            foldedTumorBoostBAFNotmissing <- diff.mbaf.vector.notmissing
            called.segments <- callsegments.matching(matching=matching.reference,segment.starts=segment.starts,segment.ends=segment.ends,a.vector=a.vector.notmissing,b.vector=b.vector.notmissing,diff.mbaf.vector=diff.mbaf.vector.notmissing,het.vector=het.vector.notmissing,normal.het.vector=normal.het.vector.notmissing,foldedTumorBoostBAF=foldedTumorBoostBAFNotmissing,chrom.vector.starts=chrom.vector.starts,total.vector.starts=total.vector.starts,min.hetero=min.hetero,trim.mean=trim.mean,alpha.equality=alpha.equality,modeOffset=modeOffset,maf=1-maf,alpha.homozygous=alpha.homozygous,min.homo.region=min.homo.region,impute.LOH=impute.LOH)
          }
#        print("After called segments")
        if(i==1)
          {
#            print("before output")
            output <- cbind(matrix(second.output[,1:5],ncol=5),called.segments$n.hetero,called.segments$mean.diff.mbaf,called.segments$mindiff.output,called.segments$maxdiff.output,called.segments$confidence95NewModes,called.segments$min.output,called.segments$max.output,called.segments$total.output)
#            print("after output")
          }
        else
          {
            output <- rbind(output,called.segments$n.hetero,called.segments$mean.diff.mbaf,called.segments$mindiff.output,called.segments$maxdiff.output,called.segments$confidence95NewModes,called.segments$min.output,called.segments$max.output,called.segments$total.output)
          }
      }
#Merge 0s
    if(nrow(output)>1)
#    if(nrow(output)>0)      
      {
        output <- as.matrix(output)
        n.segments <- nrow(output)
        merge.index <- rep(NA,n.segments)
        merge.index[1] <- 1
        merge.same <- rep(FALSE,n.segments-1)
        merge.same[(output[2:n.segments,11]==output[1:(n.segments-1),11]) & (output[2:n.segments,13]==output[1:(n.segments-1),13])] <- TRUE
#        merge.same[(output[2:n.segments,6]==output[1:(n.segments-1),6]) & (output[2:n.segments,8]==output[1:(n.segments-1),8])] <- TRUE
        for(i in 2:n.segments)
          {
            if(merge.same[i-1]) merge.index[i] <- merge.index[i-1]
            else merge.index[i] <- merge.index[i-1]+1
          }
        new.output <- matrix(NA,max(merge.index),13)
        for(i in 1:max(merge.index))
          {
            which.i <- which(merge.index==i)
            which.i.1 <- which.i[1]
            new.output[i,c(1:3,7:13)] <- output[which.i.1,c(1:3,7:13)]
            new.output[i,4] <- output[max(which.i),4]
            new.output[i,5] <- sum(as.numeric(output[which.i,5]))
            new.output[i,6] <- sum(as.numeric(output[which.i,6]))
          }
        output <- as.data.frame(new.output)
#        new.output <- matrix(NA,max(merge.index),8)
#        for(i in 1:max(merge.index))
#          {
#            which.i <- which(merge.index==i)
#            which.i.1 <- which.i[1]
#            new.output[i,c(1:3,6:8)] <- output[which.i.1,c(1:3,6:8)]
#            new.output[i,4] <- output[max(which.i),4]
#            new.output[i,5] <- sum(as.numeric(output[which.i,5]))
#          }
#        output <- as.data.frame(new.output)
      }
    else output <- as.data.frame(output)
    colnames(output) <- c("ID","chrom","loc.start","loc.end","num.mark","num.hetero","mean.diff.mbaf","mindiff.mean","maxdiff.mean","confidence95","min.mean","max.mean","total.mean")
    output
  }
