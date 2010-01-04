# Chromosome.Lengths <- c(263, 255, 214, 203, 194, 183, 171, 155, 145, 144, 144, 143, 114, 109, 106, 98, 92, 85, 67, 72, 50, 56, 164, 59)
# names(Chromosome.Lengths) <- c(as.character(1:22),"X","Y")

#New methods for allele-specific segmentation


###########################################################################/**
# @RdocFunction as.CNA
#
# @title "Creating an 'as.CNA' object"
#
# \description{
#  @get "title".
# }
# 
# @synopsis
#
# \arguments{
#   \item{genomdat.a, genomdat.b}{}
#   \item{chrom, maploc}{}
#   \item{normaldat.a, normaldat.b}{}
#   \item{genomdat.hetmatrix}{}
#   \item{normal.hetmatrix}{}
#   \item{sampleid}{}
#   \item{arraytype}{}
# }
#
# \value{
#   Returns a @list structure of class \code{as.CNA}.
# }
#
# @author
#
# \seealso{
#   @see "as.segment".
# }
#*/########################################################################### 
as.CNA <- function(genomdat.a, genomdat.b, chrom, maploc, normaldat.a, normaldat.b, genomdat.hetmatrix=NULL, normal.hetmatrix=NULL, sampleid=NULL, arraytype=NULL)
  {
    if (!is.numeric(genomdat.a)) stop("genomdat.a must be numeric")
    if (!is.numeric(genomdat.b)) stop("genomdat.b must be numeric")
    if (is.factor(chrom)) chrom <- as.character(chrom)
    if (!is.numeric(maploc)) stop("maploc must be numeric")
    if (is.vector(genomdat.a)) genomdat.a <- as.matrix(genomdat.a)
    if (is.vector(genomdat.b)) genomdat.b <- as.matrix(genomdat.b)
    n1 <- nrow(genomdat.a)
    n2 <- nrow(genomdat.b)
    p1 <- ncol(genomdat.a)
    p2 <- ncol(genomdat.b)
    n3 <- length(chrom)
    n4 <- length(maploc)
    if(n1!=n2) stop(paste("number of rows genomdat.a=",n1," number of rows genomdat.b=",n2))
    if(p1!=p2) stop(paste("number of cols genomdat.a=",p1," number of cols genomdat.b=",p2))
    if(n1!=n3) stop(paste("number of rows genomdat.a=",n1," length chrom=",n3))
    if(n1!=n4) stop(paste("number of rows genomdat.a=",n1," length chrom=",n4))
    if (is.vector(normaldat.a)) normaldat.a <- as.matrix(normaldat.a)
    if (is.vector(normaldat.b)) normaldat.b <- as.matrix(normaldat.b)
    n5 <- nrow(normaldat.a)
    n6 <- nrow(normaldat.b)
    p5 <- ncol(normaldat.a)
    p6 <- ncol(normaldat.b)
    if(n1!=n5) stop(paste("number of rows genomdat.a=",n1," number of rows normaldat.a=",n5))
    if(n1!=n6) stop(paste("number of rows genomdat.a=",n1," number of rows normaldat.b=",n6))
#        if(p1!=p5) stop(paste("number of cols genomdat.a=",p1," number of cols normaldat.a=",p5))
#        if(p1!=p6) stop(paste("number of cols genomdat.a=",p1," number of cols normaldat.b=",p6))
    if(!is.null(genomdat.hetmatrix))
      {
        if(!is.logical(genomdat.hetmatrix)) stop("genomdat.hetmatrix must be of type logical (TRUE/FALSE)")
        if (is.vector(genomdat.hetmatrix)) genomdat.hetmatrix <- as.matrix(genomdat.hetmatrix)

        n7 <- nrow(genomdat.hetmatrix)
        p7 <- ncol(genomdat.hetmatrix)
        if(n1!=n7) stop(paste("number of rows genomdat.a=",n1," number of rows genomdat.hetmatrix=",n7))
        if(p1!=p7) stop(paste("number of rows genomdat.a=",n1," number of rows genomdat.hetmatrix=",p7))
      }
    if(!is.null(normal.hetmatrix))
      {
        if(!is.logical(normal.hetmatrix)) stop("normal.hetmatrix must be of type logical (TRUE/FALSE)")
        if (is.vector(normal.hetmatrix)) normal.hetmatrix <- as.matrix(normal.hetmatrix)
        n8 <- nrow(normal.hetmatrix)
        p8 <- nrow(normal.hetmatrix)
        if(n1!=n8) stop(paste("number of rows normal.a=",n1," number of rows normal.hetmatrix=",n8))
#        if(p1!=p8) stop(paste("number of rows normal.a=",p1," number of rows normal.hetmatrix=",p8))
      }    

    if (length(sampleid) != ncol(genomdat.a))
      {
        if(!is.null(sampleid)) warning("length(sampleid) and ncol(genomdat.a) differ, names ignored\n")
        sampleid <- paste("Sample", 1:ncol(genomdat.a))
      }
    
    if (sum(is.na(chrom)|is.na(maploc))>0)
      warning("markers with missing chrom and/or maploc removed\n")
    sortindex <- order(chrom, maploc, na.last=NA)

    zzz <- vector("list",8)
    zzz[[1]] <- as.matrix(genomdat.a[sortindex,])
    zzz[[2]] <- as.matrix(genomdat.b[sortindex,])
    colnames(zzz[[1]]) <- sampleid
    colnames(zzz[[2]]) <- sampleid
    zzz[[3]] <- chrom[sortindex]
    zzz[[4]] <- maploc[sortindex]
    if(!is.null(normaldat.a) & !is.null(normaldat.b))
      {
        zzz[[5]] <- as.matrix(normaldat.a[sortindex,])
        zzz[[6]] <- as.matrix(normaldat.b[sortindex,])
        colnames(zzz[[5]]) <- sampleid
        colnames(zzz[[6]]) <- sampleid
      }
    if(!is.null(genomdat.hetmatrix))
      {
        zzz[[7]] <- as.matrix(genomdat.hetmatrix[sortindex,])
        colnames(zzz[[7]]) <- sampleid
      }
    if(!is.null(normal.hetmatrix))
      {
        zzz[[8]] <- as.matrix(normal.hetmatrix[sortindex,])
        colnames(zzz[[8]]) <- sampleid
      }
    class(zzz) <- c("as.CNA","list")
    names(zzz) <- c("genomdat.a","genomdat.b","chrom","maploc","normaldat.a","normaldat.b","genomdat.hetmatrix","normal.hetmatrix")
    zzz
  }


