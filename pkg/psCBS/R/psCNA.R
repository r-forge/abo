psCNA <- function(genomdat.a, genomdat.b, chrom, maploc, normaldat.a, normaldat.b, genomdat.hetmatrix=NULL, normaldat.hetmatrix=NULL, sampleid=NULL, arraytype=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'genomdat.a':
  if (!is.numeric(genomdat.a)) {
    stop("genomdat.a must be numeric: ", mode(genomdat.a))
  }
  if (is.vector(genomdat.a)) {
    genomdat.a <- as.matrix(genomdat.a)
  }
  nbrOfLoci <- nrow(genomdat.a);

  # Argument 'genomdat.b':
  if (!is.numeric(genomdat.b)) {
    stop("genomdat.b must be numeric: ", mode(genomdat.b))
  }
  if (is.vector(genomdat.b)) {
    genomdat.b <- as.matrix(genomdat.b)
  }

  n2 <- nrow(genomdat.b)
  if (n2 != nbrOfLoci) {
    stop("number of rows genomdat.a=", nbrOfLoci, " number of rows genomdat.b=", n2)
  }

  p1 <- ncol(genomdat.a)
  p2 <- ncol(genomdat.b)
  if (p1 != p2) {
    stop("number of cols genomdat.a=", p1, " number of cols genomdat.b=", p2)
  }


  # Argument 'chrom':
  if (is.factor(chrom)) {
    chrom <- as.character(chrom)
  }
  n3 <- length(chrom)
  if (n3 != nbrOfLoci) {
    stop("number of rows genomdat.a=", nbrOfLoci, " length chrom=", n3)
  }


  # Argument 'maploc':
  if (!is.numeric(maploc)) {
    stop("maploc must be numeric: ", mode(maploc))
  }
  n4 <- length(maploc)
  if (n4 != nbrOfLoci) {
    stop("number of rows genomdat.a=", nbrOfLoci, " length chrom=", n4)
  }


  # Argument 'normaldat.a':
  if (is.vector(normaldat.a)) {
    normaldat.a <- as.matrix(normaldat.a)
  }
  n5 <- nrow(normaldat.a)
  if (n5 != nbrOfLoci) {
    stop("number of rows genomdat.a=", nbrOfLoci, " number of rows normaldat.a=", n5)
  }

  # Argument 'normaldat.b':
  if (is.vector(normaldat.b)) {
    normaldat.b <- as.matrix(normaldat.b)
  }
  n6 <- nrow(normaldat.b)
  if (n6 != nbrOfLoci) {
    stop("number of rows genomdat.a=", nbrOfLoci, " number of rows normaldat.b=", n6)
  }

  p5 <- ncol(normaldat.a)
  p6 <- ncol(normaldat.b)


  # Argument 'genomdat.hetmatrix':
  if (!is.null(genomdat.hetmatrix)) {
    if (!is.logical(genomdat.hetmatrix)) {
      stop("genomdat.hetmatrix must be of type logical (TRUE/FALSE)")
    }
    if (is.vector(genomdat.hetmatrix)) {
      genomdat.hetmatrix <- as.matrix(genomdat.hetmatrix)
    }

    n7 <- nrow(genomdat.hetmatrix)
    p7 <- ncol(genomdat.hetmatrix)
    if (n7 != nbrOfLoci) {
      stop("number of rows genomdat.a=", nbrOfLoci, " number of rows genomdat.hetmatrix=", n7)
    }
    if (p1 != p7) {
      stop("number of rows genomdat.a=", nbrOfLoci, " number of columns genomdat.hetmatrix=", p7)
    }
  }


  # Argument 'normaldat.hetmatrix':
  if (!is.null(normaldat.hetmatrix)) {
    if (!is.logical(normaldat.hetmatrix)) {
      stop("normaldat.hetmatrix must be of type logical (TRUE/FALSE)")
    }
    if (is.vector(normaldat.hetmatrix)) {
      normaldat.hetmatrix <- as.matrix(normaldat.hetmatrix)
    }
    n8 <- nrow(normaldat.hetmatrix)
    p8 <- nrow(normaldat.hetmatrix)
    if (n8 != nbrOfLoci) {
      stop("number of rows normal.a=", nbrOfLoci, 
           " number of rows normaldat.hetmatrix=", n8)
    }
  }


  # Argument 'sampleid':
  if (length(sampleid) != ncol(genomdat.a)) {
    if (!is.null(sampleid)) {
      warning("length(sampleid) and ncol(genomdat.a) differ, names ignored\n")
    }
    sampleid <- paste("Sample", 1:ncol(genomdat.a))
  }




  if (sum(is.na(chrom) | is.na(maploc)) > 0) {
    warning("markers with missing chrom and/or maploc removed\n")
  }

  sortindex <- order(chrom, maploc, na.last=NA)
  zzz <- vector("list", length=8)
  zzz[[1]] <- as.matrix(genomdat.a[sortindex,])
  zzz[[2]] <- as.matrix(genomdat.b[sortindex,])
  colnames(zzz[[1]]) <- sampleid
  colnames(zzz[[2]]) <- sampleid
  zzz[[3]] <- chrom[sortindex]
  zzz[[4]] <- maploc[sortindex]
  if (!is.null(normaldat.a) && !is.null(normaldat.b)) {
    zzz[[5]] <- as.matrix(normaldat.a[sortindex,])
    zzz[[6]] <- as.matrix(normaldat.b[sortindex,])
    colnames(zzz[[5]]) <- NULL
    colnames(zzz[[6]]) <- NULL
    # colnames(zzz[[5]]) <- sampleid
    # colnames(zzz[[6]]) <- sampleid
  }

  if (!is.null(genomdat.hetmatrix)) {
    zzz[[7]] <- as.matrix(genomdat.hetmatrix[sortindex,])
    colnames(zzz[[7]]) <- sampleid
  }

  if (!is.null(normaldat.hetmatrix)) {
    zzz[[8]] <- as.matrix(normaldat.hetmatrix[sortindex,])
    colnames(zzz[[8]]) <- sampleid
  }

  class(zzz) <- c("psCNA", "list")
  names(zzz) <- c("genomdat.a", "genomdat.b", "chrom", "maploc", "normaldat.a", "normaldat.b", "genomdat.hetmatrix", "normaldat.hetmatrix")

  zzz
} # psCNA()




print.psCNA <- function(x, ...) {
  # Argument 'x':
  if (!inherits(x, 'psCNA')) {
    stop("First arg must be of class psCNA: ", class(x)[1]);
  }

  cat("Number of Test Samples", ncol(x[[1]]),
      "Number of Reference Samples", ncol(x[[5]]),
      "\nNumber of Probes ", nrow(x[[1]]), "\n")
} # print.psCNA()
