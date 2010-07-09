psCNA <- function(genomdat.a, genomdat.b, chrom, maploc, normaldat.a, normaldat.b, genomdat.hetmatrix=NULL, normaldat.hetmatrix=NULL, sampleid=NULL, arraytype=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'genomdat.a':
  if (!is.numeric(genomdat.a)) {
    stop("genomdat.a must be numeric: ", mode(genomdat.a));
  }
  if (is.vector(genomdat.a)) {
    genomdat.a <- as.matrix(genomdat.a);
  }
  nbrOfLoci <- nrow(genomdat.a);

  # Argument 'genomdat.b':
  if (!is.numeric(genomdat.b)) {
    stop("genomdat.b must be numeric: ", mode(genomdat.b));
  }
  if (is.vector(genomdat.b)) {
    genomdat.b <- as.matrix(genomdat.b);
  }

  if (nrow(genomdat.b) != nbrOfLoci) {
    stop("number of rows genomdat.a=", nbrOfLoci, " number of rows genomdat.b=", nrow(genomdat.b));
  }

  p1 <- ncol(genomdat.a);
  p2 <- ncol(genomdat.b);
  if (p1 != p2) {
    stop("number of cols genomdat.a=", p1, " number of cols genomdat.b=", p2);
  }


  # Argument 'chrom':
  if (is.factor(chrom)) {
    chrom <- as.character(chrom);
  }
  if (length(chrom) != nbrOfLoci) {
    stop("number of rows genomdat.a=", nbrOfLoci, " length chrom=", length(chrom));
  }


  # Argument 'maploc':
  if (!is.numeric(maploc)) {
    stop("maploc must be numeric: ", mode(maploc));
  }
  if (length(maploc) != nbrOfLoci) {
    stop("number of rows genomdat.a=", nbrOfLoci, " length chrom=", length(maploc));
  }


  # Argument 'normaldat.a':
  if (is.vector(normaldat.a)) {
    normaldat.a <- as.matrix(normaldat.a);
  }
  if (nrow(normaldat.a) != nbrOfLoci) {
    stop("number of rows genomdat.a=", nbrOfLoci, " number of rows normaldat.a=", nrow(normaldat.a));
  }

  # Argument 'normaldat.b':
  if (is.vector(normaldat.b)) {
    normaldat.b <- as.matrix(normaldat.b);
  }
  if (nrow(normaldat.b) != nbrOfLoci) {
    stop("number of rows genomdat.a=", nbrOfLoci, " number of rows normaldat.b=", nrow(normaldat.b));
  }

  p5 <- ncol(normaldat.a);
  p6 <- ncol(normaldat.b);


  # Argument 'genomdat.hetmatrix':
  if (!is.null(genomdat.hetmatrix)) {
    if (!is.logical(genomdat.hetmatrix)) {
      stop("genomdat.hetmatrix must be of type logical (TRUE/FALSE)");
    }
    if (is.vector(genomdat.hetmatrix)) {
      genomdat.hetmatrix <- as.matrix(genomdat.hetmatrix);
    }

    if (nrow(genomdat.hetmatrix) != nbrOfLoci) {
      stop("number of rows genomdat.a=", nbrOfLoci, " number of rows genomdat.hetmatrix=", nrow(genomdat.hetmatrix));
    }
    # Seems wrong? /HB 2010-07-08
    if (ncol(genomdat.hetmatrix) != p1) {
      stop("number of rows genomdat.a=", nbrOfLoci, " number of columns genomdat.hetmatrix=", ncol(genomdat.hetmatrix));
    }
  }


  # Argument 'normaldat.hetmatrix':
  if (!is.null(normaldat.hetmatrix)) {
    if (!is.logical(normaldat.hetmatrix)) {
      stop("normaldat.hetmatrix must be of type logical (TRUE/FALSE)");
    }
    if (is.vector(normaldat.hetmatrix)) {
      normaldat.hetmatrix <- as.matrix(normaldat.hetmatrix);
    }
    p8 <- nrow(normaldat.hetmatrix);
    if (nrow(normaldat.hetmatrix) != nbrOfLoci) {
      stop("number of rows normal.a=", nbrOfLoci, 
           " number of rows normaldat.hetmatrix=", nrow(normaldat.hetmatrix));
    }
  }


  # Argument 'sampleid':
  if (length(sampleid) != ncol(genomdat.a)) {
    if (!is.null(sampleid)) {
      warning("length(sampleid) and ncol(genomdat.a) differ, names ignored\n");
    }
    sampleid <- paste("Sample", seq(length=ncol(genomdat.a)));
  }




  if (sum(is.na(chrom) | is.na(maploc)) > 0) {
    warning("markers with missing chrom and/or maploc removed\n");
  }

  # Allocate empty vector
  zzz <- vector("list", length=8);
  names(zzz) <- c("genomdat.a", "genomdat.b",
                  "chrom", "maploc", 
                  "normaldat.a", "normaldat.b", 
                  "genomdat.hetmatrix", "normaldat.hetmatrix");


  sortindex <- order(chrom, maploc, na.last=NA);
  zzz$genomdat.a <- as.matrix(genomdat.a[sortindex,]);
  zzz$genomdat.b <- as.matrix(genomdat.b[sortindex,]);
  colnames(zzz$genomdat.a) <- sampleid;
  colnames(zzz$genomdat.b) <- sampleid;
  zzz$chrom <- chrom[sortindex];
  zzz$maploc <- maploc[sortindex];
  if (!is.null(normaldat.a) && !is.null(normaldat.b)) {
    zzz$normaldat.a <- as.matrix(normaldat.a[sortindex,]);
    zzz$normaldat.b <- as.matrix(normaldat.b[sortindex,]);
    colnames(zzz$normaldat.a) <- NULL;
    colnames(zzz$normaldat.b) <- NULL;
  }

  if (!is.null(genomdat.hetmatrix)) {
    zzz$genomdat.hetmatrix <- as.matrix(genomdat.hetmatrix[sortindex,]);
    colnames(zzz$genomdat.hetmatrix) <- sampleid;
  }

  if (!is.null(normaldat.hetmatrix)) {
    zzz$normaldat.hetmatrix <- as.matrix(normaldat.hetmatrix[sortindex,]);
    colnames(zzz$normaldat.hetmatrix) <- sampleid;
  }

  class(zzz) <- c("psCNA", "list");

  zzz
} # psCNA()




setMethodS3("print", "psCNA", function(x, ...) {
  cat("Number of samples", nbrOfSamples(x),
      "\nNumber of reference samples", nbrOfReferenceSamples(x),
      "\nNumber of probes ", nbrOfLoci(x), "\n")
})


setMethodS3("nbrOfLoci", "psCNA", function(x, ...) {
  nrow(x$genomdat.a);
})

setMethodS3("nbrOfSamples", "psCNA", function(x, ...) {
  ncol(x$genomdat.a);
})

setMethodS3("getSampleNames", "psCNA", function(x, ...) {
  colnames(x$genomdat.a);
})

setMethodS3("nbrOfReferenceSamples", "psCNA", function(x, ...) {
  ncol(x$normaldat.a);
})



############################################################################
# HISTORY:
# 2010-07-09
# o BACKWARD COMPATIBILITY: Now psCNA() returns a list of length 8.
# 2010-07-08
# o ROBUSTNESS: Replaced all 1:n with seq(length=n) to deal with n == 0.
# o ROBUSTNESS: Now the constructor assign list elements by names.
# o CLEANUP: No need to validate first argument for S3 methods.
# o Added getSampleNames().
# o Added nbrOfLoci().
# o Added nbrOfSamples() and nbrOfReferenceSamples().
# o Now using setMethodS3() to define methods.
############################################################################
