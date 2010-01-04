print.as.CNA <- function(x, ...)
  {
    if (!inherits(x, 'as.CNA')) stop("First arg must be of class as.CNA")
    cat("Number of Test Samples", ncol(x[[1]]),
        "Number of Reference Samples", ncol(x[[5]]),
        "\nNumber of Probes ", nrow(x[[1]]),"\n")
  }



###########################################################################
# HISTORY:
# 2010-01-03 [HB]
# o Created from as.DNAcopyMethods.R.
########################################################################### 
