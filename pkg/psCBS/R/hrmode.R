hrmode <- function(x) {
  n <- length(x)
  ox <- sort(x)
  zzz <- .Fortran("hrmode",
                  n = as.integer(n),
                  x = as.double(ox),
                  hrm = double(1),
                  PACKAGE="psCBS")
  zzz$hrm
}
