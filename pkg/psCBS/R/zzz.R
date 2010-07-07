.onAttach <- function(libname, pkgname) {
  library.dynam(pkgname, pkgname, libname, now=FALSE);
 
  pd <- utils::packageDescription(pkgname);
  packageStartupMessage(pkgname, " v", pd$Version, " (", 
    pd$Date, ") successfully loaded. See ?", pkgname, " for help."); 
}
