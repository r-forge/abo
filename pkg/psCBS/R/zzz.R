.onAttach <- function(libname, pkgname) {
  pi <- utils::packageDescription(pkgname);
  packageStartupMessage(pkgname, " v", pi$Version, " (", 
    pi$Date, ") successfully loaded. See ?", pkgname, " for help."); 
}
