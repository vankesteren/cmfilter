.onAttach <- function(libname, pkgname) {
  if (!checkOMP()) warning("OpenMP not available. For full performance, please compile for your platform with OpenMP flags enabled or (slower) pass a function (e.g., cmfilter:::prodCoef) to the decisionFunction arg in cmf().")
}