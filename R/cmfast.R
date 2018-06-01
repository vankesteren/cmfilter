#' Fast Coordinate-wise Mediation Filter backend
#'
#' @importFrom Rcpp evalCpp
#' @importFrom stats qnorm qt
#' @useDynLib cmfilter, .registration=TRUE
#' @keywords internal
cmfast <- function(x, M, y, decisionFunction = "prodcoef", nStarts, nCores, 
                   cutoff, maxIter, stableLag, p.value, pb) {
  
  # Scale the variables
  x <- as.vector(scale(x))
  M <- as.matrix(scale(M))
  y <- as.vector(scale(y))
  
  # decision function + critical value
  if (decisionFunction == "causalsteps") {
    dint <- 0
    cval <- qt(p.value/2, length(x) - 1, lower.tail = FALSE)
  } else if (decisionFunction == "prodcoef") {
    dint <- 1
    cval <- qnorm(p.value/2, lower.tail = FALSE)
  } else {
    stop("Input either 'causalsteps' or 'prodcoef' as decisionFunction.")
  }
  
  # Generate output
  return(arma_cmf(
    x         = x, 
    M         = M, 
    y         = y, 
    maxIter   = maxIter, 
    stableLag = stableLag, 
    critval   = cval,
    decfun    = dint,
    nCores    = nCores, 
    nStarts   = nStarts, 
    pb        = pb
  ))
}