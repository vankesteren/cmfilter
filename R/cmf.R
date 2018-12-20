#' Coordinate-wise Mediation Filter
#'
#' This function performs CMF on a set of potential mediators, given an input
#' and an output variable. It selects those mediators that are deemed relevant
#' by a default or a user-defined decision function, *conditional* on the other
#' mediators in the model. By doing this cyclically, and with multiple random
#' starts, the algorithm outputs an estimate of the best mediators.
#'
#' @param x exogenous variable; numeric vector or data frame with x, y, and at
#' least one M column
#' @param M potential mediators; data frame with column names
#' @param y outcome variable; numeric vector
#'
#' @param decisionFunction either a function with as inputs x, m, y, parameters,
#' and as output a TRUE (include) or FALSE (exclude) statement or a string
#' indicating the built-in decision function to use (see details)
#'
#' @param nStarts how many times to start the algorithm
#' @param nCores how many threads (cores) to use for parallel processing
#' 
#' @param cutoff a cutoff value for selection: variables are selected if they
#' display a selection rate higher than this value. Only relevant when multiple
#' starts are specified. Can also be specified post-hoc using
#' \code{\link{setCutoff}}.
#' @param maxIter the maximum number of iterations for each start
#' @param stableLag how long does the selection need to be stable before
#' deciding the algorithm has converged
#' 
#' @param pb Whether to display a progress bar (default TRUE). Only available 
#' with built-in decision functions
#' 
#' @param ... parameters passed to decisionFunction
#'
#' @details Available decision functions. These functions are implemented in 
#' C++ to speed up computation. Between brackets the additional parameter that 
#' may be passed to the function in the \code{...} argument of this function. 
#' (\code{arguments} = \code{defaultvalue}):
#'
#'  - \code{prodcoef} (\code{p.value} = 0.1): Test for the product of 
#'  coefficients, Sobel test
#'  
#'  - \code{causalsteps} (\code{p.value} = 0.1): Causal steps test min(Ta, Tb)
#'
#' @return an object of class \code{cmf}. See \code{\link{cmf-methods}}
#'
#' @export
cmf <- function(x, M, y, decisionFunction = "prodcoef",
                nStarts = 1e3, nCores = "auto", 
                cutoff = 0.5, maxIter = 25, stableLag = 5,
                pb = TRUE, ...) {

  if (nCores == "auto")  nCores <- parallel::detectCores()

  if (is.data.frame(x) && missing(M) && missing(y)) {
    d <- x
    dn <- colnames(d)
    if (all(c("x","y") %in% dn)) {
      x <- d$x
      y <- d$y
      M <- d[,!colnames(d) %in% c("x", "y")]
    } else {
      stop("Enter a data frame with x and y variables")
    }
  }
  
  # Perform CMF
  if (is.character(decisionFunction)) {
    # Perform C++ based CMF
    
    pval <- as.list(match.call())$p.value
    if (is.null(pval)) pval <- 0.1
    
    selRate <- cmfast(x, M, y, decisionFunction, nStarts, nCores, cutoff, 
                      maxIter, stableLag, pval, pb)
    selRate <- as.vector(selRate)
    names(selRate) <- colnames(M)
    res <- list(
      call = match.call(),
      selection = as.numeric(selRate > cutoff),
      selectionRate = selRate
    )
    
  } else if (is.function(decisionFunction)) {
    # Perform R-based CMF
    
    out <- cmfslow(x, M, y, decisionFunction, nStarts, nCores, 
                   cutoff, maxIter, stableLag, ...)

    if (nStarts == 1) {
      res <- list(
        call = match.call(),
        selection = out[[1]],
        selectionRate = out[[1]]
      )
    } else {
      selSums <- rowSums(data.frame(out))
      selRate <- selSums / nStarts
      res <- list(
        call = match.call(),
        selection = as.numeric(selRate > cutoff),
        selectionRate = selRate
      )
    }
  } else {
    stop("Input valid decisionFunction.")
  }
  
  
  # evaluate some call arguments
  res$call$nStarts <- nStarts
  res$call$nCores  <- nCores

  return(structure(res, class = "cmf"))
}