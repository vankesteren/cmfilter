#' R-based Coordinate-wise Mediation Filter backend
#'
#' @keywords internal
cmfslow <- function(x, M, y, decisionFunction, nStarts, nCores, cutoff, maxIter, 
                    stableLag, ...) {
  # Scale the variables
  x <- scale(x)
  M <- as.data.frame(scale(M))
  y <- scale(y)
  
  if (nCores > 1) {
    
    # Register parallel cluster
    clus <- parallel::makeCluster(nCores)
    on.exit(parallel::stopCluster(clus))
    
    # initialise the cluster with the current environment
    parallel::clusterEvalQ(clus, {library(cmfilter)})
    vars <- ls()
    parallel::clusterExport(clus, vars[vars != "clus"], parent.frame())
    
    out <- pbapply::pblapply(X = 1:nStarts, FUN = function(start) {
      
      # initialise the variables to be used
      meds <- names(M)
      len <- length(meds)
      mselSum <- msel <- integer(len)
      medsamp <- 1:len
      names(mselSum) <- names(msel) <- meds
      mselLags <- matrix(rep(1L, len*(stableLag + 1)),
                         ncol = stableLag + 1)
      
      # init with half the available degrees of freedom or half the mediators
      msel[1:len] <- generateStart(len, min(floor(nrow(M)/2), floor(len/2)))
      
      # select a random number of variables to consider
      # see James, Witten, Hastie & Tibshirani ISLR (p. 319)
      sub <- ceiling(sqrt(len))
      
      for (i in 1:maxIter) {
        medsamp <- sample(1:len, sub)
        msel <- cmfStep(x, M, y, decisionFunction, msel, medsamp, ...)
        
        # Check for convergence at lags
        mselLags <- cbind(mselLags[,-1], msel)
        converged <- all(rowSums(mselLags) %% (stableLag + 1) == 0) # Thanks Oisin!
        if (converged) break
      }
      
      return(msel)
    }, cl = clus)
    
  } else {
    
    out <- pbapply::pblapply(X = 1:nStarts, FUN = function(start) {
      
      # initialise the variables to be used
      meds <- names(M)
      len <- length(meds)
      mselSum <- msel <- integer(len)
      medsamp <- 1:len
      names(mselSum) <- names(msel) <- meds
      mselLags <- matrix(rep(1L, len*(stableLag + 1)),
                         ncol = stableLag + 1)
      
      # init with half the available degrees of freedom or half the mediators
      msel[1:len] <- generateStart(len, min(floor(nrow(M)/2), floor(len/2)))
      
      # select a random number of variables to consider
      # see James, Witten, Hastie & Tibshirani ISLR (p. 319)
      sub <- ceiling(sqrt(len))
      
      for (i in 1:maxIter) {
        medsamp <- sample(1:len, sub)
        msel <- cmfStep(x, M, y, decisionFunction, msel, medsamp, ...)
        
        # Check for convergence at lags
        mselLags <- cbind(mselLags[,-1], msel)
        converged <- all(rowSums(mselLags) %% (stableLag + 1) == 0) # Thanks Oisin!
        if (converged) break
      }
      
      return(msel)
    })
    
  }
  
  return(out)
}

#' One step the cmf function (internal)
#'
#' @param x a numeric vector, standardised exogenous variable
#' @param M a data frame with column names, potential mediators
#' @param y a numeric vector, outcome variable
#' @param decisionFunction a function with as inputs x, m, y, parameters,
#' and as output a TRUE (include) or FALSE (exclude) statement
#' @param msel binary vector of mediator selections at step i
#' @param medsamp which variables to consider
#' @param ... parameters passed to decisionFunction
#'
#' @return binary vector of mediator selections at step i+1
#'
#' @keywords internal
cmfStep <- function(x, M, y, decisionFunction, msel, medsamp, ...) {
  n <- length(x)
  currentsum <- sum(msel)
  M <- as.matrix(M)
  
  for (med in medsamp) {
    currentsum <- currentsum - msel[med]
    if (currentsum >= n) {
      # now we cannot calculate residuals
      next
    }
    
    # get mediator
    m <- M[,med]
    
    # create model matrix of included mediators
    mmsel <- msel == 1
    mmsel[med] <- FALSE
    Mx <- M[, mmsel]
    
    if (length(Mx) > 1) {
      # calculate residuals wrt this model matrix
      cp <- crossprod(Mx)
      # residual of X
      r_x <- x - Mx %*% solve(cp, crossprod(Mx, x))
      # residual of Y
      r_y <- y - Mx %*% solve(cp, crossprod(Mx, y))
    } else {
      # no mediators selected
      r_x <- x
      r_y <- y
    }
    
    # perform decision function
    if (decisionFunction(r_x, m, r_y, ...)) {
      msel[med] <- 1L
      currentsum <- currentsum + 1
    } else {
      msel[med] <- 0L
    }
  }
  msel
}


#' Generate a (random) start
#'
#' @param n total number of mediators
#' @param m number of selected mediators
#'
#' @return binary vector of mediator selections
#'
#' @keywords internal
generateStart <- function(n, m) {
  out <- integer(n)
  out[sample(1:n, m, replace = FALSE)] <- 1L
  return(out)
}
