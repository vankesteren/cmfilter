#' Coordinate-wise Mediation Filter
#'
#' This function performs CMF on a set of potential mediators, given an input
#' and an output variable. It selects those mediators that are deemed relevant
#' by a default or a user-defined decision function, *conditional* on the other
#' mediators in the model. By doing this cyclically, the algorithm converges to
#' the best fitting subset of mediators.
#'
#' @param x exogenous variable; numeric vector
#' @param M potential mediators; data frame with column names
#' @param y outcome variable; numeric vector
#'
#' @param decisionFunction either a function with as inputs x, m, y, parameters,
#' and as output a TRUE (include) or FALSE (exclude) statement or a string
#' indicating the built-in decision function to use (see details)
#'
#' @param maxIter the maximum number of iterations for each start
#' @param stableLag how long does the selection need to be stable before
#' deciding upon convergence
#'
#' @param nStarts how many times to start the algorithm, combine this with the
#' randomStart parameter
#' @param randomStart either FALSE for a start with no variables selected, or
#' the number of mediators to be randomly drawn for the initiations. TRUE means
#' half of the available degrees of freedom will be spent or half the mediators
#' will be selected; whichever is lowest.
#' @param randomOrder whether the order of the mediators is randomised before
#' the start of each iteration
#' @param cutoff a cutoff value for selection: variables are selected if they
#' display a selection rate higher than this value. Only relevant when multiple
#' starts are specified.
#'
#' @param parallel whether to run the multiple starts in parallel
#' @param nCores how many threads (cores) to use for parallel processing
#'
#' @param verbose whether to display iteration-based output and information.
#' Does not work with parallel processing.
#' @param progressBar whether to display a progress bar. Only relevant when
#' multiple starts are specified.
#'
#' @param ... parameters passed to decisionFunction
#'
#' @details Available decision functions:
#'  - corMinusPartCor p.value
#'  - prodCoef
#'
#' @return a list of class *cmf*. See \code{\link{cmf-methods}}
#'
#' @export

cmf <- function(x, M, y, decisionFunction = "corMinusPartCor",
                maxIter = 10, stableLag = 1,
                nStarts = 1, randomStart = TRUE, randomOrder = FALSE,
                cutoff = 0.5, parallel = FALSE, nCores = 2,
                verbose = FALSE, progressBar = FALSE, ...) {

  # input checking
  if (randomStart > nrow(M))
    stop("randomStart > nrow(M) /// Don't start with more mediators than cases")

  # registering decision function
  decisionFunction <- registerDecFun(decisionFunction)

  # Parallel processing: register a custom lapply function
  lapplyfun <- registerLapplyFun(parallel, progressBar, nCores)

  # cannot be verbose if parallel or progress bar
  if (parallel || progressBar) verbose <- FALSE

  # Verbose stuff: register printing function
  if (verbose) {
    w <- floor(as.numeric(options("width"))/2)
    catmsel <- registerCatFun(w, ncol(M))
  }

  # Scale the variables
  x <- scale(x)
  M <- as.data.frame(scale(M))
  y <- scale(y)

  # Generate output
  out <- lapplyfun(X = 1:nStarts, FUN = function(start) {

    # initialise the variables to be used
    meds <- names(M)
    len <- length(meds)
    mselSum <- msel <- integer(len)
    medsamp <- 1:len
    names(mselSum) <- names(msel) <- meds
    mselLags <- matrix(rep(1L, len*(stableLag + 1)),
                       ncol = stableLag + 1)

    if (randomStart > 1) {
      msel[1:len] <- generateStart(len, ceiling(randomStart))
    } else if (randomStart) {
      # init with half the available degrees of freedom or half the mediators
      msel[1:len] <- generateStart(len, min(floor(nrow(M)/2), floor(len/2)))
    }

    if (verbose) {
      cat("\nCMF Algorithm\n\n-----------\n\n")
      catmsel(msel, w)
    }

    for (i in 1:maxIter) {
      if (randomOrder) {
        medsamp <- sample(medsamp)
      }
      msel <- cmfStep(x, M, y, decisionFunction, msel, medsamp, ...)
      mselSum <- mselSum + msel

      if (verbose) catmsel(msel, w)

      # Check for convergence at lags
      mselLags <- cbind(mselLags[,-1], msel)
      converged <- all(apply(mselLags, 1, function(x) (all(x) | !any(x))))
      if (converged) {
        if (verbose) cat("\nAlgorithm converged",
                         "\n\n-----------\n\n")
        break
      }
    }

    if (i == maxIter && !converged) {
      if (nStarts == 1) warning("Maximum iteration reached, nonconvergence")
      if (verbose) cat("\nAlgorithm did not converge",
                       "\n\n-----------\n\n")
    }

    return(list(
      selection = msel,
      selSums = mselSum,
      iterations = i,
      converged = converged
    ))
  })

  if (nStarts == 1) {
    res <- list(
      call = match.call(),
      selection = out[[1]][["selection"]],
      selectionRate = out[[1]][["selSums"]] / out[[1]][["iterations"]],
      converged = out[[1]][["converged"]]
    )
  } else {
    selSum <- rowSums(sapply(out, function(l) {
      if (l[["converged"]]) {
        # upon convergence, we repeat the selection with the number of iters to
        # weigh the outcome properly
        return(l[["selection"]]*maxIter)
      } else {
        # otherwise, we use the number of occurrences of a mediator as a proxy
        # for the inclusion probability
        return(l[["selSums"]])
      }
    }))
    selRate <- selSum / (nStarts * maxIter)
    sumConverge <- sum(sapply(out, function(l) l[["converged"]]))
    res <- list(
      call = match.call(),
      selection = as.numeric(selRate > cutoff),
      selectionRate = selRate,
      converged = sumConverge / nStarts
    )
  }

  return(structure(res, class = "cmf"))
}


#' One step the cmf function (internal)
#'
#' @param x a numeric vector, standardised exogenous variable
#' @param M a data frame with column names, potential mediators
#' @param y a numeric vector, outcome variable
#' @param decisionFunction a function with as inputs x, m, y, parameters,
#' and as output a TRUE (include) or FALSE (exclude) statement
#' @param msel binary vector of mediator selections at step i
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
      xres <- x - Mx %*% solve(cp, crossprod(Mx, x))
      # residual of Y
      yres <- y - Mx %*% solve(cp, crossprod(Mx, y))
    } else {
      # no mediators selected
      xres <- x
      yres <- y
    }

    # perform decision function
    if (decisionFunction(xres, m, yres, ...)) {
      msel[med] <- 1L
      currentsum <- currentsum + 1
    } else {
      msel[med] <- 0L
    }
  }
  msel
}
