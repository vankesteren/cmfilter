#' Methods for cmf results lists
#'
#' Plotting, summarising, and printing cmf objects
#'
#' @param x cmf object
#' @param object cmf object
#' @param removeZeros whether to remove the unselected mediators from the plot
#' @param line whether to show a line at the chosen cutoff
#' @param las las argument to barplot
#' @param ylim y limits argument to barplot
#' @param topn only show the top n mediators
#' @param border the colour of the border around the bars
#' @param space the amount of space between the bars (default 0)
#' @param ... other arguments passed to barplot and summary
#'
#' @seealso \code{\link{cmf}}
#'
#' @importFrom graphics barplot abline
#'
#' @name cmf-methods
#' @method plot cmf
#'
#' @export
plot.cmf <- function(x, removeZeros = FALSE,
                     line = TRUE, las = 2, ylim = c(0, 1),
                     space = 0, border = "dark grey", ...) {
  sp <- x$selectionRate
  if (removeZeros) sp <- sp[sp != 0]
  co <- as.list(x$call)$cutoff
  if (is.null(co)) co <- 0.5
  barplot(sp, las = las, ylim = ylim, space = space,
          border = ,
          col = ifelse(sp < co, "grey", "#888888"), ...)
  if (line) abline(h = co, lty = 3)
}

#' @rdname cmf-methods
#'
#' @importFrom stats screeplot
#' @importFrom graphics axis
#'
#' @method screeplot cmf
#' @export
screeplot.cmf <- function(x, topn, border = NA, space = 0, ...) {
  if (missing(topn)) topn <- length(x$selectionRate)
  barplot(x$selectionRate[order(x$selectionRate, decreasing = TRUE)][1:topn],
          border = border, space = 0, xaxt = "n", yaxt = "n", ...)
  axis(2, pretty(x$selectionRate, n = 20), las = 1)
}


#' @rdname cmf-methods
#'
#' @method summary cmf
#' @export
summary.cmf <- function(object, ...) {
  cat("\nCMF Algorithm Results\n\n")
  cat("----------------------\n\n")
  cat("call:\n")
  print(object$call)
  cat("\n")
  ns <- sum(object$selection)
  cat("variables selected:", ns, "\n")

  rescall <- as.list(object$call)
  pars <- names(rescall)
  ns <- rescall$nStarts
  if (is.null(ns)) ns <- 1000
  cat("number of starts:", ns, "\n")
  co <- rescall$cutoff
  if (is.null(co)) co <- 0.5
  cat("cutoff probability:", co, "\n")

  apars <- pars[!pars %in% names(as.list(args(cmfilter::cmf)))]
  if (length(apars) > 0) {
    cat("\nDecision function parameters:\n")
    sapply(apars, function(par) {
      cat("  ", par, ": ",
          as.character(as.list(object$call)[[par]]), "\n",
          sep = "")
    })
  }

  cat("\n----------------------\n\n")
  m <- data.frame(SelectionRate = object$selectionRate,
                  Selected = as.logical(object$selection))

  if (nrow(m) >= 10) {
    cat("Top 10:\n")
    print(m[order(m$SelectionRate, decreasing = TRUE),][1:10,], digits = 3)
  } else {
    cat("Top ", nrow(m), ":\n", sep = "")
    print(m[order(m$SelectionRate, decreasing = TRUE),], digits = 3)
  }

  cat("\n----------------------\n\n")
}

#' @rdname cmf-methods
#'
#' @method print cmf
#' @export
print.cmf <- function(x, ...) {
  summary(x, ...)
}

#' Add samples to existing cmf results object
#'
#' This function adds additional samples to the existing results object
#'
#' @param object a cmf object
#' @param nStarts the number of starts to add (default 100)
#' @param ... not used
#'
#' @examples # generate some data
#' dat <- generateMed(a = (1:10)/20, b = (1:10)/20)
#' res <- cmf(dat)
#' # double the samples
#' res <- update(res, 1000)
#'
#' @method update cmf
#' @export
update.cmf <- function(object, nStarts = 100, ...) {
  newcall <- object$call
  newcall$nStarts <- nStarts
  
  
  # evaluate the new call in the parent environment
  newres <- eval.parent(newcall)
  
  if (is.null(object$call$nStarts)) {
    oldn <- 1e3
  } else {
    oldn <- object$call$nStarts
  }
  
  if (is.null(object$call$cutoff)) {
    co <- .5
  } else {
    co <- object$call$cutoff
  }
  
  object$selectionRate <-
    (object$selectionRate*oldn + newres$selectionRate*nStarts) / (oldn+nStarts)
  
  object$call$nStarts <- oldn+nStarts
  
  object$selection <- object$selectionRate > co
  
  object
}

#' Set the cutoff for mediator selection
#' 
#' This function sets the cutoff value on a cmf object for mediator selection. 
#' Any cutoff value between 0 and 1 is allowed, where potential mediators with 
#' empirical selection probability (selection rate) above the cutoff will be 
#' considered mediators and the others will not. The cutoff can be entered 
#' manually or, when set to "mc", be based on a monte carlo simulation. See 
#' "details"
#' 
#' @param object a cmf object
#' @param cutoff either a number between 0 and 1 or "mc" - see details
#' 
#' @return a cmf object with updated cutoff value
#' 
#' @details The monte carlo determination is based on a procedure in 
#' PCA and Factor Analysis called "Parallel Analysis". We generate data from
#' the null hypothesis with the same dimensionality as the original dataset and
#' perform the algorithm 100 times. This creates a distribution of nonmediator
#' selection rates from which we can determine the cutoff (the 99.9th 
#' percentile)
#' 
#' 
#' @export
setCutoff <- function(object, cutoff = .5) {
  if (cutoff == "mc") {
    # monte carlo determination of cutoff
    M <- eval(object$call$M)
    x <- eval(object$call$x)
    dec <- object$call$decisionFunction
    if (is.null(dec)) dec <- "prodcoef"
    if (is.function(dec)) { 
      dec <- "prodcoef"
      warning("Parallel analysis inaccurate with nonstandard decisionFunction.")
    }
    if (is.null(M)) {
      p <- ncol(x) - 2 
      n <- nrow(x)
    } else {
      p <- ncol(M)
      n <- nrow(M)
    }
    parll <- as.vector(pbsapply(1:100, function(i) {
      cmf(generateMed(n, numeric(p), numeric(p)), 
          decisionFunction = dec, 
          pb = FALSE)$selectionRate
    }))
    return(setCutoff(object, quantile(parll, 0.999)))
  }
  if (!is.numeric(cutoff) || cutoff > 1 || cutoff <= 0)
    stop("Input cutoff between 0 and 1")
  object$call$cutoff <- cutoff
  object$selection <- object$selectionRate > cutoff
  object
}

