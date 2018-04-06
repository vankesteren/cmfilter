#' Methods for cmf results lists
#'
#' @param x cmf object
#' @param object cmf object
#' @param removeZeros whether to remove the unselected mediators from the plot
#' @param line whether to show a line at the chosen cutoff
#' @param las las argument to barplot
#' @param ylim y limits argument to barplot
#' @param cutoff new cutoff for mediator selections. If left out, use
#' changepoint detection to automatically select the number of variables.
#' @param ... other arguments passed to barplot and summary
#'
#' @seealso \code{\link{cmf}}
#'
#' @importFrom graphics barplot abline
#'
#' @name cmf-methods
#'
#' @export
plot.cmf <- function(x, removeZeros = FALSE,
                      line = TRUE, las = 2, ylim = c(0, 1), ...) {
  sp <- x$selectionRate
  if (removeZeros) sp <- sp[sp != 0]
  barplot(sp, las = las, ylim = ylim, ...)
  co <- as.list(x$call)$cutoff
  if (is.null(co)) co <- 0.5
  if (line) abline(h = co, lty = 3)
}

#' @rdname cmf-methods
#'
#' @export
summary.cmf <- function(object, ...) {
  cat("\nCMF Algorithm Results\n\n")
  cat("----------------------\n\n")
  cat("call:\n")
  print(object$call)
  cat("\n")
  ns <- sum(object$selection)
  if (object$converged) {
    cat("Algorithm converged. \n")
  } else {
    cat("Algorithm did not (always) converge. Try multiple starts with cutoff.")
    cat("\n\n")
  }
  cat("variables selected:", ns, "\n")

  rescall <- as.list(object$call)
  pars <- names(rescall)
  if ("nStarts" %in% pars) {
    nStarts <- rescall$nStarts
    cat("number of starts:", nStarts, "\n")
    co <- as.list(object$call)$cutoff
    if (is.null(co)) co <- 0.5
    cat("cutoff probability:", co, "\n")
  }

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
#' @export
print.cmf <- function(x, ...) {
  summary(x, ...)
}

#' @rdname cmf-methods
#'
#' @export
setCutoff <- function(object, cutoff) {
  if (!missing(cutoff)) {
    if (!is.numeric(cutoff) || cutoff > 1 || cutoff <= 0)
      stop("Input cutoff between 0 and 1")
    object$call$cutoff <- cutoff
    object$selection <- object$selectionRate > cutoff
  } else {
    if (requireNamespace("changepoint", quietly = TRUE)) {
      ordsel <- object$selectionRate[order(object$selectionRate,
                                           decreasing = TRUE)]
      cpt <- changepoint::cpt.var(ordsel, Q = 1)
      cptcutoff <- ordsel[cpt@cpts[1]] - 1e-12
      object <- setCutoff(object, cptcutoff)
    } else {
      stop("Package \"changepoint\" needed for automatic cutoff detection")
    }
  }
  object
}
