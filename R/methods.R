#' Methods for cmf results lists
#'
#' Plotting, summarising, and printing cmf objects
#'
#' @param x cmf object
#' @param object cmf object
#' @param select optional selection vector of variables to show
#' @param line whether to show a line at the chosen cutoff
#' @param labelSelected whether to label the selected mediators in the plot, not
#' used when fewer than 20 bars are shown.
#' @param defaultColour the colour for the bars lower than the cutoff
#' @param highlightColour the colour for the bars of the selected mediators
#' @param topn only show the top n mediators
#' @param ... other arguments passed to barplot and summary
#' 
#' @examples # generate some data
#' dat <- generateMed(a = (1:10)/20, b = (1:10)/20)
#' res <- cmf(dat)
#' # screeplot of the result
#' screeplot(res)
#' # manhattan style plot the result
#' plot(res)
#' 
#'
#' @seealso \code{\link{cmf}}
#'
#' @importFrom graphics barplot abline text
#'
#' @name cmf-methods
#' @method plot cmf
#'
#' @export
plot.cmf <- function(x, select, line = TRUE, labelSelected = TRUE, 
                     defaultColour = "#00008b", highlightColour = "#e2bd36", 
                     ...) {
  
  # Compile a list of arguments
  args <- match.call(expand.dots = FALSE)$`...`
  if (is.null(args))        args        <- list()
  if (is.null(args$las))    args$las    <- 2
  if (is.null(args$ylim))   args$ylim   <- c(0, 1)
  if (is.null(args$space))  args$space  <- 0
  if (is.null(args$border)) args$border <- NA
  
  if (missing(select)) select <- seq_along(x$selectionRate)
  sp <- x$selectionRate[select]
  
  co <- as.list(x$call)$cutoff
  if (is.null(co)) co <- 0.5
  
  if (is.null(args$names.arg)) {
    if (length(select) < 20) {
      args$names.arg <- names(sp) 
    } else {
      args$names.arg <- rep("", length(select))
    }
  }
  
  args$height <- sp
  args$col    <- ifelse(sp < co, defaultColour, highlightColour)
  
  do.call(barplot, as.list(args))
  
  if (line) abline(h = co, lty = 3)
  if (length(select) >= 20 && labelSelected && length(sp[sp > co] > 0))
    text(x = which(sp > co) - .5, y = sp[sp > co], labels = names(sp[sp > co]), 
         pos = 3)
}

#' @rdname cmf-methods
#'
#' @importFrom stats screeplot
#' @importFrom graphics axis
#'
#' @method screeplot cmf
#' @export
screeplot.cmf <- function(x, topn, ...) {
  # Compile args for barplot call
  args <- match.call(expand.dots = FALSE)$`...`
  if (is.null(args))        args        <- list()
  if (is.null(args$border)) args$border <- NA
  if (is.null(args$space))  args$space  <- 0
  if (is.null(args$col))    args$col    <- "#499293"
  
  if (missing(topn)) topn <- length(x$selectionRate)
  sr <- x$selectionRate
  args$height <- sr[order(sr, decreasing = TRUE)][1:topn]
  args$xaxt   <- "n"
  args$yaxt   <- "n"
  
  do.call(barplot, as.list(args))
  if (!is.null(args$cex.axis))
    axis(2, pretty(sr, n = 20), las = 1, cex.axis = args$cex.axis)
  else
    axis(2, pretty(sr, n = 20), las = 1, cex.axis = args$cex.axis)
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
#' res <- cmf(dat, nStarts = 200)
#' # double the samples
#' res <- update(res, 500)
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
    (object$selectionRate*oldn + newres$selectionRate*nStarts) / 
    (oldn + nStarts)
  
  object$call$nStarts <- oldn + nStarts
  
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
#' @examples # generate some data
#' dat <- generateMed(a = (1:10)/20, b = (1:10)/20)
#' res <- cmf(dat)
#' # set the cutoff for this result at 0.1
#' setCutoff(res, 0.1)
#' 
#' @details The monte carlo determination is based on a procedure in 
#' PCA and Factor Analysis called "Parallel Analysis". We generate data from
#' the null hypothesis with the same dimensionality as the original dataset and
#' perform the algorithm 100 times. This creates a distribution of nonmediator
#' selection rates from which we can determine the cutoff (the 99.9th 
#' percentile)
#' 
#' @importFrom stats quantile
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
    nIter <- ceiling(1000/p)
    if (nIter == 1) {
      mockdat <- matrix(rnorm(n*(p + 2)), n)
      colnames(mockdat) <- c("x", paste0("M.", 1:p), "y")
      mockdat <- as.data.frame(mockdat)
      parll <- cmf(mockdat, decisionFunction = dec)$selectionRate
    } else {
      parll <- as.vector(pbapply::pbsapply(1:nIter, function(i) {
        mockdat <- matrix(rnorm(n*(p + 2)), n)
        colnames(mockdat) <- c("x", paste0("M.", 1:p), "y")
        mockdat <- as.data.frame(mockdat)
        cmf(mockdat, decisionFunction = dec, pb = FALSE)$selectionRate
      }))
    }
    
    return(setCutoff(object, quantile(parll, 0.999)))
  }
  if (!is.numeric(cutoff) || cutoff > 1 || cutoff <= 0)
    stop("Input cutoff between 0 and 1")
  object$call$cutoff <- cutoff
  object$selection <- object$selectionRate > cutoff
  object
}

#' Combine the results of multiple cmf objects into one
#' 
#' This function combines two cmf objects and returns one cmf object with the
#' combined results. This helps with combining results done over multiple runs,
#' for example in high-performance computing.
#' 
#' @param x a cmf object
#' @param y a cmf object
#' 
#' @examples # generate some data
#' dat <- generateMed(a = (1:10)/20, b = (1:10)/20)
#' # create two different cmf objects on this data
#' res_1 <- cmf(dat, nStarts = 500)
#' res_2 <- cmf(dat, nStarts = 500)
#' # Combine the results using the + operator
#' res_1 + res_2
#' 
#' @return a cmf object with combined results
#'  
#' @method + cmf
#' @export
`+.cmf` <- function(x, y) {
  # Find the number of iterations of each of the objects
  if (is.null(x$call$nStarts)) xn <- 1e3 else xn <- x$call$nStarts
  if (is.null(y$call$nStarts)) yn <- 1e3 else yn <- y$call$nStarts
  
  # New selection rate is weighted average of selection rates.
  cmf_out <- x
  cmf_out$call$nStarts  <- tot <- xn + yn
  cmf_out$selectionRate <- (xn * x$selectionRate + yn * y$selectionRate) / tot
  
  # Recalculate cutoff and return
  if (is.null(x$call$cutoff)) co <- .5 else co <- x$call$cutoff
  setCutoff(cmf_out, co)
}
